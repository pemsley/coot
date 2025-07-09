
/*********************************************************
 *
 *  Copyright (C) 2014 by Vitaliy Vitsentiy
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *********************************************************/


#ifndef __ctpl_thread_pool_H__
#define __ctpl_thread_pool_H__

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <exception>
#include <future>
#include <mutex>
#include <boost/lockfree/queue.hpp>
#include <iostream> // Added for debugging output

#ifndef _ctplThreadPoolLength_
#define _ctplThreadPoolLength_  100
#endif


// thread pool to run user's functors with signature
//      ret func(int id, other_params)
// where id is the index of the thread that runs the functor
// ret is some return type


namespace ctpl {

    class thread_pool {

    public:

        thread_pool() : q(_ctplThreadPoolLength_) {
            this->init();
            std::cerr << "ctpl::thread_pool(): Default constructor called." << std::endl;
        }

        // Constructor with nThreads
        thread_pool(int nThreads, int queueSize = _ctplThreadPoolLength_) : q(queueSize) {
            this->init();
            std::cerr << "ctpl::thread_pool(" << nThreads << " threads, queueSize=" << queueSize << "): Constructor called." << std::endl;
            try {
                this->resize(nThreads);
                std::cerr << "ctpl::thread_pool: Resize successful." << std::endl;
            } catch (const std::system_error& e) {
                std::cerr << "Error creating threads during pool initialization (std::system_error): " << e.what() << std::endl;
                throw; // Re-throw to allow the calling code to handle it
            } catch (const std::exception& e) {
                std::cerr << "An unexpected error occurred during pool initialization (std::exception): " << e.what() << std::endl;
                throw;
            } catch (...) {
                std::cerr << "An unknown error occurred during pool initialization." << std::endl;
                throw;
            }
        }

        // the destructor waits for all the functions in the queue to be finished
        ~thread_pool() {
            std::cerr << "ctpl::thread_pool: Destructor called, stopping pool." << std::endl;
            this->stop(true);
            std::cerr << "ctpl::thread_pool: Pool stopped." << std::endl;
        }

        // get the number of running threads in the pool
        int size() { return static_cast<int>(this->threads.size()); }

        // number of idle threads
        int n_idle() { return this->nWaiting; }
        std::thread & get_thread(int i) { return *this->threads[i]; }

        // change the number of threads in the pool
        // should be called from one thread, otherwise be careful to not interleave, also with this->stop()
        // nThreads must be >= 0
        void resize(int nThreads) {
            std::cerr << "ctpl::thread_pool::resize(" << nThreads << "): called. Current size: " << this->threads.size() << std::endl;
            if (!this->isStop && !this->isDone) {
                int oldNThreads = static_cast<int>(this->threads.size());
                if (oldNThreads <= nThreads) {  // if the number of threads is increased
                    std::cerr << "ctpl::thread_pool::resize: Increasing threads from " << oldNThreads << " to " << nThreads << std::endl;
                    this->threads.resize(nThreads);
                    this->flags.resize(nThreads);

                    for (int i = oldNThreads; i < nThreads; ++i) {
                        this->flags[i] = std::make_shared<std::atomic<bool>>(false);
                        std::cerr << "ctpl::thread_pool::resize: Setting up thread " << i << std::endl;
                        this->set_thread(i);
                        std::cerr << "ctpl::thread_pool::resize: Thread " << i << " set up." << std::endl;
                    }
                }
                else {  // the number of threads is decreased
                    std::cerr << "ctpl::thread_pool::resize: Decreasing threads from " << oldNThreads << " to " << nThreads << std::endl;
                    for (int i = oldNThreads - 1; i >= nThreads; --i) {
                        std::cerr << "ctpl::thread_pool::resize: Detaching thread " << i << std::endl;
                        *this->flags[i] = true;  // this thread will finish
                        // Check if the unique_ptr is valid before detaching, though detach on empty thread is safe
                        if (this->threads[i] && this->threads[i]->joinable()) { // Added joinable check as well
                            this->threads[i]->detach();
                        } else {
                            std::cerr << "WARNING: Attempted to detach invalid or non-joinable thread at index " << i << std::endl;
                        }
                    }
                    {
                        std::unique_lock<std::mutex> lock(this->mutex);
                        this->cv.notify_all();
                    }
                    // This resize operation will trigger destructors of unique_ptrs being removed.
                    // If any std::thread object was malformed from an earlier failure, its destructor here might crash.
                    std::cerr << "ctpl::thread_pool::resize: Resizing threads vector to " << nThreads << std::endl;
                    this->threads.resize(nThreads);  
                    this->flags.resize(nThreads);  
                }
            }
            std::cerr << "ctpl::thread_pool::resize: finished." << std::endl;
        }

        // empty the queue
        void clear_queue() {
            std::function<void(int id)> * _f;
            while (this->q.pop(_f))
                delete _f;  // empty the queue
            std::cerr << "ctpl::thread_pool::clear_queue: Queue cleared." << std::endl;
        }

        // pops a functional wrapper to the original function
        std::function<void(int)> pop() {
            std::function<void(int id)> * _f = nullptr;
            this->q.pop(_f);
            std::unique_ptr<std::function<void(int id)>> func(_f);  // at return, delete the function even if an exception occurred
            
            std::function<void(int)> f;
            if (_f)
                f = *_f;
            return f;
        }


        // wait for all computing threads to finish and stop all threads
        // may be called asyncronously to not pause the calling thread while waiting
        // if isWait == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
        void stop(bool isWait = false) {
            std::cerr << "ctpl::thread_pool::stop(isWait=" << (isWait ? "true" : "false") << "): called." << std::endl;
            if (!isWait) {
                if (this->isStop) {
                    std::cerr << "ctpl::thread_pool::stop: Already marked as stopping." << std::endl;
                    return;
                }
                this->isStop = true;
                for (int i = 0, n = this->size(); i < n; ++i) {
                    std::cerr << "ctpl::thread_pool::stop: Commanding thread " << i << " to stop." << std::endl;
                    *this->flags[i] = true;  // command the threads to stop
                }
                this->clear_queue();  // empty the queue
            }
            else {
                if (this->isDone || this->isStop) {
                    std::cerr << "ctpl::thread_pool::stop: Already marked as done or stopping." << std::endl;
                    return;
                }
                this->isDone = true;  // give the waiting threads a command to finish
                std::cerr << "ctpl::thread_pool::stop: Marked as done." << std::endl;
            }
            {
                std::unique_lock<std::mutex> lock(this->mutex);
                this->cv.notify_all();  // stop all waiting threads
                std::cerr << "ctpl::thread_pool::stop: Notified all waiting threads." << std::endl;
            }
            for (int i = 0; i < static_cast<int>(this->threads.size()); ++i) {  // wait for the computing threads to finish
                if (this->threads[i] && this->threads[i]->joinable()) { // Check if valid and joinable
                    std::cerr << "ctpl::thread_pool::stop: Joining thread " << i << std::endl;
                    this->threads[i]->join();
                } else {
                    std::cerr << "ctpl::thread_pool::stop: Skipping non-joinable or invalid thread " << i << std::endl;
                }
            }
            // if there were no threads in the pool but some functors in the queue, the functors are not deleted by the threads
            // therefore delete them here
            this->clear_queue();
            this->threads.clear();
            this->flags.clear();
            std::cerr << "ctpl::thread_pool::stop: Finished cleanup." << std::endl;
        }

        template<typename F, typename... Rest>
        auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))> {
            auto pck = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
                std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
            );

            auto _f = new std::function<void(int id)>([pck](int id) {
                (*pck)(id);
            });
            this->q.push(_f);

            std::unique_lock<std::mutex> lock(this->mutex);
            this->cv.notify_one();

            return pck->get_future();
        }

        // run the user's function that excepts argument int - id of the running thread. returned value is templatized
        // operator returns std::future, where the user can get the result and rethrow the catched exceptins
        template<typename F>
        auto push(F && f) ->std::future<decltype(f(0))> {
            auto pck = std::make_shared<std::packaged_task<decltype(f(0))(int)>>(std::forward<F>(f));

            auto _f = new std::function<void(int id)>([pck](int id) {
                (*pck)(id);
            });
            this->q.push(_f);

            std::unique_lock<std::mutex> lock(this->mutex);
            this->cv.notify_one();

            return pck->get_future();
        }


    private:

        // deleted
        thread_pool(const thread_pool &);// = delete;
        thread_pool(thread_pool &&);// = delete;
        thread_pool & operator=(const thread_pool &);// = delete;
        thread_pool & operator=(thread_pool &&);// = delete;

        void set_thread(int i) {
            std::shared_ptr<std::atomic<bool>> flag(this->flags[i]);  // a copy of the shared ptr to the flag
            auto f = [this, i, flag/* a copy of the shared ptr to the flag */]() {
                std::atomic<bool> & _flag = *flag;
                std::function<void(int id)> * _f;
                bool isPop = this->q.pop(_f);
                while (true) {
                    while (isPop) {  // if there is anything in the queue
                        std::unique_ptr<std::function<void(int id)>> func(_f);  // at return, delete the function even if an exception occurred
                        (*_f)(i);

                        if (_flag)
                            return;  // the thread is wanted to stop, return even if the queue is not empty yet
                        else
                            isPop = this->q.pop(_f);
                    }

                    // the queue is empty here, wait for the next command
                    std::unique_lock<std::mutex> lock(this->mutex);
                    ++this->nWaiting;
                    this->cv.wait(lock, [this, &_f, &isPop, &_flag](){ isPop = this->q.pop(_f); return isPop || this->isDone || _flag; });
                    --this->nWaiting;

                    if (!isPop)
                        return;  // if the queue is empty and this->isDone == true or *flag then return
                }
            };
            // Potential crash point: if new std::thread(f) fails to construct properly
            // Original: this->threads[i].reset(new std::thread(f));
            try {
                std::cerr << "ctpl::set_thread(" << i << "): Attempting to create new std::thread." << std::endl;
                this->threads[i] = std::make_unique<std::thread>(f); // Using make_unique is slightly safer
                std::cerr << "ctpl::set_thread(" << i << "): std::thread created successfully." << std::endl;
            } catch (const std::system_error& e) {
                std::cerr << "ERROR: Failed to create thread " << i << " (std::system_error): " << e.what() << std::endl;
                // Leave unique_ptr as nullptr if creation failed.
                this->threads[i].reset();
                throw; // Re-throw to propagate to the constructor's catch block
            } catch (const std::exception& e) {
                std::cerr << "ERROR: Unexpected exception creating thread " << i << " (std::exception): " << e.what() << std::endl;
                this->threads[i].reset();
                throw;
            } catch (...) {
                std::cerr << "ERROR: Unknown exception creating thread " << i << std::endl;
                this->threads[i].reset();
                throw;
            }
        }

        void init() {
            this->nWaiting = 0;
            this->isStop = false;
            this->isDone = false;
            // No need to clear vectors here, they are managed by constructors/destructors
        }

        std::vector<std::unique_ptr<std::thread>> threads;
        std::vector<std::shared_ptr<std::atomic<bool>>> flags;
        mutable boost::lockfree::queue<std::function<void(int id)> *> q;
        std::atomic<bool> isDone;
        std::atomic<bool> isStop;
        std::atomic<int> nWaiting;  // how many threads are waiting

        std::mutex mutex;
        std::condition_variable cv;
    };

}


#endif // __ctpl_thread_pool_H__

