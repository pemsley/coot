/*
 * api/test-thread-limits.cc
 *
 * Does libcootapi respect set_max_number_of_threads()?
 *
 * This program counts every thread that the process creates (by interposing
 * pthread_create()) and, for a set of API calls, reports
 *
 *   pool     threads alive before the call (the container's thread pool)
 *   extra    most threads alive during the call, over and above "pool"
 *   created  number of threads created by the call
 *
 * for several values of N = set_max_number_of_threads(N). A call passes if
 * "extra" is not more than N: it did not start more threads than it was
 * allowed to. The pool itself should have N threads.
 *
 * Some functions that the API does not use (yet) are tested in the same way,
 * because they are in the libraries that the API is built from.
 *
 * The last test writes to the log from several threads at the same time, as
 * the threads of partition_map_by_chain() do. If the logger is not thread-safe
 * that test loses messages or crashes.
 *
 * This matters most for WebAssembly, where threads come from a fixed-size
 * pool of workers, and asking for more threads than are in the pool hangs
 * the program - but a thread limit that is not a limit is a bug anywhere.
 *
 * Usage:
 *
 *   MOORHEN_TEST_DATA_DIR=/path/to/test-data test-thread-limits [N ...]
 *
 * using the same data as test-molecules-container, of which these are needed:
 * moorhen-tutorial-structure-number-1.pdb, moorhen-tutorial-map-number-1.mtz
 * and LZA.cif.
 *
 * The exit status is the number of failed checks.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <pthread.h>

#include "utils/coot-utils.hh"
#include "utils/logging.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-map-heavy.hh"
#include "ideal/model-bond-deltas.hh"
#include "molecules-container.hh"

extern logging logger;

// ------------------------------------------------------------------------
//                     count the threads of this process
// ------------------------------------------------------------------------

#if defined(__EMSCRIPTEN__)
extern "C" int emscripten_builtin_pthread_create(pthread_t *, const pthread_attr_t *,
                                                 void *(*)(void *), void *);
#define HAVE_THREAD_CENSUS 1
#elif defined(__linux__)
#include <dlfcn.h>
#define HAVE_THREAD_CENSUS 1
#endif

namespace thread_census {

   std::atomic<int> n_alive(0);      // created and not yet finished
   std::atomic<int> n_alive_peak(0);
   std::atomic<int> n_created(0);

   struct start_info_t {
      void *(*start_routine)(void *);
      void *arg;
   };

   // the destructor runs when the thread function returns and when the thread
   // calls pthread_exit().
   struct departure_t {
      ~departure_t() { --n_alive; }
   };

   void *trampoline(void *p) {
      start_info_t si = *static_cast<start_info_t *>(p);
      delete static_cast<start_info_t *>(p);
      departure_t departure;
      return si.start_routine(si.arg);
   }

   // a thread is counted from the moment that it is asked for, so that the
   // peak does not depend on how quickly the new thread starts running.
   void arrival() {
      ++n_created;
      int n = ++n_alive;
      int p = n_alive_peak.load();
      while (n > p && ! n_alive_peak.compare_exchange_weak(p, n)) {}
   }

   // wait for the number of threads to stop changing (threads that were
   // detached, e.g. when a thread pool is shrunk, finish in their own time).
   int settle() {
      int n = n_alive.load();
      int n_same = 0;
      for (int i=0; i<400 && n_same<10; i++) {
         std::this_thread::sleep_for(std::chrono::milliseconds(5));
         int n_now = n_alive.load();
         if (n_now == n) n_same++; else n_same = 0;
         n = n_now;
      }
      return n;
   }
}

#ifdef HAVE_THREAD_CENSUS
extern "C" int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
                              void *(*start_routine)(void *), void *arg) {

#if defined(__EMSCRIPTEN__)
   auto real_pthread_create = emscripten_builtin_pthread_create;
#else
   typedef int (*pthread_create_t)(pthread_t *, const pthread_attr_t *, void *(*)(void *), void *);
   static pthread_create_t real_pthread_create =
      reinterpret_cast<pthread_create_t>(dlsym(RTLD_NEXT, "pthread_create"));
#endif

   thread_census::start_info_t *si = new thread_census::start_info_t{start_routine, arg};
   thread_census::arrival();
   int status = real_pthread_create(thread, attr, thread_census::trampoline, si);
   if (status != 0) {
      --thread_census::n_alive;
      delete si;
   }
   return status;
}
#endif

// ------------------------------------------------------------------------
//                                 the tests
// ------------------------------------------------------------------------

std::string
reference_data(const std::string &file) {
   char *env = getenv("MOORHEN_TEST_DATA_DIR");
   if (env)
      return coot::util::append_dir_file(env, file);
   return file;
}

class check_t {
public:
   unsigned int n_max;
   std::string name;
   int pool;
   int extra;
   int created;
   bool pass;
   std::string note;
};

std::vector<check_t> checks;

void report(const check_t &c) {
   std::cout << "THREADS:: N=" << std::setw(2) << c.n_max << " " << std::left << std::setw(46) << c.name << std::right
             << " pool " << std::setw(2) << c.pool
             << " extra " << std::setw(3) << c.extra
             << " created " << std::setw(4) << c.created
             << (c.pass ? "  pass" : "  FAIL");
   if (! c.note.empty()) std::cout << "  (" << c.note << ")";
   std::cout << std::endl;
}

// run f() and see how many threads it used.
void probe(unsigned int n_max, const std::string &name, const std::function<std::string()> &f) {

   check_t c;
   c.n_max = n_max;
   c.name  = name;
   c.pool  = thread_census::settle();
   thread_census::n_alive_peak = c.pool;
   int n_created_start = thread_census::n_created;
   c.note    = f();
   c.extra   = thread_census::n_alive_peak - c.pool;
   c.created = thread_census::n_created - n_created_start;
   c.pass    = (c.extra <= static_cast<int>(n_max));
   checks.push_back(c);
   report(c);
}

// the thread pool of a new container should not be bigger than the limit that
// has already been set.
void test_constructor(unsigned int n_max) {

   coot::set_max_number_of_threads(n_max);
   int n_before = thread_census::settle();
   {
      molecules_container_t mc(false);
      check_t c;
      c.n_max   = n_max;
      c.name    = "molecules_container_t constructor";
      c.pool    = thread_census::settle() - n_before;
      c.extra   = 0;
      c.created = c.pool;
      c.pass    = (c.pool <= static_cast<int>(n_max));
      c.note    = "pool is the size of the new thread pool";
      checks.push_back(c);
      report(c);
   }
   thread_census::settle();
}

void test_api_calls(unsigned int n_max) {

   molecules_container_t mc(false);
   mc.set_show_timings(false);
   mc.set_refinement_is_verbose(false);
   mc.set_make_backups(false);
   mc.set_max_number_of_threads(n_max);

   // is the thread pool the size that we asked for?
   {
      check_t c;
      c.n_max   = n_max;
      c.name    = "set_max_number_of_threads()";
      c.pool    = thread_census::settle();
      c.extra   = 0;
      c.created = 0;
      c.pass    = (c.pool == static_cast<int>(n_max));
      c.note    = "pool should be N";
      checks.push_back(c);
      report(c);
   }

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   if (! mc.is_valid_model_molecule(imol) || ! mc.is_valid_map_molecule(imol_map)) {
      std::cout << "ERROR:: failed to read the test data - is MOORHEN_TEST_DATA_DIR set?" << std::endl;
      exit(99);
   }
   mc.set_imol_refinement_map(imol_map);
   mc.fill_rotamer_probability_tables();

   auto ok = [] () { return std::string(); };

   // ------------------------------ maps

   mc.set_map_is_contoured_with_thread_pool(false);
   probe(n_max, "get_map_contours_mesh() (std::thread)",
         [&] { mc.get_map_contours_mesh(imol_map, 55, 10, 10, 20, 0.4); return ok(); });
   mc.set_map_is_contoured_with_thread_pool(true);
   probe(n_max, "get_map_contours_mesh() (thread pool)",
         [&] { mc.get_map_contours_mesh(imol_map, 55, 10, 10, 20, 0.4); return ok(); });

   // ------------------------------ validation

   probe(n_max, "density_correlation_analysis()",
         [&] { mc.density_correlation_analysis(imol, imol_map); return ok(); });
   probe(n_max, "mmrrcc()",
         [&] { mc.mmrrcc(imol, "A", imol_map); return ok(); });
   probe(n_max, "get_rotamer_dodecs_instanced()",
         [&] { mc.get_rotamer_dodecs_instanced(imol); return ok(); });
   probe(n_max, "all_molecule_contact_dots()",
         [&] { mc.all_molecule_contact_dots(imol, 1); return ok(); });
   probe(n_max, "get_molecular_representation_mesh() surface",
         [&] { mc.get_molecular_representation_mesh(imol, "//A/10-40", "Chain", "MolecularSurface", 2); return ok(); });

   // ------------------------------ refinement

   probe(n_max, "refine_residues_using_atom_cid()",
         [&] { mc.refine_residues_using_atom_cid(imol, "//A/100", "TRIPLE", 100); return ok(); });
   probe(n_max, "refine_residue_range()",
         [&] { mc.refine_residue_range(imol, "A", 30, 40, 100); return ok(); });
   probe(n_max, "minimize_energy()",
         [&] { mc.minimize_energy(imol, "//A/50-54", 50, false, 20.0, false, 0.0, true); return ok(); });

   // ------------------------------ ligands

   if (mc.import_cif_dictionary(reference_data("LZA.cif"), -999999) == 1) {
      int imol_ligand = mc.get_monomer_and_position_at("LZA", -999999, 0, 0, 0);
      probe(n_max, "get_ligand_validation_vs_dictionary()",
            [&] { mc.get_ligand_validation_vs_dictionary(imol_ligand, "//", true); return ok(); });
      probe(n_max, "fit_ligand_right_here() 9 conformers",
            [&] {
               std::vector<int> fits = mc.fit_ligand_right_here(imol, imol_map, imol_ligand,
                                                                61.15, 45.47, 36.40, 1.0, true, 9, 2);
               return fits.empty() ? std::string("no ligand was fitted - is this the right map?") : std::string();
            });
   } else {
      std::cout << "WARNING:: failed to read LZA.cif - the ligand tests were not run" << std::endl;
   }

   // ------------------------------ a map for each chain
   //
   // Split chain A so that we have a model of several chains.

   int imol_chains = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   const std::string new_chain_ids = "CDEFGHIJK";
   for (unsigned int i=0; i<new_chain_ids.size(); i++)
      mc.change_chain_id(imol_chains, "A", std::string(1, new_chain_ids[i]), true, 30 * i, 30 * i + 29);
   probe(n_max, "partition_map_by_chain() (11 chains)",
         [&] {
            std::vector<int> maps = mc.partition_map_by_chain(imol_map, imol_chains);
            return std::to_string(maps.size()) + " maps";
         });
}

// A limit of 0 means "use the default". It must not make a thread pool of 0 threads,
// because nothing that is given to such a pool is ever run.
void test_limit_of_zero() {

   int n_before = thread_census::settle();
   {
      molecules_container_t mc(false);
      mc.set_max_number_of_threads(0);
      check_t c;
      c.n_max   = 0;
      c.name    = "set_max_number_of_threads(0)";
      c.pool    = thread_census::settle() - n_before;
      c.extra   = 0;
      c.created = 0;
      c.pass    = (c.pool >= 1);
      c.note    = "pool must not be empty";
      checks.push_back(c);
      report(c);
   }
   thread_census::settle();
}

// Functions that the API does not use, called directly. There is no container
// here, so there is no thread pool: "pool" should be 0.
void test_library_calls(unsigned int n_max) {

   coot::set_max_number_of_threads(n_max);

   clipper::Xmap<float> xmap;
   bool map_is_filled = coot::util::map_fill_from_mtz(&xmap, reference_data("moorhen-tutorial-map-number-1.mtz"),
                                                      "FWT", "PHWT", "W", false);
   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::ERROR_CODE err = mol->ReadCoorFile(reference_data("moorhen-tutorial-structure-number-1.pdb").c_str());
   if (! map_is_filled || err != mmdb::Error_NoError) {
      std::cout << "ERROR:: failed to read the test data - is MOORHEN_TEST_DATA_DIR set?" << std::endl;
      exit(99);
   }

   auto ok = [] () { return std::string(); };

   int sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "A", 30, "", 60, "", "*", "*", "*", "*");
   probe(n_max, "util::make_nxmap_from_xmap()",
         [&] { clipper::NXmap<float> nxmap = coot::util::make_nxmap_from_xmap(xmap, mol, sel_hnd, 3.0); return ok(); });
   mol->DeleteSelection(sel_hnd);

   probe(n_max, "util::soi_variance::proc()",
         [&] {
            coot::util::soi_variance sv(xmap);
            sv.proc(0.5);
            std::remove("soi.map"); // it writes this
            return ok();
         });

   probe(n_max, "util::compare_structure_factors()",
         [&] {
            coot::util::compare_structure_factors(xmap, xmap);
            std::remove("compare-sfs.table"); // it writes this
            return ok();
         });

   coot::protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();
   probe(n_max, "model_bond_deltas::resolve()",
         [&] {
            coot::model_bond_deltas mbd(mol, coot::protein_geometry::IMOL_ENC_ANY, &geom);
            mbd.resolve();
            return ok();
         });

   delete mol;
}

// The logger is used from threads (for example those of partition_map_by_chain()),
// so it has to cope with being called from several threads at the same time.
// If it can't, then this either crashes or loses messages.
bool test_logger_from_many_threads() {

   const unsigned int n_threads = 8;
   const unsigned int n_messages_per_thread = 5000;
   std::size_t n_before = logger.get_log_history_from(0).size();
   std::vector<std::thread> threads;
   for (unsigned int i=0; i<n_threads; i++)
      threads.push_back(std::thread([n_messages_per_thread] (unsigned int i_thread) {
         for (unsigned int j=0; j<n_messages_per_thread; j++)
            logger.log(log_t::DEBUG, "test-thread-limits message from thread", static_cast<int>(i_thread));
      }, i));
   for (auto &thread : threads) thread.join();
   std::size_t n_after = logger.get_log_history_from(0).size();
   std::size_t n_expected = n_before + n_threads * n_messages_per_thread;
   bool pass = (n_after == n_expected);
   std::cout << "THREADS:: logger.log() from " << n_threads << " threads at once: " << (n_after - n_before)
             << " of " << (n_threads * n_messages_per_thread) << " messages were stored"
             << (pass ? "  pass" : "  FAIL") << std::endl;
   return pass;
}

int main(int argc, char **argv) {

#ifndef HAVE_THREAD_CENSUS
   std::cout << "test-thread-limits: threads can't be counted on this platform - nothing was tested" << std::endl;
   return 0;
#endif

   std::vector<unsigned int> limits;
   for (int i=1; i<argc; i++) {
      int n = atoi(argv[i]);
      if (n > 0) limits.push_back(n);
   }
   if (limits.empty())
      limits = {1, 2, 4};

   for (unsigned int n_max : limits) {
      test_constructor(n_max);
      test_api_calls(n_max);
      thread_census::settle();
      test_library_calls(n_max);
      thread_census::settle();
   }
   test_limit_of_zero();

   unsigned int n_fail = 0;
   std::cout << "\n------------------------- test-thread-limits summary -------------------------\n";
   for (const auto &c : checks) {
      report(c);
      if (! c.pass) n_fail++;
   }
   std::cout << "THREADS:: " << checks.size() << " checks, " << n_fail << " failed" << std::endl;

   // this is the last test because, if the logger is not thread-safe, it can crash.
   if (! test_logger_from_many_threads())
      n_fail++;

   return n_fail;
}
