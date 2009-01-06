#   remote_control.py
#   Copyright (C) 2008  Bernhard Lohkamp, University of York
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

global coot_listener_socket
coot_listener_socket = False

# Open a coot listener socket, return nothing much, but do set! 
# %coot-listener-socket.
# 
# Hmmm... why make a side-effect like that?  Why not return
# %coot-listener-socket so that the caller can set it?  There may be
# a reason...
# 
# And the reason is that I can then call
# coot-listener-idle-function-proc without having to use a c++
# variable.
#
# BL says:
# this almost does what it should!
# it breaks the scripting window however! Or for that matter all python
# scripting which is not comming from the socket! Baeh! May
# need to thread it after all, but how?!?
def open_coot_listener_socket(port_number, host_name):

    import socket
    global coot_listener_socket

    print "in open_coot_listener_socket port: %s host %s" %(port_number, host_name)

    soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    host_adress = "127.0.0.1"

    soc.connect((host_adress, port_number))

    print "Coot listener socket ready! from remote_control.py"
    soc.send("Coot listener socket ready!\n")
    coot_listener_socket = soc
    if soc:
        set_coot_listener_socket_state_internal(1)

    # threads?
    import thread
    status = thread.start_new_thread(coot_listener_idle_function_proc,())
    print "status trhead ", status

    #import threading
    #class MyThread(threading.Thread):
    #    def run(self):
    #        coot_listener_idle_function_proc()
    #MyThread().start()
    #print "thread count ", threading.activeCount()

def coot_listener_idle_function_procbbbbbb():

    #global coot_listener_socket
    msg = 'start'
    while msg != 'exit':
        data = coot_listener_socket.recv(1024)
        print 'Received Message: ', data
        #msg = raw_input('Enter a Message: ')
        coot_listener_socket.send('return value')
    coot_listener_socket.close()
#    import time
#    if (coot_listener_socket):
#        time.sleep(1000)
#    else:
#        print "coot_listener_idle_function_proc bad sock: ", coot_listener_socket


def coot_listener_error_handler(key, args=""):

    print "coot_listener_error_handler handling error in %s with args %s" %(key, args)


def coot_listener_idle_function_proc():
    global coot_listener_socket
    import time
    #print " in idle func"
    if (coot_listener_socket):
        import thread
        while (1):
            continue_qm = listen_coot_listener_socket(coot_listener_socket)
            if (not continue_qm):
                set_coot_listener_socket_state_internal(0)
                coot_listener_socket.shutdown(2)
                coot_listener_socket.close()
                coot_listener_socket = False
                print "server gone - listener thread ends", coot_listener_socket
                return
            else:
                time.sleep(0.01)
                return
    else:
        print "coot_listener_idle_func_proc bad sock: ", coot_listener_socket

def listen_coot_listener_socket(soc):

    end_transmission_string = "end"  # just 'end' for test

    def evaluate_character_list(string):
        print "received in evaluate",  string
        if (string == end_transmission_string):
            print "finish socket"
            soc.send("Closing session")
            return False
        try:
            ret = eval(string)
        except SyntaxError:
            if (string == "\n"): print "CR"
            try:
                exec string in globals()
                ret = None
            except:
                ret = "INFO:: cannot eval or exec given string: " + string
        except:
            ret = "Info:: input error"
        # can only send strings back in this way
        ret = "return value: " + str(ret)
        soc.send(ret)
        return True

    # main body
    #
    data = False
    soc.setblocking(0)
    try:
        #print "waiting to receive"
        data = soc.recv(1024)
        print "received", data
    except:
        pass
    if (not data):
        #print "nothing on the line..."
        return True
    else:
        res = evaluate_character_list(data)
        return res


#open_coot_listener_socket(50007, "bla")
#import time
#time.sleep(20)

#!/usr/bin/python
import socket
import sys

if len(sys.argv) < 3:
    print 'Usage: %s [hostname] [portnumber]' % sys.argv[0]
    sys.exit(1)

hostname = sys.argv[1]
port = int(sys.argv[2])
msg = 'start'

def make_server_for_remote_control():
    #Setup a standard internet socket.
    #The sockopt call lets this server use the given port even if
    #it was recently used by another server
    sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
    sock.bind((hostname,port))
    sock.listen(1)
    print 'Waiting for a Request'

    #Handle a client request
    request,clientAddress = sock.accept()
    print 'Request received from: ', clientAddress
    data = request.recv(1024)
    print 'Received connection Msg: ', data
    while msg != 'exit':
        msg = raw_input('Enter a Message: ')
        request.send(msg)
        data = request.recv(1024)
        print 'Received Msg: ', data
    request.shutdown(2) #Stop the client from reading or writing anything.
    sock.close()
