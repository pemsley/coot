#   remote_control.py
#   Copyright (C) 2008, 2009  Bernhard Lohkamp, University of York
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
def open_coot_listener_socket(port_number, host_adress = "127.0.0.1"):

    import socket
    global coot_listener_socket

    print("in open_coot_listener_socket port: %s host %s" %(port_number, host_adress))

    soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    soc.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    #host_adress = "127.0.0.1"

    soc.connect((host_adress, port_number))

    print("Coot listener socket ready!")
    soc.send("Coot listener socket ready!\n")
    coot_listener_socket = soc
    if soc:
        set_coot_listener_socket_state_internal(1)

    # threads? Needed?
    status = run_python_thread(coot_listener_idle_function_proc,())

# yet another go to make a coot port reader work.  This time, we use
# a gtk-timer to read stuff from the socket.  
# 
# The gtk-timer function must return True to be called again.  When we
# want to close the socket reader, simply make the function return False.
# 
def open_coot_listener_socket_with_timeout(port_number, host_adress = "127.0.0.1"):

    try:
        import gobject
    except:
        print("BL WARNING:: no gobject available, so no socket. Sorry!")
        return
    
    import socket
    global coot_listener_socket

    print("in open_coot_listener_socket port: %s host %s" %(port_number, host_adress))

    soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    #host_adress = "127.0.0.1"

    try:
        soc.connect((host_adress, port_number))
    except:
        print("BL INFO:: cannot connect to socket on host %s with port %s" %(host_adress, port_number))
        return

    print("Coot listener socket ready!")
    soc.send("Coot listener socket ready!\n")
    coot_listener_socket = soc

    gobject.timeout_add(1000, coot_socket_timeout_func)

# based on coot_listener_idle_func_proc
#
def coot_socket_timeout_func():
    
    global coot_listener_socket

    global total_socket_data
    if (coot_listener_socket):
        while (1):
            continue_qm = listen_coot_listener_socket(coot_listener_socket)
            if (not continue_qm):
                set_coot_listener_socket_state_internal(0)
                try:
                    coot_listener_socket.shutdown(2)
                except:
                    try:
                        coot_listener_socket.shutdown(1)
                    except:
                        print("BL INFO:: problems shutting down the controling socket, ")
                        print("maybe already down")
                coot_listener_socket.close()
                coot_listener_socket = False
                print("server gone - listener thread ends", coot_listener_socket)
                return False
            else:
                # keep listening
                #time.sleep(0.01)
                return True
    else:
        print("coot_socket_timeout_func bad sock: ", coot_listener_socket)
        return False


def coot_listener_idle_function_procbbbbbb():

    #global coot_listener_socket
    msg = 'start'
    while msg != 'exit':
        data = coot_listener_socket.recv(1024)
        print('Received Message: ', data)
        #msg = raw_input('Enter a Message: ')
        coot_listener_socket.send('return value')
    coot_listener_socket.close()
#    import time
#    if (coot_listener_socket):
#        time.sleep(1000)
#    else:
#        print "coot_listener_idle_function_proc bad sock: ", coot_listener_socket


def coot_listener_error_handler(key, args=""):

    print("coot_listener_error_handler handling error in %s with args %s" %(key, args))


def coot_listener_idle_function_proc():
    global coot_listener_socket
    import time
    #print " in idle func"
    if (coot_listener_socket):
        while (1):
            continue_qm = listen_coot_listener_socket(coot_listener_socket)
            if (not continue_qm):
                set_coot_listener_socket_state_internal(0)
                try:
                    coot_listener_socket.shutdown(2)
                except:
                    try:
                        coot_listener_socket.shutdown(1)
                    except:
                        print("BL INFO:: problem shutting down the controling socket, ")
                        print("maybe already down")
                coot_listener_socket.close()
                coot_listener_socket = False
                print("server gone - listener thread ends", coot_listener_socket)
                return
            else:
                time.sleep(0.01)
                return
    else:
        print("coot_listener_idle_func_proc bad sock: ", coot_listener_socket)

global total_socket_data
total_socket_data = []

def listen_coot_listener_socket(soc):

    global total_socket_data
    close_connection_string = "# close"
    end_connection_string = "# end"    # just 'end' for test

    def evaluate_character_list(string_list):
        print("received in evaluate", string_list)
        # new eval?!
        for line in string_list.split("\n"):
            if (line == close_connection_string):
                print("finish socket")
                #soc.send("Closing session")
                return False
            try:
                ret = eval(line)
            except SyntaxError:
                if (line == "\n"): print("CR")
                try:
                    exec(line, globals())
                    ret = None
                except:
                    ret = "INFO:: cannot eval or exec given string: " + string
            except:
                    ret = "Info:: input error"
            # can only send strings back in this way
            ret = "return value: " + str(ret)
            try:
                soc.send(ret)
            except:
                pass
        return True

    # main body
    #
    data = False
    soc.setblocking(0)
    def check_aliveness():
        # first check if serve is still there?
        try:
            soc.sendall("")
        except:
            print("BL INFO:: appear that serve is down, closing down connection")
            return False
        return True

    try:
        while True:
            data = soc.recv(1024)

            import time
            time.sleep(0.1)
            if (end_connection_string in data):
                total_socket_data.append(data[:data.rfind(end_connection_string)])
                break

            if (not data == ""):
                total_socket_data.append(data)
            else:
                return check_aliveness()
                
            if (len(total_socket_data) > 1):
                # check if end_connection_string was split
                last_pair = total_socket_data[-2] + total_socket_data[-1]
                if (end_connection_string in last_pair):
                    total_socket_data[-2] = last_pair[:last_pair.find(end_connection_string)]
                    total_socket_data.pop()
                    break
        total_socket_data = ''.join(total_socket_data)
        if (not data):
            #print "nothing on the line..."
            return True
        else:
            ret = evaluate_character_list(total_socket_data)
            total_socket_data = [] # reset total data
            return ret
    except:
        return True
        


#open_coot_listener_socket(50007, "bla")

