#!/usr/bin/python
import socket
import sys

if len(sys.argv) < 3:
    print 'Usage: %s [hostname] [portnumber]' % sys.argv[0]
    sys.exit(1)

try:
    hostname = sys.argv[1]
    port = int(sys.argv[2])
except:
    print 'Usage: %s [hostname] [portnumber]' % sys.argv[0]
    sys.exit(1)

def make_server_for_remote_control():
    global sock
    msg = 'start'
    print "BL DEBUG:: msg, host, port", msg, hostname,port
    #Setup a standard internet socket.
    #The sockopt call lets this server use the given port even if
    #it was recently used by another server
    sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
    sock.bind((hostname,port))
    sock.listen(1)
    print 'Waiting for a Request'

    #Handle a client request
    request, clientAddress = sock.accept()
    print 'Request received from: ', clientAddress
    data = request.recv(1024)
    print 'Received connection Msg: ', data
    while msg != 'exit':
        msg = raw_input('Enter a Message: ')
        request.send(msg + "# end")
        data = request.recv(1024)
        print 'Received Msg: ', data
    request.shutdown(2) #Stop the client from reading or writing anything.
    sock.close()
    
make_server_for_remote_control()
