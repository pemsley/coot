# create_server.py

import socket
import sys

hostname = sys.argv[1]
port = int(sys.argv[2])

def connection_proc(port, hostname):

    close_transmission_string = "# close\n# end\n"
    end_transmission_string = "# end\n"

    def end_transmission():
        request.send(end_transmission_string)

    def close_transmission():
        request.send(close_transmission_string)

    # main body

    print("BL DEBUG:: host, port", hostname,port)
    #Setup a standard internet socket.
    #The sockopt call lets this server use the given port even if
    #it was recently used by another server
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
    sock.bind((hostname,port))
    sock.listen(1)

    #Handle a client request
    request, clientAddress = sock.accept()
    print("################ Got a connection!")

    request.send("\n")

    print("BL DEBUG:: read ACT")
    request.send("read_pdb('monomer-ACT.pdb')\n")
    end_transmission()

    import time
    time.sleep(5)
    # more!?
    print("BL DEBUG:: read PIN")
    request.send("imol = coot.read_pdb('monomer-PIN.pdb')\n")
    request.send("set_rotation_centre(5, 5, 5)\n")
    request.send("move_molecule_to_screen_centre(imol)\n")
    end_transmission()
    # close transmission
    print("BL DEBUG:: send close")
    close_transmission()
    sock.close()

    #data = request.recv(1024)
    #print 'Received connection Msg: ', data
    #while msg != 'exit':
    #    msg = raw_input('Enter a Message: ')
    #    request.send(msg)
    #    data = request.recv(1024)
    #    print 'Received Msg: ', data
    #request.shutdown(2) #Stop the client from reading or writing anything.

connection_proc(port, hostname)
