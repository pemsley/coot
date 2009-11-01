

# We are given the binary sys type from the command line, and also
# whether or not this is a pre-release ("pre-release" is passed on
# the command line).
# 
# We need to construct the right url for the binary and go and check
# if it is there.


# Return a boolean for whether or not pre-release was passed as a
# command line arg.
def extract_pre_release_flag(cla):
    return "pre-release" in cla

# Return the argument after interesting-string, or null if no such
# arg.
def extract_next_arg_from_command_line(interesting_string, cla):
    pos = cla.index(interesting_string)
    if (len(cla) > pos+1):
        return cla[pos+1]
    else:
        return False
  
# Return the binary type as a string (return False for this part on
# string_not_found).  The binary string is any string that is
# proceeded with an arg of "binary"
def extract_binary_type(cla):
  return extract_next_arg_from_command_line("binary", cla)

# Return what we presume to be the first argument in the command line
# for the program that called this program.  We do this so that we
# can find the bin dir and bin/.. dir and test whether we can write
# into it.
# FIXME not sure what needs to be done for python and esp Windows yet!
# 
def extract_coot_command_line(cla):
  return extract_next_arg_from_command_line("command-line", cla)


def get_url_as_string(my_url):

    import urllib
    try:
        s = urllib.urlopen(my_url).read()
    except:
        s = False # not sure if that's the right way
    return s

import sys, os
command_line = sys.argv[1:]
pre_release_flag = extract_pre_release_flag(command_line)
binary_type      = extract_binary_type(command_line)
is_windows       = os.name == 'nt'

if is_windows:
    web_name = "lohkamp"
else:
    web_name = "emsley"
web_name = "emsley"

#   print "pre-release:", pre_release_flag
#   print "binary-type:", binary_type

# for testing
#binary_type = "Linux-i386-fedora-10-python-gtk2"

if binary_type:

    if pre_release_flag:
        check_url = "http://www.ysbl.york.ac.uk/~"  + \
                    web_name + \
                    "/software/binaries/nightlies/pre-release/" + \
                    "type-binary-" + \
                    binary_type  +\
                    "-latest.txt"
    else:
        check_url = "http://www.ysbl.york.ac.uk/~" + \
                    web_name + \
                    "/software/binaries/release/" + \
                    "type-binary-" + \
                    binary_type + \
                    "-latest.txt"

    # print "BL DEBUG:: Now get:", check_url
    s = get_url_as_string(check_url)

    print "LATEST-BUILD-ON-SERVER:", s

print command_line
