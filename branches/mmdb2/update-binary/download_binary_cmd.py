
# Return a boolean for whether or not pre-release was passed as a
# command line arg.
def extract_pre_release_flag(cla):
    return "pre-release" in cla

# Return the argument after interesting-string, or False if no such
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

# Return the binary type as a string (return #f for this part on
# string-not-found).  The binary string is any string that is
# proceeded with an arg of "binary"
def extract_version_type(cla):
    return extract_next_arg_from_command_line("version", cla)

import sys, os
command_line = sys.argv[1:]
vers             = extract_version_type(command_line)
binary_type      = extract_binary_type(command_line)
pre_release_flag = extract_pre_release_flag(command_line)
is_windows       = os.name == 'nt'

if is_windows:
    web_name = "lohkamp"
else:
    web_name = "emsley"

print "vers:", vers
print "binary-type:", binary_type
print "pre-release:", pre_release_flag

# for testing
#binary_type = "Linux-i386-fedora-10-python-gtk2"

if (binary_type and vers):

    if pre_release_flag:
        binary_url = "http://www.ysbl.york.ac.uk/~"  + \
                    web_name + \
                    "/software/binaries/nightlies/pre-release/" + \
                    vers + \
                    "-binary-" + \
                    binary_type  + \
                    ".tar.gz"
    else:
        binary_url = "http://www.ysbl.york.ac.uk/~" + \
                    web_name + \
                    "/software/binaries/release/" + \
                    vers + \
                    "-binary-" + \
                    binary_type + \
                    "tar.gz"

    print "download_bin_cmd: get url:", binary_url
