# hello.py 
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2003, 2004, 2005, 2006 Paul Emsley, University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Primarily for Indian Names.
#
# Say we are given str: (list "M." "D." "Albert" "Dorkins").  We want
# to return ""Albert" not "M.")  We reject names that are one
# character long and names that are 2 characters long that end in
# ".".  So, "M." is rejected, but "Ma" is not.
#
# An excercise for the committed is to also reject run-together
# dotted initials such as "D.K.".  I was sufficiently committed.
# BL says: I'll try to too! Although I dont think this is a valuable
# (Win)Coot function, but fun though...
#
def first_non_trivial_name(str_list):

    name = "Person with No Name"
    if len(str_list) == 0 : pass
    else:
      for ls in str_list:
          if (not ls.endswith(".") and len(ls)>1):
             name = ls
             break
    return name

import time
import string
import os, sys

if (sys.platform == 'darwin'):
    user = os.getenv('USER')
else:
    user = os.getenv('USERNAME')

hour = int(time.strftime("%H", time.localtime()))
if hour < 12: time_str = "Morning"
elif hour < 18: time_str = "Afternoon"
else : time_str = "Evening"

name_strings = string.split(str(user))

# BL says: sorry but in windows we dont have a proper LC_* or language
# setting, so no swapping of names for Japanese, Koreans etc.
# and I cant test (and be asked to do it in general, maybe later...

personal_name = string.capitalize(first_non_trivial_name(name_strings))

hello_str =  "Good %(a)s %(b)s, Welcome to Coot. %(c)s" \
            %{"a":time_str, "b":personal_name, "c":coot_version()}
# alternative
# hello_str1 =  "%s %s %s%s" % ("Good",time_str, user, ", Welcome to Coot.")

print hello_str
set_display_intro_string(hello_str)

