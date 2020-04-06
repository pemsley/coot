# hello.py
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2003, 2004, 2005, 2006 Paul Emsley, University of York
# Copyright 2016 by Medical Research Council
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

import time
import string
import os, sys
import getpass

def coot_says_hello():

   def first_non_trivial_name(str_list):

       name = "Person with No Name"
       if len(str_list) == 0 : pass
       else:
          for ls in str_list:
             if (not ls.endswith(".") and len(ls)>1):
                name = ls
                break
       return name

   hour = int(time.strftime("%H", time.localtime()))
   if hour < 12: time_str = "Morning"
   elif hour < 18: time_str = "Afternoon"
   else : time_str = "Evening"
   try:

      import pwd
      user = getpass.getuser()
      name_string = pwd.getpwnam(user).pw_gecos
      name_strings = name_string.split()

   except:
      name_string = getpass.getuser()
      name_strings = name_string.split()

   # reverse name_strings if locale is japanese
   l1 = os.getenv("LANG")
   l2 = os.getenv("LANGUAGE")
   if l1 == "ja": name_strings.reverse() # perhaps beginswith
   if l2 == "ja": name_strings.reverse()
   personal_name = string.capitalize(first_non_trivial_name(name_strings))
   hello_str = "Good %s %s, Welcome to Coot version %s" %(time_str, personal_name, coot_version())
   print(hello_str)
   set_display_intro_string(hello_str)

# check for main
#if __name__ == "__main__":

coot_says_hello()
