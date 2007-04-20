# Copyright 2004, 2005 by The University of York
# Copyright 2005 by Bernhard Lohkamp
 
#;;; This program is free software; you can redistribute it and/or modify
#;;; it under the terms of the GNU General Public License as published by
#;;; the Free Software Foundation; either version 2 of the License, or (at
#;;; your option) any later version.
 
#;;; This program is distributed in the hope that it will be useful, but
#;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;; General Public License for more details.
 
#;;; You should have received a copy of the GNU General Public License
#;;; along with this program; if not, write to the Free Software
#;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;; direction is either 'forwards or 'backwards
#;; 
#;; start-resno is higher than stop-resno if we are building backwards
#;; 
#;; (fit-gap 0 "A"  23 26)   ; we'll build forwards
#;; (fit-gap 0 "A"  26 23)   ; we'll build backwards
#;;

# translated by BL
import time
import os

user = os.getenv('USERNAME')

hour = int(time.strftime("%H", time.localtime()))
if hour < 12: time_str = "Morning"
elif hour < 18: time_str = "Afternoon"
else : time_str = "Evening"
hello_str =  "Good %(a)s %(b)s, Welcome to Coot." % {"a":time_str, "b":user}
# alternative
# hello_str1 =  "%s %s %s%s" % ("Good",time_str, user, ", Welcome to Coot.")

print hello_str
# print hello_str1
set_display_intro_string(hello_str)
