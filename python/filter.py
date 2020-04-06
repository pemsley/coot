# filter.py
# Copyright 2006, 2007 by Bernhard Lohkamp
# Copyright 2006, 2007 by Paul Emsley, The University of York
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


# Basic python function, filter the objects in list ls by function
# fn.  e.g. filter(lambda x: not x%2, [0,1,2,3]) -> [0,2]
# 
# BL says:: python has this function build in (now), so use python's
# version
#
def filter_coot(fn,ls):

  ret_ls = []
  if len(ls) > 0:
   for i in range(len(ls)):
     print(i, ls[i])
     if fn(ls[i]):
        print(i, ls[i])
        ret_ls.append(ls[i])
     else: pass
  else:
    ret_ls = ls
  return ret_ls


