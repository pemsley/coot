# group_settings.py
# Copyright 2004, 2006 by Paul Emsley, The University of York
# Copyright 2007 by Bernhard Lohkamp, The University of York
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
 

# Add here settings for your group - rather than expecting each user
# to add them to their own .coot files.
# 

# York setting for molprobity. 
# 
#probe_command = "/y/people/lohkamp/coot/Linux/bin/probe.2.12.070727.linuxi386"
#reduce_command = "/y/people/lohkamp/coot/Linux/bin/reduce.3.10.070814.linuxi386"
# try to get the abspath of probe and reduce automatically
# alternatively, specify yourself here.
global probe_command
global reduce_command

probe_command = find_exe("probe", "PATH", "CBIN", "CCP4_BIN")
if not probe_command:
    probe_command = "probe" # useless (?) fallback

reduce_command = find_exe("reduce", "PATH", "CBIN", "CCP4_BIN")
if not reduce_command:
    reduce_command = "reduce" # useless (?) fallback

    
