#
# Copyright 2007 The University of York
# Author: Paul Emsley
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

if test -e ltmain.sh ; then
    # this should not be in the repo, I guess
    rm ltmain.sh
fi

echo glibtoolize --copy
glibtoolize --copy

echo aclocal -I macros
aclocal -I macros

echo autoconf
autoconf

echo automake --add-missing --copy
automake --add-missing --copy

