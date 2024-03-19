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
sys=`uname`

OS=`uname`
hostname=`hostname`
if [ -z "$HOST" ]; then
 HOST=$hostname
fi

if test -e ltmain.sh ; then
    # this should not be in the repo, I guess
    rm ltmain.sh
fi


# Bernie wants libtoolize
echo libtoolize --copy
libtoolize --copy


# we can only add the directories if they exist, otherwise aclocal barfs (sigh)

aclocal_extras=

# these need to correspond with the build directories
dir_list="$HOME/autobuild/$OS-$HOST-gtk4/share/aclocal $HOME/autobuild/$OS-$HOST-gtk3/share/aclocal"

for dir in $dir_list
    do
       if [ -e $dir ] ; then
	  aclocal_extras="$aclocal_extras -I $dir"
       fi
    done

case "$sys" in
    MINGW32_NT-* )
	echo We have WIN
	aclocal_extras="-I /usr/local/share/aclocal -I /mingw/share/aclocal" ;;
esac
echo aclocal -I macros $aclocal_extras
aclocal -I macros $aclocal_extras

echo autoconf
autoconf

echo automake --add-missing --copy
automake --add-missing --copy

