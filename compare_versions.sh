#!/bin/bash
tmp_v1=`grep AC_INIT configure.ac`
tmp_v2=`grep '<property name="version">' ui/coot-gtk4.ui`
ver1=$(echo "$tmp_v1" | tr -cd '[:digit:].')
ver2=$(echo "$tmp_v2" | tr -cd '[:digit:].')
echo version configure $ver1
echo version ui file   $ver2
if [ $ver1 == $ver2 ] ; then
	echo All good
	echo $ver1 > coot-version
else
	echo damn, versions dont match
fi
