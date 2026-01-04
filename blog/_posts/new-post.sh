
if [ -z "$1" ] ; then
   echo Usage: $0 post-name-stub
   exit 1
fi

d=$(date +%Y-%m-%d)
echo $d

f=${d}-$1.md

cat << ! > $f
---
layout: post
title:  "x"
!

echo date: $(date) >> $f
# echo categories: x >> $f
echo --- >> $f


nvim $f

echo git add $f

