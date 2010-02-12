
if [ $# -lt 1 ] ; then
   echo $0 arg
else 
   grep -nH $1 */*.cc */*.c */*.hh */*.h */*.scm */*.py
fi

