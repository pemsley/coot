
if [ $# != 2 ] ; then
   echo must provide guile version and 2 file arguments
   echo we got: $*
   exit 2
fi

pre="$1"
post="$2"

echo "#ifdef USE_PYTHON " > "$post"
cat "$pre" >> "$post" 
echo "#endif // USE_PYTHON " >> "$post"
