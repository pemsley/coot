
if [ $# != 2 ] ; then
   echo must provide 2 file arguments
   exit 2
fi

pre="$1"
post="$2"

sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/' $pre > $post
