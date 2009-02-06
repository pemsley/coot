
if [ $# != 3 ] ; then
   echo must provide guile version and 2 file arguments
   echo we got: $*
   exit 2
fi

pre="$2"
post="$3"

guile_version=$1

case $guile_version in

   1.6*) 
   echo sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre $post
   sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre > $post
    ;; 
   
   1.8*)
   echo sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/' -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre $post
   sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/' -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' -e '/.libguile.h./{x;s/.*/#include <cstdio>/;G;}' $pre > $post
   ;; 

esac



