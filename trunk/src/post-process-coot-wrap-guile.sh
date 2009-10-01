
if [ $# != 4 ] ; then
   echo must provide guile-version gtk2-flag and 2 file arguments
   echo we got: $*
   exit 2
fi

pre="$3"
post="$4"

guile_version=$1
gtk2=$2 

case $guile_version in

   1.6*) 
   echo sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre $post
        sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre > $post
    ;; 
   
   1.8*)
   echo sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/' -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' -e '/.libguile.h./{x;s/.*/#include <cstdio>/;G;}' $pre ..to.. $post
        sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/' -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' -e '/.libguile.h./{x;s/.*/#include <cstdio>/;G;}' $pre > $pre.tmp
   echo adding gtk2 and guile ifdefs
   echo "#ifdef USE_GUILE"  > $post
   if [ "$gtk2" = gtk2 ] ; then
      echo "#ifdef COOT_USE_GTK2_INTERFACE"  >> $post
   else 
      echo "#ifndef COOT_USE_GTK2_INTERFACE"  >> $post
   fi    
   cat $pre.tmp >> $post
   echo "#endif // COOT_USE_GTK2_INTERFACE" >> $post
   echo "#endif // USE_GUILE" >> $post
   ;; 

esac



