


if [ -z "$BASH_VERSINFO" ] ; then
    echo wrong shell
    exit 1
fi

function add_ifdefs { 

pre_tmp=$1
post=$2
    
echo adding gtk2 and guile ifdefs from $pre_tmp $post

echo "#ifdef USE_PYTHON"      > $post
echo "#include \"Python.h\"" >> $post
echo "#endif"                >> $post

echo "#ifdef USE_GUILE"      >> $post
echo "#include <cstddef> "   >> $post

cat $pre_tmp >> $post
echo "#endif // USE_GUILE" >> $post

}


#####  start ################

if [ $# != 4 ] ; then
   echo must provide path-to-guile-config gtk2-flag and 2 file arguments
   echo we got: $*
   exit 2
fi


pre="$3"
post="$4"

# sometimes "" gets passed as $1
guile_config=$1
gtk2=$2 

echo ::::::::::::::::::::: post-process-coot-wrap-guile 1: $1 
echo ::::::::::::::::::::: post-process-coot-wrap-guile 2: $2 
echo ::::::::::::::::::::: post-process-coot-wrap-guile 3: $3 
echo ::::::::::::::::::::: post-process-coot-wrap-guile 4: $4 

if [ -z "$guile_config" ] ; then
   echo WARNING:: $0 guile_config arg was missing
else
   guile_version=`$guile_config --version 2>&1`
fi

# these can be blank (see below)
echo post-process-coot-wrap-guile: guile_config:  $guile_config
echo post-process-coot-wrap-guile: guile_version: $guile_version

# if guile_config was blank (as can be the case when we compile
# without guile, then we want a new blank file for $post.
# Otherwise, do the filtering and addtion of the ifdefs
#
if [ -z "$guile_config" ] ; then 
   if [ -e "$post" ] ; then 
      rm -f "$post"
   fi
   touch "$post"

else 


   case $guile_version in

      *1.6*) 
      echo sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre $post
           sed -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' $pre > $pre.tmp
      add_ifdefs $pre.tmp $post
      ;; 
   
      *1.8*)
      # SCM_MUST_MALLOC to be replaced by scm_gc_malloc is more complicated.  scm_gc_malloc takes 2 args and 
      # SCM_MUST_MALLOC takes one.  Leave it to be fixed in SWIG.
      echo sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/'  \
               -e 's/SCM_STRINGP/scm_is_string/'              \
               -e 's/SCM_STRING_LENGTH/scm_c_string_length/'  \
               -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' \
               -e 's/static char .gswig_const_COOT_PYTHON_DIR/static const char *gswig_const_COOT_PYTHON_DIR/' \
               -e '/.libguile.h./{x;s/.*/#include <cstdio>/;G;}' \
               $pre ..to.. $pre.tmp
           sed -e 's/SCM_STRING_CHARS/scm_to_locale_string/'      \
               -e 's/SCM_STRINGP/scm_is_string/'                  \
               -e 's/SCM_STRING_LENGTH/scm_c_string_length/'      \
               -e 's/static char .gswig_const_COOT_SCHEME_DIR/static const char *gswig_const_COOT_SCHEME_DIR/' \
               -e 's/static char .gswig_const_COOT_PYTHON_DIR/static const char *gswig_const_COOT_PYTHON_DIR/' \
               -e '/.libguile.h./{x;s/.*/#include <cstdio>/;G;}' \
               $pre > $pre.tmp
      add_ifdefs $pre.tmp $post
      ;; 
   esac

fi




