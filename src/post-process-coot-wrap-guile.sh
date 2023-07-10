

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


if [ $# != 2 ] ; then
    echo bad args
    exit 2
fi

pre=coot_wrap_guile_pre.cc
post=coot_wrap_guile.cc

add_ifdefs $pre $post
