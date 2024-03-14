
dir=$1
if [ -z "$dir" ] ; then
    echo Usage: $0 dir-name
    exit
fi

pass_overs="1h1s.pdb autogen.sh cd44A.pdb config.sub config.guess configure.ac ctpl_boost.h Makefile.am Makefile.in surface.h surface.cpp Simple.frag Simple.vert nautilus_lib.pdb cootaneer-llk-2.40.dat protein.db"

for source in $dir/*.* ;
do
    # echo source $source
    fix_this=true
    for pass_over in $pass_overs
    do
        if [ $source = $dir/$pass_over ] ; then
            fix_this=false
        fi
    done
    extension="${source##*.}"
    if [ "$extension" = "new" ] ; then
        fix_this=false
    fi
    if [ $fix_this = true ] ; then
        python3 license-header-replace.py $source
        if [ -e $source.new ] ; then
            # diff $source $source.new
            mv $source.new $source
        fi
    else
        echo dont fix this $source
    fi
done
