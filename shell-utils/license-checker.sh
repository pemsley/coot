
pass_overs="save save/*/* */save/*/*"

if [ 1 = 2 ] ; then
    for po in $pass_overs ;
    do
       echo $po
    done
fi

for file in * */* */*/* ;
do
    check_this=true
    for po in $pass_overs ;
    do
        if [ "$file" = $po ] ; then
            check_this=false
        fi
    done
    if [ $check_this = true ] ; then
       echo check this $file
       :
    else
       echo pass $file
    fi
done
