
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
        extension="${file##*.}"
        if [ "$extension" = ".png" ] ; then
            check_this=false
        fi
        if [ "$extension" = ".tab" ] ; then
            check_this=false
        fi
        if [ -d $file ] ; then
           check_this=false
        fi
        if [ $check_this = true ] ; then
            echo "" ----- $file
            head -20 $file | grep -i license
            head -20 $file | grep -i copyright
        fi
    else
       echo pass $file
    fi
done
