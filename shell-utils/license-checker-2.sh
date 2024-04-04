
for file in $(git ls-files) ;
do
    echo ---- $file
    head -20 $file | grep -i license
    head -20 $file | grep -i copyright
done
