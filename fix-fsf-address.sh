
declare -a old_addresses_arr=("675 Mass Ave,.*USA" "59 Temple Place.*USA" "Cambridge,.*USA")

# also I needed these because "USA" was on the next line
# declare -a old_addresses_arr=("59 Temple Place.*1307")

for addr in "${old_addresses_arr[@]}" ;
do
   echo looking for $addr
   # only look in directies at depth 0,1,2 (not save directories)
   files=$(find . -type f -depth -3 -exec grep -nIHl "$addr" {} \;)
   echo files: $files
   for file in $files ; 
   do
      ext="${file##*.}"
      if [ "$ext" != "sed" ] ; then
         echo $file
         # for xxx.*USA pattern
         # sed -i.sed "s/$addr/51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA/" $file
         # for xxx.*1307 pattern
         sed -i.sed "s/$addr/51 Franklin Street, Fifth Floor, Boston, MA 02110-1307/" $file
      fi
   done
done
