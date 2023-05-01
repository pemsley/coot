sort rel-todo | uniq -c \
   | grep -v Done \
   | grep -v ' Fixed'\
   | grep -v ' Fixed\.'\
   | grep -v '\* Punt'\
   | grep -v "Bus error" \
   | grep -v "Thread debuggin" \
   | grep -v "In unknown file" \
   | grep -v "Can.t understand" \
   | grep -v "Can.t reproduce" \
   | grep -v " #0 " \
   | grep -v "\.\.\." \
   | grep -v "Core was generated" \
   | grep -v " It does now" \
   | grep -v "Program terminated with signal 11," \
   | awk 'NF>1 && $1 != "1" && $2 != ":::"'
