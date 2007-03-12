

if [ -z "$1" ] ; then
  echo usage $0 filename
  exit 1
fi

# ok to proceed with adding

mv "$1" "$1".tmp

echo "#ifndef BEGIN_C_DECLS"                 > "$1".tmp.a 
echo "#ifdef __cplusplus"                   >> "$1".tmp.a 
echo "#define BEGIN_C_DECLS extern \"C\" {" >> "$1".tmp.a 
echo "#define END_C_DECLS }"                >> "$1".tmp.a 
echo "#else"                                >> "$1".tmp.a 
echo "#define BEGIN_C_DECLS"                >> "$1".tmp.a 
echo "#define END_C_DECLS"                  >> "$1".tmp.a 
echo "#endif"                               >> "$1".tmp.a 
echo "#endif"                               >> "$1".tmp.a
echo ""              >> "$1".tmp.a
echo "BEGIN_C_DECLS" >> "$1".tmp.a 
echo ""              >> "$1".tmp.a


echo ""             > "$1".tmp.b
echo "END_C_DECLS" >> "$1".tmp.b
echo ""            >> "$1".tmp.b

cat "$1".tmp.a  "$1".tmp "$1".tmp.b > "$1"

rm  "$1".tmp.a  "$1".tmp "$1".tmp.b
