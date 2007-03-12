
if [ -z "$1" ] ; then
  echo usage $0 filename
  exit 1
fi


sed -e 's?#  define gettext(String) (String)?/*#  define gettext(String) (String) */?' \
    -e 's?#  define dgettext(Domain,Message) (Message)?/* #  define dgettext(Domain,Message) (Message) */?' \
    $1 > $1.tmp

mv $1.tmp $1

