

AC_DEFUN([AM_COOT_SYS_BUILD_TYPE],
[
AC_MSG_CHECKING([build type])


OS=`uname`
systype=unknown 
if test "$OS" = "Darwin" ; then
   osversion=`sw_vers -productVersion`
   # uname -a gives processor type in last field on Darwin
   processor=`uname -a | awk '{print $(NF)}'`
   systype=MacOSX-${osversion}-${processor}
fi

if test "$OS" = "MINGW32_NT-5.1" ; then
   systype=`uname -m`
fi

if test $OS = Linux ; then 

  architecture=`uname -i`
  # uname -i and uname -p (strangely) return unknown on my ubuntu
  if test $architecture = unknown ; then
     architecture=`uname -m`
  fi

  which rpm > /dev/null
  have_rpm=$?
  if  test $have_rpm = 0 ; then 
      for i in fedora redhat centos ; do
	dist=`rpm -q --qf '%{name}' ${i}-release`
	if test $? = 0 ; then
	  dist_name=`echo ${dist} | sed s/\-release//g`
	  dist_ver=`rpm -q --qf '%{version}' ${i}-release`
	  break
	else
	  dist_name='unknown'
	fi
      done
   else
      dist_name='unknown'
   fi

  case ${dist_name} in
    redhat )
    case ${dist_ver} in
      [0-9] | [0-9].[0-9]* )
        systype=${architecture}-redhat-${dist_ver}
      ;;
      * )
        systype=${architecture}-rhel-`echo ${dist_ver} | sed s/[A-Za-z]//g`
      ;;
    esac
    ;;
    fedora | centos )
      systype=${architecture}-${dist_name}-${dist_ver}
    ;;
    * )
      # echo going the star path in dist_name case
      if test -r /etc/issue; then

dnl     $1 and $2 were being substituted to "" (nothing) and I
dnl     couldn't see how to protect them.  So use cut, tr and head
dnl     instead.
dnl         dist_name=`awk 'NR==1{print tolower($1)}' /etc/issue`
dnl         dist_ver=`awk 'NR==1{print tolower($2)}' /etc/issue`
	dist_name=`cut -d" " -f 1 /etc/issue | tr A-Z a-z | head -1`
	dist_ver=` cut -d" " -f 2 /etc/issue | tr A-Z a-z | head -1`
dnl         echo dist_name $dist_name
dnl         echo dist_ver $dist_ver
        systype=${architecture}-${dist_name}-${dist_ver}
      else
        systype=${architecture}-unknown-Linux
      fi
    ;;
  esac
fi

# echo :::::: architecture has been set to $architecture
# echo :::::: systype has been set to $systype

if test "$coot_python" = true ; then
   # echo setting python_tag to -python
   python_tag=-python
else 
   # echo setting python_tag to blank
   python_tag=
fi
  
if test $coot_gtk2 = TRUE ; then
   gtk2=-gtk2
else
   gtk2=
fi

COOT_SYS_BUILD_TYPE=${OS}-${systype}${python_tag}${gtk2}

AC_MSG_RESULT([$COOT_SYS_BUILD_TYPE])
AC_SUBST(COOT_SYS_BUILD_TYPE)
])
