

AC_DEFUN([AM_COOT_SYS_BUILD_TYPE],
[
AC_MSG_CHECKING([build type])

echo ..........................here\! ...............

OS=`uname`
systype=unknown 
if test "$OS" = "Darwin" ; then
   osversion=`sw_vers -productVersion`
   # uname -a gives processor type in last field on Darwin
   processor=`uname -a | awk '{print $(NF)}'`
   systype=MacOSX-${osversion}-${processor}
fi

if test $OS = Linux ; then 
  which rpm > /dev/null
  have_rpm=$?
  if  [ $have_rpm = 0 ] ; then 
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
        if [ $arch = x86_64 ] ; then 
	   if [ $dist_ver = 4WS ] ; then
	      echo RedHat 4 Linux x86_64 detected. need to update libtool
              update_libtool=1
	   fi
        fi
      ;;
    esac
    ;;
    fedora | centos )
      systype=${architecture}-${dist_name}-${dist_ver}
    ;;
    * )
      if test -r /etc/issue; then
        dist_name=`awk 'NR==1{print tolower($1)}' /etc/issue`
        dist_ver=`awk 'NR==1{print tolower($2)}' /etc/issue`
        systype=${architecture}-${dist_name}-${dist_ver}
      else
        systype=${architecture}-unknown-Linux
      fi
    ;;
  esac
fi

if test "$coot_python" = true ; then
   python_tag=-python
else 
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
