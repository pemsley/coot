

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

if test "$OS" = Linux ; then 

  architecture=$(uname -i)
  # uname -i and uname -p (strangely) return unknown on my ubuntu
  if test $architecture = unknown ; then
     architecture=$(uname -m)
  fi

  if rpm -a >/dev/null 2>&1 ; then
      for i in fedora-release redhat-release redhat-release-workstation centos-release sl-release openSUSE-release ; do
        dist=$(rpm -q --qf '%{name}' ${i})
        if test $? = 0 ; then
          dist_name=$(echo ${dist} | sed -e 's/\-release//g' -e 's/\-workstation//g')
          dist_ver=`rpm -q --qf '%{version}' ${i}`
          break
        else
          dist_name='unknown'
        fi
      done
   else
      dist_name='unknown'
   fi

  case ${dist_name} in

    redhat* )
    case ${dist_ver} in
      [0-9] | [0-9].[0-9]* )
        systype=${architecture}-redhat-${dist_ver}
      ;;
      * )
        systype=${architecture}-rhel-$(echo ${dist_ver} | sed s/[A-Za-z]//g)
        if test "$architecture" = x86_64 ; then 
	   if test $dist_ver = 4WS  ; then
	      echo RedHat 4 Linux x86_64 detected. need to update libtool
              update_libtool=1
              # stupid la file of libcur puts /usr/lib in the link path (early).  This causes link 
              # problems on 64 bit RHEL4.  So lets fix curl-config and libcurl.la
              post_process_libcurl=1
	   fi
        fi
      ;;
    esac
    ;;

    fedora | centos | openSUSE )
      systype=${architecture}-${dist_name}-${dist_ver}
    ;;


    sl ) 
      systype=${architecture}-scientific-linux-${dist_ver}
      ;; 


    * )
      if test -r /etc/issue; then

dnl     $1 and $2 were being substituted to "" (nothing) and I
dnl     couldn't see how to protect them.  So use cut, tr and head
dnl     instead.
dnl         dist_name=`awk 'NR==1{print tolower($1)}' /etc/issue`
dnl         dist_ver=`awk 'NR==1{print tolower($2)}' /etc/issue`

        dist_name=`cut -d" " -f 1 /etc/issue | tr A-Z a-z | head -1`
        dist_ver=` cut -d" " -f 2 /etc/issue | tr A-Z a-z | head -1`

        if test "$dist_name" = "debian" ; then
	   dist_ver=$(cut -d" " -f 3 /etc/issue | tr A-Z a-z | tr / - | head -1)
        fi
        systype=${architecture}-${dist_name}-${dist_ver}
      else
        systype=${architecture}-unknown-Linux
      fi
    ;;
  esac
fi

dnl echo :::::: architecture has been set to $architecture
dnl echo :::::: systype has been set to $systype

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
