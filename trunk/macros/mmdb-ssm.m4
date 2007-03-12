
AC_DEFUN([AM_WITH_MMDBSSM],
[AC_PROVIDE([AM_USE_MMDBSSM])


AC_ARG_WITH(ssmlib-prefix, 
	AC_HELP_STRING( [--with-ssmlib-prefix=PFX], [Prefix where SSMLib has been installed] ),
	[ with_ssmlib_prefix="$withval" ],
 with_ssmlib_prefix="")

AC_MSG_CHECKING([for SSMLib])

if test x$with_ssmlib_prefix != x; then

   MMDBSSM_CXXFLAGS="-DHAVE_SSMLIB"
   MMDBSSM_LIBS="-L$with_ssmlib_prefix/lib -lssm"

ac_mmdb_dirs='
.
include
include/ssm
include/mmdb
lib
src
lib/src
lib/src/mmdb'

   for ac_dir in $ac_mmdb_dirs; do
      if test -r "$with_ssmlib_prefix/$ac_dir/ssm_superpose.h"; then
         ac_MMDBSSM_CXXFLAGS="-I$with_ssmlib_prefix/$ac_dir"
         break
      fi
   done

  MMDBSSM_CXXFLAGS="$MMDBSSM_CXXFLAGS $ac_MMDBSSM_CXXFLAGS"
  
else 

   MMDBSSM_CXXFLAGS=""
   MMDBSSM_LIBS=""
   with_ssmlib_prefix=no 

fi

AC_MSG_RESULT([$with_ssmlib_prefix])

AC_SUBST(MMDBSSM_CXXFLAGS)
AC_SUBST(MMDBSSM_LIBS)
])
