
mv interface.h interface.h.tmp

cat << ! > interface.h
#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif
#endif

BEGIN_C_DECLS
!

cat interface.h.tmp >> interface.h

echo END_C_DECLS >> interface.h
