
cootmonomersdir = $(pkgdatadir)/lib/data/monomers

dist_cootmonomers_DATA = \
a/A.cif                 d/DA.cif                g/GLY.cif               m/MET.cif               t/THR.cif \
a/ALA.cif               d/DC.cif                g/GOL.cif               m/MSE.cif               t/TRP.cif \
a/ARG.cif               d/DG.cif                h/HIS.cif               n/NA.cif                t/TYR.cif \
a/ASN.cif               d/DT.cif                h/HOH.cif               p/PHE.cif               u/U.cif   \
a/ASP.cif               e/EDO.cif               i/ILE.cif               p/PO4.cif               v/VAL.cif \
c/C.cif                 g/G.cif                 l/LEU.cif               p/PRO.cif               m/MG.cif  \
c/CIT.cif               g/GLN.cif               l/LYS.cif               s/SER.cif \
c/CYS.cif               g/GLU.cif               list/mon_lib_list.cif   s/SO4.cif

install-data-local:
	if test -d $(srcdir); then \
	   mkdir $(DESTDIR)/$(pkgdatadir)/lib/data/monomers; \
	   ls -d $(srcdir)/* ; \
           for cif_dir in $(srcdir)/* ; \
               do \
               l_cif_dir=$$(basename $$cif_dir) ; \
	       $(MKDIR_P) $(DESTDIR)/$(pkgdatadir)/lib/data/monomers/$$l_cif_dir ; \
	       for cif_file in $$cif_dir/*.cif ; do \
                   if test -e $$cif_file ; then \
		       cp $$cif_file $(DESTDIR)/$(pkgdatadir)/lib/data/monomers/$$l_cif_dir ; \
		   fi \
                   done ; \
               done \
        fi
