
cootmonomersdir = $(pkgdatadir)/lib/data/monomers

cootmonomers_files = ener_lib.cif \
1/1PE.cif               c/CU.cif                h/HOH.cif               p/PEG.cif \
5/5MC.cif               c/CYS.cif               i/ILE.cif               p/PG4.cif \
a/A.cif                 d/DA.cif                i/IMD.cif               p/PGE.cif \
a/ACE.cif               d/DC.cif                i/IOD.cif               p/PHE.cif \
a/ACT.cif               d/DG.cif                i/IPA.cif               p/PLP.cif \
a/ACY.cif               d/DMS.cif               k/K.cif                 p/PO4.cif \
a/ADP.cif               d/DT.cif                k/KCX.cif               p/PRO.cif \
a/ALA.cif               e/EDO.cif               l/LEU.cif               p/PSU.cif \
a/AMP.cif               e/EPE.cif               l/LYS.cif               p/PTR.cif \
a/ANP.cif               f/FAD.cif               list/mon_lib_list.cif   s/SAH.cif \
a/ARG.cif               f/FE.cif                m/MAN.cif               s/SAM.cif \
a/ASN.cif               f/FE2.cif               m/MES.cif               s/SEP.cif \
a/ASP.cif               f/FES.cif               m/MET.cif               s/SER.cif \
a/ATP.cif               f/FMN.cif               m/MG.cif                s/SF4.cif \
b/BGC.cif               f/FMT.cif               m/MLI.cif               s/SO4.cif \
b/BMA.cif               f/FUC.cif               m/MN.cif                t/THR.cif \
b/BME.cif               g/G.cif                 m/MPD.cif               t/TPO.cif \
b/BR.cif                g/GAL.cif               m/MSE.cif               t/TRP.cif \
c/C.cif                 g/GDP.cif               n/NA.cif                t/TRS.cif \
c/CA.cif                g/GLC.cif               n/NAD.cif               t/TYR.cif \
c/CD.cif                g/GLN.cif               n/NAG.cif               u/U.cif   \
c/CIT.cif               g/GLU.cif               n/NAP.cif               u/UNK.cif \
c/CL.cif                g/GLY.cif               n/NDP.cif               u/UNX.cif \
c/CLR.cif               g/GOL.cif               n/NH2.cif               v/VAL.cif \
c/CME.cif               g/GTP.cif               n/NI.cif                z/ZN.cif  \
c/CO.cif                h/HEC.cif               n/NO3.cif \
c/COA.cif               h/HEM.cif               o/OMG.cif \
c/CSO.cif               h/HIS.cif               p/PCA.cif \
b/BTB.cif               e/EBE.cif               p/PE8.cif               t/TAM.cif

EXTRA_DIST = $(cootmonomers_files)

install-data-local:
	@for file in $(cootmonomers_files); do \
        cif_dir=$$(dirname $$file) ; \
	mkdir -p "$(DESTDIR)$(cootmonomersdir)/$$cif_dir" ; \
	install -m 644 "$(srcdir)/$$file" "$(DESTDIR)$(cootmonomersdir)/$$file"; \
	done
