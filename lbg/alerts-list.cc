
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// alerts from SilOS-it's biscu-it.
//

#include "lbg.hh"

typedef std::pair<std::string, std::string> p;

std::vector<std::pair<std::string, std::string> >
lbg_info_t::alert_smarts() const {

   // missing
   //
   // [C,c](=O)[C,c](=O) diketo group
   // SS disulphide
   // four-member lactone C1(=O)OCC1
   // hydrazine N[NH2]
   // imine C=[N!R]
   // N oxide [NX2,nX3][OX1] N oxide
   // N-C-halo NC[F,Cl,Br,I]
   // N-nitroso [#7]-N-O
   // perfluorinated chain [CX4](F)(F)[CX4](F)F
   // peroxide OO
   // stilbene c1ccccc1C=Cc2ccccc2
   // SC=O thioester
   

   unsigned int n_smarts = 116;
   std::pair<std::string, std::string> v[] = { 
      p("*1[O,S,N]*1", "three-membered hetrocycle"),
      p("[S,C](=[O,S])[F,Br,Cl,I]", ""),
      p("[CX4][Cl,Br,I]", "alkyl halide"),
      p("[C,c]S(=O)(=O)O[C,c]", "sulfonic acid"),
      p("[$([CH]),$(CC)]#CC(=O)[C,c]", "Michael acceptor"),
      p("[$([CH]),$(CC)]#CC(=O)O[C,c]", "Michael acceptor"),
      p("n[OH]", "N-hydroxyl pyridine"),
      p("[$([CH]),$(CC)]#CS(=O)(=O)[C,c]", "Michael acceptor"),
      p("C=C(C=O)C=O", ""),
      p("n1c([F,Cl,Br,I])cccc1", "2-halo pyridine"),
      p("[CH1](=O)", "aldehyde"),
      p("[O,o][O,o]", ""),
      p("[C;!R]=[N;!R]", ""),
      p("[N!R]=[N!R]", "diazo group"),
      p("[#6](=O)[#6](=O)", ""),
      p("[S,s][S,s]", ""),
      p("[N,n][NH2]", ""),
      p("C(=O)N[NH2]", "acyl hydrazine"),
      p("[C,c]=S", "thiocarbonyl group"),
      p("[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]", "isolated alkene"),
      p("C1(=[O,N])C=CC(=[O,N])C=C1", "chinone"),
      p("C1(=[O,N])C(=[O,N])C=CC=C1", "chinone"),
      p("a21aa3a(aa1aaaa2)aaaa3", "polycyclic aromatic hydrocarbon"),
      p("a31a(a2a(aa1)aaaa2)aaaa3", "polycyclic aromatic hydrocarbon"),
      p("a1aa2a3a(a1)A=AA=A3=AA=A2", "polycyclic aromatic hydrocarbon"),
      p("c1cc([NH2])ccc1", "aniline"),
      p("[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]", "heavy metal"),
      p("I", "iodine"),
      p("OS(=O)(=O)[O-]", "sulphate"),
      p("[N+](=O)[O-]", "nitro group"),
      p("C(=O)N[OH]", "hydroxamic acid"),
      p("C1NC(=O)NC(=O)1", "hydantoin"),
      p("[SH]", "thiol"),
      p("[S-]", "thiol"),
      p("c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]", "halogenated ring"),
      p("c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]", "halogenated ring"),
      p("[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1", "cycloheptane"),
      p("[CR1]1[CR1][CR1]cc[CR1][CR1]1", "cycloheptane"),
      p("[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2]1", "cycloheptane"),
      p("[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1", "cyclooctane"),
      p("[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1", "cyclooctane"),
      p("[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "azepane"),
      p("[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "azocane"),
      p("C#C", "carbon triple bond"),
      p("[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", "crown ether"),
      p("[$([N+R]),$([n+R]),$([N+]=C)][O-]", ""),
      p("[C,c]=N[OH]", "oxime"),
      p("[C,c]=NOC=O", "oxime"),
      p("[C,c](=O)[CX4,CR0X3,O][C,c](=O)", "beta-keto/anhydride"),
      p("c1ccc2c(c1)ccc(=O)o2", "cumarine"),
      p("[O+,o+,S+,s+]", "charged oxygen or sulfur"),
      p("N=C=O", "isocyanate"),
      p("[NX3,NX4][F,Cl,Br,I]", "N-halo"),
      p("c1ccccc1OC(=O)[#6]", "phenol ester"),
      p("[CR0]=[CR0][CR0]=[CR0]", "polyene"),
      p("[C+,c+,C-,c-]", "carbo cation/anion"),
      p("N=[N+]=[N-]", "azido group"),
      p("C12C(NC(N1)=O)CSC2", "biotin analogue"),
      p("c1c([OH])c([OH,NH2,NH])ccc1", "catechol"),
      p("P", "phosphor"),
      p("[N,O,S]C#N", "cyanate/aminonitrile/thiocyanate"),
      p("C=C=O", "ketone"),
      p("[Si][F,Cl,Br,I]", "silicon halogen"),
      p("[SX2]O", "sulfur-oxygen single bond"),
      p("[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)", "triphenyl methyl-silyl"),
      p("O1CCCCC1OC2CCC3CCCCC3C2", "saponine derivative"),
      p("N=[CR0][N,n,O,S]", "imine"),
      p("[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2", "benzidine"),
      p("C=[C!r]C#N", ""),
      p("[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1", "diaminobenzene"),
      p("[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1", "diaminobenzene"),
      p("[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", "diaminobenzene"),
      p("[OH]c1ccc([OH,NH2,NH])cc1", "hydroquinone"),
      p("c1ccccc1OC(=O)O", "phenol carbonate"),
      p("[SX2H0][N]", "sulfur-nitrogen single bond"),
      p("c12ccccc1(SC(S)=N2)", "thiobenzothiazole"),
      p("c12ccccc1(SC(=S)N2)", "thiobenzothiazole"),
      p("c1nnnn1C=O", "amidotetrazole"),
      p("s1c(S)nnc1NC=O", "N-acyl-2-amino-5-mercapto-1,3,4-thiadiazole"),
      p("S1C=CSC1=S", "methylidene-1,3-dithiole"),
      p("C(=O)Onnn", "ester of HOBT"),
      p("OS(=O)(=O)C(F)(F)F", "triflate"),
      p("N#CC[OH]", "cyanhydrin"),
      p("N#CC(=O)", "acyl cyanide"),
      p("S(=O)(=O)C#N", ""),
      p("N[CH2]C#N", "cyanamide"),
      p("C1(=O)NCC1", ""),
      p("S(=O)(=O)[O-,OH]", "sulfonic acid"),
      p("NC[F,Cl,Br,I]", ""),
      p("C=[C!r]O", ""),
      p("[NX2+0]=[O+0]", ""),
      p("[OR0,NR0][OR0,NR0]", "oxygen-nitrogen single bond"),
      p("C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", "more than 2 ester groups"),
      p("[CX2R0][NX3R0]", "enamine"),
      p("c1ccccc1[C;!R]=[C;!R]c2ccccc2", ""),
      p("[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]", "het-C-het not in ring"),
      p("[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]", "quaternary nitrogen"),
      p("[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", "quaternary nitrogen"),
      p("[*]=[N+]=[*]", "quaternary nitrogen"),
      p("[SX3](=O)[O-,OH]", "sulfinic acid"),
      p("N#N", "azo group"),
      p("F.F.F.F", ""),
      p("[R0;D2][R0;D2][R0;D2][R0;D2]", "aliphatic long chain"),
      p("[cR,CR]~C(=O)NC(=O)~[cR,CR]", "phthalimide"),
      p("C=!@CC=[O,S]", "Michael acceptor"),
      p("[#6,#8,#16][C,c](=O)O[C,c]", ""),
      p("c[C;R0](=[O,S])[C,c]", ""),
      p("c[SX2][C;!R]", ""),
      p("C=C=C", ""),
      p("c1nc([F,Cl,Br,I,S])ncc1", ""),
      p("c1ncnc([F,Cl,Br,I,S])c1", ""),
      p("c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])", ""),
      p("[C,c]S(=O)(=O)c1ccc(cc1)F", ""),
      p("[15N]", "radioactive"),
      p("[13C]", "radioactive"),
      p("[18O]", "radioactive"),
      p("[34S]", "radioactive")
   };

   std::vector<std::pair<std::string, std::string> > vv(n_smarts);
   for (unsigned int i=0; i<n_smarts; i++)
      vv[i] = v[i];
   return vv;
}

#endif

