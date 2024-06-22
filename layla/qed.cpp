#include "qed.hpp"
#ifdef __EMSCRIPTEN__
#include "../lhasa/glog_replacement.hpp"
#else
#include <glib.h>
#endif
#include <cmath>
#include <type_traits>

// #include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/ChemTransforms/ChemTransforms.h>
// #include <rdkit/GraphMol/FileParsers/MolSupplier.h>
// #include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Descriptors/MolDescriptors.h>
#include <rdkit/GraphMol/Descriptors/Crippen.h>
#include <rdkit/GraphMol/Descriptors/MolSurf.h>
#include <rdkit/GraphMol/Descriptors/Lipinski.h>

namespace coot::layla::RDKit {

namespace impl {

    auto make_acceptors() -> std::vector<std::unique_ptr<const ::RDKit::ROMol>> {
        try {
            const unsigned int smarts_array_len = 11;
            const char* smarts_array[smarts_array_len] = {
                "[oH0;X2]", "[OH1;X2;v2]", "[OH0;X2;v2]", 
                "[OH0;X1;v2]", "[O-;X1]", "[SH0;X2;v2]", "[SH0;X1;v2]",
                "[S-;X1]", "[nH0;X2]", "[NH0;X1;v3]", "[$([N;+0;X3;v3]);!$(N[C,S]=O)]"
            };
            std::vector<std::unique_ptr<const ::RDKit::ROMol>> ret;
            for(unsigned int i = 0; i < smarts_array_len; i++) {
                ret.push_back(std::unique_ptr<const ::RDKit::ROMol>(::RDKit::SmartsToMol(smarts_array[i])));
            }
            return ret;
        }
        catch(const std::exception& e) {
            g_warning("QED make_acceptors(): failed to initialize static const: %s", e.what());
            return std::vector<std::unique_ptr<const ::RDKit::ROMol>>();
        }
    }
    auto make_structural_alerts() -> std::vector<std::unique_ptr<const ::RDKit::ROMol>> {
        try {
            const unsigned int smarts_array_len = 116;
            const char* smarts_array[] = {
                "*1[O,S,N]*1", "[S,C](=[O,S])[F,Br,Cl,I]", "[CX4][Cl,Br,I]", "[#6]S(=O)(=O)O[#6]",
                "[$([CH]),$(CC)]#CC(=O)[#6]", "[$([CH]),$(CC)]#CC(=O)O[#6]", "n[OH]",
                "[$([CH]),$(CC)]#CS(=O)(=O)[#6]", "C=C(C=O)C=O", "n1c([F,Cl,Br,I])cccc1", "[CH1](=O)", "[#8][#8]",
                "[C;!R]=[N;!R]", "[N!R]=[N!R]", "[#6](=O)[#6](=O)", "[#16][#16]", "[#7][NH2]", "C(=O)N[NH2]",
                "[#6]=S", "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
                "C1(=[O,N])C=CC(=[O,N])C=C1", "C1(=[O,N])C(=[O,N])C=CC=C1", "a21aa3a(aa1aaaa2)aaaa3",
                "a31a(a2a(aa1)aaaa2)aaaa3", "a1aa2a3a(a1)A=AA=A3=AA=A2", "c1cc([NH2])ccc1",
                "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]", "I", "OS(=O)(=O)[O-]", "[N+](=O)[O-]", "C(=O)N[OH]",
                "C1NC(=O)NC(=O)1", "[SH]", "[S-]", "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
                "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]", "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
                "[CR1]1[CR1][CR1]cc[CR1][CR1]1", "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
                "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
                "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "C#C",
                "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
                "[#6]=N[OH]", "[#6]=NOC=O", "[#6](=O)[CX4,CR0X3,O][#6](=O)", "c1ccc2c(c1)ccc(=O)o2",
                "[O+,o+,S+,s+]", "N=C=O", "[NX3,NX4][F,Cl,Br,I]", "c1ccccc1OC(=O)[#6]", "[CR0]=[CR0][CR0]=[CR0]",
                "[C+,c+,C-,c-]", "N=[N+]=[N-]", "C12C(NC(N1)=O)CSC2", "c1c([OH])c([OH,NH2,NH])ccc1", "P",
                "[N,O,S]C#N", "C=C=O", "[Si][F,Cl,Br,I]", "[SX2]O", "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
                "O1CCCCC1OC2CCC3CCCCC3C2", "N=[CR0][N,n,O,S]",
                "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
                "C=[C!r]C#N", "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
                "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
                "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", "[OH]c1ccc([OH,NH2,NH])cc1",
                "c1ccccc1OC(=O)O", "[SX2H0][N]", "c12ccccc1(SC(S)=N2)", "c12ccccc1(SC(=S)N2)", "c1nnnn1C=O",
                "s1c(S)nnc1NC=O", "S1C=CSC1=S", "C(=O)Onnn", "OS(=O)(=O)C(F)(F)F", "N#CC[OH]", "N#CC(=O)",
                "S(=O)(=O)C#N", "N[CH2]C#N", "C1(=O)NCC1", "S(=O)(=O)[O-,OH]", "NC[F,Cl,Br,I]", "C=[C!r]O",
                "[NX2+0]=[O+0]", "[OR0,NR0][OR0,NR0]", "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", "[CX2R0][NX3R0]",
                "c1ccccc1[C;!R]=[C;!R]c2ccccc2", "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
                "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
                "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", "[*]=[N+]=[*]", "[SX3](=O)[O-,OH]",
                "N#N", "F.F.F.F", "[R0;D2][R0;D2][R0;D2][R0;D2]", "[cR,CR]~C(=O)NC(=O)~[cR,CR]", "C=!@CC=[O,S]",
                "[#6,#8,#16][#6](=O)O[#6]", "c[C;R0](=[O,S])[#6]", "c[SX2][C;!R]", "C=C=C",
                "c1nc([F,Cl,Br,I,S])ncc1", "c1ncnc([F,Cl,Br,I,S])c1", "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
                "[#6]S(=O)(=O)c1ccc(cc1)F", "[15N]", "[13C]", "[18O]", "[34S]"
            };
            std::vector<std::unique_ptr<const ::RDKit::ROMol>> ret;
            for(unsigned int i = 0; i < smarts_array_len; i++) {
                ret.push_back(std::unique_ptr<const ::RDKit::ROMol>(::RDKit::SmartsToMol(smarts_array[i])));
            }
            return ret;
        }
        catch(const std::exception& e) {
            g_warning("QED make_structural_alerts(): failed to initialize static const: %s", e.what());
            return std::vector<std::unique_ptr<const ::RDKit::ROMol>>();
        }
    }

    inline double QEDproperties_sum(const QED::QEDproperties& props) noexcept {
        return props.MW + props.ALOGP + props.HBA + props.HBD + props.PSA + props.ROTB + props.AROM + props.ALERTS;
    }
  
}

const QED::QEDproperties QED::WEIGHT_MAX  = QEDproperties({0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00});
const QED::QEDproperties QED::WEIGHT_MEAN = QEDproperties({0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95});
const QED::QEDproperties QED::WEIGHT_NONE = QEDproperties({1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00});

auto make_aliphatic_rings() ->  std::unique_ptr<const ::RDKit::ROMol> {
    auto ret = std::unique_ptr<const ::RDKit::ROMol>(nullptr);
    try {
        ret = std::unique_ptr<const ::RDKit::ROMol>(::RDKit::SmartsToMol("[$([A;R][!a])]"));
    } catch(const std::exception& e) {
        g_warning("QED make_aliphatic_rings(): failed to initialize static const: %s", e.what());
    }
    return ret;
}

const std::unique_ptr<const ::RDKit::ROMol> QED::AliphaticRings = make_aliphatic_rings();

const std::vector<std::unique_ptr<const ::RDKit::ROMol>> QED::Acceptors = impl::make_acceptors();
const std::vector<std::unique_ptr<const ::RDKit::ROMol>> QED::StructuralAlerts = impl::make_structural_alerts();
const std::vector<QED::ADSparameter> QED::adsParameters = {
// {"MW",
  ADSparameter({2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677,
               65.37051707, 104.9805561}),
// {"ALOGP",
  ADSparameter({3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154,
               0.576295591, 131.3186604}),
// {"HBA",
  ADSparameter({2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953,
               1.300669958, 148.7763046}),
// {"HBD",
  ADSparameter({1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843,
               0.920922555, 258.1632616}),
// {"PSA",
  ADSparameter({1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824,
               28.51324732, 104.5686167}),
// {"ROTB",
  ADSparameter({0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166,
               2.758063707, 105.4420403}),
// {"AROM",
  ADSparameter({3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384,
               0.375760881, 312.3372610}),
// {"ALERTS",
  ADSparameter({0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477,
               0.875193782, 417.7253140})
};

double QED::ads(double x, const ADSparameter& p) noexcept {
    double exp1 = 1 + exp(-1 * (x - p.C + p.D / 2) / p.E);
    double exp2 = 1 + exp(-1 * (x - p.C - p.D / 2) / p.F);
    double dx = p.A + p.B / exp1 * (1 - 1 / exp2);
    return dx / p.DMAX;
}

QED::QEDproperties QED::properties(const ::RDKit::ROMol& mol_raw) {
    auto mol = std::unique_ptr<const ::RDKit::ROMol>(::RDKit::MolOps::removeHs(mol_raw));

    auto has_substruct_match = [](const ::RDKit::ROMol& mol, const ::RDKit::ROMol& pattern){
        ::RDKit::MatchVectType _res;
        return ::RDKit::SubstructMatch(mol, pattern, _res);
    };

    auto get_hba = [&mol, has_substruct_match](){
        unsigned int ret = 0;
        if(Acceptors.empty()) {
            g_warning("QED: Acceptors is empty. Number of hydrogen bonds acceptors will be incorrect.");
        }
        for(const auto& pattern: Acceptors) {
            if(has_substruct_match(*mol, *pattern)) {
                std::vector<::RDKit::MatchVectType> _matches;
                auto num_of_matches = ::RDKit::SubstructMatch(*mol, *pattern, _matches);
                ret += num_of_matches;
            }
        }
        return ret;
    };
    auto get_arom = [&mol]() -> unsigned int {
        if(!AliphaticRings) {
            g_warning("QED: AromaticRings is null. Number of aromatic rings will be incorrect.");
            return 0;
        }
        auto rmol = std::unique_ptr<const ::RDKit::ROMol>(::RDKit::deleteSubstructs(*mol, *AliphaticRings));
        /// Replaces Chem.GetSSSR
        std::vector<std::vector<int>> rings;
        ::RDKit::MolOps::findSSSR(*rmol, rings);
        return rings.size();
    };
    auto get_alerts = [&mol, has_substruct_match](){
        unsigned int ret = 0;
        if(StructuralAlerts.empty()) {
            g_warning("QED: StructuralAlerts is empty. Number of alerts will be incorrect.");
        }
        for(const auto& alert: StructuralAlerts) {
            if(has_substruct_match(*mol, *alert)) {
                ret += 1;
            }
        }
        return ret;
    };

    auto qedProperties = QEDproperties({
        ::RDKit::Descriptors::calcAMW(*mol),// MW=rdmd._CalcMolWt(mol),
        ::RDKit::Descriptors::calcClogP(*mol), // ALOGP=Crippen.MolLogP(mol),
        static_cast<double>(get_hba()),// HBA=sum(
        // len(mol.GetSubstructMatches(pattern)) for pattern in Acceptors
        // if mol.HasSubstructMatch(pattern)),
        static_cast<double>(::RDKit::Descriptors::calcNumHBD(*mol)),// HBD=rdmd.CalcNumHBD(mol),
        ::RDKit::Descriptors::calcTPSA(*mol),// PSA=MolSurf.TPSA(mol),
        static_cast<double>(::RDKit::Descriptors::calcNumRotatableBonds(*mol, ::RDKit::Descriptors::Strict)),// ROTB=rdmd.CalcNumRotatableBonds(mol, rdmd.NumRotatableBondsOptions.Strict),
        static_cast<double>(get_arom()),// AROM=len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), AliphaticRings))),
        static_cast<double>(get_alerts())// ALERTS=sum(1 for alert in StructuralAlerts if mol.HasSubstructMatch(alert)),
    });
    // The replacement
    // AROM=Lipinski.NumAromaticRings(mol),
    // is not identical. The expression above tends to count more rings
    // N1C2=CC=CC=C2SC3=C1C=CC4=C3C=CC=C4
    // OC1=C(O)C=C2C(=C1)OC3=CC(=O)C(=CC3=C2C4=CC=CC=C4)O
    // CC(C)C1=CC2=C(C)C=CC2=C(C)C=C1  uses 2, should be 0 ?
    return qedProperties;
}


double QED::qed(const ::RDKit::ROMol& mol, std::optional<QEDproperties> qedPropertiesOpt, QEDproperties w) {
    typedef std::underlying_type<QEDPropName>::type utype;
    QEDproperties qedProperties = qedPropertiesOpt.has_value() ? qedPropertiesOpt.value() : properties(mol);

    double t;
    t  = qedProperties.MW *     log(ads(qedProperties.MW,     adsParameters[static_cast<utype>(QEDPropName::MW)]));
    t += qedProperties.ALOGP *  log(ads(qedProperties.ALOGP,  adsParameters[static_cast<utype>(QEDPropName::ALOGP)]));
    t += qedProperties.HBA *    log(ads(qedProperties.HBA,    adsParameters[static_cast<utype>(QEDPropName::HBA)]));
    t += qedProperties.HBD *    log(ads(qedProperties.HBD,    adsParameters[static_cast<utype>(QEDPropName::HBD)]));
    t += qedProperties.PSA *    log(ads(qedProperties.PSA,    adsParameters[static_cast<utype>(QEDPropName::PSA)]));
    t += qedProperties.ROTB *   log(ads(qedProperties.ROTB,   adsParameters[static_cast<utype>(QEDPropName::ROTB)]));
    t += qedProperties.AROM *   log(ads(qedProperties.AROM,   adsParameters[static_cast<utype>(QEDPropName::AROM)]));
    t += qedProperties.ALERTS * log(ads(qedProperties.ALERTS, adsParameters[static_cast<utype>(QEDPropName::ALERTS)]));

    return exp(t / impl::QEDproperties_sum(w));
}




} // namespace
