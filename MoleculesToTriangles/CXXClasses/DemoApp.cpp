#include <GLFW/glfw3.h>

#include "SceneSetup.h"
#include "RendererGLSL.hpp"
#include "MyMolecule.h"
#include "CompoundSelection.h"
#include "MolecularRepresentationInstance.h"
#include "Camera.h"
#include "Light.h"
#include "CameraPort.h"
#include "ColorScheme.h"

int main(int argc, char** argv)
{
    if (argc<2) {
        std::cout << "Usage:" << argv[0] << " inputFile.pdb\n";
        exit(0);
    }
    CameraPort *cameraPort = new CameraPort();

    int RC;
    mmdb::InitMatType();
    auto mmdb = new mmdb::Manager();
    mmdb->SetFlag( mmdb::MMDBF_PrintCIFWarnings );
    RC = mmdb->ReadCoorFile (argv[1]);
    if (RC) {
        std::cout << "error could not read file " << argv[1] << std::endl;
    }

    auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "grey");
    auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "grey");
    auto dark_helix_cr  = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "forestgreen");
    auto stand_cr       = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "firebrick");
    auto brown_dna_cr   = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "brown");

    auto ss_cs = ColorScheme::colorBySecondaryScheme();

    auto        ramp_cs = ColorScheme::colorRampChainsScheme();
    auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
    auto ele_cs         = ColorScheme::colorByElementScheme();
    auto chains_cs      = ColorScheme::colorChainsScheme();
    auto handles = chains_cs->prepareForMMDB(mmdb);
    auto handles_2 = ramp_cs->prepareForMMDB(mmdb);

    ss_cs->addRule(init_cr);
    ss_cs->addRule(none_cr);
    ss_cs->addRule(dark_helix_cr);
    ss_cs->addRule(stand_cr);
    ramp_cs->addRule(brown_dna_cr);

    cameraPort->initialiseMoleculesToTriangles();
    // cameraPort->addRepresentationInstance(mmdb, ColorScheme::colorByElementScheme(), "ALL", "Ribbon");
    // cameraPort->addRepresentationInstance(mmdb, ss_cs, "//A", "MolecularSurface");
    cameraPort->addRepresentationInstance(mmdb, ss_cs, "//A", "Ribbon"); 
    //       "Ribbon", "Sticks", "Calpha", "Spheres", "Cylinders", "MolecularSurface", "VdWSurface", "AccessibleSurface", "HydrogenBonds"

    cameraPort->runLoop("Message");
    
    return 0;
}

