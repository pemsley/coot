
#ifndef USE_RDKIT_HH
#define USE_RDKIT_HH

#include <vector>
#include <algorithm>
#include <string>

#ifndef WINDOWS_MINGW
#define ENABLE_NLS // fixes dcgettext() header problems on including
		   // libintl.h (via RDKitBase.h etc (including boost
		   // stuff).
#endif

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
/* use MolDescriptors.h rather than Lipinski.h as in rdkit Doc */
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/GraphMol.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <DistGeom/DistGeomUtils.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>

// Fix RDKIT_VERSION if not exist
#ifndef RDKIT_VERSION
/* Version check macro
   Can be used like #if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2018, 3, 1)) */
#define RDKIT_VERSION_CHECK(year, month, rev) ((year*1000)+(month*10)+(rev))

/* RDKIT_VERSION is (year*1000) + (month*10) + (rev) */
#define RDKIT_VERSION RDKIT_VERSION_CHECK(2017, 9, 3)  // a guess
#endif

#endif // USE_RDKIT_HH
