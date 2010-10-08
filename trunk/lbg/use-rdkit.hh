
#ifndef USE_RDKIT_HH
#define USE_RDKIT_HH

#include <vector>
#include <algorithm>
#include <string>

#define ENABLE_NLS // fixes dcgettext() header problems on including
		   // libintl.h (via RDKitBase.h etc (including boost
		   // stuff).

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/GraphMol.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <DistGeom/DistGeomUtils.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>


#endif // USE_RDKIT_HH
