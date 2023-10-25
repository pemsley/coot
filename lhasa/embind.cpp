#include "embind.hpp"
#include "../layla/ligand_editor_canvas.hpp"
#include "../layla/utils.hpp"

EMSCRIPTEN_BINDINGS(lhasa) {
  function("remove_non_polar_hydrogens", &coot::layla::remove_non_polar_hydrogens);
  class_<CootLigandEditorCanvas>("LhasaCanvas")
    // remove the semicolon when adding functions
    .constructor<>();
    // .function("demo", &LhasaDemo::demo);
}