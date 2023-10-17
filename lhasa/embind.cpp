#include "embind.hpp"
#include "../layla/ligand_editor_canvas.hpp"

EMSCRIPTEN_BINDINGS(lhasa) {
//   function("lerp", &lerp);
  class_<CootLigandEditorCanvas>("LhasaCanvas")
  // remove the semicolon when adding functions
    .constructor<>();
    // .function("demo", &LhasaDemo::demo);
}