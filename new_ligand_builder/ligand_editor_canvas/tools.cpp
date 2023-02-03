#include "tools.hpp"

using namespace coot::ligand_editor_canvas;

ActiveTool::ActiveTool() noexcept {
    this->variant = ActiveTool::Variant::None;
}

ActiveTool::Variant ActiveTool::get_variant() const noexcept {
    return this->variant;
}