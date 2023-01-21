#!/bin/bash

# Touching makefile_marker triggers build.rs which generates the header file
pushd ligand_editor_canvas && ([ -e "ligand_editor_canvas.hpp" ] && cargo build -r || (touch makefile_marker && cargo build -r)); popd