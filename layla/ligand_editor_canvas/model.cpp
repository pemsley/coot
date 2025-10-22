/* layla/ligand_editor_canvas/model.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "model.hpp"
#include "render.hpp"
#include "../qed.hpp"
#include <exception>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MolOps.h>
#include <cmath>
#include <boost/range/iterator_range.hpp>
#include <string>
#include <tuple>
#include <utility>
#include <cctype>
#include <complex>
#ifdef __EMSCRIPTEN__
#include "../../lhasa/glog_replacement.hpp"
#endif

using namespace coot::ligand_editor_canvas;

const float CanvasMolecule::ATOM_HITBOX_RADIUS = 12.f;
const float CanvasMolecule::BOND_LINE_SEPARATION = 0.3f;
const float CanvasMolecule::BOND_DISTANCE_BOUNDARY = 10.f;
const float CanvasMolecule::BASE_SCALE_FACTOR = 30.f;
// 33.75 degrees
const float CanvasMolecule::VERTICAL_SUPERATOM_ANGLE_THRESHOLD = 11.f/48.f * M_PI;

const char* coot::ligand_editor_canvas::display_mode_to_string(DisplayMode mode) noexcept {
    switch (mode) {
        case DisplayMode::AtomIndices:{
            return "Atom Indices";
        }
        case DisplayMode::AtomNames:{
            return "Atom Names";
        }
        default:
        case DisplayMode::Standard:{
            return "Standard";
        }
    }
}
std::optional<DisplayMode> coot::ligand_editor_canvas::display_mode_from_string(const char* value_raw) noexcept {
    std::string value(value_raw);
    if(value == "Standard") {
        return DisplayMode::Standard;
    } else if (value == "Atom Indices") {
        return DisplayMode::AtomIndices;
    } else if (value == "Atom Names") {
        return DisplayMode::AtomNames;
    } else {
        return std::nullopt;
    }
}

CanvasMolecule::MaybeAtomOrBond CanvasMolecule::resolve_click(int x, int y, float canvas_scale) const noexcept {
    float scale = BASE_SCALE_FACTOR * canvas_scale;
    auto x_offset = this->x_canvas_translation * scale;
    auto y_offset = this->y_canvas_translation * scale;
    // atoms first 
    for(const auto& atom: this->atoms) {
        float atom_x = atom.x * scale + x_offset;
        float atom_y = atom.y * scale + y_offset;
        // Circle equation. Checks if click coords are within atom's radius
        if (ATOM_HITBOX_RADIUS * ATOM_HITBOX_RADIUS >= std::pow(atom_x - x,2.f) + std::pow(atom_y - y,2.f)) {
            return atom;
        }
    }
    // then bonds
    for(const auto& bond: this->bonds) {
        // 1. Find the point lying in the middle of the segment representing the bond
        float first_atom_x = bond->first_atom_x * scale + x_offset;
        float first_atom_y = bond->first_atom_y * scale + y_offset;
        float second_atom_x = bond->second_atom_x * scale + x_offset;
        float second_atom_y = bond->second_atom_y * scale + y_offset;
        float bond_center_x = (first_atom_x + second_atom_x) / 2.0f;
        float bond_center_y = (first_atom_y + second_atom_y) / 2.0f;
        // 2. Compute its' distance to the edges and let it be the radius of a circle
        //    inside of which the click location has to be found (as a pre-condition)
        float bounding_circle_radius_squared = std::pow(first_atom_x - bond_center_x,2.f) + std::pow(first_atom_y - bond_center_y,2.f);
        // If our distance from the bond's center point is greater than the radius
        // then this means that we're outside of the circle of interest and so the current bond can be skipped.
        // Operating on squared values saves us from the expensive square root operation.
        if (bounding_circle_radius_squared < std::pow(bond_center_x - x,2.f) + std::pow(bond_center_y - y,2.f)) {
            continue;
        }
        // 3. Use the formula for distance to a point from a line to determine if
        //    the click location qualifies as matching.

        float distance = std::fabs((second_atom_x - first_atom_x)*(first_atom_y - y) - (first_atom_x - x)*(second_atom_y - first_atom_y)) 
            / std::sqrt(std::pow(first_atom_x - second_atom_x,2.f) + std::pow(first_atom_y - second_atom_y,2.f));
        // We've got a match
        if (BOND_DISTANCE_BOUNDARY >= distance) {
            return *bond;
        }
    }
    // nothing matches the click
    return std::nullopt;
}

void CanvasMolecule::apply_canvas_translation(int delta_x, int delta_y, float canvas_scale) noexcept {
    float scale = BASE_SCALE_FACTOR * canvas_scale;
    this->x_canvas_translation += (float) delta_x / scale;
    this->y_canvas_translation += (float) delta_y / scale;
}

std::pair<float,float> CanvasMolecule::get_on_screen_coords(float x, float y, const std::pair<int, int>& viewport_offset, float canvas_scale) const noexcept {
    float scale = BASE_SCALE_FACTOR * canvas_scale;
    auto x_offset = this->x_canvas_translation * scale - viewport_offset.first;
    auto y_offset = this->y_canvas_translation * scale - viewport_offset.second;
    return std::make_pair(x * scale + x_offset, y * scale + y_offset);
}

std::optional<std::pair<float,float>> CanvasMolecule::get_on_screen_coords_of_atom(unsigned int atom_idx, const std::pair<int, int>& viewport_offset, float canvas_scale) const noexcept {
    if(this->atoms.size() <= atom_idx) {
        return std::nullopt;
    }
    const Atom& a = this->atoms[atom_idx];
    return this->get_on_screen_coords(a.x, a.y, viewport_offset, canvas_scale);
}

graphene_rect_t CanvasMolecule::get_on_screen_bounding_rect(const std::pair<int, int>& viewport_offset, float canvas_scale) const noexcept {
    float scale = BASE_SCALE_FACTOR * canvas_scale;
    auto x_offset = this->x_canvas_translation * scale;
    auto y_offset = this->y_canvas_translation * scale;
    graphene_rect_t ret;
    ret.origin.x = this->bounding_atom_coords.first.x * scale + x_offset;
    ret.origin.y = this->bounding_atom_coords.first.y * scale + y_offset;
    ret.size.width = (this->bounding_atom_coords.second.x - this->bounding_atom_coords.first.x) * scale + x_offset;
    ret.size.height = (this->bounding_atom_coords.second.y - this->bounding_atom_coords.first.y) * scale + y_offset;
    return ret;
}

std::optional<CanvasMolecule::QEDInfo> CanvasMolecule::get_qed_info() const noexcept {
    return this->qed_info;
}


void CanvasMolecule::perform_flip(FlipMode flip_mode) {
    for(auto& atom: this->cached_atom_coordinate_map.value()) {
        if(flip_mode == FlipMode::Horizontal) {
            atom.second.x *= -1;
        } else {
            atom.second.y *= -1;
        }
    }
}

void CanvasMolecule::rotate_by_angle(double radians) {
    std::complex<double> rotation_mult(std::cos(radians),std::sin(radians));
    for(auto& atom: this->cached_atom_coordinate_map.value()) {
        std::complex<double> atom_cn(atom.second.x,atom.second.y);
        atom_cn *= rotation_mult;
        atom.second.x = atom_cn.real();
        atom.second.y = atom_cn.imag();
    }
}

std::tuple<float,float,float> CanvasMolecule::hightlight_to_rgb(CanvasMolecule::HighlightType htype) noexcept {
    switch (htype) {
        case CanvasMolecule::HighlightType::Edition: {
            return std::make_tuple(1.0, 0.5, 1.0);
        }
        case CanvasMolecule::HighlightType::Error: {
            return std::make_tuple(1.0, 0.0, 0.0);
        }
        case CanvasMolecule::HighlightType::Selection: {
            return std::make_tuple(0.0, 0.75, 1.0);
        }
        default:
        case CanvasMolecule::HighlightType::Hover: {
            return std::make_tuple(0.0, 1.0, 0.5);
        }
    }
}

std::optional<CanvasMolecule::HighlightType> CanvasMolecule::determine_dominant_highlight(CanvasMolecule::highlight_t hcode) noexcept {
    if (hcode == 0) {
        return std::nullopt;
    }
    auto has_highlight = [hcode](HighlightType h){
        return hcode & static_cast<highlight_t>(h);
    };

    if(has_highlight(HighlightType::Hover)) {
        return HighlightType::Hover;
    } else if(has_highlight(HighlightType::Edition)) {
        return HighlightType::Edition;
    } else if(has_highlight(HighlightType::Selection)) {
        return HighlightType::Selection;
    } else if(has_highlight(HighlightType::Error)) {
        return HighlightType::Error;
    }
    return std::nullopt;
}

std::tuple<float,float,float> CanvasMolecule::atom_color_to_rgb(CanvasMolecule::AtomColor color) noexcept {
    switch (color) {
        case AtomColor::Green:{
            return std::make_tuple(0.0,0.75,0.0);
        }
        case AtomColor::Blue:{
            return std::make_tuple(0.0,0.0,1.0);
        }
        case AtomColor::Red:{
            return std::make_tuple(1.0,0.0,0.0);
        }
        case AtomColor::Brown:{
            return std::make_tuple(0.5,0.5,0.0);
        }
        case AtomColor::DarkRed:{
            return std::make_tuple(0.5,0.0,0.0);
        }
        case AtomColor::Orange:{
            return std::make_tuple(1.0,0.5,0.0);
        }
        case AtomColor::DarkBlue:{
            return std::make_tuple(0.0,0.0,0.5);
        }
        case AtomColor::Black:
        default: {
            return std::make_tuple(0.0,0.0,0.0);
        }
    }
}

std::string CanvasMolecule::atom_color_to_html(CanvasMolecule::AtomColor color) noexcept {
    switch (color) {
        case AtomColor::Green:{
            return "#00C000";
        }
        case AtomColor::Blue:{
            return "#0000FF";
        }
        case AtomColor::Red:{
            return "#FF0000";
        }
        case AtomColor::Brown:{
            return "#808000";
        }
        case AtomColor::DarkRed:{
            return "#800000";
        }
        case AtomColor::Orange:{
            return "#FF8000";
        }
        case AtomColor::DarkBlue:{
            return "#000080";
        }
        case AtomColor::Black:
        default: {
            return "#000000";
        }
    }
}

CanvasMolecule::AtomColor CanvasMolecule::atom_color_from_rdkit(const RDKit::Atom * atom) noexcept {
    auto atomic_number = atom->getAtomicNum();
    switch(atomic_number) {
        // Nitrogen
        case 7: {
            return AtomColor::Blue;
        }
        // Oxygen
        case 8: {
            return AtomColor::Red;
        }
        // Phosphorus
        case 15: {
            return AtomColor::Orange;
        }
        // Sulphur
        case 16: {
            return AtomColor::Brown;
        }
        // Flourine
        case 9:
        // Chlorine
        case 17: {
            return AtomColor::Green;
        }
        // Bromine
        case 35: {
            return AtomColor::DarkRed;
        }
        // Iodine
        case 53: {
            return AtomColor::DarkBlue;
        }
        // Carbon
        case 6:
        default: {
            return AtomColor::Black;
        }
    }
}

CanvasMolecule::Atom::Appendix::Appendix() noexcept 
    :charge(0),
    reversed(false),
    vertical(false) {
    
}

std::pair<float,float> CanvasMolecule::Bond::get_vector() const noexcept {
    float bond_vector_x = second_atom_x - first_atom_x;
    float bond_vector_y = second_atom_y - first_atom_y;
    return std::make_pair(bond_vector_x,bond_vector_y);
}

std::pair<float,float> CanvasMolecule::Bond::get_versor() const noexcept {
    auto [bond_vector_x, bond_vector_y] = this->get_vector();
    float bond_vector_len = std::sqrt(std::pow(bond_vector_x,2.f) + std::pow(bond_vector_y,2.f));
    if (bond_vector_len == 0) {
        return std::make_pair(0.f, 0.f);
    } else {
        return std::make_pair(bond_vector_x/bond_vector_len,bond_vector_y/bond_vector_len);
    }
}

std::pair<float,float> CanvasMolecule::Bond::get_perpendicular_versor() const noexcept {
    auto [x, y] = this->get_versor();
    return std::make_pair(-y,x);
}

float CanvasMolecule::Bond::get_length() const noexcept {
    auto [bond_vector_x, bond_vector_y] = this->get_vector();
    return std::sqrt(std::pow(bond_vector_x,2.f) + std::pow(bond_vector_y,2.f));
}

void CanvasMolecule::draw(impl::Renderer& ren, DisplayMode display_mode, const std::pair<int, int>& viewport_offset, float canvas_scale) const noexcept {
    impl::MoleculeRenderContext renctx(*this, ren, display_mode, viewport_offset, canvas_scale);
    renctx.draw_atoms();
    renctx.draw_bonds();
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol, bool allow_invalid_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->cached_atom_coordinate_map = std::nullopt;
    this->bounding_atom_coords = std::make_pair(RDGeom::Point2D(0,0),RDGeom::Point2D(0,0));
    this->lower_from_rdkit(!allow_invalid_mol);
    this->x_canvas_translation = 0;
    this->y_canvas_translation = 0;
}

void CanvasMolecule::update_source_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = rdkit_mol;
}

CanvasMolecule::BondType CanvasMolecule::bond_type_from_rdkit(RDKit::Bond::BondType rdkit_bond) {
    switch (rdkit_bond) {
        case RDKit::Bond::SINGLE: {
            return BondType::Single;
        }
        case RDKit::Bond::DOUBLE: {
            return BondType::Double;
        }
        case RDKit::Bond::TRIPLE: {
            return BondType::Triple;
        }
        case RDKit::Bond::AROMATIC:{
            throw std::runtime_error("An aromatic bond remained after kekulization!");
        }
        case RDKit::Bond::UNSPECIFIED:
        case RDKit::Bond::QUADRUPLE:
        case RDKit::Bond::QUINTUPLE:
        case RDKit::Bond::HEXTUPLE:
        case RDKit::Bond::ONEANDAHALF:
        case RDKit::Bond::TWOANDAHALF:
        case RDKit::Bond::THREEANDAHALF:
        case RDKit::Bond::FOURANDAHALF:
        case RDKit::Bond::FIVEANDAHALF:
        case RDKit::Bond::IONIC:
        case RDKit::Bond::HYDROGEN:
        case RDKit::Bond::THREECENTER:
        case RDKit::Bond::DATIVEONE:
        case RDKit::Bond::DATIVE:
        case RDKit::Bond::DATIVEL:
        case RDKit::Bond::DATIVER:
        case RDKit::Bond::OTHER:
        case RDKit::Bond::ZERO:
        default: {
            g_warning("Unhandled RDKit bond type: %i", rdkit_bond);
            return BondType::Single;
        }
    }
}

CanvasMolecule::BondGeometry CanvasMolecule::bond_geometry_from_rdkit(RDKit::Bond::BondDir dir) noexcept {
    //Not handled here: EITHERDOUBLE
    switch (dir) {
        default:{
            g_warning("Unhandled RDKit bond geometry: %i! Falling back to flat.", dir);
        }
        case RDKit::Bond::NONE:{
            return BondGeometry::Flat;
        }
        case RDKit::Bond::UNKNOWN:{
            return BondGeometry::Unspecified;
        }
        case RDKit::Bond::BEGINWEDGE:{
            return BondGeometry::WedgeTowardsSecond;
        }
        case RDKit::Bond::BEGINDASH:{
            return BondGeometry::DashedTowardsSecond;
        }
        case RDKit::Bond::ENDDOWNRIGHT:{
            // todo: make sure that this makes sense
            return BondGeometry::WedgeTowardsFirst;
        }
        case RDKit::Bond::ENDUPRIGHT:{
            // todo: make sure that this makes sense
            return BondGeometry::DashedTowardsFirst;
        }
    }
}

CanvasMolecule::BondGeometry CanvasMolecule::cycle_bond_geometry(BondGeometry geom) noexcept {
    switch (geom) {
        default: 
        case BondGeometry::Flat: {
            return BondGeometry::WedgeTowardsFirst;
        }
        case BondGeometry::Unspecified: {
            return BondGeometry::Flat;
        }
        case BondGeometry::WedgeTowardsFirst: {
            return BondGeometry::WedgeTowardsSecond;
        }
        case BondGeometry::WedgeTowardsSecond: {
            return BondGeometry::DashedTowardsFirst;
        }
        case BondGeometry::DashedTowardsFirst: {
            return BondGeometry::DashedTowardsSecond;
        }
        case BondGeometry::DashedTowardsSecond: {
            return BondGeometry::Unspecified;
        }
    }
}

RDKit::Bond::BondDir CanvasMolecule::bond_geometry_to_rdkit(BondGeometry geom) noexcept {
    switch (geom) {
        default: 
        case BondGeometry::Flat: {
            return RDKit::Bond::BondDir::NONE;
        }
        case BondGeometry::Unspecified: {
            return RDKit::Bond::BondDir::UNKNOWN;
        }
        case BondGeometry::WedgeTowardsFirst: {
            return RDKit::Bond::BondDir::ENDDOWNRIGHT;
        }
        case BondGeometry::WedgeTowardsSecond: {
            return RDKit::Bond::BondDir::BEGINWEDGE;
        }
        case BondGeometry::DashedTowardsFirst: {
            return RDKit::Bond::BondDir::ENDUPRIGHT;
        }
        case BondGeometry::DashedTowardsSecond: {
            return RDKit::Bond::BondDir::BEGINDASH;
        }
    }
}


RDKit::Bond::BondType CanvasMolecule::bond_type_to_rdkit(CanvasMolecule::BondType ty) noexcept {
    switch (ty) {
        default:
        case BondType::Single:{
            return RDKit::Bond::SINGLE;
        }
        case BondType::Double:{
            return RDKit::Bond::DOUBLE;
        }
        case BondType::Triple:{
            return RDKit::Bond::TRIPLE;
        }
    }
}


RDGeom::INT_POINT2D_MAP CanvasMolecule::compute_molecule_geometry() const {
    // The following code is heavily based on RDKit documentation.

    const RDGeom::INT_POINT2D_MAP* previous_coordinate_map = nullptr;
    
    // Pruned coordinate map will only be instantiated if the last coordinate map
    // contains atoms which are no longer in the molecule (and thus need to be removed).
    std::unique_ptr<RDGeom::INT_POINT2D_MAP> pruned_previous_coordinate_map = nullptr;

    if (this->cached_atom_coordinate_map.has_value()) {
        // g_debug("Computing 2D coords using a reference");
        const RDGeom::INT_POINT2D_MAP& prev_coord_map_ref = this->cached_atom_coordinate_map.value();
        previous_coordinate_map = &prev_coord_map_ref;
        // We need to make sure that each atom in the last_atom_coordinate_map still exists in the molecule.
        // If it doesn't, RDDepict::compute2DCoords() throws an exception. We cannot let this happen.
        for(const auto& [atom_idx,_point]: prev_coord_map_ref) {
            try {
                auto* _atom_ptr = this->rdkit_molecule->getAtomWithIdx(atom_idx);
                // All good.
            } catch(...) { // Atom does not exist
                g_info("Atom with id=%u missing in the molecule, but found in 2D coordinate reference. It will be omitted.",atom_idx);
                // If there's any atom which is no longer there, we need to copy the whole coordinate map
                // into pruned_previous_coordinate_map and then remove the atom there.
                if(! pruned_previous_coordinate_map) {
                    // Copy the original coordinate map
                    pruned_previous_coordinate_map = std::make_unique<RDGeom::INT_POINT2D_MAP>(this->cached_atom_coordinate_map.value());
                    // Also we need to update the previous_coordinate_map pointer to refer to our pruned map.
                    previous_coordinate_map = pruned_previous_coordinate_map.get();
                }
                pruned_previous_coordinate_map->erase(pruned_previous_coordinate_map->find(atom_idx));
            }
        }
    } else {
        g_info("Computing fresh 2D coords (without previous reference).");
    }

    try {
        RDDepict::compute2DCoords(*this->rdkit_molecule,previous_coordinate_map,true,true);
    } catch(std::exception& e) {
        throw std::runtime_error(std::string("Failed to compute 2D coords with RDKit! ")+e.what());
    }


    RDKit::MatchVectType matchVect;
    if(!RDKit::SubstructMatch(*this->rdkit_molecule, *this->rdkit_molecule, matchVect)) {
        throw std::runtime_error("SubstractMatch failed.");
    }

    // Maps atom indices to 2D points
    RDGeom::INT_POINT2D_MAP coordinate_map;

    RDKit::Conformer& conf = this->rdkit_molecule->getConformer();
    for(auto mv: matchVect) {
        RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
        RDGeom::Point2D pt2( pt3.x , pt3.y );
        coordinate_map[mv.second] = pt2;
    }
    // what is going on here?
    // That doesn't seem to change much
    // RDDepict::compute2DCoords( *this->rdkit_molecule, &coordinate_map, true, true);

    if(coordinate_map.empty()) {
        throw std::runtime_error("RDKit coordinate mapping is empty");
    }

    return coordinate_map;
}

void CanvasMolecule::process_alignment_in_rings() {
    const auto& rings = this->rdkit_molecule->getRingInfo();
    // g_debug("Number of rings: %lu", rings->atomRings().size());

    for(const auto& ring: rings->atomRings()) {

        float ring_center_x = 0.f;
        float ring_center_y = 0.f;
        for(int atom_idx: ring) {
            ring_center_x += this->atoms.at(atom_idx).x;
            ring_center_y += this->atoms.at(atom_idx).y;
        }
        ring_center_x /= ring.size();
        ring_center_y /= ring.size();
        
        int i = 0;
        int j = 1;

        // Go over every bond
        int ring_size = ring.size();
        while (i != ring_size) {
            unsigned int atom_one_idx = ring[i];
            unsigned int atom_two_idx = ring[j];
            auto bonds_of_atom_one = this->bond_map.find(atom_one_idx);
            if(bonds_of_atom_one == this->bond_map.end()) {
                throw std::runtime_error("Critical internal error: Could not find a bond while processing rings.");
            }
            // Find iterator pointing to the bond
            auto bond = std::find_if(bonds_of_atom_one->second.begin(),bonds_of_atom_one->second.end(),[=](const auto bond){
                return (bond->first_atom_idx == atom_one_idx && bond->second_atom_idx == atom_two_idx) 
                    || (bond->first_atom_idx == atom_two_idx && bond->second_atom_idx == atom_one_idx);
            });
            if (bond == bonds_of_atom_one->second.end()) {
                throw std::runtime_error("Critical internal error: Could not find a bond while processing rings.");
            }
            auto bond_ptr = *bond;
            if (bond_ptr->type == BondType::Double) {
                float x_offset_from_center = (bond_ptr->first_atom_x + bond_ptr->second_atom_x) / 2.f - ring_center_x;
                // negative y on screen is actually "higher" so we need to flip the sign
                float y_offset_from_center = ring_center_y - (bond_ptr->first_atom_y + bond_ptr->second_atom_y) / 2.f;
                bool sign_of_x_offset_from_center = x_offset_from_center > 0.f;
                bool sign_of_y_offset_from_center = y_offset_from_center > 0.f;
                bool x_requirement = (bond_ptr->second_atom_x > bond_ptr->first_atom_x) == sign_of_y_offset_from_center;
                // negative y on screen is actually "higher" so we need to flip the sign
                bool y_requirement = (bond_ptr->second_atom_y <= bond_ptr->first_atom_y) != sign_of_x_offset_from_center;
                bool bond_direction = x_requirement && y_requirement;
                // g_debug(
                //     "Bond: %i->%i DeltaX: %f DeltaY: %f CX: %f CY: %f XO: %f SignXO: %i YO: %f SignYO: %i ReqX: %i ReqY: %i DIR: %i",
                //     bond->first_atom_idx,
                //     bond->second_atom_idx,
                //     bond->second_atom_x - bond->first_atom_x,
                //     bond->first_atom_y - bond->second_atom_y,
                //     ring_center_x,ring_center_y,
                //     x_offset_from_center,
                //     sign_of_x_offset_from_center,
                //     y_offset_from_center,
                //     sign_of_y_offset_from_center,
                //     x_requirement,
                //     y_requirement,bond_direction
                // );
                bond_ptr->bond_drawing_direction = bond_direction ? DoubleBondDrawingDirection::Primary : DoubleBondDrawingDirection::Secondary;
            }
            // Lastly, process appendices' alignment relative to ring center
            auto& atom = this->atoms.at(atom_one_idx);
            if (atom.appendix.has_value()) {
                bool should_reverse = atom.x < ring_center_x;
                if(!atom.appendix->vertical) {
                    atom.appendix->reversed = should_reverse;
                }
            }
            i++;
            j++;
            // Process the last bond
            if (j == ring_size) {
                // Loop j to point to the first atom in the ring
                j = 0;
            }
        }
    }
}

void CanvasMolecule::shorten_double_bonds() {
    typedef std::pair<const Bond*,float> bond_ptr_and_angle;
    for (auto& bond: this->bonds) {
        if(bond->type != BondType::Double) {
            continue;
        }
        if(bond->bond_drawing_direction.has_value()) {
            if(bond->bond_drawing_direction == DoubleBondDrawingDirection::Centered) {
                continue;
            }
        }
        auto find_angle_between_bonds = [&](const Bond* other_bond, bool flip){
            // We can find the angle between two bonds
            // by computing cosinus arcus (reverse cosine)
            // of ( dot product / (length a) * (length b) )
            auto [bond_vec_x,bond_vec_y] = bond->get_vector();
            auto [other_bond_vec_x,other_bond_vec_y] = other_bond->get_vector();
            if(flip) {
                other_bond_vec_x *= -1.f;
                other_bond_vec_y *= -1.f;
            }
            auto dot_product = (bond_vec_x*other_bond_vec_x) + (bond_vec_y*other_bond_vec_y);
            auto bond_length = bond->get_length();
            auto other_bond_length = other_bond->get_length();
            auto result = std::acos(dot_product/(bond_length*other_bond_length));
            return result;
        };
        // 1. Find the adjacent bond(s)
        auto find_adjacent_bonds = [this,&bond,&find_angle_between_bonds]() -> std::pair<std::vector<bond_ptr_and_angle>,std::vector<bond_ptr_and_angle>> {
            // Adjacent bonds touching the first atom
            std::vector<bond_ptr_and_angle> first_bonds;
            // Adjacent bonds touching the second atom
            std::vector<bond_ptr_and_angle> second_bonds;

            auto first_bonds_iter = this->bond_map.find(bond->first_atom_idx);
            if(first_bonds_iter != this->bond_map.end()) {
                // Going over bonds of the first atom in the currently evaluated bond
                for(const auto &i: first_bonds_iter->second) {
                    if(i->first_atom_idx == bond->first_atom_idx) {
                        if(i->second_atom_idx == bond->second_atom_idx) {
                            // We're looking at the 'bond' itself. We must skip it.
                            continue;
                        }
                        first_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), false)));
                    } else if(i->second_atom_idx == bond->first_atom_idx) {
                        // i's second atom is bond's first, so we need to flip the sign of the bond vectors
                        // so that we can correctly compute the angle between them.
                        first_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), true)));
                    } else {
                        throw std::runtime_error("Internal error: bond_map is inconsistent!");
                    }
                }
            }
            auto second_bonds_iter = this->bond_map.find(bond->second_atom_idx);
            if(second_bonds_iter != this->bond_map.end()) {
                // Going over bonds of the second atom in the currently evaluated bond
                for(const auto &i: second_bonds_iter->second) {
                    if(i->first_atom_idx == bond->second_atom_idx) {
                        // i's first atom is bond's second, so we need to flip the sign of the bond vectors
                        // so that we can correctly compute the angle between them.
                        second_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), true)));
                    } else if(i->second_atom_idx == bond->second_atom_idx) {
                        if(i->first_atom_idx == bond->first_atom_idx) {
                            // We're looking at the 'bond' itself. We must skip it.
                            continue;
                        }
                        second_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), false)));
                    } else {
                        throw std::runtime_error("Internal error: bond_map is inconsistent!");
                    }
                }
            }
            return std::make_pair(first_bonds,second_bonds);
        };
        auto compute_shortening_proportion = [&](const Bond* other_bond, float angle_between_bonds){
            // 3. Do a little trigonometry to find the length to be shortened
            auto absolute_shortened_length = BOND_LINE_SEPARATION / std::tan(angle_between_bonds/2.f);
            // 4. Find the proportion of the shortening
            auto bond_length = bond->get_length();
            return absolute_shortened_length / bond_length;
        };
        auto [first_bonds, second_bonds] = find_adjacent_bonds();

        auto element_with_smallest_angle_between_bonds = [&](const std::vector<bond_ptr_and_angle>& adjacent_bonds) -> bond_ptr_and_angle {
            auto min = std::min_element(adjacent_bonds.cbegin(),adjacent_bonds.cend(),[&](const auto& lhs, const auto& rhs){
                auto [lbond, langle] = lhs;
                auto [rbond, rangle] = rhs;
                return langle < rangle;
            });
            return std::make_pair(min->first,min->second);
        };
        if(!first_bonds.empty()) {
            auto [adjbond, angle] = element_with_smallest_angle_between_bonds(first_bonds);
            bond->first_shortening_proportion = compute_shortening_proportion(adjbond, angle);
        }
        if(!second_bonds.empty()) {
            auto [adjbond, angle] = element_with_smallest_angle_between_bonds(second_bonds);
            bond->second_shortening_proportion = compute_shortening_proportion(adjbond, angle);
        }
    }
}

void CanvasMolecule::build_internal_molecule_representation(const RDGeom::INT_POINT2D_MAP &coordinate_map) {
    // First, clear what we have

    this->atoms.clear();
    this->bonds.clear();
    this->bond_map.clear();
    this->bounding_atom_coords = std::make_pair(RDGeom::Point2D(0,0),RDGeom::Point2D(0,0));

    /// Used to avoid duplicating bonds
    std::set<unsigned int> processed_atoms_indices;
    
    // 1. Process atoms and compute bonds
    for (const auto& [atom_idx_i, plane_point] : coordinate_map) {
        unsigned int atom_idx = atom_idx_i; // just changing signedness
        const auto* rdkit_atom = this->rdkit_molecule->getAtomWithIdx(atom_idx);
        auto canvas_atom = CanvasMolecule::Atom();
        canvas_atom.color = atom_color_from_rdkit(rdkit_atom);
        if(rdkit_atom->hasProp("name")) {
            std::string atom_name;
            rdkit_atom->getProp("name", atom_name);
            canvas_atom.name = atom_name;
        }
        canvas_atom.highlight = 0;
        canvas_atom.idx = atom_idx;
        canvas_atom.symbol = rdkit_atom->getSymbol();
        canvas_atom.x = plane_point.x;
        canvas_atom.y = plane_point.y;

        if(canvas_atom.x < bounding_atom_coords.first.x) {
            bounding_atom_coords.first.x = canvas_atom.x;
        }
        if(canvas_atom.x > bounding_atom_coords.second.x) {
            bounding_atom_coords.second.x = canvas_atom.x;
        }
        if(canvas_atom.y < bounding_atom_coords.first.y) {
            bounding_atom_coords.first.y = canvas_atom.y;
        }
        if(canvas_atom.y > bounding_atom_coords.second.y) {
            bounding_atom_coords.second.y = canvas_atom.y;
        }

        auto surrounding_hydrogen_count = rdkit_atom->getTotalNumHs(false);
        auto surrounding_non_hydrogen_count = 0;
        auto charge = rdkit_atom->getFormalCharge();
        if (charge != 0) {
            Atom::Appendix ap;
            ap.charge = charge;
            canvas_atom.appendix = ap;
        }

        // Bond pointers to be stored in the `bond_map`
        std::vector<std::shared_ptr<Bond>> bonds_to_be_cached;
        // Used to determine if the 'appendix' should be 'reversed' or 'vertical'
        std::optional<std::vector<std::pair<float,float>>> coordinates_of_bonded_atoms;

        for (const auto& bond: boost::make_iterator_range(this->rdkit_molecule->getAtomBonds(rdkit_atom))) {
            // Based on `getAtomBonds` documentation.
            // Seems weird but we have to do it that way.
            const auto* bond_ptr = (*this->rdkit_molecule)[bond];
            auto  first_atom_idx = bond_ptr->getBeginAtomIdx();
            auto second_atom_idx = bond_ptr->getEndAtomIdx();
            auto the_other_atom_idx = (first_atom_idx == atom_idx) ? second_atom_idx : first_atom_idx;
            const auto* the_other_atom =  this->rdkit_molecule->getAtomWithIdx(the_other_atom_idx);
            if(the_other_atom->getSymbol() != "H") {
                surrounding_non_hydrogen_count++;
                if(coordinates_of_bonded_atoms.has_value()) {
                    coordinates_of_bonded_atoms->push_back(std::make_pair(coordinate_map.at(the_other_atom_idx).x, coordinate_map.at(the_other_atom_idx).y));
                } else {
                    coordinates_of_bonded_atoms = std::vector<std::pair<float,float>>({std::make_pair(coordinate_map.at(the_other_atom_idx).x, coordinate_map.at(the_other_atom_idx).y)});
                }
            } 
            // else {
            //     g_warning("Skipping explicit hydrogen bound to atom with idx=%u!",canvas_atom.idx);
            //     continue;
            // }

            // We don't want to have duplicate bonds of atoms that we have already processed
            // so we skip them.
            if(processed_atoms_indices.find(first_atom_idx) != processed_atoms_indices.end() 
            || processed_atoms_indices.find(second_atom_idx) != processed_atoms_indices.end()) {
                continue;
            }

            auto canvas_bond = CanvasMolecule::Bond();

            canvas_bond.first_atom_idx = first_atom_idx;
            canvas_bond.first_atom_x = coordinate_map.at(first_atom_idx).x;
            canvas_bond.first_atom_y = coordinate_map.at(first_atom_idx).y;
            
            canvas_bond.second_atom_idx = second_atom_idx;
            canvas_bond.second_atom_x = coordinate_map.at(second_atom_idx).x;
            canvas_bond.second_atom_y = coordinate_map.at(second_atom_idx).y;

            canvas_bond.highlight = 0;
            canvas_bond.type = bond_type_from_rdkit(bond_ptr->getBondType());
            canvas_bond.geometry = bond_geometry_from_rdkit(bond_ptr->getBondDir());

            auto canvas_bond_ptr = std::make_shared<Bond>(canvas_bond);
            this->bonds.push_back(canvas_bond_ptr);

            bonds_to_be_cached.push_back(canvas_bond_ptr);
            auto cached_bonds_for_other_atom = this->bond_map.find(the_other_atom_idx);
            if(cached_bonds_for_other_atom == this->bond_map.end()) {
                std::vector<std::shared_ptr<Bond>> vec;
                vec.push_back(canvas_bond_ptr);
                this->bond_map.emplace(std::pair(the_other_atom_idx,std::move(vec)));
            } else {
                cached_bonds_for_other_atom->second.push_back(canvas_bond_ptr);
            }
        }

        auto cached_bonds_for_this_atom = this->bond_map.find(atom_idx);
        if (cached_bonds_for_this_atom == this->bond_map.end()) {
            this->bond_map.emplace(std::pair(atom_idx,std::move(bonds_to_be_cached)));
        } else {
            // In the for-loop above,
            // We skip bonds with atoms that have already been processed.
            // Therefore just appending the whole vector 
            // should never result in a duplicate.
            std::move(bonds_to_be_cached.begin(),bonds_to_be_cached.end(),std::back_inserter(cached_bonds_for_this_atom->second));
        }

        bool terminus = (surrounding_non_hydrogen_count < 2);
        bool in_an_a_corner = (surrounding_non_hydrogen_count == 2);
        // If this is something which has hydrogens to be drawn
        if (canvas_atom.symbol != "H" && (canvas_atom.symbol != "C" || terminus) && surrounding_hydrogen_count > 0) {
            Atom::Appendix ap = canvas_atom.appendix.value_or(Atom::Appendix());
            ap.superatoms = "H";
            if(surrounding_hydrogen_count > 1) {
                ap.superatoms += std::to_string(surrounding_hydrogen_count);
            }

            if(coordinates_of_bonded_atoms.has_value()) {
                const auto& coords = coordinates_of_bonded_atoms.value();
                if(terminus) {
                    if(coords.size() != 1) {
                        throw std::runtime_error("Internal error: A terminus should have exactly one non-hydrogen neighbor.");
                    }
                    auto [x_coordinate_of_bonded_atom, _y_coordinate_of_bonded_atom] = coords.front();
                    float x_diff = x_coordinate_of_bonded_atom - canvas_atom.x;
                    ap.reversed = x_diff > 0.2;
                }
                else if(in_an_a_corner) {
                    if(coords.size() != 2) {
                        throw std::runtime_error("Internal error: An atom in a corner should have exactly two non-hydrogen neighbors.");
                    }
                    auto [x1, y1] = coords.front();
                    auto [x2, y2] = coords.back();
                    float x_diff_1 = x1 - canvas_atom.x;
                    float y_diff_1 = y1 - canvas_atom.y;
                    float x_diff_2 = x2 - canvas_atom.x;
                    float y_diff_2 = y2 - canvas_atom.y;

                    // Normalize vectors
                    float len_1 = std::sqrt(x_diff_1 * x_diff_1 + y_diff_1 * y_diff_1);
                    float len_2 = std::sqrt(x_diff_2 * x_diff_2 + y_diff_2 * y_diff_2);
                    x_diff_1 /= len_1;
                    y_diff_1 /= len_1;
                    x_diff_2 /= len_2;
                    y_diff_2 /= len_2;

                    float mid_x = x_diff_1 + x_diff_2;
                    float mid_y = y_diff_1 + y_diff_2;

                    float alpha = std::atan(mid_y / mid_x);
                    ap.vertical = std::abs(alpha) > M_PI_2 - CanvasMolecule::VERTICAL_SUPERATOM_ANGLE_THRESHOLD;

                    if(ap.vertical) {
                        if(mid_y > 0) {
                            ap.reversed = true;
                        }
                    } else {
                        ap.reversed = mid_x > 0.2;
                    }
                    // g_info("Mid vect=(%f,%f) Alpha=%f Vertical=%i Reversed=%i", mid_x, mid_y, alpha / (M_PI * 2) * 360, ap.vertical, ap.reversed);
                }
            }
            canvas_atom.appendix = ap;
        }

        auto setup_potential_centered_double_bond = [&] (const unsigned int &atom_idx) {

            auto bonds_of_this_atom = this->bond_map.find(atom_idx);
            if (bonds_of_this_atom == this->bond_map.end()) {
                return;
            }
            if (!terminus || bonds_of_this_atom->second.empty()) {
                return;
            }
            // This means that we only have one bond
            Bond* bond = bonds_of_this_atom->second.front().get();
            if (bond->type != BondType::Double) {
                return;
            }
            unsigned int the_other_atom_idx = (bond->first_atom_idx == atom_idx) ? bond->second_atom_idx : bond->first_atom_idx;
            if (this->rdkit_molecule->getAtomWithIdx(the_other_atom_idx)->getAtomicNum() == 6) {
                // This should always be a valid iterator at this point
                auto bonds_of_the_other_atom = this->bond_map.find(the_other_atom_idx);
                if(bonds_of_the_other_atom->second.size() != 3 && bonds_of_the_other_atom->second.size() != 1) {
                    return;
                }
            }
            bond->bond_drawing_direction = DoubleBondDrawingDirection::Centered;
        };
        setup_potential_centered_double_bond(atom_idx);

        this->atoms.push_back(std::move(canvas_atom));
        
        // Mark the atom as processed
        processed_atoms_indices.insert(atom_idx);
    }
    std::sort(this->atoms.begin(),this->atoms.end(),[](const auto& lhs, const auto& rhs){
        return lhs.idx < rhs.idx;
    });
    // Make sure that double bonds are aligned properly
    this->process_alignment_in_rings();
    this->shorten_double_bonds();
}

void CanvasMolecule::lower_from_rdkit(bool sanitize_after, bool with_qed) {

    // 2. Do the lowering

    // 2.0 Kekulize
    if(sanitize_after) {
        RDKit::MolOps::Kekulize(*this->rdkit_molecule);
    } else {
        try {
            RDKit::MolOps::Kekulize(*this->rdkit_molecule);
        } catch(std::exception& e) {
            g_warning("Could not kekulize molecule: %s", e.what());
        }
    }

    /// 2.1 Compute geometry
    auto geometry = this->compute_molecule_geometry();

    // 2.2 Build internal repr
    this->build_internal_molecule_representation(geometry);
    this->cached_atom_coordinate_map = std::move(geometry);

    // 2.3 Reverse kekulization on the original molecule after lowering.
    if (sanitize_after) {
        RDKit::MolOps::sanitizeMol(*this->rdkit_molecule);
    }

    // QED update
    if (with_qed) {
        if (sanitize_after) {
            this->update_qed_info();
        } else {
            try {
                this->update_qed_info();
            } catch(std::exception& e) {
                g_warning("Could not update QED info: %s", e.what());
            }
        }
        
    }
    // Process problematic areas
    this->process_problematic_areas(!sanitize_after);
}

void CanvasMolecule::process_problematic_areas(bool allow_invalid_molecules) {
    this->clear_highlights(HighlightType::Error);
    if(!allow_invalid_molecules) {
        return;
    }
    try {
        auto problems = RDKit::MolOps::detectChemistryProblems(*this->rdkit_molecule);
        for(const auto& eptr: problems) {
            auto* raw_eptr = eptr.get();
            auto* atom_sanitize_exception = dynamic_cast<RDKit::AtomSanitizeException*>(raw_eptr);
            if(atom_sanitize_exception) {
                
                this->add_atom_highlight(atom_sanitize_exception->getAtomIdx(), HighlightType::Error);
            }
        }
    } catch(std::exception& e) {
        g_warning("Could not process problematic areas: %s", e.what());
    }
}

void CanvasMolecule::update_qed_info() {

    using QED = layla::RDKit::QED;

    coot::layla::RDKit::QED::QEDproperties raw_props = QED::properties(*this->rdkit_molecule);
    auto qed_score_and_ads = QED::qed(*this->rdkit_molecule, raw_props);
    QEDInfo new_info;

    new_info.alogp                             = raw_props.ALOGP;
    new_info.molecular_polar_surface_area      = raw_props.PSA;
    new_info.molecular_weight                  = raw_props.MW;
    new_info.number_of_alerts                  = raw_props.ALERTS;
    new_info.number_of_aromatic_rings          = raw_props.AROM;
    new_info.number_of_hydrogen_bond_acceptors = raw_props.HBA;
    new_info.number_of_hydrogen_bond_donors    = raw_props.HBD;
    new_info.number_of_rotatable_bonds         = raw_props.ROTB;
    new_info.qed_score = qed_score_and_ads.qed_score;
    new_info.ads_mw = qed_score_and_ads.ads_mw;
    new_info.ads_alogp = qed_score_and_ads.ads_alogp;
    new_info.ads_hba = qed_score_and_ads.ads_hba;
    new_info.ads_hbd = qed_score_and_ads.ads_hbd;
    new_info.ads_psa = qed_score_and_ads.ads_psa;
    new_info.ads_rotb = qed_score_and_ads.ads_rotb;
    new_info.ads_arom = qed_score_and_ads.ads_arom;
    new_info.ads_alert = qed_score_and_ads.ads_alerts;

    g_debug("Updated QED: ALOGP=%f PSA=%f MW=%f ALERTS=%u AROM=%u HBA=%u HBD=%u ROTB=%u QED=%f",
        new_info.alogp,
        new_info.molecular_polar_surface_area,
        new_info.molecular_weight,
        new_info.number_of_alerts,
        new_info.number_of_aromatic_rings,
        new_info.number_of_hydrogen_bond_acceptors,
        new_info.number_of_hydrogen_bond_donors,
        new_info.number_of_rotatable_bonds,
        new_info.qed_score
    );

    this->qed_info = new_info;
}

void CanvasMolecule::add_atom_highlight(int atom_idx, HighlightType htype) {
    //g_debug("Highlighted atom with idx=%i",atom_idx);
    this->atoms.at(atom_idx).highlight |= static_cast<highlight_t>(htype);
}

void CanvasMolecule::add_bond_highlight(unsigned int atom_a, unsigned int atom_b, HighlightType htype) {

    //g_debug("Highlighted bond between atoms with indices %i and %i",atom_a,atom_b);
    auto bonds_for_atom_a = this->bond_map.find(atom_a);
    if (bonds_for_atom_a == this->bond_map.end()) {
        throw std::runtime_error("Bond doesn't exist");
    }
    auto target = std::find_if(bonds_for_atom_a->second.begin(), bonds_for_atom_a->second.end(), [=](const auto& bond) {
          return ((bond->second_atom_idx == atom_b) && (bond->first_atom_idx == atom_a))
              || ((bond->second_atom_idx == atom_a) && (bond->first_atom_idx == atom_b));
    });
    if (target == bonds_for_atom_a->second.end()) {
        throw std::runtime_error("Bond doesn't exist");
    }
    (*target)->highlight |= static_cast<highlight_t>(htype);
}

void CanvasMolecule::clear_highlights(HighlightType htype) {
    for(auto& bond: this->bonds) {
        auto& highlight = bond->highlight;
        highlight &= ~static_cast<highlight_t>(htype);
    }
    for(auto& atom: this->atoms) {
        auto& highlight = atom.highlight;
        highlight &= ~static_cast<highlight_t>(htype);
    }
}

void CanvasMolecule::add_highlight_to_all_bonds(HighlightType htype) {
    for(auto& bond: this->bonds) {
        auto& highlight = bond->highlight;
        highlight |= static_cast<highlight_t>(htype);
    }
}

void CanvasMolecule::clear_cached_atom_coordinate_map() {
    this->cached_atom_coordinate_map = std::nullopt;
}

void CanvasMolecule::update_cached_atom_coordinate_map_after_atom_removal(unsigned int removed_atom_idx) {
    if (this->cached_atom_coordinate_map.has_value()) {
        auto& coordinate_map = this->cached_atom_coordinate_map.value();
        // There's no point in working with an empty cache.
        // Might as well delete it.
        if (coordinate_map.empty()) {
            this->cached_atom_coordinate_map = std::nullopt;
            return;
        }
        auto to_be_removed = coordinate_map.find(removed_atom_idx);
        if (to_be_removed == coordinate_map.end()) {
            g_warning("Atom to be removed (idx=%u) does not exist in the cached coordinate map!",removed_atom_idx);
            return;
        }
        // We can now remove the atom from the cache.
        coordinate_map.erase(to_be_removed);
        // If that was the only atom in the cache, let's get rid of the cache.
        if (coordinate_map.empty()) {
            this->cached_atom_coordinate_map = std::nullopt;
            return;
        }
        // After the atom has been removed, all of the indices which are
        // greater than its' index have to be decremented.
        //
        // Now, I don't want to mess around with iterator invalidation.
        // So the easiest course of action is to just copy away all
        // of the affected elements, change them, remove them from the original map
        // and put them back there again.
        std::vector<std::pair<int,RDGeom::Point2D>> altered_elements;

        // This is the highest index AFTER the atom has been removed from the cache.
        // The only reason we get it is that we can make a smart allocation for our vector.
        auto highest_idx = coordinate_map.crbegin()->first;
        // Allocate all the space that we might need
        altered_elements.reserve(highest_idx - removed_atom_idx + 1);
        // Get the affected elements, copy them over to the vector, with decreased indices
        for(auto i = coordinate_map.upper_bound(removed_atom_idx);i != coordinate_map.end(); i++) {
            altered_elements.push_back(std::make_pair(i->first - 1,i->second));
        }
        // Now, remove the original elements from the map
        coordinate_map.erase(coordinate_map.upper_bound(removed_atom_idx),coordinate_map.end());
        // Then put them back again after they've been modified
        for(auto x: std::move(altered_elements)) {
            coordinate_map.emplace(x.first,x.second);
        }
    }
}
