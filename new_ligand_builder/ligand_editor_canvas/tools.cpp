#include "tools.hpp"
#include "core.hpp"
#include "model.hpp"
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include "../ligand_editor_canvas.hpp"

using namespace coot::ligand_editor_canvas;

ActiveTool::ActiveTool() noexcept {
    this->variant = ActiveTool::Variant::None;
}

BondModifier::BondModifier(BondModifierMode mode) noexcept {
    this->mode = mode;
    this->molecule_idx_and_first_atom_of_new_bond = std::nullopt;
    this->is_in_drag = false;
}

bool BondModifier::is_creating_bond() const noexcept {
    return this->is_in_drag;
}

void BondModifier::begin_creating_bond(unsigned int molecule_idx, unsigned int atom_idx) noexcept {
    this->is_in_drag = true;
    this->molecule_idx_and_first_atom_of_new_bond = std::make_pair(molecule_idx,atom_idx);
}

void BondModifier::finish_creating_bond() noexcept {
    this->is_in_drag = false;
    this->molecule_idx_and_first_atom_of_new_bond = std::nullopt;
}

std::optional<std::pair<unsigned int,unsigned int>> BondModifier::get_molecule_idx_and_first_atom_of_new_bond() const noexcept {
    return this->molecule_idx_and_first_atom_of_new_bond;
}


TransformManager::TransformManager() noexcept {
    this->state = IdleState();
}

bool TransformManager::is_active() const noexcept {
    return !std::holds_alternative<IdleState>(this->state);
}

void TransformManager::begin_transform(int x, int y, Mode mode) noexcept {
    switch(mode) {
        case Mode::Rotation: {
            auto rot = RotationState();
            rot.original_rotation_pos = std::make_pair(x,y);
            rot.current_rotation_pos = std::make_pair(x, y);
            rot.last_absolute_angle = 0.f;
            this->state = rot;
            break;
        }
        case Mode::Translation: {
            auto tr = TranslationState();
            tr.prev_move_pos = std::make_pair(x,y);
            tr.current_move_pos = std::make_pair(x, y);
            this->state = tr;
            break;
        }
    }
}

double TransformManager::RotationState::get_current_absolute_angle(bool snap_to_angle) const {
    auto [x1,y1] = this->original_rotation_pos;
    auto [x2,y2] = this->current_rotation_pos;
    auto diff_x = x2 - x1;
    auto diff_y = y1 - y2;
    auto diff_original = ((double)(diff_x + diff_y)) / 125.0;
    if(snap_to_angle) {
        const double snap = 15.f / 180 * M_PI;
        int divided = diff_original / snap;
        if(divided != 0)  {
            diff_original = divided * snap;
        } else {
            diff_original = 0;
        }
    }
    return diff_original;
}

void TransformManager::update_current_cursor_pos(int x, int y, bool snap) noexcept {
    RotationState* rot = std::get_if<RotationState>(&this->state);
    if(rot) {
        rot->last_absolute_angle = rot->get_current_absolute_angle(snap);
        rot->current_rotation_pos = std::make_pair(x, y);
        return;
    }
    TranslationState* tr = std::get_if<TranslationState>(&this->state);
    if(tr) {
        tr->prev_move_pos = tr->current_move_pos;
        tr->current_move_pos = std::make_pair(x, y);
        return;
    }
}

void TransformManager::end_transform() noexcept {
    this->state = IdleState();
    this->canvas_mol_idx = std::nullopt;
}

void TransformManager::set_canvas_molecule_index(unsigned int idx) noexcept {
    this->canvas_mol_idx = idx;
}

void TransformManager::apply_current_transform_state(impl::WidgetCoreData* widget_data, bool snap_to_angle, bool about_to_end) const {
    if(!this->canvas_mol_idx.has_value()) {
        return;
    }
    auto& mol = widget_data->molecules->at(this->canvas_mol_idx.value());
    const RotationState* rot = std::get_if<RotationState>(&this->state);
    if(rot) {
        auto angle = rot->get_current_angle_diff(snap_to_angle);
        auto abs_angle = rot->get_current_absolute_angle(snap_to_angle) / M_PI * 180;
        mol.rotate_by_angle(angle);
        mol.lower_from_rdkit(!widget_data->allow_invalid_molecules);
        std::string msg;
        if (about_to_end) {
            msg = "Molecule rotated by: " + std::to_string(abs_angle) + " degrees.";
        } else {
            msg = "Current molecule rotation: " + std::to_string(abs_angle) + " degrees.";
        }
        widget_data->update_status(msg.c_str());
        return;
    }
    const TranslationState* tr = std::get_if<TranslationState>(&this->state);
    if(tr) {
        auto [offset_x,offset_y] = tr->get_current_offset();
        mol.apply_canvas_translation(offset_x, offset_y);
        return;
    }
}

std::pair<int,int> TransformManager::TranslationState::get_current_offset() const {
    auto [x1,y1] = this->prev_move_pos;
    auto [x2,y2] = this->current_move_pos;
    return std::make_pair(x2 - x1, y2 - y1);
}

double TransformManager::RotationState::get_current_angle_diff(bool snap_to_angle) const {
    auto original = this->get_current_absolute_angle(snap_to_angle);
    auto diff = original - this->last_absolute_angle;
    return diff;
}

StructureInsertion::StructureInsertion(StructureInsertion::Structure st) noexcept 
:structure(st) {

}

StructureInsertion::Structure StructureInsertion::get_structure() const noexcept {
    return this->structure;
}


CanvasMolecule::BondType BondModifier::get_target_bond_type() const noexcept {
    switch (this->mode) {
        default:
        case BondModifierMode::Single:{
            return CanvasMolecule::BondType::Single;
        }
        case BondModifierMode::Double:{
            return CanvasMolecule::BondType::Double;
        }
        case BondModifierMode::Triple:{
            return CanvasMolecule::BondType::Triple;
        }
    }
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->variant = ActiveTool::Variant::ElementInsertion;
    this->element_insertion = insertion;
}

ActiveTool::ActiveTool(BondModifier modifier) noexcept {
    this->variant = ActiveTool::Variant::BondModifier;
    this->bond_modifier = modifier;
}

ActiveTool::ActiveTool(DeleteTool deltool) noexcept {
    this->variant = ActiveTool::Variant::Delete;
    this->delete_tool = deltool;
}

ActiveTool::ActiveTool(ChargeModifier chargemod) noexcept {
    this->variant = ActiveTool::Variant::ChargeModifier;
    this->charge_modifier = chargemod;
}

ActiveTool::ActiveTool(MoveTool mov) noexcept {
    this->variant = ActiveTool::Variant::MoveTool;
    this->move_tool = mov;
}

ActiveTool::ActiveTool(StructureInsertion insertion) noexcept {
    this->variant = ActiveTool::Variant::StructureInsertion;
    this->structure_insertion = insertion;
}

ActiveTool::ActiveTool(GeometryModifier modifier) noexcept {
    this->variant = ActiveTool::Variant::GeometryModifier;
    this->geometry_modifier = modifier;
}

ActiveTool::ActiveTool(FormatTool fmt) noexcept {
    this->variant = ActiveTool::Variant::Format;
    this->format_tool = fmt;
}

ActiveTool::ActiveTool(RotateTool rot) noexcept {
    this->variant = ActiveTool::Variant::RotateTool;
    this->rotate_tool = rot;
}

ActiveTool::ActiveTool(FlipTool flip) noexcept {
    this->variant = ActiveTool::Variant::FlipTool;
    this->flip_tool = flip;
}

ActiveTool::ActiveTool(RemoveHydrogensTool rh) noexcept {
    this->variant = ActiveTool::Variant::RemoveHydrogens;
    this->remove_hydrogens_tool = rh;
}

ActiveTool::Variant ActiveTool::get_variant() const noexcept {
    return this->variant;
}

void ActiveTool::set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept {
    this->widget_data = static_cast<impl::WidgetCoreData*>(owning_widget);
}

void ActiveTool::check_variant(ActiveTool::Variant expected) const {
    if (expected != this->variant) {
        throw std::runtime_error("Unexpected ActiveTool variant.");
    }
}



void ActiveTool::insert_atom(int x, int y) {
    check_variant(Variant::ElementInsertion);
    auto& element_insertion = this->element_insertion;
    unsigned int atomic_number = element_insertion.get_atomic_number();
    auto el_name = RDKit::PeriodicTable::getTable()->getElementSymbol(atomic_number);
    g_debug("Inserting element '%u' (%s) at %i %i.",atomic_number,el_name.c_str(),x,y);
    //1. Find what we've clicked at
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            this->widget_data->begin_edition();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                g_debug("Resolved insertion destination atom: idx=%i, symbol=%s",atom.idx,atom.symbol.c_str());
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* new_atom = new RDKit::Atom(std::string(el_name));
                rdkit_mol->replaceAtom(atom.idx, new_atom);
                this->widget_data->update_status("Atom has been replaced.");
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                g_warning("TODO: Implement handling insertion at bonds (if any should happen)");
            }
            this->widget_data->finalize_edition();
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not insert atom: ") + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::alter_bond(int x, int y) {
    check_variant(Variant::BondModifier);
    BondModifier& mod = this->bond_modifier;
    
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            this->widget_data->begin_edition();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                mod.begin_creating_bond(molecule_idx, atom.idx);
            } else {
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* rdkit_bond = rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                auto old_type = rdkit_bond->getBondType();
                rdkit_bond->setBondType(CanvasMolecule::bond_type_to_rdkit(mod.get_target_bond_type()));
                try {
                    this->sanitize_molecule(*rdkit_mol.get());
                }catch(std::exception& e) {
                    // rollback
                    g_warning("Rolling back invalid molecule change that makes it unable to sanitize with the following error: %s",e.what());
                    rdkit_bond->setBondType(old_type);
                    this->sanitize_molecule(*rdkit_mol.get());
                    // rethrow
                    throw std::runtime_error(std::string("Invalid bond modification: ") + e.what());
                }
                this->widget_data->update_status("Bond has been altered.");
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
                this->widget_data->finalize_edition();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter bond: ") + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

bool ActiveTool::is_creating_bond() const {
    check_variant(Variant::BondModifier);
    const BondModifier& mod = this->bond_modifier;
    return mod.is_creating_bond();
}

void ActiveTool::finish_creating_bond(int x, int y) {
    check_variant(Variant::BondModifier);
    BondModifier& mod = this->bond_modifier;
    auto click_result = this->widget_data->resolve_click(x, y);
    auto [original_molecule_idx, first_atom_idx] = mod.get_molecule_idx_and_first_atom_of_new_bond().value();
    mod.finish_creating_bond();

    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto second_atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                if(original_molecule_idx != molecule_idx) {
                    this->widget_data->update_status("Cannot create bond between different molecules!");
                    this->widget_data->rollback_current_edition();
                    return;
                }
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                if(first_atom_idx == second_atom.idx) {
                    auto* new_atom = new RDKit::Atom(6);
                    auto new_atom_idx = rdkit_mol->addAtom(new_atom,false,true);
                    rdkit_mol->addBond(new_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(mod.get_target_bond_type()));
                    g_info("New atom added: idx=%i",new_atom_idx);
                    this->widget_data->update_status("New carbon atom added.");
                } else {
                    rdkit_mol->addBond(first_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(mod.get_target_bond_type()));
                    this->widget_data->update_status("Created new bond between atoms.");
                }
                this->sanitize_molecule(*rdkit_mol.get());
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
                this->widget_data->finalize_edition();
            } else {
                this->widget_data->update_status("Can't link bond to a bond!");
                this->widget_data->rollback_current_edition();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter/create bond: ") + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        std::string msg = "The new bond goes nowhere.";
        this->widget_data->update_status(msg.c_str());
        this->widget_data->rollback_current_edition();
    }
}


void ActiveTool::alter_geometry(int x, int y) {
    check_variant(Variant::GeometryModifier);
    GeometryModifier& mod = this->geometry_modifier;
    
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                g_warning("The GeometryModifier tool does not operate on atoms. Nothing to do.");
                //auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
            } else {
                this->widget_data->begin_edition();
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* rdkit_bond = rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);

                auto bond_geometry = CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir());
                auto new_bond_geometry = CanvasMolecule::cycle_bond_geometry(bond_geometry);
                g_debug("Target bond geometry: %u",static_cast<unsigned int>(new_bond_geometry));
                rdkit_bond->setBondDir(CanvasMolecule::bond_geometry_to_rdkit(new_bond_geometry));

                this->widget_data->update_status("Geometry of bond has been altered.");
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
                g_debug("Final bond geometry: %u",static_cast<unsigned int>(CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir())));
                this->widget_data->finalize_edition();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter bond geometry: ") + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::alter_charge(int x, int y) {
    check_variant(Variant::ChargeModifier);
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                this->widget_data->begin_edition();
                
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                // Do we need this here?
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                auto* rdkit_atom = rdkit_mol->getAtomWithIdx(atom.idx);
                int old_charge = rdkit_atom->getFormalCharge();

                const RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
                auto valence_list = table->getValenceList(rdkit_atom->getAtomicNum());
                auto valence_list_to_string = [&](){
                    std::string valence_list_string = "[";
                    auto vl_it = valence_list.begin();
                    while(vl_it != valence_list.end()) {
                        valence_list_string += std::to_string(*vl_it);
                        vl_it++;
                        if(vl_it != valence_list.end()) {
                            valence_list_string += ", ";
                        }
                    }
                    valence_list_string += "]";
                    return valence_list_string;
                };

                auto valence_list_string = valence_list_to_string();
                g_debug(
                    "Valence list from RDKit: %s Num of outer shell electrons: %i", 
                    valence_list_string.c_str(),
                    table->getNouterElecs(rdkit_atom->getAtomicNum())
                );

                g_warning("todo: Fix-up computing plausible charges for atoms in the charge tool.");
                // We won't have to clear the list when it'll have the right contents (something to be fixed)
                valence_list.clear();
                // valence_list.push_back(0);
                for(int i = -4; i != 5; i++) {
                    valence_list.push_back(i);
                }

                auto it = std::find(valence_list.begin(),valence_list.end(),old_charge);
                if(it != valence_list.end()) {
                    it++;
                }
                if(it == valence_list.end()) {
                    it = valence_list.begin();
                }

                int new_charge = *it;

                valence_list_string = valence_list_to_string();
                g_info("Old formal charge: %i New formal charge: %i List: %s",old_charge,new_charge,valence_list_string.c_str());
                rdkit_atom->setFormalCharge(new_charge);

                this->widget_data->update_status("Charge of atom has been altered.");

                this->sanitize_molecule(*rdkit_mol.get());
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);

                this->widget_data->finalize_edition();
            } else { // a bond
                // auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                g_warning("The ChargeModifier tool does not operate on bonds. Nothing to do.");
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter charge: ") + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::delete_at(int x, int y) {
    check_variant(Variant::Delete);
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            this->widget_data->begin_edition();
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            RDKit::MolOps::Kekulize(*rdkit_mol.get());
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                rdkit_mol->removeAtom(atom.idx);
                canvas_mol.update_cached_atom_coordinate_map_after_atom_removal(atom.idx);
                this->widget_data->update_status("Atom has been deleted.");
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                rdkit_mol->removeBond(bond.first_atom_idx, bond.second_atom_idx);
                this->widget_data->update_status("Bond has been deleted.");
            }
            this->sanitize_molecule(*rdkit_mol.get());
            canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
            this->widget_data->finalize_edition();
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string status_msg = "Coud not delete atom: "; status_msg += e.what();
            this->widget_data->update_status(status_msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::insert_structure(int x, int y) {
    check_variant(Variant::StructureInsertion);
    auto click_result = this->widget_data->resolve_click(x, y);
    using Structure = StructureInsertion::Structure;
    auto structure_kind = this->structure_insertion.get_structure();

    RDKit::RWMol* rdkit_mol_ptr = nullptr;

    auto append_carbon = [&](unsigned int target_idx, RDKit::Bond::BondType bond_type = RDKit::Bond::SINGLE) -> unsigned int {
        RDKit::Atom* new_carbon = new RDKit::Atom(6);
        auto new_carbon_idx = rdkit_mol_ptr->addAtom(new_carbon,false,true);
        rdkit_mol_ptr->addBond(target_idx,new_carbon_idx,bond_type);
        return new_carbon_idx;
    };

    auto append_carbon_chain = [&](unsigned int chain_start_idx, std::size_t atom_count) -> unsigned int {
        unsigned int current_atom = chain_start_idx;
        for(std::size_t i = 0; i < atom_count; i++) {
            current_atom = append_carbon(current_atom);
        }
        return current_atom;
    };

    auto append_structure_to_atom = [&](unsigned int atom_idx){
        auto first_carbon = append_carbon(atom_idx);
        switch (structure_kind) {
            case Structure::CycloPropaneRing: {
                auto third_carbon = append_carbon_chain(first_carbon,2);
                rdkit_mol_ptr->addBond(first_carbon,third_carbon,RDKit::Bond::SINGLE);
                break;
            }
            case Structure::CycloButaneRing: {
                auto last_carbon = append_carbon_chain(first_carbon,3);
                rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
                break;
            }
            case Structure::CycloPentaneRing: {
                auto last_carbon = append_carbon_chain(first_carbon,4);
                rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
                break;
            }
            case Structure::CycloHexaneRing: {
                auto last_carbon = append_carbon_chain(first_carbon,5);
                rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
                break;
            }
            case Structure::BenzeneRing: {
                auto second_carbon = append_carbon(first_carbon);
                auto third_carbon = append_carbon(second_carbon,RDKit::Bond::DOUBLE);
                auto fourth_carbon = append_carbon(third_carbon);
                auto fifth_carbon = append_carbon(fourth_carbon,RDKit::Bond::DOUBLE);
                auto sixth_carbon = append_carbon(fifth_carbon);
                rdkit_mol_ptr->addBond(first_carbon,sixth_carbon,RDKit::Bond::DOUBLE);
                break;
            }
            case Structure::CycloHeptaneRing: {
                auto last_carbon = append_carbon_chain(first_carbon,6);
                rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
                break;
            }
            case Structure::CycloOctaneRing: {
                auto last_carbon = append_carbon_chain(first_carbon,7);
                rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
                break;
            }
        }
    };

    try {
        if(click_result.has_value()) {
            this->widget_data->begin_edition();
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
            rdkit_mol_ptr = rdkit_mol.get();

            RDKit::MolOps::Kekulize(*rdkit_mol.get());

            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                append_structure_to_atom(atom.idx);
                this->widget_data->update_status("Carbon ring has been appended to the atom.");
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                auto last_carbon = bond.second_atom_idx;
                auto first_carbon = bond.first_atom_idx;
                switch (structure_kind) {
                    case Structure::CycloPropaneRing: {
                        auto second_carbon = append_carbon(first_carbon);
                        rdkit_mol->addBond(second_carbon,last_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                    case Structure::CycloButaneRing: {
                        auto fourth_carbon = append_carbon_chain(first_carbon,2);
                        rdkit_mol->addBond(last_carbon,fourth_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                    case Structure::CycloPentaneRing: {
                        auto fifth_carbon = append_carbon_chain(first_carbon,3);
                        rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                    case Structure::CycloHexaneRing: {
                        auto sixth_carbon = append_carbon_chain(first_carbon,4);
                        rdkit_mol->addBond(last_carbon,sixth_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                    case Structure::BenzeneRing: {
                        auto second_carbon = append_carbon(first_carbon,RDKit::Bond::DOUBLE);
                        auto third_carbon = append_carbon(second_carbon);
                        auto fourth_carbon = append_carbon(third_carbon,RDKit::Bond::DOUBLE);
                        auto fifth_carbon = append_carbon(fourth_carbon);
                        rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::DOUBLE);
                        break;
                    }
                    case Structure::CycloHeptaneRing: {
                        auto seventh_carbon = append_carbon_chain(first_carbon,5);
                        rdkit_mol->addBond(last_carbon,seventh_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                    case Structure::CycloOctaneRing: {
                        auto eigth_carbon = append_carbon_chain(first_carbon,6);
                        rdkit_mol->addBond(last_carbon,eigth_carbon,RDKit::Bond::SINGLE);
                        break;
                    }
                }

                this->widget_data->update_status("Carbon ring has been added, adjacent to the bond.");
            }
            this->sanitize_molecule(*rdkit_mol.get());
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
            this->widget_data->finalize_edition();
            
        } else {
            // Nothing has been clicked on.
            g_debug("The click could not be resolved to any atom or bond.");
            if(this->widget_data->rdkit_molecules->empty()) {
                g_debug("There are no molecules. Structure insertion will therefore create a new one.");
                auto* widget_ptr = static_cast<impl::CootLigandEditorCanvasPriv*>(this->widget_data);
                auto rdkit_mol = std::make_shared<RDKit::RWMol>();
                rdkit_mol->addAtom(new RDKit::Atom(6),false,true);
                rdkit_mol_ptr = rdkit_mol.get();
                append_structure_to_atom(0);
                rdkit_mol->removeAtom((unsigned int) 0);
                // This function calls "begin_edition" and "finalize_edition", 
                // so we can't call "begin_edition" here above.
                coot_ligand_editor_append_molecule(COOT_COOT_LIGAND_EDITOR_CANVAS(widget_ptr), rdkit_mol);
                this->widget_data->update_status("New molecule created from carbon ring.");
            }
        }
    } catch(std::exception& e) {
        g_warning("An error occured: %s",e.what());
        std::string status_msg = "Coud not insert structure: "; status_msg += e.what();
        if (this->widget_data->is_in_edition()) {
            this->widget_data->rollback_current_edition();
        } else { // Needed if a new molecule was created on an empty canvas
            this->widget_data->undo_edition();
        }
        this->widget_data->update_status(status_msg.c_str());
    }
}

void ActiveTool::format_at(int x, int y) {
    check_variant(Variant::Format);
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            this->widget_data->begin_edition();
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            canvas_mol.clear_cached_atom_coordinate_map();
            canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
            this->widget_data->finalize_edition();
            this->widget_data->update_status("Molecule has been formatted.");
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string status_msg = "Coud not format molecule: "; status_msg += e.what();
            this->widget_data->update_status(status_msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::begin_transform(int x, int y, TransformManager::Mode mode) {
    auto mol_opt = this->widget_data->resolve_click(x, y);
    if(mol_opt.has_value()) {
        auto [atom_or_bond,mol_id] = mol_opt.value();
        transform_manager.begin_transform(x, y, mode);
        transform_manager.set_canvas_molecule_index(mol_id);
        this->widget_data->begin_edition();
    }
}

bool ActiveTool::is_in_transform() const noexcept {
    return this->transform_manager.is_active();
}

void ActiveTool::update_transform_cursor_pos(int x, int y, bool snap_to_angle) noexcept {
    if(!transform_manager.is_active()) {
        return;
    }
    transform_manager.update_current_cursor_pos(x, y, snap_to_angle);
    transform_manager.apply_current_transform_state(this->widget_data, snap_to_angle, false);
}

void ActiveTool::end_transform(bool snap_to_angle) {
    if(!transform_manager.is_active()) {
        return;
    }
    transform_manager.apply_current_transform_state(this->widget_data, snap_to_angle, true);
    transform_manager.end_transform();
    this->widget_data->finalize_edition();
}

void ActiveTool::sanitize_molecule(RDKit::RWMol& mol) const {
    if (!this->widget_data->allow_invalid_molecules) {
        RDKit::MolOps::sanitizeMol(mol);
    }
}

void ActiveTool::flip(int x, int y) {
    check_variant(Variant::FlipTool);
    auto& flip = this->flip_tool;
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            this->widget_data->begin_edition();
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            canvas_mol.perform_flip(flip.get_mode());
            canvas_mol.lower_from_rdkit(!this->widget_data->allow_invalid_molecules);
            this->widget_data->finalize_edition();
            this->widget_data->update_status("Molecule has been flipped.");
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string status_msg = "Coud not flip molecule: "; status_msg += e.what();
            this->widget_data->update_status(status_msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

/// Valid for Variant::BondModifier.
std::optional<std::pair<unsigned int,unsigned int>> ActiveTool::get_molecule_idx_and_first_atom_of_new_bond() const {
    check_variant(Variant::BondModifier);
    auto& mod = this->bond_modifier;
    return mod.get_molecule_idx_and_first_atom_of_new_bond();
}

ElementInsertion::ElementInsertion(ElementInsertion::Element el) noexcept {
    this->element = el;
}

ElementInsertion::ElementInsertion(const char* symbol) {
    RDKit::PeriodicTable* t = RDKit::PeriodicTable::getTable();
    unsigned int atomic_number = t->getAtomicNumber(symbol);
    this->element = atomic_number;
}

unsigned int ElementInsertion::get_atomic_number() const noexcept {
    if(std::holds_alternative<Element>(this->element)) {
        switch (std::get<Element>(this->element)) {
            default:
            case ElementInsertion::Element::C:{
                return 6;
            }
            case ElementInsertion::Element::N:{
                return 7;
            }
            case ElementInsertion::Element::O:{
                return 8;
            }
            case ElementInsertion::Element::S:{
                return 16;
            }
            case ElementInsertion::Element::P:{
                return 15;
            }
            case ElementInsertion::Element::H:{
                return 1;
            }
            case ElementInsertion::Element::F:{
                return 9;
            }
            case ElementInsertion::Element::Cl:{
                return 17;
            }
            case ElementInsertion::Element::Br:{
                return 35;
            }
            case ElementInsertion::Element::I:{
                return 53;
            }
        }
    } else { // raw atomic number
        return std::get<unsigned int>(this->element);
    }
}

FlipMode FlipTool::get_mode() const noexcept {
    return this->mode;
}

FlipTool::FlipTool(FlipMode mode) noexcept {
    this->mode = mode;
}