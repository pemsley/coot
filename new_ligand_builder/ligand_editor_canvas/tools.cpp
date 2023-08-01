#include "tools.hpp"
#include "core.hpp"
#include "model.hpp"
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include "../ligand_editor_canvas.hpp"

using namespace coot::ligand_editor_canvas;


void Tool::on_click(impl::WidgetCoreData& widget_data, int x, int y) {
    // nothing by default
}

void Tool::on_blank_space_click(impl::WidgetCoreData& widget_data, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
}

void Tool::on_release(impl::WidgetCoreData& widget_data, int x, int y) {
    // nothing by default
}

void Tool::after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&) {
    // nothing by default
}

bool Tool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&) {
    // nothing by default
    return true;
}

void Tool::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&, CanvasMolecule::Bond&) {
    // nothing by default
    g_debug("The tool does not operate on bonds.");
}

void Tool::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&, CanvasMolecule::Atom&) {
    // nothing by default
    g_debug("The tool does not operate on atoms.");
}


std::string Tool::get_exception_message_prefix() const noexcept {
    return "An error occured: ";
}

Tool::~Tool() {

}

void ActiveTool::on_click(int x, int y) {
    if(!this->tool) {
        return;
    }

    this->tool->on_click(*this->widget_data, x, y);
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            if(!this->tool->on_molecule_click(*this->widget_data, molecule_idx, rdkit_mol, canvas_mol)) {
                return;
            }
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                this->tool->on_atom_click(*this->widget_data, molecule_idx, rdkit_mol, canvas_mol, atom);
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                this->tool->on_bond_click(*this->widget_data, molecule_idx, rdkit_mol, canvas_mol, bond);
            }
            this->tool->after_molecule_click(*this->widget_data, molecule_idx, rdkit_mol, canvas_mol);
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = this->tool->get_exception_message_prefix() + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        this->tool->on_blank_space_click(*this->widget_data, x, y);
    }
}

void ActiveTool::on_release(int x, int y) {
    if(!this->tool) {
        return;
    }
    this->tool->on_release(*this->widget_data, x, y);
}

void TransformTool::on_click(impl::WidgetCoreData& widget_data, int x, int y) {
    auto mol_opt = widget_data.resolve_click(x, y);
    if(mol_opt.has_value()) {
        auto [atom_or_bond,mol_id] = mol_opt.value();
        this->transform_manager->begin_transform(x, y, this->mode);
        this->transform_manager->set_canvas_molecule_index(mol_id);
        widget_data.begin_edition();
    }
}

bool TransformTool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&) {
    return false;
}

TransformTool::TransformTool(TransformManager::Mode mode) noexcept {
    this->mode = mode;
    this->transform_manager = nullptr;
}

void TransformTool::set_transform_manager(TransformManager* mgr) noexcept {
    this->transform_manager = mgr;
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

ActiveTool::ActiveTool() noexcept {
    // nothing
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->tool = std::make_unique<ElementInsertion>(std::move(insertion));
}

ActiveTool::ActiveTool(BondModifier modifier) noexcept {
    this->tool = std::make_unique<BondModifier>(std::move(modifier));
}

ActiveTool::ActiveTool(DeleteTool deltool) noexcept {
    this->tool = std::make_unique<DeleteTool>(std::move(deltool));
}

ActiveTool::ActiveTool(ChargeModifier chargemod) noexcept {
    this->tool = std::make_unique<ChargeModifier>(std::move(chargemod));
}

ActiveTool::ActiveTool(TransformTool mov) noexcept {
    mov.set_transform_manager(&this->transform_manager);
    this->tool = std::make_unique<TransformTool>(std::move(mov));
}

ActiveTool::ActiveTool(StructureInsertion insertion) noexcept {
    this->tool = std::make_unique<StructureInsertion>(std::move(insertion));
}

ActiveTool::ActiveTool(GeometryModifier modifier) noexcept {
    this->tool = std::make_unique<GeometryModifier>(std::move(modifier));
}

ActiveTool::ActiveTool(FormatTool fmt) noexcept {
    this->tool = std::make_unique<FormatTool>(std::move(fmt));
}

ActiveTool::ActiveTool(FlipTool flip) noexcept {
    this->tool = std::make_unique<FlipTool>(std::move(flip));
}

ActiveTool::ActiveTool(RemoveHydrogensTool rh) noexcept {
    this->tool = std::make_unique<RemoveHydrogensTool>(std::move(rh));
}

void ActiveTool::set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept {
    this->widget_data = static_cast<impl::WidgetCoreData*>(owning_widget);
}

bool ElementInsertion::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    // maybe more useful to override in the future
    return true;
}

void ElementInsertion::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) {
    g_warning("TODO: Implement handling insertion at bonds (if any should happen)");
}

void ElementInsertion::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) {
    unsigned int atomic_number = this->get_atomic_number();
    auto el_name = RDKit::PeriodicTable::getTable()->getElementSymbol(atomic_number);
    g_debug("Appending element '%u' (%s) to destination atom: idx=%i, symbol=%s.",atomic_number,el_name.c_str(),atom.idx,atom.symbol.c_str());
    widget_data.begin_edition();
    auto* new_atom = new RDKit::Atom(std::string(el_name));
    rdkit_mol->replaceAtom(atom.idx, new_atom);
    widget_data.update_status("Atom has been replaced.");
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
}

void ElementInsertion::after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    // maybe more useful to override in the future
}

std::string ElementInsertion::get_exception_message_prefix() const noexcept {
    return "Could not insert atom: ";
}

bool BondModifier::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    widget_data.begin_edition();
    return true;
}

void BondModifier::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) {
    auto* rdkit_bond = rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);
    RDKit::MolOps::Kekulize(*rdkit_mol.get());
    auto old_type = rdkit_bond->getBondType();
    rdkit_bond->setBondType(CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
    try {
        this->sanitize_molecule(widget_data, *rdkit_mol.get());
    }catch(std::exception& e) {
        // rollback
        g_warning("Rolling back invalid molecule change that makes it unable to sanitize with the following error: %s",e.what());
        rdkit_bond->setBondType(old_type);
        this->sanitize_molecule(widget_data, *rdkit_mol.get());
        // rethrow
        throw std::runtime_error(std::string("Invalid bond modification: ") + e.what());
    }
    widget_data.update_status("Bond has been altered.");
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
}

void BondModifier::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) {
    this->begin_creating_bond(mol_idx, atom.idx);
}

std::string BondModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter/create bond: ";
}

bool ActiveTool::is_creating_bond() const noexcept {
    const BondModifier* mod = dynamic_cast<const BondModifier*>(this->tool.get());
    if (mod) {
        return mod->is_creating_bond();
    }
    return false;
}

void BondModifier::on_release(impl::WidgetCoreData& widget_data, int x, int y) {
    if(!this->is_creating_bond()) {
        return;
    }
    auto click_result = widget_data.resolve_click(x, y);
    auto [original_molecule_idx, first_atom_idx] = this->get_molecule_idx_and_first_atom_of_new_bond().value();
    this->finish_creating_bond();
    widget_data.currently_created_bond = std::nullopt;

    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto second_atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                if(original_molecule_idx != molecule_idx) {
                    widget_data.update_status("Cannot create bond between different molecules!");
                    widget_data.rollback_current_edition();
                    return;
                }
                auto& rdkit_mol = widget_data.rdkit_molecules->at(molecule_idx);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                if(first_atom_idx == second_atom.idx) {
                    auto* new_atom = new RDKit::Atom(6);
                    auto new_atom_idx = rdkit_mol->addAtom(new_atom,false,true);
                    rdkit_mol->addBond(new_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
                    g_info("New atom added: idx=%i",new_atom_idx);
                    widget_data.update_status("New carbon atom added.");
                } else {
                    rdkit_mol->addBond(first_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
                    widget_data.update_status("Created new bond between atoms.");
                }
                this->sanitize_molecule(widget_data, *rdkit_mol.get());
                auto& canvas_mol = widget_data.molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
                widget_data.finalize_edition();
            } else {
                widget_data.update_status("Can't link bond to a bond!");
                widget_data.rollback_current_edition();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter/create bond: ") + e.what();
            widget_data.update_status(msg.c_str());
            widget_data.rollback_current_edition();
        }
    } else {
        std::string msg = "The new bond goes nowhere.";
        widget_data.update_status(msg.c_str());
        widget_data.rollback_current_edition();
    }
}

void GeometryModifier::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) {
    widget_data.begin_edition();
    auto* rdkit_bond = rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);

    auto bond_geometry = CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir());
    auto new_bond_geometry = CanvasMolecule::cycle_bond_geometry(bond_geometry);
    g_debug("Target bond geometry: %u",static_cast<unsigned int>(new_bond_geometry));
    rdkit_bond->setBondDir(CanvasMolecule::bond_geometry_to_rdkit(new_bond_geometry));

    widget_data.update_status("Geometry of bond has been altered.");
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    g_debug("Final bond geometry: %u",static_cast<unsigned int>(CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir())));
    widget_data.finalize_edition();
}

std::string GeometryModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter bond geometry: ";
}

void ChargeModifier::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) {
    widget_data.begin_edition();
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

    widget_data.update_status("Charge of atom has been altered.");

    Tool::sanitize_molecule(widget_data, *rdkit_mol.get());
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);

    widget_data.finalize_edition();
}

std::string ChargeModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter charge: ";
}

bool DeleteTool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    widget_data.begin_edition();
    RDKit::MolOps::Kekulize(*rdkit_mol.get());
    return true;
}

void DeleteTool::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) {
    rdkit_mol->removeBond(bond.first_atom_idx, bond.second_atom_idx);
    widget_data.update_status("Bond has been deleted.");
}

void DeleteTool::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) {
    rdkit_mol->removeAtom(atom.idx);
    canvas_mol.update_cached_atom_coordinate_map_after_atom_removal(atom.idx);
    widget_data.update_status("Atom has been deleted.");
}

void DeleteTool::after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    Tool::sanitize_molecule(widget_data, *rdkit_mol.get());
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
}

std::string DeleteTool::get_exception_message_prefix() const noexcept {
    return "Could not delete atom/bond: ";
}


unsigned int StructureInsertion::append_carbon(RDKit::RWMol* rdkit_mol_ptr, unsigned int target_idx, RDKit::Bond::BondType bond_type) {
    RDKit::Atom* new_carbon = new RDKit::Atom(6);
    auto new_carbon_idx = rdkit_mol_ptr->addAtom(new_carbon,false,true);
    rdkit_mol_ptr->addBond(target_idx,new_carbon_idx,bond_type);
    return new_carbon_idx;  
}

unsigned int StructureInsertion::append_carbon_chain(RDKit::RWMol* rdkit_mol_ptr, unsigned int chain_start_idx, std::size_t atom_count) {
    unsigned int current_atom = chain_start_idx;
    for(std::size_t i = 0; i < atom_count; i++) {
        current_atom = append_carbon(rdkit_mol_ptr, current_atom);
    }
    return current_atom;
}

void StructureInsertion::append_structure_to_atom(RDKit::RWMol* rdkit_mol_ptr, unsigned int atom_idx) const {
    auto first_carbon = append_carbon(rdkit_mol_ptr, atom_idx);
    switch (this->get_structure()) {
        case Structure::CycloPropaneRing: {
            auto third_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,2);
            rdkit_mol_ptr->addBond(first_carbon,third_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloButaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,3);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloPentaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,4);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloHexaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,5);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::BenzeneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon);
            auto third_carbon = append_carbon(rdkit_mol_ptr, second_carbon, RDKit::Bond::DOUBLE);
            auto fourth_carbon = append_carbon(rdkit_mol_ptr, third_carbon);
            auto fifth_carbon = append_carbon(rdkit_mol_ptr, fourth_carbon, RDKit::Bond::DOUBLE);
            auto sixth_carbon = append_carbon(rdkit_mol_ptr, fifth_carbon);
            rdkit_mol_ptr->addBond(first_carbon,sixth_carbon,RDKit::Bond::DOUBLE);
            break;
        }
        case Structure::CycloHeptaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,6);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloOctaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,7);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
    }
}


bool StructureInsertion::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    widget_data.begin_edition();
    RDKit::MolOps::Kekulize(*rdkit_mol.get());
    return true;
}

void StructureInsertion::on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) {
    auto last_carbon = bond.second_atom_idx;
    auto first_carbon = bond.first_atom_idx;
    auto structure_kind = this->get_structure();
    auto* rdkit_mol_ptr = rdkit_mol.get();
    switch (structure_kind) {
        case Structure::CycloPropaneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon);
            rdkit_mol->addBond(second_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloButaneRing: {
            auto fourth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,2);
            rdkit_mol->addBond(last_carbon,fourth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloPentaneRing: {
            auto fifth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,3);
            rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloHexaneRing: {
            auto sixth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,4);
            rdkit_mol->addBond(last_carbon,sixth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::BenzeneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon,RDKit::Bond::DOUBLE);
            auto third_carbon = append_carbon(rdkit_mol_ptr, second_carbon);
            auto fourth_carbon = append_carbon(rdkit_mol_ptr, third_carbon,RDKit::Bond::DOUBLE);
            auto fifth_carbon = append_carbon(rdkit_mol_ptr, fourth_carbon);
            rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::DOUBLE);
            break;
        }
        case Structure::CycloHeptaneRing: {
            auto seventh_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,5);
            rdkit_mol->addBond(last_carbon,seventh_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloOctaneRing: {
            auto eigth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,6);
            rdkit_mol->addBond(last_carbon,eigth_carbon,RDKit::Bond::SINGLE);
            break;
        }
    }
    widget_data.update_status("Carbon ring has been added, adjacent to the bond.");
}

void StructureInsertion::on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) {
    this->append_structure_to_atom(rdkit_mol.get(),atom.idx);
    widget_data.update_status("Carbon ring has been appended to the atom.");
}

void StructureInsertion::after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    this->sanitize_molecule(widget_data, *rdkit_mol);
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
}

std::string StructureInsertion::get_exception_message_prefix() const noexcept {
    return "Could not insert structure: ";
}

void StructureInsertion::on_blank_space_click(impl::WidgetCoreData& widget_data, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
    if(widget_data.rdkit_molecules->empty()) {
        g_debug("There are no molecules. Structure insertion will therefore create a new one.");
        auto* widget_ptr = static_cast<impl::CootLigandEditorCanvasPriv*>(&widget_data);
        auto rdkit_mol = std::make_shared<RDKit::RWMol>();
        rdkit_mol->addAtom(new RDKit::Atom(6),false,true);
        append_structure_to_atom(rdkit_mol.get(),0);
        rdkit_mol->removeAtom((unsigned int) 0);
        // This function calls "begin_edition" and "finalize_edition", 
        // so we can't call "begin_edition" here above.
        coot_ligand_editor_append_molecule(COOT_COOT_LIGAND_EDITOR_CANVAS(widget_ptr), rdkit_mol);
        widget_data.update_status("New molecule created from carbon ring.");
        // todo: make sure that this is crash-safe vs edit/undo
    }
}

bool RemoveHydrogensTool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    g_warning("todo: Implement RemoveHydrogensTool");
    return false;
}

std::string RemoveHydrogensTool::get_exception_message_prefix() const noexcept {
    return "Could not remove hydrogens from molecule: ";
}


bool FlipTool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    widget_data.begin_edition();
    canvas_mol.perform_flip(this->mode);
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
    widget_data.update_status("Molecule has been flipped.");
    return false;
}

std::string FlipTool::get_exception_message_prefix() const noexcept {
    return "Could not flip molecule: ";
}


bool FormatTool::on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) {
    widget_data.begin_edition();
    canvas_mol.clear_cached_atom_coordinate_map();
    canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
    widget_data.finalize_edition();
    widget_data.update_status("Molecule has been formatted.");
    return false;
}

std::string FormatTool::get_exception_message_prefix() const noexcept {
    return "Could not format molecule: ";
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

void Tool::sanitize_molecule(impl::WidgetCoreData& widget_data, RDKit::RWMol& mol) {
    if (!widget_data.allow_invalid_molecules) {
        RDKit::MolOps::sanitizeMol(mol);
    }
}

std::optional<std::pair<unsigned int,unsigned int>> ActiveTool::get_molecule_idx_and_first_atom_of_new_bond() const noexcept {
    const BondModifier* mod = dynamic_cast<const BondModifier*>(this->tool.get());
    if(mod) {
        return mod->get_molecule_idx_and_first_atom_of_new_bond();
    } else {
        return std::nullopt;
    }
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

FlipTool::FlipTool(FlipMode mode) noexcept {
    this->mode = mode;
}