// mmdb-shim — layer 4 of 4: out-of-line definitions.
//
// Everything here needs one or more complete classes from the layers above: the
// g() resolvers (detached wrapper -> _local, else parent's gemmi vector), the
// UDData put/get helpers, the Atom/Residue/Chain/Model/Manager methods that touch
// siblings or the Manager, the detached-construction constructors, the CID/atom
// selection matchers (namespace detail), and the free helper functions.
#pragma once

#include "_shim_manager.hh"

namespace mmdb {

   // ---- g() resolvers ----
   // A wrapper with no parent is "detached" (Coot's `new mmdb::Atom` idiom: build
   // standalone, set fields, then Add*() into a parent). While detached, g()
   // resolves to a wrapper-owned local gemmi object; Add*() copies that local into
   // the parent's gemmi vector and rebinds (sets parent + index). Index-based
   // resolution makes the vector push/reallocation harmless for siblings.
   inline gemmi::Model &Model::g() const { return mgr ? mgr->st.models[mi] : const_cast<Model *>(this)->_local; }
   inline gemmi::Chain &Chain::g() const { return model ? model->g().chains[ci] : const_cast<Chain *>(this)->_local; }
   inline gemmi::Residue &Residue::g() const { return chain ? chain->g().residues[ri] : const_cast<Residue *>(this)->_local; }
   inline gemmi::Atom &Atom::g() const {
      if (res && res->chain && res->chain->model && res->chain->model->mgr) {
         auto &mgr = *res->chain->model->mgr;
         int mi = res->chain->model->mi;
         if (mi >= 0 && mi < (int)mgr.st.models.size()) {
            return res->g().atoms[ai];
         }
      }
      return const_cast<Atom *>(this)->_local;
   }

   // ---- UDData helpers ----
   inline Manager::UDReg *_ud_desc(Manager *mgr, UDR_TYPE myType, int handle, int kind,
                                   int &err) {
      if (!mgr || handle < 1 || handle > (int)mgr->ud_regs.size()) {
         err = UDDATA_WrongHandle;
         return nullptr;
      }
      Manager::UDReg &d = mgr->ud_regs[handle - 1];  // handles are 1-based (see _regUD)
      if (d.type != myType || d.kind != kind) {
         err = UDDATA_WrongUDRType;
         return nullptr;
      }
      err = UDDATA_Ok;
      return &d;
   }
   inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, int v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 0, e);
      if (!d) return e;
      if ((int)s._udi.size() <= d->slot) s._udi.resize(d->slot + 1, 0);
      s._udi[d->slot] = v;
      return UDDATA_Ok;
   }
   inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, realtype v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 1, e);
      if (!d) return e;
      if ((int)s._udr.size() <= d->slot) s._udr.resize(d->slot + 1, 0.0);
      s._udr[d->slot] = v;
      return UDDATA_Ok;
   }
   inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, cpstr v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 2, e);
      if (!d) return e;
      if ((int)s._uds.size() <= d->slot) s._uds.resize(d->slot + 1);
      s._uds[d->slot] = v ? v : "";
      return UDDATA_Ok;
   }
   inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, int &v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 0, e);
      if (!d) return e;
      if ((int)s._udi.size() <= d->slot) return UDDATA_NoData;
      v = s._udi[d->slot];
      return UDDATA_Ok;
   }
   inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, realtype &v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 1, e);
      if (!d) return e;
      if ((int)s._udr.size() <= d->slot) return UDDATA_NoData;
      v = s._udr[d->slot];
      return UDDATA_Ok;
   }
   inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, pstr &v) {
      int e;
      auto *d = _ud_desc(mgr, t, h, 2, e);
      if (!d) return e;
      if ((int)s._uds.size() <= d->slot) return UDDATA_NoData;
      v = (pstr)s._uds[d->slot].c_str();
      return UDDATA_Ok;  // borrowed
   }

   // ---- Atom out-of-line ----
   inline pstr Atom::GetAtomName() const {
      // MMDB returns the PDB-column-aligned 4-char atom name (e.g. " N  ", " CA ",
      // " CG2"); gemmi stores the trimmed name ("N"/"CA"/"CG2"). Reproduce MMDB
      // alignment via gemmi's padded_name() (left-pad by element) + right-pad to 4
      // — the exact rule gemmi's own mmdb.hpp bridge uses. Coot's atom_spec_t names
      // are these 4-char strings, so returning the trimmed name breaks every lookup.
      std::string padded = g().padded_name();
      if (padded.size() < 4) padded.resize(4, ' ');
      std::snprintf(_name_buf, sizeof(_name_buf), "%s", padded.c_str());
      return _name_buf;
   }
   // Coot passes MMDB-aligned 4-char names (" CA "); gemmi stores trimmed names
   // ("CA") and re-pads on output (GetAtomName / PDB write) — store trimmed so
   // gemmi's own formatting stays correct.
   inline void Atom::SetAtomName(const AtomName aName) { g().name = aName ? gemmi::trim_str(aName) : ""; }
   inline pstr Atom::GetElementName() {
      // MMDB returns the PDB-column-aligned element: 2 chars, right-justified,
      // UPPERCASE (" H", " C", "NA", "FE"). gemmi's Element::name() is unpadded and
      // mixed-case ("H"/"C"/"Na"/"Fe"), so Coot's element tests — e.g.
      // get_number_of_hydrogen_atoms() comparing `ele == " H"` — never match. Align
      // to MMDB: uppercase then right-pad into a 2-wide field.
      std::string e = g().element.name();
      for (char &c : e) c = std::toupper((unsigned char)c);
      if (e.size() < 2) e.insert(e.begin(), 2 - e.size(), ' ');
      std::snprintf(_elem_buf, sizeof(_elem_buf), "%s", e.c_str());
      return _elem_buf;
   }
   inline void Atom::SetElementName(const Element elName) { g().element = gemmi::Element(elName); }
   inline pstr Atom::GetChainID() { return res ? res->GetChainID() : mmdb_empty_pstr(); }
   inline int Atom::GetSeqNum() { return res ? res->GetSeqNum() : 0; }
   inline Chain *Atom::GetChain() { return res ? res->GetChain() : nullptr; }
   inline Model *Atom::GetModel() { return (res && res->chain) ? res->chain->model : nullptr; }
   inline pstr Atom::GetLabelCompID() { return res ? res->GetLabelCompID() : nullptr; }
   inline pstr Atom::GetLabelAsymID() { return res ? res->GetLabelAsymID() : nullptr; }
   inline int Atom::GetLabelSeqID() { return res ? res->GetLabelSeqID() : 0; }
   inline int Atom::GetLabelEntityID() { return res ? res->GetLabelEntityID() : 0; }
   inline int Atom::GetResidueNo() { return res ? res->GetResidueNo() : 0; }
   inline int Atom::GetSSEType() { return res ? res->SSE : SSE_None; }
   inline bool Atom::isSolvent() { return res ? res->isSolvent() : false; }
   inline bool Atom::isNTerminus() { return res ? res->isNTerminus() : false; }
   inline bool Atom::isCTerminus() { return res ? res->isCTerminus() : false; }
   inline pstr Atom::GetInsCode() { return res ? res->GetInsCode() : mmdb_empty_pstr(); }
   inline pstr Atom::GetResName() { return res ? res->GetResName() : mmdb_empty_pstr(); }
   inline int Atom::GetModelNum() { return res ? res->GetModelNum() : 0; }
   inline int Atom::GetIndex() { return ai; }
   inline void Atom::SetCoordinates(realtype xx, realtype yy, realtype zz,
                                    realtype occ, realtype tF) {
      auto &a = g();
      a.pos = gemmi::Position(xx, yy, zz);
      a.occ = (float)occ;
      a.b_iso = (float)tF;
   }

   // ---- Residue out-of-line ----
   inline pstr Residue::GetResName() {
      std::snprintf(_resname_buf, sizeof(_resname_buf), "%s", g().name.c_str());
      return _resname_buf;
   }
   inline int &Residue::GetSeqNum() { return g().seqid.num.value; }
   inline pstr Residue::GetInsCode() {
      _inscode_buf[0] = g().seqid.icode == ' ' ? '\0' : g().seqid.icode;
      _inscode_buf[1] = '\0';
      return _inscode_buf;
   }
   inline pstr Residue::GetChainID() { return chain ? chain->GetChainID() : mmdb_empty_pstr(); }
   inline int Residue::GetModelNum() { return (chain && chain->model) ? chain->model->GetSerNum() : 0; }
   inline PAtom Residue::GetAtom(const AtomName aname, const Element elname, const AltLoc aloc) {
      // MMDB matches the PDB-column-aligned 4-char name (real Coot calls
      // `GetAtom(" CA ")`), but gemmi stores names trimmed ("CA"). Trim the query so
      // both padded and unpadded lookups resolve — mirrors the selection matchers
      // (detail::inList / SelectAtoms), which already trim both sides. Element/altLoc
      // still disambiguate when supplied (e.g. carbon-alpha " CA " vs calcium "CA  ",
      // which share a trimmed name).
      const std::string want = aname ? gemmi::trim_str(aname) : std::string();
      for (Atom *a : atoms) {
         if (std::string(a->g().name) != want) continue;
         if (elname && *elname && a->g().element.name() != std::string(elname)) continue;
         if (aloc && *aloc && a->g().altloc != aloc[0]) continue;
         return a;
      }
      return nullptr;
   }
   inline PAtom Residue::AddAtom(Manager &m, gemmi::Atom a) {
      g().atoms.push_back(std::move(a));
      Atom *aw = m.newAtom();
      aw->mgr = &m;
      aw->res = this;
      aw->ai = (int)atoms.size();
      atoms.push_back(aw);
      return aw;
   }
   // Adopt a Coot-`new`d detached atom (the `new mmdb::Atom; …; res->AddAtom(at)`
   // idiom). Copy its local gemmi into this residue, rebind the wrapper, and — when
   // this residue is manager-bound — transfer ownership to the manager so `~Manager`
   // frees it (MMDB semantics); `~Atom` drops it back out on an early Coot `delete`.
   inline int Residue::AddAtom(PAtom atm) {
      g().atoms.push_back(atm->_local);
      atm->res = this;
      atm->mgr = mgr;
      atm->ai = (int)atoms.size();
      if (atm->_altloc_buf[0]) g().atoms[atm->ai].altloc = atm->_altloc_buf[0];
      atoms.push_back(atm);
      if (mgr) mgr->_atom_allocs.insert(atm);
      _sync_atom();
      return 0;
   }
   inline void Residue::DeleteAtom(int pos) {
      if (pos < 0 || pos >= (int)atoms.size()) return;
      g().atoms.erase(g().atoms.begin() + pos);
      atoms[pos]->alive = false;
      atoms[pos]->ai = -1;
      atoms.erase(atoms.begin() + pos);
      for (int k = pos; k < (int)atoms.size(); ++k) atoms[k]->ai = k;
   }

   // Unlink one atom wrapper from this residue without freeing it: erase from the
   // wrapper table AND the parallel gemmi atom vector (kept in lockstep), then
   // reindex trailing atoms. Used by ~Atom when Coot `delete`s a live atom.
   inline void Residue::_detach_atom(Atom *a) {
      auto it = std::find(atoms.begin(), atoms.end(), a);
      if (it == atoms.end()) return;
      int pos = (int)(it - atoms.begin());
      if (pos < (int)g().atoms.size()) g().atoms.erase(g().atoms.begin() + pos);
      atoms.erase(atoms.begin() + pos);
      for (int k = pos; k < (int)atoms.size(); ++k) atoms[k]->ai = k;
   }

   // Coot's `delete atom;` removes an atom from the hierarchy. Detach it from the
   // parent residue, the manager/model flat lists, and any selections so no stale
   // pointer survives; then drop it from the ownership set. A no-op during teardown
   // (memory is freed wholesale by ~Manager) and for never-adopted detached atoms.
   inline Atom::~Atom() {
      if (!mgr || mgr->_bulk_free) return;
      Manager *m = mgr;
      if (res) {
         Model *mod = (res->chain ? res->chain->model : nullptr);
         res->_detach_atom(this);
         if (mod) {
            auto &ma = mod->all_atoms;
            ma.erase(std::remove(ma.begin(), ma.end(), this), ma.end());
         }
      }
      auto &aa = m->all_atoms;
      aa.erase(std::remove(aa.begin(), aa.end(), this), aa.end());
      for (int h = 1; h <= (int)m->selections.size(); ++h) {
         if (isInSelection(h)) {
            auto &sa = m->selections[h - 1].atoms;
            sa.erase(std::remove(sa.begin(), sa.end(), this), sa.end());
         }
      }
      m->_atom_allocs.erase(this);
   }

   // Coot's `delete residue_p;` (and Chain::DeleteResidue) removes a residue. Free its
   // atoms (MMDB: deleting a residue deletes its atoms), then DEFER the structural
   // removal: null this residue's slot in the parent chain but keep the gemmi
   // placeholder so siblings' `ri` stays valid; _compact_residues drops both. A no-op
   // during Manager teardown (memory is freed wholesale by ~Manager).
   inline Residue::~Residue() {
      if (!mgr || mgr->_bulk_free) return;
      Manager *m = mgr;
      // free my atoms: null each atom's res first so ~Atom doesn't mutate this->atoms
      // mid-iteration (~Atom still cleans all_atoms/model/selections/registry).
      std::vector<Atom *> ats = atoms;
      atoms.clear();
      for (Atom *a : ats) {
         a->res = nullptr;
         delete a;
      }
      // tombstone my slot in the parent chain (keep the gemmi placeholder residue)
      if (chain) {
         auto &rv = chain->residues;
         for (size_t i = 0; i < rv.size(); ++i)
            if (rv[i] == this) {
               rv[i] = nullptr;
               break;
            }
      }
      m->_res_allocs.erase(this);
   }

   inline void Chain::_compact_residues() {
      // Drop tombstoned (null) wrapper slots and their gemmi placeholder residues in
      // lock-step (wrapper[i] <-> g().residues[i]); then reindex ri to the new order.
      bool any_null = false;
      for (Residue *r : residues)
         if (!r) {
            any_null = true;
            break;
         }
      if (!any_null) return;
      std::vector<Residue *> keptw;
      std::vector<gemmi::Residue> keptg;
      gemmi::Chain &gc = g();
      keptw.reserve(residues.size());
      keptg.reserve(gc.residues.size());
      for (size_t i = 0; i < residues.size(); ++i) {
         if (residues[i]) {
            keptw.push_back(residues[i]);
            if (i < gc.residues.size()) keptg.push_back(std::move(gc.residues[i]));
         }
      }
      residues.swap(keptw);
      gc.residues.swap(keptg);
      for (int k = 0; k < (int)residues.size(); ++k) residues[k]->ri = k;
   }

   // ---- Chain out-of-line ----
   inline bool Chain::isAminoacidChain() {
      for (Residue *r : residues)
         if (r && r->isAminoacid()) return true;
      return false;
   }
   inline bool Chain::isNucleotideChain() {
      for (Residue *r : residues)
         if (r && r->isNucleotide()) return true;
      return false;
   }
   inline bool Chain::isSolventChain() {
      if (residues.empty()) return false;
      for (Residue *r : residues)
         if (r && !r->isSolvent()) return false;
      return true;
   }
   inline pstr Chain::GetChainID() {
      std::snprintf(_chainid_buf, sizeof(_chainid_buf), "%s", g().name.c_str());
      return _chainid_buf;
   }
   inline PResidue Chain::AddResidue(Manager &m, gemmi::Residue r) {
      g().residues.push_back(std::move(r));
      Residue *rw = m.newRes();
      rw->mgr = &m;
      rw->chain = this;
      rw->ri = (int)residues.size();
      for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
         Atom *aw = m.newAtom();
         aw->mgr = &m;
         aw->res = rw;
         aw->ai = ai;
         rw->atoms.push_back(aw);
      }
      residues.push_back(rw);
      return rw;
   }
   inline PResidue Chain::InsResidue(Manager &m, int pos, gemmi::Residue r) {
      _compact_residues();  // don't insert/reindex across deferred-delete tombstones
      g().residues.insert(g().residues.begin() + pos, std::move(r));
      Residue *rw = m.newRes();
      rw->mgr = &m;
      rw->chain = this;
      rw->ri = pos;
      residues.insert(residues.begin() + pos, rw);
      for (int k = pos + 1; k < (int)residues.size(); ++k) residues[k]->ri = k;
      for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
         Atom *aw = m.newAtom();
         aw->mgr = &m;
         aw->res = rw;
         aw->ai = ai;
         rw->atoms.push_back(aw);
      }
      return rw;
   }

   // ---- Model out-of-line ----
   inline PChain Model::GetChain(const ChainID chID) {
      for (Chain *c : chains)
         if (c->g().name == chID) return c;
      return nullptr;
   }

   // ---- Manager out-of-line ----
   inline void Manager::build_from_gemmi() {
      models.clear();
      all_atoms.clear();
      for (int mi = 0; mi < (int)st.models.size(); ++mi) {
         Model *mw = newModel();
         mw->mgr = this;
         mw->mi = mi;
         auto &gm = st.models[mi];
         for (int ci = 0; ci < (int)gm.chains.size(); ++ci) {
            Chain *cw = newChain();
            cw->mgr = this;
            cw->model = mw;
            cw->ci = ci;
            auto &gc = gm.chains[ci];
            for (int ri = 0; ri < (int)gc.residues.size(); ++ri) {
               Residue *rw = newRes();
               rw->mgr = this;
               rw->chain = cw;
               rw->ri = ri;
               auto &gr = gc.residues[ri];
               for (int ai = 0; ai < (int)gr.atoms.size(); ++ai) {
                  Atom *aw = newAtom();
                  aw->mgr = this;
                  aw->res = rw;
                  aw->ai = ai;
                  aw->WhatIsSet = ASET_Coordinates | ASET_Occupancy | ASET_tempFactor;
                  const gemmi::SMat33<float> &an = gr.atoms[ai].aniso;
                  if (an.u11 != 0.f || an.u22 != 0.f || an.u33 != 0.f) aw->WhatIsSet |= ASET_Anis_tFac;
                  rw->atoms.push_back(aw);
                  all_atoms.push_back(aw);
                  mw->all_atoms.push_back(aw);
               }
               rw->_sync_atom();
               rw->_load_id();
               cw->residues.push_back(rw);
            }
            mw->chains.push_back(cw);
         }
         models.push_back(mw);
      }
      _load_metadata();
   }

   // Map gemmi's structure-level metadata (connections / cispeps / helices /
   // sheets) onto the MMDB per-Model record containers. gemmi is the reader; the
   // shim just re-shapes. Connections/helices/sheets are not model-scoped in gemmi,
   // so they go on model 1 (MMDB's usual home); cispeps honour their model_num.
   inline void Manager::_load_metadata() {
      link_pool.clear();
      linkr_pool.clear();
      cispep_pool.clear();
      helix_pool.clear();
      sheet_pool.clear();
      strand_pool.clear();
      strandarr_pool.clear();
      author_pool.clear();

      // PDB title AUTHOR records (gemmi meta.authors)
      title.author.data.clear();
      for (const std::string &au : st.meta.authors) {
         author_pool.emplace_back();
         std::snprintf(author_pool.back().Line, sizeof(author_pool.back().Line), "%s", au.c_str());
         title.author.data.push_back(&author_pool.back());
      }
      if (models.empty()) return;

      auto fill_ends = [](const gemmi::AtomAddress &a, ChainID &cid, ResName &rn,
                          int &seq, InsCode &ic, AtomName *an, AltLoc *al) {
         std::snprintf(cid, sizeof(ChainID), "%s", a.chain_name.c_str());
         std::snprintf(rn, sizeof(ResName), "%s", a.res_id.name.c_str());
         seq = a.res_id.seqid.num.value;
         ic[0] = (a.res_id.seqid.icode && a.res_id.seqid.icode != ' ') ? a.res_id.seqid.icode : '\0';
         ic[1] = '\0';
         if (an) std::snprintf(*an, sizeof(AtomName), "%s", a.atom_name.c_str());
         if (al) {
            (*al)[0] = a.altloc ? a.altloc : '\0';
            (*al)[1] = '\0';
         }
      };

      // --- LINK records (gemmi Connection) -> model 1 ---
      Model *m1 = models[0];
      for (const gemmi::Connection &cn : st.connections) {
         link_pool.emplace_back();
         Link &l = link_pool.back();
         fill_ends(cn.partner1, l.chainID1, l.resName1, l.seqNum1, l.insCode1, &l.atName1, &l.aloc1);
         fill_ends(cn.partner2, l.chainID2, l.resName2, l.seqNum2, l.insCode2, &l.atName2, &l.aloc2);
         l.dist = cn.reported_distance;
         m1->_links.push_back(&l);
         // a connection carrying a Refmac link id is also a LINKR record
         if (!cn.link_id.empty()) {
            linkr_pool.emplace_back();
            LinkR &lr = linkr_pool.back();
            std::snprintf(lr.linkRID, sizeof(lr.linkRID), "%s", cn.link_id.c_str());
            AtomName an;
            AltLoc al;
            fill_ends(cn.partner1, lr.chainID1, lr.resName1, lr.seqNum1, lr.insCode1, &an, &al);
            std::snprintf(lr.atName1, sizeof(AtomName), "%s", an);
            std::snprintf(lr.aloc1, sizeof(AltLoc), "%s", al);
            fill_ends(cn.partner2, lr.chainID2, lr.resName2, lr.seqNum2, lr.insCode2, &an, &al);
            std::snprintf(lr.atName2, sizeof(AtomName), "%s", an);
            std::snprintf(lr.aloc2, sizeof(AltLoc), "%s", al);
            lr.dist = cn.reported_distance;
            m1->_linkrs.push_back(&lr);
         }
      }

      // --- CISPEP records (gemmi CisPep) -> model by model_num (default 1) ---
      for (const gemmi::CisPep &cp : st.cispeps) {
         int mnum = cp.model_num > 0 ? cp.model_num : 1;
         Model *mw = GetModel(mnum);
         if (!mw) mw = m1;
         cispep_pool.emplace_back();
         CisPep &c = cispep_pool.back();
         InsCode ic1, ic2;
         int s1, s2;
         fill_ends(cp.partner_c, c.chainID1, c.pep1, s1, ic1, nullptr, nullptr);
         fill_ends(cp.partner_n, c.chainID2, c.pep2, s2, ic2, nullptr, nullptr);
         c.seqNum1 = s1;
         std::snprintf(c.icode1, sizeof(InsCode), "%s", ic1);
         c.seqNum2 = s2;
         std::snprintf(c.icode2, sizeof(InsCode), "%s", ic2);
         c.modNum = mnum;
         if (!std::isnan(cp.reported_angle)) c.measure = cp.reported_angle;
         mw->_cispeps.push_back(&c);
      }

      // --- HELIX records (gemmi Helix) -> model 1 ---
      for (const gemmi::Helix &gh : st.helices) {
         helix_pool.emplace_back();
         Helix &h = helix_pool.back();
         AtomName an;
         AltLoc al;
         fill_ends(gh.start, h.initChainID, h.initResName, h.initSeqNum, h.initICode, &an, &al);
         fill_ends(gh.end, h.endChainID, h.endResName, h.endSeqNum, h.endICode, &an, &al);
         h.helixClass = (int)gh.pdb_helix_class;
         h.length = gh.length;
         h.serNum = (int)helix_pool.size();
         m1->helices.AddData(&h);
      }

      // --- SHEET / STRAND records (gemmi Sheet) -> model 1 ---
      if (!st.sheets.empty()) {
         m1->sheets.nSheets = (int)st.sheets.size();
         m1->_sheet_ptrs.assign(st.sheets.size(), nullptr);  // backs Sheets::sheet (Sheet**)
         for (size_t is = 0; is < st.sheets.size(); ++is) {
            const gemmi::Sheet &gs = st.sheets[is];
            sheet_pool.emplace_back();
            Sheet &sh = sheet_pool.back();
            std::snprintf(sh.sheetID, sizeof(sh.sheetID), "%s", gs.name.c_str());
            sh.nStrands = (int)gs.strands.size();
            strandarr_pool.emplace_back();
            std::vector<PStrand> &sarr = strandarr_pool.back();
            sarr.reserve(gs.strands.size());
            for (const gemmi::Sheet::Strand &gst : gs.strands) {
               strand_pool.emplace_back();
               Strand &str = strand_pool.back();
               AtomName an;
               AltLoc al;
               fill_ends(gst.start, str.initChainID, str.initResName, str.initSeqNum, str.initICode, &an, &al);
               fill_ends(gst.end, str.endChainID, str.endResName, str.endSeqNum, str.endICode, &an, &al);
               std::snprintf(str.sheetID, sizeof(str.sheetID), "%s", gs.name.c_str());
               str.strandNo = (int)sarr.size() + 1;
               str.sense = gst.sense;
               sarr.push_back(&str);
            }
            sh.strand = sarr.data();
            m1->_sheet_ptrs[is] = &sh;
         }
         m1->sheets.sheet = m1->_sheet_ptrs.data();
      }
   }

   // ---- Manager::PutAtom (hierarchy insertion) ----
   inline int Manager::PutAtom(int index, PAtom A, int serNum) {
      if (!A) return 0;
      Residue *src = A->res;
      // ensure a model exists (Coot calls PutAtom on a fresh, empty Manager)
      Model *mw = models.empty() ? nullptr : models[0];
      if (!mw) {
         st.models.emplace_back(1);
         mw = newModel();
         mw->mgr = this;
         mw->mi = 0;
         models.push_back(mw);
      }
      // find or create the chain implied by the source atom's chain
      std::string cid = (src && src->chain) ? src->chain->g().name : std::string("A");
      Chain *cw = mw->GetChain(cid.c_str());
      if (!cw) cw = mw->CreateChain(cid.c_str());
      // find or create the residue implied by (seqNum, insCode)
      int seq = src ? src->g().seqid.num.value : 0;
      char ic = src ? src->g().seqid.icode : ' ';
      char icn = ic ? ic : ' ';
      Residue *rw = nullptr;
      for (Residue *r : cw->residues) {
         gemmi::Residue &gr = r->g();
         if (gr.seqid.num.value == seq && (gr.seqid.icode ? gr.seqid.icode : ' ') == icn) {
            rw = r;
            break;
         }
      }
      if (!rw) {
         gemmi::Residue gr;
         gr.name = src ? src->g().name : std::string("UNK");
         gr.seqid.num = seq;
         gr.seqid.icode = icn;
         rw = cw->AddResidue(*this, gr);
         rw->_load_id();
      }
      // append a copy of the atom's gemmi backing + register it in the flat tables
      Atom *aw = rw->AddAtom(*this, A->g());
      aw->WhatIsSet = A->WhatIsSet;
      aw->Het = A->Het;
      std::memcpy(aw->segID, A->segID, sizeof aw->segID);
      aw->g().serial = serNum ? serNum : (index > 0 ? index : (int)all_atoms.size() + 1);
      rw->_sync_atom();
      all_atoms.push_back(aw);
      mw->all_atoms.push_back(aw);
      return (int)all_atoms.size();  // 1-based position (GetAtomI(pos) returns aw)
   }

   // ---- selection matching ----
   namespace detail {
      inline std::string trimws(const std::string &s) {
         size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
         return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
      }
      // Whitespace-insensitive membership test. MMDB atom names are stored space-
      // padded ("_CA_", "_O__"), while CID/selection queries are unpadded ("CA",
      // "O"); real MMDB matches them regardless of padding, so trim both sides.
      // Harmless for chain IDs / residue / element names (already unpadded).
      inline bool inList(cpstr list, const std::string &v) {
         if (!list || !*list || std::strcmp(list, "*") == 0) return true;
         // MMDB negation: a leading '!' inverts the match (e.g. "!HOH" = any residue
         // that is not water). Coot's Select() uses this for chain/residue/element/
         // atom-name filters; without it every residue is (wrongly) excluded.
         if (list[0] == '!') return !inList(list + 1, v);
         std::string vt = trimws(v);
         const char *p = list;
         while (*p) {
            const char *c = std::strchr(p, ',');
            std::string tok(p, c ? (size_t)(c - p) : std::strlen(p));
            if (trimws(tok) == vt) return true;
            if (!c) break;
            p = c + 1;
         }
         return false;
      }
      inline bool altMatch(cpstr list, char alt) {
         if (!list || std::strcmp(list, "*") == 0) return true;
         std::string a = alt ? std::string(1, alt) : std::string();
         if (!*list) return a.empty();  // "" -> only blank altLoc
         return inList(list, a);
      }
   }  // namespace detail

   inline void Manager::Select(int selHnd, SELECTION_TYPE sType, int iModel,
                               cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames,
                               cpstr ANames, cpstr Elements, cpstr altLocs, SELECTION_KEY selKey) {
      Selection &sel = selections[selHnd - 1];
      if (sel.type == STYPE_UNDEFINED) sel.type = sType;
      std::vector<Atom *> oldA = sel.atoms;
      std::vector<Residue *> oldR = sel.residues;
      std::vector<Chain *> oldC = sel.chains;

      std::vector<Atom *> mAtoms;
      std::vector<Residue *> mResidues;
      std::vector<Chain *> mChains;
      for (Model *mw : models) {
         if (iModel > 0 && mw->GetSerNum() != iModel) continue;
         for (Chain *cw : mw->chains) {
            if (!detail::inList(Chains, cw->g().name)) continue;
            bool anyResidue = false;
            for (Residue *rw : cw->residues) {
               int sn = rw->g().seqid.num.value;
               char ric = rw->g().seqid.icode ? rw->g().seqid.icode : ' ';
               // (seqNum, insCode) range: an explicit insCode only constrains the
               // boundary residue; blank/"*" includes every insCode at that seqNum.
               if (ResNo1 != ANY_RES) {
                  if (sn < ResNo1) continue;
                  if (sn == ResNo1 && Ins1 && Ins1[0] && std::strcmp(Ins1, "*") && ric < Ins1[0]) continue;
               }
               if (ResNo2 != ANY_RES) {
                  if (sn > ResNo2) continue;
                  if (sn == ResNo2 && Ins2 && Ins2[0] && std::strcmp(Ins2, "*") && ric > Ins2[0]) continue;
               }
               if (!detail::inList(RNames, rw->g().name)) continue;
               bool anyAtom = false;
               for (Atom *aw : rw->atoms) {
                  if (!detail::inList(ANames, aw->g().name)) continue;
                  if (!detail::inList(Elements, aw->g().element.name())) continue;
                  if (!detail::altMatch(altLocs, aw->g().altloc)) continue;
                  anyAtom = true;
                  if (sType == STYPE_ATOM) mAtoms.push_back(aw);
               }
               if (anyAtom) anyResidue = true;
               if (anyAtom && sType == STYPE_RESIDUE) mResidues.push_back(rw);
            }
            // STYPE_CHAIN: a chain matching the chain filter (and, if given, having a
            // residue that passes the residue/atom filters) is selected whole.
            if (sType == STYPE_CHAIN && anyResidue) mChains.push_back(cw);
         }
      }
      auto combine = [&](auto &cur, auto &matched) {
         using Vec = typename std::decay<decltype(cur)>::type;
         std::set<typename Vec::value_type> curset(cur.begin(), cur.end());
         std::set<typename Vec::value_type> mset(matched.begin(), matched.end());
         if (selKey == SKEY_NEW) {
            cur = matched;
         } else if (selKey == SKEY_OR) {
            for (auto *x : matched)
               if (!curset.count(x)) cur.push_back(x);
         } else if (selKey == SKEY_AND) {
            Vec o;
            for (auto *x : cur)
               if (mset.count(x)) o.push_back(x);
            cur = o;
         } else if (selKey == SKEY_XOR) {
            Vec o;
            for (auto *x : cur)
               if (!mset.count(x)) o.push_back(x);
            for (auto *x : matched)
               if (!curset.count(x)) o.push_back(x);
            cur = o;
         } else if (selKey == SKEY_CLR) {
            Vec o;
            for (auto *x : cur)
               if (!mset.count(x)) o.push_back(x);
            cur = o;
         }
      };
      if (sType == STYPE_ATOM)
         combine(sel.atoms, mAtoms);
      else if (sType == STYPE_RESIDUE)
         combine(sel.residues, mResidues);
      else if (sType == STYPE_CHAIN)
         combine(sel.chains, mChains);
      for (Atom *a : oldA) a->_setInSel(selHnd, false);
      for (Atom *a : sel.atoms) a->_setInSel(selHnd, true);
      for (Residue *r : oldR) r->_setInSel(selHnd, false);
      for (Residue *r : sel.residues) r->_setInSel(selHnd, true);
      for (Chain *c : oldC) c->_setInSel(selHnd, false);
      for (Chain *c : sel.chains) c->_setInSel(selHnd, true);
   }

   // select-from-selection: combine selHnd2's contents into selHnd1
   inline void Manager::Select(int selHnd1, SELECTION_TYPE sType, int selHnd2,
                               SELECTION_KEY sKey) {
      Selection &s1 = selections[selHnd1 - 1];
      Selection &s2 = selections[selHnd2 - 1];
      if (s1.type == STYPE_UNDEFINED) s1.type = sType;
      std::vector<Atom *> oldA = s1.atoms;
      std::vector<Residue *> oldR = s1.residues;
      auto combine = [&](auto &cur, auto &m) {
         using Vec = typename std::decay<decltype(cur)>::type;
         std::set<typename Vec::value_type> curset(cur.begin(), cur.end());
         std::set<typename Vec::value_type> mset(m.begin(), m.end());
         if (sKey == SKEY_NEW)
            cur = m;
         else if (sKey == SKEY_OR) {
            for (auto *x : m)
               if (!curset.count(x)) cur.push_back(x);
         } else if (sKey == SKEY_AND) {
            Vec o;
            for (auto *x : cur)
               if (mset.count(x)) o.push_back(x);
            cur = o;
         } else if (sKey == SKEY_XOR) {
            Vec o;
            for (auto *x : cur)
               if (!mset.count(x)) o.push_back(x);
            for (auto *x : m)
               if (!curset.count(x)) o.push_back(x);
            cur = o;
         } else if (sKey == SKEY_CLR) {
            Vec o;
            for (auto *x : cur)
               if (!mset.count(x)) o.push_back(x);
            cur = o;
         }
      };
      if (sType == STYPE_ATOM)
         combine(s1.atoms, s2.atoms);
      else if (sType == STYPE_RESIDUE)
         combine(s1.residues, s2.residues);
      for (Atom *a : oldA) a->_setInSel(selHnd1, false);
      for (Atom *a : s1.atoms) a->_setInSel(selHnd1, true);
      for (Residue *r : oldR) r->_setInSel(selHnd1, false);
      for (Residue *r : s1.residues) r->_setInSel(selHnd1, true);
   }

   inline void Manager::SelectAtom(int selHnd, PAtom atom, SELECTION_KEY sKey, bool) {
      Selection &sel = selections[selHnd - 1];
      if (sel.type == STYPE_UNDEFINED) sel.type = STYPE_ATOM;
      if (sKey == SKEY_NEW) {
         for (Atom *a : sel.atoms) a->_setInSel(selHnd, false);
         sel.atoms.clear();
      }
      if (atom && !atom->isInSelection(selHnd)) {
         sel.atoms.push_back(atom);
         atom->_setInSel(selHnd, true);
      }
   }

   // Pragmatic CID parser: "/model/chain/seqNum1[.ins1]-seqNum2[.ins2]/atom"
   // (best-effort; strips (resname)/[element]/:altloc suffixes; parses insertion
   // codes after '.'). Not the full MMDB CID grammar but covers Coot's usage.
   inline void Manager::Select(int selHnd, SELECTION_TYPE sType, cpstr CID,
                               SELECTION_KEY sKey) {
      std::string s = CID ? CID : "";
      std::vector<std::string> t;
      size_t p = (!s.empty() && s[0] == '/') ? 1 : 0;
      while (p <= s.size()) {
         size_t q = s.find('/', p);
         t.push_back(s.substr(p, q == std::string::npos ? std::string::npos : q - p));
         if (q == std::string::npos) break;
         p = q + 1;
      }
      auto tok = [&](size_t i) { return i < t.size() ? t[i] : std::string(); };
      auto strip = [](std::string v, const char *seps) {
         size_t c = v.find_first_of(seps);
         return c == std::string::npos ? v : v.substr(0, c);
      };
      // Assign tokens to model/chain/residue/atom. A leading '/' (or any '/') means
      // the model field is present at tok(0). A slash-less CID has NO model/chain
      // prefix: MMDB reads a bare numeric token as a residue seqNum ("262" = residue
      // 262 in every chain), and a bare non-numeric token as a chain id ("A").
      std::string model_s, chain_s, res_s, atom_s;
      if (s.find('/') != std::string::npos) {
         model_s = tok(0);
         chain_s = tok(1);
         res_s = tok(2);
         atom_s = tok(3);
      } else {
         const std::string only = tok(0);
         if (!only.empty() && (std::isdigit((unsigned char)only[0]) || only[0] == '-'))
            res_s = only;
         else
            chain_s = only;
      }
      int iModel = 0;
      if (!model_s.empty() && model_s != "*" && model_s != "0") iModel = atoi(model_s.c_str());
      std::string chains = chain_s.empty() ? "*" : chain_s;
      int r1 = ANY_RES, r2 = ANY_RES;
      std::string ins1 = "*", ins2 = "*";
      // split "num[.ins]" into number + insertion code
      auto parse_resid = [](const std::string &v, int &num, std::string &ins) {
         size_t dot = v.find('.');
         num = atoi(v.substr(0, dot).c_str());
         ins = (dot == std::string::npos) ? std::string() : v.substr(dot + 1);
      };
      std::string rr = strip(res_s, "(");  // drop (resname)
      if (!rr.empty() && rr != "*") {
         size_t dash = rr.find('-', rr[0] == '-' ? 1 : 0);
         if (dash == std::string::npos) {
            parse_resid(rr, r1, ins1);
            r2 = r1;
            ins2 = ins1;
         } else {
            parse_resid(rr.substr(0, dash), r1, ins1);
            parse_resid(rr.substr(dash + 1), r2, ins2);
         }
      }
      std::string anames = strip(strip(atom_s, "["), ":");  // drop [element]/:altloc
      if (anames.empty()) anames = "*";
      Select(selHnd, sType, iModel, chains.c_str(), r1, ins1.c_str(), r2, ins2.c_str(), "*",
             anames.c_str(), "*", "*", sKey);
   }

   inline int Manager::GetNumberOfAtoms(cpstr CID) {
      int h = NewSelection();
      Select(h, STYPE_ATOM, CID, SKEY_NEW);
      int n = (int)selections[h - 1].atoms.size();
      DeleteSelection(h);
      return n;
   }

   inline void Manager::GetAtomStatistics(int selHnd, RAtomStat AS) {
      AS = AtomStat();
      std::vector<Atom *> &atoms = selections[selHnd - 1].atoms;
      AS.nAtoms = (int)atoms.size();
      if (atoms.empty()) return;
      double sx = 0, sy = 0, sz = 0;
      AS.xmin = AS.xmax = atoms[0]->x();
      AS.ymin = AS.ymax = atoms[0]->y();
      AS.zmin = AS.zmax = atoms[0]->z();
      for (Atom *a : atoms) {
         double X = a->x(), Y = a->y(), Z = a->z();
         sx += X;
         sy += Y;
         sz += Z;
         AS.xmin = X < AS.xmin ? X : AS.xmin;
         AS.xmax = X > AS.xmax ? X : AS.xmax;
         AS.ymin = Y < AS.ymin ? Y : AS.ymin;
         AS.ymax = Y > AS.ymax ? Y : AS.ymax;
         AS.zmin = Z < AS.zmin ? Z : AS.zmin;
         AS.zmax = Z > AS.zmax ? Z : AS.zmax;
      }
      AS.xm = sx / atoms.size();
      AS.ym = sy / atoms.size();
      AS.zm = sz / atoms.size();
   }

   inline void Manager::SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x,
                                     realtype y, realtype z, realtype r, SELECTION_KEY sKey) {
      Selection &sel = selections[selHnd - 1];
      if (sel.type == STYPE_UNDEFINED) sel.type = sType;
      std::vector<Atom *> oldA = sel.atoms;
      std::vector<Residue *> oldR = sel.residues;
      gemmi::Position c(x, y, z);
      double r2 = r * r;
      std::vector<Atom *> mAtoms;
      std::vector<Residue *> mResidues;
      for (Model *mw : models)
         for (Chain *cw : mw->chains)
            for (Residue *rw : cw->residues) {
               bool any = false;
               for (Atom *aw : rw->atoms)
                  if (aw->g().pos.dist_sq(c) <= r2) {
                     any = true;
                     if (sType == STYPE_ATOM) mAtoms.push_back(aw);
                  }
               if (any && sType == STYPE_RESIDUE) mResidues.push_back(rw);
            }
      auto combine = [&](auto &cur, auto &m) {
         std::set<typename std::decay<decltype(cur)>::type::value_type> cs(cur.begin(), cur.end());
         if (sKey == SKEY_NEW)
            cur = m;
         else if (sKey == SKEY_OR) {
            for (auto *p : m)
               if (!cs.count(p)) cur.push_back(p);
         }
      };
      if (sType == STYPE_ATOM)
         combine(sel.atoms, mAtoms);
      else if (sType == STYPE_RESIDUE)
         combine(sel.residues, mResidues);
      for (Atom *a : oldA) a->_setInSel(selHnd, false);
      for (Atom *a : sel.atoms) a->_setInSel(selHnd, true);
      for (Residue *r : oldR) r->_setInSel(selHnd, false);
      for (Residue *r : sel.residues) r->_setInSel(selHnd, true);
   }

   // SeekContacts (both overloads) is defined in mmdb-shim/src/contacts.cc using
   // gemmi::NeighborSearch — keeps the heavy neighbor.hpp out of Coot's many TUs.

   // ---- detached-construction constructors + subtree ops (need complete types) ----
   inline Atom::Atom(Residue *r) {
      if (r) r->AddAtom(this);
   }
   inline Residue::Residue(Chain *c) {
      if (c) c->AddResidue(this);
   }
   inline Chain::Chain(Model *m, const ChainID id) {
      if (m) m->AddChain(this);
      SetChainID(id);
   }

   // peptide-bond distance threshold for backbone C-N (a real bond is ~1.33 A).
   inline bool Residue::isNTerminus() {
      if (!chain || ri <= 0) return true;  // first (or detached) residue
      Residue *prev = chain->residues[ri - 1];
      if (!prev) return true;  // previous slot is a deferred-delete tombstone
      const gemmi::Atom *N = g().get_n();
      const gemmi::Atom *prevC = prev->g().get_c();
      if (!N || !prevC) return true;         // missing backbone -> terminus
      return N->pos.dist(prevC->pos) > 1.7;  // not bonded to previous C
   }
   inline bool Residue::isCTerminus() {
      if (!chain || ri < 0 || ri >= (int)chain->residues.size() - 1) return true;  // last/detached
      Residue *next = chain->residues[ri + 1];
      if (!next) return true;  // next slot is a deferred-delete tombstone
      const gemmi::Atom *C = g().get_c();
      const gemmi::Atom *nextN = next->g().get_n();
      if (!C || !nextN) return true;
      return C->pos.dist(nextN->pos) > 1.7;  // not bonded to next N
   }
   inline Model *Residue::GetModel() { return chain ? chain->model : nullptr; }

   inline void Chain::Copy(PChain src) {
      Manager *pool = mgr ? mgr : src->mgr;
      g() = src->g();  // deep gemmi copy (residues + atoms)
      residues.clear();
      if (!pool) return;
      gemmi::Chain &gc = g();
      for (int r = 0; r < (int)gc.residues.size(); ++r) {
         Residue *rw = pool->newRes();
         rw->mgr = mgr;
         rw->chain = this;
         rw->ri = r;
         for (int a = 0; a < (int)gc.residues[r].atoms.size(); ++a) {
            Atom *aw = pool->newAtom();
            aw->mgr = mgr;
            aw->res = rw;
            aw->ai = a;
            rw->atoms.push_back(aw);
         }
         rw->_sync_atom();
         rw->_load_id();
         residues.push_back(rw);
      }
   }

   inline void Model::Copy(PModel src) {
      Manager *pool = mgr ? mgr : src->mgr;
      g() = src->g();
      chains.clear();
      if (!pool) return;
      gemmi::Model &gm = g();
      for (int c = 0; c < (int)gm.chains.size(); ++c) {
         Chain *cw = pool->newChain();
         cw->mgr = mgr;
         cw->model = this;
         cw->ci = c;
         for (int r = 0; r < (int)gm.chains[c].residues.size(); ++r) {
            Residue *rw = pool->newRes();
            rw->mgr = mgr;
            rw->chain = cw;
            rw->ri = r;
            for (int a = 0; a < (int)gm.chains[c].residues[r].atoms.size(); ++a) {
               Atom *aw = pool->newAtom();
               aw->mgr = mgr;
               aw->res = rw;
               aw->ai = a;
               rw->atoms.push_back(aw);
            }
            rw->_sync_atom();
            rw->_load_id();
            cw->residues.push_back(rw);
         }
         chains.push_back(cw);
      }
   }

   inline PChain Model::CreateChain(const ChainID id) {
      Chain *c = mgr ? mgr->newChain() : new Chain();
      c->mgr = mgr;
      c->model = this;
      c->ci = (int)chains.size();
      g().chains.emplace_back(id ? id : "");
      chains.push_back(c);
      return c;
   }

   inline pstr Atom::GetAtomID(pstr S) {
      if (S) std::snprintf(S, 100, "/%d/%s/%d(%s)/%s", GetModelNum(), GetChainID(),
                           res ? res->GetSeqNum() : 0, GetResName(), GetAtomName());
      return S;
   }

   // one-letter residue code (mmdb_tables.h) via gemmi's tabulated residues
   inline void Get1LetterCode(cpstr res3, pstr res1) {
      if (!res1) return;
      char c = gemmi::find_tabulated_residue(res3 ? res3 : "").one_letter_code;
      res1[0] = c ? (char)std::toupper((unsigned char)c) : 'X';
      res1[1] = '\0';
   }
   inline void Get1LetterCode(cpstr res3, char &res1) {
      char b[2];
      Get1LetterCode(res3, b);
      res1 = b[0];
   }

   // sort a contact array by distance (mmdb_coormngr.h SortContacts) — sortkey ignored
   inline void SortContacts(PContact contacts, int nContacts, int /*sortkey*/) {
      if (contacts && nContacts > 1)
         std::sort(contacts, contacts + nContacts,
                   [](const Contact &a, const Contact &b) { return a.dist < b.dist; });
   }

   // centroid of an atom array (mmdb_coormngr.h GetMassCenter)
   inline void GetMassCenter(PPAtom A, int nA, realtype &xc, realtype &yc, realtype &zc) {
      double sx = 0, sy = 0, sz = 0;
      int n = 0;
      for (int i = 0; i < nA; ++i)
         if (A[i]) {
            sx += A[i]->x();
            sy += A[i]->y();
            sz += A[i]->z();
            ++n;
         }
      if (n) {
         xc = sx / n;
         yc = sy / n;
         zc = sz / n;
      } else {
         xc = yc = zc = 0;
      }
   }

}  // namespace mmdb
