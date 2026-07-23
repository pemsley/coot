// mmdb-shim — layer 2 of 4: the coordinate hierarchy (Atom / Residue / Chain /
// Model).
//
// These classes ARE the stable wrapper nodes: each holds Manager* + parent* + a
// cached sibling index, resolves to the live gemmi object via g(), and exposes its
// children as a canonical vector<T*> that doubles as the identity cache AND the
// PPAtom/PPResidue table MMDB hands back. A wrapper with no parent is "detached"
// and backs onto a wrapper-owned _local gemmi object until Add*() adopts it.
//
// Only class definitions live here; method bodies that need a complete Manager (or
// sibling class) are declared here and defined out-of-line in _shim_inline.hh.
#pragma once

#include "_shim_types.hh"

namespace mmdb {

   // ===========================================================================
   class Atom : public UDStore {
   public:
      Manager *mgr = nullptr;
      Residue *res = nullptr;  // parent; null => detached (use _local)
      int ai = 0;              // cached index within parent residue's atoms
      bool alive = true;
      gemmi::Atom _local;        // backing store while detached (see g() resolvers)
      int Het = 0;               // heteroatom flag (MMDB public field; Coot sets it)
      int Ter = 0;               // chain-terminator flag (gemmi has none -> always 0)
      word WhatIsSet = 0;        // ASET_* mask; ASET_Anis_tFac set on load if aniso present
      AtomName label_atom_id{};  // mmcif label_atom_id (shim-owned; Coot sets on build)

      Atom() = default;
      explicit Atom(Residue *r);  // construct + add to residue (out-of-line)
      // MMDB owns atoms via `new`/`delete`: Coot writes `delete atom;` to remove an
      // atom from the hierarchy (coot-molecule.cc et al.). So atoms are individually
      // heap-allocated (Manager::newAtom) and tracked in Manager::_atom_allocs; this
      // destructor detaches from the parent residue + flat lists + selections when a
      // live atom is deleted, and is a no-op during Manager teardown (_bulk_free).
      ~Atom();  // out-of-line (needs complete Manager/Residue)

      gemmi::Atom &g() const;  // resolve to live gemmi (defined after Manager)

      // --- rewritten field accessors (pure B) ---
      // Scalar fields -> reference-returning accessors, so a uniform `->field`->
      // `->field()` rewrite covers both reads and writes. (occ/b_iso/charge are
      // narrower than realtype in gemmi, so those refs are float/schar-typed — the
      // rare take-address-of-realtype sites surface at Coot build time.)
      // non-const (writable ref) + const (by value) overloads, so reads work on a
      // `const mmdb::Atom` and writes work through `->x() = v` on a non-const one.
      realtype &x() { return g().pos.x; }
      realtype x() const { return g().pos.x; }
      realtype &y() { return g().pos.y; }
      realtype y() const { return g().pos.y; }
      realtype &z() { return g().pos.z; }
      realtype z() const { return g().pos.z; }
      float &occupancy() { return g().occ; }
      float occupancy() const { return g().occ; }
      float &tempFactor() { return g().b_iso; }
      float tempFactor() const { return g().b_iso; }
      signed char &charge() { return g().charge; }
      signed char charge() const { return g().charge; }
      int &serNum() { return g().serial; }
      int serNum() const { return g().serial; }
      // altLoc is a char[] (C-string) in MMDB; gemmi stores a single char. Return a
      // buffer-backed C-string ("" when unset) so strcmp/strcpy-style code works.
      // The non-const overload returns a WRITABLE buffer so `strncpy(at->altLoc(),..)`
      // compiles; the buffer's first char is pushed back into gemmi by Residue::AddAtom
      // (the buffer is refreshed from gemmi on entry, so reads stay correct).
      pstr altLoc() {
         _altloc_buf[0] = g().altloc;
         _altloc_buf[1] = '\0';
         return _altloc_buf;
      }
      const char *altLoc() const {
         _altloc_buf[0] = g().altloc;
         _altloc_buf[1] = '\0';
         return _altloc_buf;
      }
      void set_occupancy(realtype v) { g().occ = (float)v; }
      void set_tempFactor(realtype v) { g().b_iso = (float)v; }
      void set_altLoc(char c) { g().altloc = c; }
      void SetCharge(realtype ch) { g().charge = (signed char)ch; }
      // coordinate/occupancy/B ESDs (MMDB public fields) — gemmi has none, so shim-
      // owned; reference-returning so the rewritten `->sigX` covers reads and writes.
      float &sigX() { return _sigx; }
      float &sigY() { return _sigy; }
      float &sigZ() { return _sigz; }
      float &sigOcc() { return _sigocc; }
      float &sigTemp() { return _sigtemp; }
      bool isMetal() const { return gemmi::Element(g().element).is_metal(); }
      // anisotropic B tensor — gemmi's SMat33<float> aniso. Reference-returning so the
      // rewritten `->u11` covers both reads and writes. The mutable accessor marks the
      // tensor present (ASET_Anis_tFac) so a write (e.g. SHELX import) sets the flag as
      // real MMDB does. Const reads never set it; a non-const read over-approximates,
      // which is harmless — the PDB/mmCIF writer emits ANISOU on the actual values.
      float &u11() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u11;
      }
      float &u22() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u22;
      }
      float &u33() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u33;
      }
      float &u12() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u12;
      }
      float &u13() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u13;
      }
      float &u23() {
         WhatIsSet |= ASET_Anis_tFac;
         return g().aniso.u23;
      }
      float u11() const { return g().aniso.u11; }
      float u22() const { return g().aniso.u22; }
      float u33() const { return g().aniso.u33; }
      float u12() const { return g().aniso.u12; }
      float u13() const { return g().aniso.u13; }
      float u23() const { return g().aniso.u23; }
      // bonds — not modelled yet (gemmi connections); report none.
      int GetNBonds() { return 0; }
      void GetBonds(PAtomBond &atomBond, int &n) {
         atomBond = nullptr;
         n = 0;
      }
      int AddBond(PAtom /*a*/, int /*order*/, int /*nAdd*/ = 1) { return 0; }
      SegID segID{};  // shim-owned (gemmi has no segID); MMDB public char[] field

      // --- method surface (hot subset; rest stubbed) ---
      pstr GetAtomName() const;  // aligned name, MMDB semantics (const: called on const Atom)
      void SetAtomName(const AtomName aName);
      pstr GetElementName();
      void SetElementName(const Element elName);
      pstr GetChainID();
      int GetSeqNum();
      pstr GetInsCode();
      pstr GetResName();
      Residue *&GetResidue() { return res; }  // ref: rewritten `->residue` is assignable
      void SetResidue(Residue *r) { res = r; }
      Chain *GetChain();  // out-of-line (needs complete Residue/Chain)
      Model *GetModel();  // out-of-line
      int GetModelNum();
      // residue-delegating accessors (bound by the Python API); out-of-line.
      pstr GetLabelCompID();
      pstr GetLabelAsymID();
      int GetLabelSeqID();
      int GetLabelEntityID();
      int GetResidueNo();
      int GetSSEType();
      bool isSolvent();
      bool isNTerminus();
      bool isCTerminus();
      bool isTer() const { return false; }  // gemmi has no TER atoms; see notes
      void SetCoordinates(realtype xx, realtype yy, realtype zz,
                          realtype occ, realtype tF);
      int GetIndex();
      void MakeTer() { Ter = 1; }  // mark as chain terminator
      pstr GetAtomID(pstr S);      // "/mdl/chain/seq(res).ins/name[elem]:alt" (out-of-line)
      int GetUDData(int h, pstr &v) { return ud_get(mgr, UDR_ATOM, *this, h, v); }
      // copy another atom's data into this one (mmdb Atom::Copy — no hierarchy refs)
      void Copy(PAtom a) {
         g() = a->g();
         Het = a->Het;
         WhatIsSet = a->WhatIsSet;
         std::memcpy(segID, a->segID, sizeof segID);
      }
      // apply a 4x4 (rot+trans) or 3x3+vec to the coordinates (mmdb Atom::Transform)
      void Transform(const mat44 &tm) {
         gemmi::Position &p = g().pos;
         double x = p.x, y = p.y, z = p.z;
         p.x = tm[0][0] * x + tm[0][1] * y + tm[0][2] * z + tm[0][3];
         p.y = tm[1][0] * x + tm[1][1] * y + tm[1][2] * z + tm[1][3];
         p.z = tm[2][0] * x + tm[2][1] * y + tm[2][2] * z + tm[2][3];
      }
      void Transform(const mat33 &tm, vect3 &v) {
         gemmi::Position &p = g().pos;
         double x = p.x, y = p.y, z = p.z;
         p.x = tm[0][0] * x + tm[0][1] * y + tm[0][2] * z + v[0];
         p.y = tm[1][0] * x + tm[1][1] * y + tm[1][2] * z + v[1];
         p.z = tm[2][0] * x + tm[2][1] * y + tm[2][2] * z + v[2];
      }
      // UDData
      int PutUDData(int h, int v) { return ud_put(mgr, UDR_ATOM, *this, h, v); }
      int PutUDData(int h, realtype v) { return ud_put(mgr, UDR_ATOM, *this, h, v); }
      int PutUDData(int h, cpstr v) { return ud_put(mgr, UDR_ATOM, *this, h, v); }
      int GetUDData(int h, int &v) { return ud_get(mgr, UDR_ATOM, *this, h, v); }
      int GetUDData(int h, realtype &v) { return ud_get(mgr, UDR_ATOM, *this, h, v); }

   private:
      friend class Residue;  // AddAtom pushes the strncpy'd altLoc buffer to gemmi
      mutable AtomName _name_buf{};
      Element _elem_buf{};
      mutable char _altloc_buf[4]{};
      float _sigx = 0, _sigy = 0, _sigz = 0, _sigocc = 0, _sigtemp = 0;
   };

   // ===========================================================================
   class Residue : public UDStore {
   public:
      Manager *mgr = nullptr;
      Chain *chain = nullptr;  // parent; null => detached (use _local)
      int ri = 0;
      bool alive = true;
      gemmi::Residue _local;      // backing store while detached
      std::vector<Atom *> atoms;  // canonical child wrappers == PPAtom table
      PPAtom atom = nullptr;      // MMDB public atom-table field; kept = atoms.data()
      int nAtoms = 0;             // MMDB public field; kept = atoms.size()
      void _sync_atom() {
         atom = atoms.data();
         nAtoms = (int)atoms.size();
      }
      // mmcif label_* (shim-owned; Coot sets when building dictionary residues)
      ResName label_comp_id{};
      ChainID label_asym_id{};
      int label_seq_id = 0, label_entity_id = 0;
      pstr GetLabelCompID() { return label_comp_id; }
      pstr GetLabelAsymID() { return label_asym_id; }
      int GetLabelSeqID() { return label_seq_id; }
      int GetLabelEntityID() { return label_entity_id; }
      int GetResidueNo() { return ri; }  // 0-based index within its chain
      int GetNofAltLocations() {         // distinct non-blank altLocs
         std::set<char> a;
         for (Atom *at : atoms) {
            char c = at->g().altloc;
            if (c && c != ' ') a.insert(c);
         }
         return a.empty() ? 1 : (int)a.size();
      }
      // sugar / modified-residue classification via gemmi's tabulated residues.
      bool isSugar() {
         gemmi::ResidueKind k = gemmi::find_tabulated_residue(g().name).kind;
         return k == gemmi::ResidueKind::PYR || k == gemmi::ResidueKind::KET;
      }
      // MMDB isModRes reflects PDB MODRES records (a non-standard, modified form of a
      // standard residue). gemmi has no per-residue MODRES flag on the model tree, so
      // approximate: an amino/nucleic residue whose one-letter code is lower-case
      // (gemmi marks non-standard monomers that way). Water/ligands are excluded.
      bool isModRes() {
         const gemmi::ResidueInfo ri = gemmi::find_tabulated_residue(g().name);
         return ri.found() && !ri.is_standard() &&
                (ri.is_amino_acid() || ri.is_nucleic_acid());
      }

      Residue() = default;
      explicit Residue(Chain *c);  // construct + add to chain (out-of-line)
      // MMDB owns residues via `new`/`delete` (Coot writes `delete residue_p;`). Like
      // ~Atom, this frees the residue's atoms, then DEFERS structural removal: it nulls
      // this residue's slot in the parent chain (a tombstone) but leaves the gemmi
      // placeholder residue in place so siblings' `ri` stays valid — _compact_residues
      // drops both later. No-op during Manager teardown (_bulk_free).
      ~Residue();  // out-of-line (needs complete Manager/Chain)

      // MMDB public char-array fields. Coot reads `residue->name` and writes
      // `strncpy(residue->insCode,..)`. Kept as the interface: synced gemmi->buffer on
      // load (_load_id, in build_from_gemmi) and buffer->gemmi at the adopt point
      // (_store_id, in Chain::Add/InsResidue). SetResName/SetResID keep both in step.
      ResName name{};
      InsCode insCode{};
      void _load_id() {
         std::snprintf(name, sizeof name, "%s", g().name.c_str());
         insCode[0] = g().seqid.icode && g().seqid.icode != ' ' ? g().seqid.icode : '\0';
         insCode[1] = '\0';
      }
      void _store_id() {
         g().name = name;
         g().seqid.icode = insCode[0] ? insCode[0] : ' ';
      }

      gemmi::Residue &g() const;

      int GetNumberOfAtoms() { return (int)atoms.size(); }
      int GetNumberOfAtoms(bool /*countTers*/) { return (int)atoms.size(); }
      PAtom GetAtom(int atomNo) {
         return (atomNo >= 0 && atomNo < (int)atoms.size()) ? atoms[atomNo] : nullptr;
      }
      PAtom GetAtom(const AtomName aname, const Element elname = nullptr,
                    const AltLoc aloc = nullptr);
      void GetAtomTable(PPAtom &atomTable, int &n) {
         atomTable = atoms.data();
         n = (int)atoms.size();
      }
      PAtom AddAtom(Manager &m, gemmi::Atom a);  // append: O(1)
      // Adopt a detached atom (Coot's `new mmdb::Atom` idiom). Copies the atom's
      // local gemmi into this residue's gemmi (detached or bound, via g()) and
      // rebinds the wrapper. Pushes the strncpy'd altLoc buffer back into gemmi.
      int AddAtom(PAtom atm);  // out-of-line: needs complete Manager (_atom_allocs)
      void DeleteAtom(int pos);
      void _detach_atom(Atom *a);  // unlink (no free); used by ~Atom on Coot `delete atom`
      void TrimAtomTable() {}  // compact after deletions — shim keeps them in sync

      pstr GetResName();
      void SetResName(const ResName n) {
         g().name = n ? n : "";
         std::snprintf(name, sizeof name, "%s", n ? n : "");
      }
      void SetResID(const ResName resName, int seqNo, const InsCode ic) {
         g().name = resName ? resName : "";
         g().seqid.num.value = seqNo;
         g().seqid.icode = (ic && ic[0]) ? ic[0] : ' ';
         std::snprintf(name, sizeof name, "%s", resName ? resName : "");
         insCode[0] = (ic && ic[0]) ? ic[0] : '\0';
         insCode[1] = '\0';
      }
      int &GetSeqNum();  // writable (rewrite maps `->seqNum` reads and writes)
      pstr GetInsCode();
      pstr GetChainID();
      int GetModelNum();
      int &GetIndex() { return ri; }  // ref: rewritten `->index` is assignable
      Chain *GetChain() { return chain; }
      Model *GetModel();  // out-of-line (Chain incomplete here)
      // terminus tests — peptide-bond-aware: N-terminus if no preceding residue's C is
      // within bonding distance of this N, C-terminus if this C bonds no following N
      // (out-of-line: need Chain + backbone atom geometry).
      bool isNTerminus();
      bool isCTerminus();
      pstr GetResidueID(pstr S) {  // "seqnum(name):inscode"
         if (S) std::snprintf(S, 100, "%d(%s):%s", GetSeqNum(), name, insCode);
         return S;
      }
      Residue *next = nullptr;  // MMDB has this; wired lazily if needed
      int SSE = SSE_None;       // secondary-structure element (shim-owned public field)
      bool isAminoacid() { return gemmi::find_tabulated_residue(g().name).is_amino_acid(); }
      bool isNucleotide() { return gemmi::find_tabulated_residue(g().name).is_nucleic_acid(); }
      bool isDNARNA() { return isNucleotide(); }
      bool isSolvent() { return gemmi::find_tabulated_residue(g().name).is_water(); }
      // UDData
      int PutUDData(int h, int v) { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
      int PutUDData(int h, realtype v) { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
      int PutUDData(int h, cpstr v) { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
      int GetUDData(int h, int &v) { return ud_get(mgr, UDR_RESIDUE, *this, h, v); }
      int GetUDData(int h, realtype &v) { return ud_get(mgr, UDR_RESIDUE, *this, h, v); }

   private:
      ResName _resname_buf{};
      InsCode _inscode_buf{};
   };

   // ===========================================================================
   class Chain : public UDStore {
   public:
      Manager *mgr = nullptr;
      Model *model = nullptr;  // parent; null => detached (use _local)
      int ci = 0;
      bool alive = true;
      gemmi::Chain _local;  // backing store while detached
      std::vector<Residue *> residues;

      gemmi::Chain &g() const;

      int GetNumberOfResidues() { return (int)residues.size(); }
      PResidue GetResidue(int resNo) {
         return (resNo >= 0 && resNo < (int)residues.size()) ? residues[resNo] : nullptr;
      }
      // find by (seqNum, insCode) — MMDB's 2-arg overload. Skips deferred-delete
      // tombstone (null) slots.
      PResidue GetResidue(int seqNum, const InsCode insCode) {
         char ic = (insCode && insCode[0]) ? insCode[0] : ' ';
         for (Residue *r : residues) {
            if (!r) continue;  // tombstone
            gemmi::Residue &gr = r->g();
            char ric = gr.seqid.icode ? gr.seqid.icode : ' ';
            if (gr.seqid.num.value == seqNum && ric == ic) return r;
         }
         return nullptr;
      }
      void GetResidueTable(PPResidue &t, int &n) {
         t = residues.data();
         n = (int)residues.size();
      }
      // MMDB DeleteResidue is DEFERRED: it frees the residue and leaves a NULL slot
      // (nResidues unchanged) until TrimResidueTable/FinishStructEdit. `delete residue`
      // (via ~Residue) does the same. Coot relies on this (e.g. change_chain_id iterates
      // to the original count while deleting). We keep the gemmi placeholder residue too
      // so surviving residues' `ri` stays aligned; _compact_residues() drops both later.
      void DeleteResidue(int resNo) {
         if (resNo < 0 || resNo >= (int)residues.size()) return;
         if (residues[resNo]) delete residues[resNo];  // ~Residue nulls the slot
      }
      void DeleteResidue(int seqNum, const InsCode ic) {  // by (seqNum, insCode)
         PResidue r = GetResidue(seqNum, ic);
         if (r) delete r;  // ~Residue nulls its slot; keeps the gemmi placeholder
      }
      void TrimResidueTable() { _compact_residues(); }
      // Drop tombstoned (null) residue slots and their gemmi placeholder residues in
      // lock-step, then reindex ri. Out-of-line: needs a complete Residue.
      void _compact_residues();
      pstr GetChainID();
      pstr GetChainID(pstr buf) {
         if (buf) std::snprintf(buf, sizeof(ChainID), "%s", g().name.c_str());
         return buf;
      }
      Manager *GetCoordHierarchy() { return mgr; }  // parent manager
      void SetChainID(const ChainID id) { g().name = id ? id : ""; }
      Chain() = default;
      Chain(Model *m, const ChainID id);  // construct + add to model (out-of-line)
      void Copy(PChain src);              // deep-copy subtree (out-of-line: needs Manager)
      // Reorder residues (and their gemmi backing) ascending by (seqNum, insCode),
      // MMDB's default. Keeps the wrapper vector and gemmi vector in lock-step and
      // re-indexes ri. sortKey variants beyond ascending-by-number are uncommon in
      // Coot and treated as the default.
      void SortResidues(int /*sortKey*/ = 0) {
         _compact_residues();  // never sort across deferred-delete tombstones
         int n = (int)residues.size();
         if (n < 2) return;
         std::vector<int> ord(n);
         for (int i = 0; i < n; ++i) ord[i] = i;
         gemmi::Chain &gc = g();
         std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
            const gemmi::Residue &ra = gc.residues[a], &rb = gc.residues[b];
            if (ra.seqid.num.value != rb.seqid.num.value) return ra.seqid.num.value < rb.seqid.num.value;
            char ia = ra.seqid.icode ? ra.seqid.icode : ' ', ib = rb.seqid.icode ? rb.seqid.icode : ' ';
            return ia < ib;
         });
         std::vector<gemmi::Residue> gnew;
         gnew.reserve(n);
         std::vector<Residue *> wnew;
         wnew.reserve(n);
         for (int k = 0; k < n; ++k) {
            gnew.push_back(std::move(gc.residues[ord[k]]));
            wnew.push_back(residues[ord[k]]);
         }
         gc.residues = std::move(gnew);
         residues = std::move(wnew);
         for (int k = 0; k < n; ++k) residues[k]->ri = k;
      }
      bool isAminoacidChain();  // defined out-of-line (needs Residue predicates)
      bool isNucleotideChain();
      bool isSolventChain();
      PResidue AddResidue(Manager &m, gemmi::Residue r);  // append
      PResidue InsResidue(Manager &m, int pos, gemmi::Residue r);
      // Adopt a detached residue (its atom wrappers already point at it, so they
      // ride along once its gemmi is copied in and the wrapper is rebound).
      int AddResidue(PResidue res) {
         res->_store_id();                  // push name/insCode buffers into gemmi
         g().residues.push_back(res->g());  // res detached -> its _local (with atoms)
         res->chain = this;
         res->mgr = mgr;
         res->ri = (int)residues.size();
         residues.push_back(res);
         return 0;
      }
      int InsResidue(PResidue res, int pos) {
         _compact_residues();  // don't insert/reindex across tombstones
         if (pos < 0) pos = 0;
         if (pos > (int)residues.size()) pos = (int)residues.size();
         res->_store_id();
         g().residues.insert(g().residues.begin() + pos, res->g());
         res->chain = this;
         res->mgr = mgr;
         res->ri = pos;
         residues.insert(residues.begin() + pos, res);
         for (int k = pos + 1; k < (int)residues.size(); ++k) residues[k]->ri = k;
         return 0;
      }

   private:
      ChainID _chainid_buf{};
   };

   // ===========================================================================
   class Model : public UDStore {
   public:
      Manager *mgr = nullptr;  // null => detached (use _local)
      int mi = 0;              // 0-based internal; GetModel is 1-based externally
      gemmi::Model _local{1};  // backing store while detached (gemmi Model num is int)
      std::vector<Chain *> chains;

      gemmi::Model &g() const;

      int GetNumberOfChains() { return (int)chains.size(); }
      PChain GetChain(int chainNo) {
         return (chainNo >= 0 && chainNo < (int)chains.size()) ? chains[chainNo] : nullptr;
      }
      PChain GetChain(const ChainID chID);
      // Adopt a detached chain (Coot's `new mmdb::Chain` idiom): copy its local
      // gemmi (with any residues/atoms) into this model and rebind, cascading mgr
      // to the sub-tree that was built while detached (mgr was null).
      int AddChain(PChain chn) {
         g().chains.push_back(chn->g());
         chn->model = this;
         chn->mgr = mgr;
         chn->ci = (int)chains.size();
         chains.push_back(chn);
         for (Residue *r : chn->residues) {
            r->mgr = mgr;
            for (Atom *a : r->atoms) a->mgr = mgr;
         }
         return 0;
      }
      int GetSerNum() { return mi + 1; }
      // delete chain at index: erase gemmi + wrapper, reindex the tail
      void DeleteChain(int chainNo) {
         if (chainNo < 0 || chainNo >= (int)chains.size()) return;
         g().chains.erase(g().chains.begin() + chainNo);
         chains.erase(chains.begin() + chainNo);
         for (int k = chainNo; k < (int)chains.size(); ++k) chains[k]->ci = k;
      }
      void DeleteChain(const ChainID chainID) {
         for (int i = 0; i < (int)chains.size(); ++i)
            if (chains[i]->g().name == (chainID ? chainID : "")) {
               DeleteChain(i);
               return;
            }
      }
      void GetChainTable(PPChain &t, int &n) {
         t = chains.data();
         n = (int)chains.size();
      }
      std::vector<Atom *> all_atoms;  // flat, filled by build_from_gemmi
      PPAtom GetAllAtoms() { return all_atoms.data(); }
      int GetNumberOfAtoms() { return (int)all_atoms.size(); }
      int GetNumberOfAtoms(bool /*countTers*/) { return (int)all_atoms.size(); }
      // Secondary-structure assignment: mocked. gemmi's DSSP has its SS prediction
      // disabled upstream ("commented out ... wasn't correct anyway"), so there is no
      // gemmi-backed SS to forward to. Return the non-OK code so callers treat SS as
      // unavailable rather than trusting a bogus assignment. (residue SSE stays None.)
      int CalcSecStructure(bool /*flag*/) { return SSERC_noResidues; }
      // LINK records — gemmi-loaded (Manager::_load_metadata) plus Coot-created ones
      // (AddLink) stored here; GetLink is 1-based like MMDB.
      std::vector<Link *> _links;
      int GetNumberOfLinks() { return (int)_links.size(); }
      PLink GetLink(int i) { return (i >= 1 && i <= (int)_links.size()) ? _links[i - 1] : nullptr; }
      void AddLink(PLink link) {
         if (link) _links.push_back(link);
      }
      // Refmac LINKR records — gemmi Connections that carry a link_id (_load_metadata).
      std::vector<LinkR *> _linkrs;
      int GetNumberOfLinkRs() { return (int)_linkrs.size(); }
      PLinkR GetLinkR(int i) { return (i >= 1 && i <= (int)_linkrs.size()) ? _linkrs[i - 1] : nullptr; }
      void AddLinkR(PLinkR lr) {
         if (lr) _linkrs.push_back(lr);
      }
      std::vector<CisPep *> _cispeps;
      int GetNumberOfCisPeps() { return (int)_cispeps.size(); }
      PCisPep GetCisPep(int i) { return (i >= 1 && i <= (int)_cispeps.size()) ? _cispeps[i - 1] : nullptr; }
      void AddCisPep(PCisPep cp) {
         if (cp) _cispeps.push_back(cp);
      }
      void RemoveCisPeps() { _cispeps.clear(); }
      // secondary structure. Records live in `helices`/`sheets` below, populated
      // either from gemmi on load (build_from_gemmi) or by Coot's own SS computation
      // via the access_model subclass (which reaches these public members directly).
      // 1-based indexing to match MMDB.
      int GetNumberOfHelices() { return (int)helices.data.size(); }
      PHelix GetHelix(int i) { return (i >= 1 && i <= (int)helices.data.size()) ? helices.data[i - 1] : nullptr; }
      int GetNumberOfSheets() { return sheets.nSheets; }
      PSheet GetSheet(int i) { return (i >= 1 && i <= sheets.nSheets && sheets.sheet) ? sheets.sheet[i - 1] : nullptr; }
      Sheets sheets;                    // SS records (gemmi-backed on load; access_model fills)
      Helices helices;                  // "     "     "
      std::vector<PSheet> _sheet_ptrs;  // backing array for sheets.sheet (gemmi load)
      PSheets GetSheets() { return &sheets; }
      int GetModelID() { return mi + 1; }
      pstr GetModelID(pstr buf) {
         if (buf) std::snprintf(buf, 16, "%d", mi + 1);
         return buf;
      }
      int CalcSecStructure(int /*flag*/, int /*selHnd*/) { return SSERC_noResidues; }  // mocked; see bool overload
      void Copy(PModel src);                                                           // deep-copy subtree (out-of-line)
      Manager *GetCoordHierarchy() { return mgr; }                                     // parent manager
      int GetNumberOfResidues() {
         int n = 0;
         for (Chain *c : chains) n += c->GetNumberOfResidues();
         return n;
      }
      LinkContainer _linkc;
      PLinkContainer GetLinks() {
         _linkc.data.assign(_links.begin(), _links.end());
         return &_linkc;
      }
      void RemoveLinks() { _links.clear(); }
      // Reorder chains (and gemmi backing) by chain ID. sortKey selects ascending
      // (default) or descending; other MMDB sort keys collapse to ID order.
      void SortChains(int sortKey = 0) {
         int n = (int)chains.size();
         if (n < 2) return;
         bool desc = (sortKey == SORT_CHAIN_ChainID_Desc);
         std::vector<int> ord(n);
         for (int i = 0; i < n; ++i) ord[i] = i;
         gemmi::Model &gm = g();
         std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
            return desc ? (gm.chains[a].name > gm.chains[b].name)
                        : (gm.chains[a].name < gm.chains[b].name);
         });
         std::vector<gemmi::Chain> gnew;
         gnew.reserve(n);
         std::vector<Chain *> wnew;
         wnew.reserve(n);
         for (int k = 0; k < n; ++k) {
            gnew.push_back(std::move(gm.chains[ord[k]]));
            wnew.push_back(chains[ord[k]]);
         }
         gm.chains = std::move(gnew);
         chains = std::move(wnew);
         for (int k = 0; k < n; ++k) chains[k]->ci = k;
      }
      PChain CreateChain(const ChainID id);  // add empty chain (out-of-line: needs Manager)
      int GetNumberOfStrands(int sheetNo) {
         PSheet s = GetSheet(sheetNo);
         return s ? s->nStrands : 0;
      }
      PStrand GetStrand(int sheetNo, int strandNo) {
         PSheet s = GetSheet(sheetNo);
         return (s && strandNo >= 1 && strandNo <= s->nStrands && s->strand) ? s->strand[strandNo - 1] : nullptr;
      }
   };

}  // namespace mmdb
