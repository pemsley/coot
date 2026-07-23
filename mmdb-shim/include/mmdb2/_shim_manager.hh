// mmdb-shim — layer 3 of 4: the Manager (MMDB Root/CoorManager/SelManager rolled
// into one).
//
// Manager owns the live gemmi::Structure plus the wrapper trees built over it, the
// individually heap-allocated Atom/Residue nodes (so Coot's `delete atom;` idiom
// works), the stable-address pools for gemmi-derived metadata records, the
// handle-based selection engine, and the UDData registry. I/O, selection matching,
// PutAtom and the gemmi metadata mapping are declared here and defined out-of-line
// in _shim_inline.hh / src/*.cc.
#pragma once

#include "_shim_hierarchy.hh"

namespace mmdb {

   // ===========================================================================
   class Manager {
   public:
      gemmi::Structure st;
      // Atoms are heap-allocated individually (not pooled) so Coot's MMDB idiom
      // `delete atom;` frees exactly one node. The manager owns every atom it hands
      // out and frees the survivors at teardown; `~Atom` removes itself from this set
      // when Coot deletes it early. `_bulk_free` tells `~Atom` to skip detach work
      // while the manager is tearing everything down.
      std::set<Atom *> _atom_allocs;
      std::set<Residue *> _res_allocs;  // residues heap-allocated too (Coot `delete residue_p`)
      bool _bulk_free = false;
      ~Manager() {
         _bulk_free = true;
         for (Atom *a : _atom_allocs) delete a;
         _atom_allocs.clear();
         for (Residue *r : _res_allocs) delete r;
         _res_allocs.clear();
      }
      // stable-address pools (Chain/Model still pooled — see _atom_allocs / _res_allocs note)
      std::deque<Chain> chain_pool;
      std::deque<Model> model_pool;
      std::vector<Model *> models;
      // stable-address pools for gemmi-derived metadata records (LINK / CISPEP /
      // HELIX / SHEET). Filled by build_from_gemmi -> _load_metadata(); owned here so
      // the Model containers can hold bare pointers into them.
      std::deque<Link> link_pool;
      std::deque<LinkR> linkr_pool;
      std::deque<CisPep> cispep_pool;
      std::deque<Helix> helix_pool;
      std::deque<Sheet> sheet_pool;
      std::deque<Strand> strand_pool;
      std::deque<std::vector<PStrand>> strandarr_pool;  // backing for Sheet::strand (Strand**)
      std::deque<Author> author_pool;                   // backing for title.author records
      void _load_metadata();                            // out-of-line: needs complete gemmi metadata types

      Atom *newAtom() {
         Atom *a = new Atom();
         a->mgr = this;
         _atom_allocs.insert(a);
         return a;
      }
      Residue *newRes() {
         Residue *r = new Residue();
         r->mgr = this;
         _res_allocs.insert(r);
         return r;
      }
      Chain *newChain() {
         chain_pool.emplace_back();
         return &chain_pool.back();
      }
      Model *newModel() {
         model_pool.emplace_back();
         return &model_pool.back();
      }

      int GetNumberOfModels() { return (int)models.size(); }
      PModel GetModel(int modelNo) {  // MMDB: 1 <= modelNo <= nModels
         int i = modelNo - 1;
         return (i >= 0 && i < (int)models.size()) ? models[i] : nullptr;
      }
      // per-model chain access (modelNo is 1-based, chainNo 0-based) — mmdb_coormngr.h
      int GetNumberOfChains(int modelNo) {
         PModel m = GetModel(modelNo);
         return m ? m->GetNumberOfChains() : 0;
      }
      PChain GetChain(int modelNo, int chainNo) {
         PModel m = GetModel(modelNo);
         return m ? m->GetChain(chainNo) : nullptr;
      }
      // Re-index/renumber after edits. Sibling indices are kept in sync as the shim
      // mutates (so PDBCLEAN_INDEX is implicit); PDBCLEAN_SERIAL renumbers atom serials
      // 1..N in hierarchy order. Other clean flags are not needed by the shim.
      word PDBCleanup(word CleanKey) {
         // INDEX cleanup compacts deferred residue deletions (drops tombstones) — do it
         // before renumbering so serials/indices count only surviving atoms.
         if (CleanKey & PDBCLEAN_INDEX) {
            for (Model *mw : models)
               for (Chain *cw : mw->chains) cw->_compact_residues();
            _rebuild_all_atoms();
         }
         if (CleanKey & (PDBCLEAN_SERIAL | PDBCLEAN_INDEX)) {
            int s = 1;
            for (Atom *a : all_atoms) a->g().serial = s++;
         }
         return 0;
      }

      // PDB title records — Coot reaches `title` via an access_mol subclass; the
      // TITLE string comes from gemmi (_struct.title), authors are filled on load.
      Title title;
      pstr GetStructureTitle(pstr T) {
         if (T) std::strcpy(T, st.get_info("_struct.title").c_str());  // caller allocates (MMDB contract)
         return T;
      }

      // Orthogonal symmetry transformation for operator Nop (0-based) + cell shifts,
      // via gemmi's space group + unit cell (shared helper). Returns 0 on success,
      // nonzero if there is no usable space group / the operator is out of range.
      int GetTMatrix(mat44 &TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) {
         return gemmi_sym_tmatrix(st.cell, st.spacegroup_hm, TMatrix, Nop,
                                  cellshift_a, cellshift_b, cellshift_c);
      }

      void build_from_gemmi();

      // adopt a detached model (Coot: `new mmdb::Model` -> AddChain… -> AddModel).
      // Copy its local gemmi into st, rebind, cascade mgr through the sub-tree.
      int AddModel(PModel mw) {
         st.models.push_back(mw->g());
         mw->mgr = this;
         mw->mi = (int)models.size();
         models.push_back(mw);
         for (Chain *cw : mw->chains) {
            cw->mgr = this;
            for (Residue *rw : cw->residues) {
               rw->mgr = this;
               for (Atom *aw : rw->atoms) {
                  aw->mgr = this;
                  all_atoms.push_back(aw);
                  mw->all_atoms.push_back(aw);
               }
            }
         }
         return 0;
      }

      // clone another manager's structure (mmdb Manager::Copy(PManager, COPY_MASK)).
      // Copies the whole gemmi Structure and rebuilds all wrappers — clean & correct.
      void Copy(PManager m, int /*CopyMask*/) {
         if (m) {
            st = m->st;
            build_from_gemmi();
         }
      }

      // ---- crystal cell & symmetry (gemmi UnitCell / SpaceGroup) ----
      std::string _sg_buf, _symop_buf;
      void GetCell(realtype &a, realtype &b, realtype &c, realtype &al, realtype &be,
                   realtype &ga, realtype &vol, int &orthcode) {
         const gemmi::UnitCell &u = st.cell;
         a = u.a;
         b = u.b;
         c = u.c;
         al = u.alpha;
         be = u.beta;
         ga = u.gamma;
         vol = u.volume;
         orthcode = 1;
      }
      void GetCell(realtype &a, realtype &b, realtype &c, realtype &al, realtype &be,
                   realtype &ga, realtype &vol) {
         int oc;
         GetCell(a, b, c, al, be, ga, vol, oc);
      }
      void SetCell(realtype a, realtype b, realtype c, realtype al, realtype be,
                   realtype ga, int /*OrthCode*/ = 1) { st.cell.set(a, b, c, al, be, ga); }
      void Orth2Frac(realtype x, realtype y, realtype z, realtype &u, realtype &v, realtype &w) {
         gemmi::Fractional f = st.cell.fractionalize(gemmi::Position(x, y, z));
         u = f.x;
         v = f.y;
         w = f.z;
      }
      void Frac2Orth(realtype u, realtype v, realtype w, realtype &x, realtype &y, realtype &z) {
         gemmi::Position p = st.cell.orthogonalize(gemmi::Fractional(u, v, w));
         x = p.x;
         y = p.y;
         z = p.z;
      }
      pstr GetSpaceGroup() {
         _sg_buf = st.spacegroup_hm;
         return (pstr)_sg_buf.c_str();
      }
      pstr GetSpaceGroupFix() { return GetSpaceGroup(); }
      int SetSpaceGroup(cpstr sg) {
         st.spacegroup_hm = sg ? sg : "";
         return 0;
      }
      int GetNumberOfSymOps() {
         const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(st.spacegroup_hm);
         return sg ? (int)sg->operations().order() : 0;
      }
      pstr GetSymOp(int Nop) {
         const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(st.spacegroup_hm);
         if (!sg) return nullptr;
         int i = 0;
         for (gemmi::Op op : sg->operations()) {
            if (i++ == Nop) {
               _symop_buf = op.triplet();
               return (pstr)_symop_buf.c_str();
            }
         }
         return nullptr;
      }

      // ---- selection ----
      struct Selection {
         SELECTION_TYPE type = STYPE_UNDEFINED;
         std::vector<Atom *> atoms;
         std::vector<Residue *> residues;
         std::vector<Chain *> chains;
      };
      std::vector<Selection> selections;  // handle is 1-based index

      int NewSelection() {
         selections.emplace_back();
         return (int)selections.size();
      }
      void DeleteSelection(int selHnd) {
         if (selHnd < 1 || selHnd > (int)selections.size()) return;
         Selection &s = selections[selHnd - 1];
         for (Atom *a : s.atoms) a->_setInSel(selHnd, false);
         for (Residue *r : s.residues) r->_setInSel(selHnd, false);
         for (Chain *c : s.chains) c->_setInSel(selHnd, false);
         s = Selection();
      }
      void GetSelIndex(int selHnd, PPAtom &SelAtom, int &n) {
         Selection &s = selections[selHnd - 1];
         SelAtom = s.atoms.data();
         n = (int)s.atoms.size();
      }
      void GetSelIndex(int selHnd, PPResidue &SelRes, int &n) {
         Selection &s = selections[selHnd - 1];
         SelRes = s.residues.data();
         n = (int)s.residues.size();
      }
      void GetSelIndex(int selHnd, PPChain &SelChain, int &n) {
         Selection &s = selections[selHnd - 1];
         SelChain = s.chains.data();
         n = (int)s.chains.size();
      }
      // select atoms by serial-number range (iSer1..iSer2; 0,0 => all).
      void SelectAtoms(int selHnd, int iSer1, int iSer2, SELECTION_KEY key) {
         if (selHnd < 1 || selHnd > (int)selections.size()) return;
         Selection &s = selections[selHnd - 1];
         std::vector<Atom *> pick;
         for (Atom *a : all_atoms) {
            int sn = a->g().serial;
            if ((iSer1 == 0 && iSer2 == 0) || (sn >= iSer1 && sn <= iSer2)) pick.push_back(a);
         }
         if (key == SKEY_OR) {
            for (Atom *a : pick)
               if (!a->isInSelection(selHnd)) s.atoms.push_back(a);
         } else {
            for (Atom *a : s.atoms) a->_setInSel(selHnd, false);
            s.atoms = pick;
         }
         s.type = STYPE_ATOM;
         for (Atom *a : s.atoms) a->_setInSel(selHnd, true);
      }

      // full spatial+CID atom selection (mmdb_selmngr.h) — sphere around (x,y,z) with
      // chain/resname/atomname/element filters ("!X" = exclusion, "*" = any).
      void SelectAtoms(int selHnd, int /*iModel*/, cpstr Chains, int ResNo1, cpstr /*Ins1*/,
                       int ResNo2, cpstr /*Ins2*/, cpstr RNames, cpstr ANames, cpstr Elements,
                       cpstr /*altLocs*/, cpstr /*segIDs*/, cpstr /*charges*/,
                       realtype /*occ1*/, realtype /*occ2*/, realtype x, realtype y, realtype z,
                       realtype radius, SELECTION_KEY key) {
         if (selHnd < 1 || selHnd > (int)selections.size()) return;
         Selection &s = selections[selHnd - 1];
         // self-contained comma-list matcher ("*"=any, "!X"=exclude); `detail::` is
         // declared after Manager, so don't depend on it in this inline body.
         auto trimws = [](const std::string &s) -> std::string {
            size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
            return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
         };
         auto inlist = [&trimws](cpstr list, const std::string &v) -> bool {
            if (!list || !*list || std::strcmp(list, "*") == 0) return true;
            std::string vt = trimws(v);
            for (const char *p = list; *p;) {
               const char *c = std::strchr(p, ',');
               std::string tok(p, c ? (size_t)(c - p) : std::strlen(p));
               if (trimws(tok) == vt) return true;
               if (!c) break;
               p = c + 1;
            }
            return false;
         };
         auto match = [&](cpstr list, const std::string &v) -> bool {
            if (!list || !*list || std::strcmp(list, "*") == 0) return true;
            if (list[0] == '!') return !inlist(list + 1, v);
            return inlist(list, v);
         };
         gemmi::Position pt(x, y, z);
         double r2 = radius * radius;
         std::vector<Atom *> pick;
         for (Atom *a : all_atoms) {
            if (radius > 0 && a->g().pos.dist_sq(pt) > r2) continue;
            Residue *r = a->res;
            int sn = r->GetSeqNum();
            if (ResNo1 != ANY_RES && sn < ResNo1) continue;
            if (ResNo2 != ANY_RES && sn > ResNo2) continue;
            if (!match(Chains, r->chain->g().name)) continue;
            if (!match(RNames, std::string(r->GetResName()))) continue;
            if (!match(ANames, std::string(a->GetAtomName()))) continue;
            if (!match(Elements, gemmi::Element(a->g().element).name())) continue;
            pick.push_back(a);
         }
         if (key == SKEY_OR) {
            for (Atom *a : pick)
               if (!a->isInSelection(selHnd)) s.atoms.push_back(a);
         } else {
            for (Atom *a : s.atoms) a->_setInSel(selHnd, false);
            s.atoms = pick;
         }
         s.type = STYPE_ATOM;
         for (Atom *a : s.atoms) a->_setInSel(selHnd, true);
      }

      // --- misc hierarchy/bond/UDData ops used by Coot ---
      void RemoveBonds() {}  // gemmi has no persistent bond table
      // Partial-hierarchy delete (mmdb Manager::Delete). Coot's use is
      // Delete(MMDBFCM_SC) to drop secondary-structure/connectivity records before
      // writing; also honour Coord (atoms) and Cryst (cell/SG) for completeness.
      void Delete(int DelKey) {
         bool all = DelKey == MMDBFCM_All;
         if (all || (DelKey & MMDBFCM_SC)) {
            for (Model *m : models) {
               m->_links.clear();
               m->_linkrs.clear();
               m->_cispeps.clear();
               m->helices.data.clear();
               m->sheets.nSheets = 0;
               m->sheets.sheet = nullptr;
               m->_sheet_ptrs.clear();
            }
            link_pool.clear();
            linkr_pool.clear();
            cispep_pool.clear();
            helix_pool.clear();
            sheet_pool.clear();
            strand_pool.clear();
            strandarr_pool.clear();
         }
         if (all || (DelKey & MMDBFCM_Cryst)) {
            st.cell = gemmi::UnitCell();
            st.spacegroup_hm.clear();
         }
         if (all || (DelKey & MMDBFCM_Coord)) {
            st.models.clear();
            build_from_gemmi();
         }
      }
      void DeleteAllModels() {
         st.models.clear();
         build_from_gemmi();
      }  // clears the hierarchy
      void DeleteModel(int modelNo) {  // 1-based; erase model + rebuild wrappers
         int i = modelNo - 1;
         if (i >= 0 && i < (int)st.models.size()) {
            st.models.erase(st.models.begin() + i);
            build_from_gemmi();
         }
      }
      pstr GetInputBuffer(pstr buf, int &count) {
         count = 0;
         if (buf) buf[0] = '\0';
         return buf;
      }
      // Insert (a copy of) an atom into the hierarchy (mmdb Manager::PutAtom). MMDB
      // keeps a flat atom array with a parallel hierarchy rebuilt by FinishStructEdit;
      // the shim's storage IS the hierarchy, so PutAtom finds/creates the chain and
      // residue implied by the atom's source residue and appends a copy there. Only
      // append (index<=0 or top) is supported — the semantics Coot relies on
      // (create_mmdbmanager_from_atom_selection_straight). Returns the atom's 1-based
      // position (so GetAtomI(pos) returns it). Defined out-of-line (needs Add*).
      int PutAtom(int index, PAtom atom, int serNum = 0);
      // hierarchy-level UDData (UDR_HIERARCHY) — Manager owns its own UDStore.
      UDStore _ud;
      int PutUDData(int h, int v) { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
      int PutUDData(int h, realtype v) { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
      int PutUDData(int h, cpstr v) { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
      int GetUDData(int h, int &v) { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
      int GetUDData(int h, realtype &v) { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
      int GetUDData(int h, pstr &v) { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
      // primary CID-range selection (STYPE via Select; SelectAtoms forwards as STYPE_ATOM)
      void Select(int selHnd, SELECTION_TYPE sType, int iModel, cpstr Chains,
                  int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames,
                  cpstr ANames, cpstr Elements, cpstr altLocs, SELECTION_KEY selKey = SKEY_OR);
      void SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1,
                       int ResNo2, cpstr Ins2, cpstr RNames, cpstr ANames,
                       cpstr Elements, cpstr altLocs, SELECTION_KEY selKey = SKEY_OR) {
         Select(selHnd, STYPE_ATOM, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2,
                RNames, ANames, Elements, altLocs, selKey);
      }
      void SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x, realtype y,
                        realtype z, realtype r, SELECTION_KEY sKey = SKEY_OR);
      // select-from-selection: combine selHnd2's contents into selHnd1 per sKey
      void Select(int selHnd1, SELECTION_TYPE sType, int selHnd2, SELECTION_KEY sKey);
      // atoms within [d1,d2] of any atom in the given set (defined in contacts.cc)
      void SelectNeighbours(int selHnd, SELECTION_TYPE sType, PPAtom atoms, int nAtoms,
                            realtype d1, realtype d2, SELECTION_KEY sKey = SKEY_OR);
      void SetFlag(int /*flags*/) {}  // no-op: read/write behaviour is fixed
      void SetFlag(cpstr /*flags*/) {}
      int PutPDBString(cpstr /*card*/) { return Error_NoError; }  // no-op
      // No persistent bond table. Verified safe: Coot's only caller (make_bonds in
      // coot-utils/bonded-atoms.cc) ignores the mmdb bond table and recomputes bonds
      // itself from geometry, so a no-op here matches observed Coot behaviour.
      int MakeBonds(bool /*calc*/) { return 0; }

      // flat atom access (across the whole hierarchy)
      std::vector<Atom *> all_atoms;
      int GetNumberOfAtoms() { return (int)all_atoms.size(); }
      int GetNumberOfAtoms(bool /*countTers*/) { return (int)all_atoms.size(); }
      int GetNumberOfAtoms(cpstr CID);  // count atoms matching CID (defined below)
      // MMDB GetAtomI is 1-based: returns Atom[index-1].
      PAtom GetAtomI(int i) { return (i >= 1 && i <= (int)all_atoms.size()) ? all_atoms[i - 1] : nullptr; }
      void GetAtomTable(PPAtom &t, int &n) {
         t = all_atoms.data();
         n = (int)all_atoms.size();
      }
      void GetModelTable(PPModel &t, int &n) {
         t = models.data();
         n = (int)models.size();
      }
      void GetAtomStatistics(int selHnd, RAtomStat AS);  // defined below
      int MakeSelIndex(int selHnd) {
         return (selHnd >= 1 && selHnd <= (int)selections.size())
                    ? (int)selections[selHnd - 1].atoms.size()
                    : 0;
      }
      void SelectAtom(int selHnd, PAtom atom, SELECTION_KEY sKey, bool makeIndex = true);
      // CID-string selection, e.g. "/1/A/10-20/CA"
      void Select(int selHnd, SELECTION_TYPE sType, cpstr CID, SELECTION_KEY sKey);

      // ---- contacts (gemmi NeighborSearch; TMatrix path uses a uniform grid) ----
      // TMatrix is MMDB's optional symmetry transform applied to the 2nd set: when
      // given, contacts.cc transforms that set and searches against it (symmetry
      // mates); when null, gemmi NeighborSearch over the untransformed model is used.
      void SeekContacts(PPAtom A1, int n1, PPAtom A2, int n2, realtype d1,
                        realtype d2, int seqDist, PContact &contact, int &ncontacts,
                        int maxlen = 0, pmat44 TMatrix = nullptr, long group = 0);
      void SeekContacts(PPAtom A, int n, realtype d1, realtype d2, int seqDist,
                        PContact &contact, int &ncontacts, int maxlen = 0,
                        pmat44 TMatrix = nullptr, long group = 0);
      // single-atom vs selection (forwards to the array overload with a 1-elem array)
      void SeekContacts(PAtom a, PPAtom A2, int n2, realtype d1, realtype d2, int seqDist,
                        PContact &contact, int &ncontacts, int maxlen = 0,
                        pmat44 TMatrix = nullptr, long group = 0) {
         PAtom a1[1] = {a};
         SeekContacts(a1, 1, A2, n2, d1, d2, seqDist, contact, ncontacts, maxlen, TMatrix, group);
      }

      // Compact deferred residue deletions across the whole hierarchy: MMDB defers
      // DeleteResidue (tombstone the slot, keep the count) until FinishStructEdit, so
      // here we drop the null tombstone slots + their gemmi placeholder residues and
      // rebuild the flat atom lists. (Atoms already stay in sync eagerly.)
      void _rebuild_all_atoms() {
         all_atoms.clear();
         for (Model *mw : models) {
            mw->mgr = this;
            mw->all_atoms.clear();
            for (Chain *cw : mw->chains) {
               cw->mgr = this;
               for (Residue *rw : cw->residues)
                  if (rw) {
                     rw->mgr = this;
                     for (Atom *aw : rw->atoms) {
                        // Rebind ownership pointers: residues added via AddResidue/
                        // InsResidue (e.g. add_terminal_residue) carry atoms whose mgr
                        // still points at the deep-copy temporary (or is null). Atom
                        // UDData routes through Atom::mgr, so without this the new atoms
                        // fail Put/GetUDData with WrongHandle — which drops their bonds
                        // (the atom-index UDD never lands) and any UD colouring.
                        aw->mgr = this;
                        aw->res = rw;
                        all_atoms.push_back(aw);
                        mw->all_atoms.push_back(aw);
                     }
                  }
            }
         }
      }
      int FinishStructEdit() {
         for (Model *mw : models)
            for (Chain *cw : mw->chains) cw->_compact_residues();
         _rebuild_all_atoms();
         return 0;
      }

      // ---- UDData registry ----
      struct UDReg {
         UDR_TYPE type;
         int kind;
         std::string name;
         int slot;
      };  // kind:0=int,1=real,2=str
      std::vector<UDReg> ud_regs;
      int ud_counts[5][3] = {{0}};  // [UDR_TYPE][kind] -> next slot

      int RegisterUDInteger(UDR_TYPE t, cpstr name) { return _regUD(t, 0, name); }
      int RegisterUDReal(UDR_TYPE t, cpstr name) { return _regUD(t, 1, name); }
      int RegisterUDString(UDR_TYPE t, cpstr name) { return _regUD(t, 2, name); }
      int GetUDDHandle(UDR_TYPE t, cpstr name) {
         for (int i = 0; i < (int)ud_regs.size(); ++i)
            if (ud_regs[i].type == t && ud_regs[i].name == name) return i + 1;
         return 0;  // MMDB: 0 == "not registered" — Coot relies on `if (handle == 0) Register…`
      }

   private:
      // MMDB UDData handles are 1-based (0 is reserved for "not registered", see
      // GetUDDHandle). Return a 1-based handle; _ud_desc maps back with handle-1.
      int _regUD(UDR_TYPE t, int kind, cpstr name) {
         ud_regs.push_back({t, kind, name ? name : "", ud_counts[t][kind]++});
         return (int)ud_regs.size();
      }

   public:
      // ---- I/O (defined in mmdb-shim/src/io.cc; keeps heavy gemmi write/read
      //      headers out of the ~229 Coot TUs that include mmdb_manager.h) ----
      ERROR_CODE ReadPDBASCII(cpstr fname);
      ERROR_CODE ReadCoorFile(cpstr fname);  // auto-detects PDB / mmCIF
      ERROR_CODE WritePDBASCII(cpstr fname);
      ERROR_CODE WriteCIFASCII(cpstr fname);
   };

}  // namespace mmdb
