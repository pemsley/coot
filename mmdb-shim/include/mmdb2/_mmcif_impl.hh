// -*- mode: c++; -*-
//
// mmdb-shim: mmdb::mmcif::* implemented as a thin veneer over gemmi::cif.
// See coot-shim-prefer-gemmi memory: NEVER hand-write CIF parsing — gemmi does
// the real work; this only presents the MMDB-shaped API Coot's geometry/ uses.
//
// Copyright 2026 by Medical Research Council Laboratory of Molecular Biology
//
// Included at the end of _shim_impl.hh (so pstr/cpstr/realtype already exist).
#ifndef COOT_MMDB_SHIM_MMCIF_IMPL_HH
#define COOT_MMDB_SHIM_MMCIF_IMPL_HH

#include <gemmi/cifdoc.hpp>
#include <gemmi/cif.hpp>       // gemmi::cif::read_file
#include <gemmi/to_cif.hpp>    // gemmi::cif::write_cif_to_stream

#include <deque>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace mmdb {
namespace mmcif {

// ---- return codes (mirror mmdb_mmcif_.h) --------------------------------
enum {
  CIFRC_Loop           =  2,
  CIFRC_Structure      =  1,
  CIFRC_Ok             =  0,
  CIFRC_StructureNoTag = -1,
  CIFRC_LoopNoTag      = -2,
  CIFRC_NoCategory     = -3,
  CIFRC_WrongFormat    = -4,
  CIFRC_NoTag          = -5,
  CIFRC_NotAStructure  = -6,
  CIFRC_NotALoop       = -7,
  CIFRC_WrongIndex     = -8,
  CIFRC_NoField        = -9,
  CIFRC_Created        = -12,
  CIFRC_CantOpenFile   = -13,
  CIFRC_NoDataLine     = -14,
  CIFRC_NoData         = -15
};

// ---- file flags ----------------------------------------------------------
enum {
  CIFFL_PrintWarnings     = 0x00000001,
  CIFFL_StopOnWarnings    = 0x00000002,
  CIFFL_SuggestCategories = 0x00000004,
  CIFFL_SuggestTags       = 0x00000008
};

enum MMCIF_ITEM {
  MMCIF_None = 0, MMCIF_Struct = 1, MMCIF_Loop = 2,
  MMCIF_Data = 3, MMCIF_Category = 4
};

class Loop;
class Struct;
class Category;
class Data;
class File;
typedef Loop     *PLoop;
typedef Struct   *PStruct;
typedef Category *PCategory;
typedef Data     *PData;
typedef File     *PFile;

namespace detail {
  // MMDB category names have no trailing dot; gemmi wants "_cat." — normalise
  // to WITH-dot internally so full tag = cat + subtag.
  inline std::string with_dot(const char *cat) {
    std::string s(cat ? cat : "");
    if (s.empty() || s.back() != '.') s += '.';
    return s;
  }
  inline std::string strip_dot(const std::string &s) {
    if (!s.empty() && s.back() == '.') return s.substr(0, s.size() - 1);
    return s;
  }
}

// =========================================================================
//  Category — just a named handle (Coot uses GetCategoryName / GetCategoryID)
// =========================================================================
class Category {
 public:
  std::string cat;                  // WITH trailing dot
  MMCIF_ITEM  kind = MMCIF_Category;
  std::deque<std::string> sret;
  Category() {}
  pstr GetCategoryName() {
    sret.push_back(detail::strip_dot(cat));
    return (pstr) sret.back().c_str();
  }
  MMCIF_ITEM GetCategoryID() { return kind; }
};

// =========================================================================
//  Loop
// =========================================================================
class Loop {
 public:
  Data       *owner = nullptr;      // null for a bare `new Loop`
  std::string cat;                  // WITH trailing dot
  gemmi::cif::Loop *direct = nullptr;  // FindLoop binds the gemmi loop directly
  bool        write_mode = false;
  // write buffer (row-major); sub-tags only (no category prefix)
  std::vector<std::string>              wtags;
  std::vector<std::vector<std::string>> wrows;
  // storage for borrowed pstr returns (Coot never frees these)
  std::deque<std::string> sret;

  Loop() {}

  gemmi::cif::Loop *gloop() const;  // read loop, or nullptr (defined after Data)

  int  GetLoopLength();
  int  GetNofTags();
  pstr GetTag(int tagNo);
  pstr GetField(int rowNo, int tagNo);

  pstr GetString  (cpstr TName, int nrow, int &RC);
  int  GetReal    (realtype &R, cpstr TName, int nrow, bool Remove = false);
  int  GetInteger (int &I, cpstr TName, int nrow, bool Remove = false);

  void AddLoopTag (cpstr T, bool Remove = true) { (void) Remove; wcol(T, true); }
  void PutString  (cpstr S, cpstr T, int nrow) { wput(T, nrow, S ? S : "."); }
  void PutInteger (int  I, cpstr T, int nrow)  { wput(T, nrow, std::to_string(I)); }
  void PutReal    (realtype R, cpstr T, int nrow, int prec = 8) {
    char b[64]; std::snprintf(b, sizeof b, "%.*f", prec, (double) R); wput(T, nrow, b);
  }
  void PutReal    (realtype R, cpstr T, int nrow, cpstr /*format*/) { PutReal(R, T, nrow, 8); }

  // write helpers
  int  wcol(cpstr T, bool create);
  void wput(cpstr T, int nrow, const std::string &val);
  void flush(gemmi::cif::Block &b);
};

inline int Loop::wcol(cpstr T, bool create) {
  for (size_t i = 0; i < wtags.size(); ++i)
    if (wtags[i] == T) return (int) i;
  if (!create) return -1;
  wtags.push_back(T);
  for (auto &row : wrows) row.resize(wtags.size());
  return (int) wtags.size() - 1;
}

inline void Loop::wput(cpstr T, int nrow, const std::string &val) {
  write_mode = true;
  int col = wcol(T, true);
  if (nrow < 0) nrow = 0;
  while ((int) wrows.size() <= nrow) wrows.emplace_back(wtags.size());
  wrows[nrow][col] = val;
}

inline void Loop::flush(gemmi::cif::Block &b) {
  if (wtags.empty()) return;
  gemmi::cif::Loop &gl = b.init_mmcif_loop(cat, wtags);  // tags become cat+subtag
  gl.values.clear();
  gl.values.reserve(wrows.size() * wtags.size());
  for (auto &row : wrows)
    for (size_t c = 0; c < wtags.size(); ++c) {
      const std::string &v = c < row.size() ? row[c] : std::string();
      gl.values.push_back(v.empty() ? "." : gemmi::cif::quote(v));
    }
}

inline int Loop::GetLoopLength() {
  gemmi::cif::Loop *g = gloop();
  return g ? (int) g->length() : (int) wrows.size();
}
inline int Loop::GetNofTags() {
  gemmi::cif::Loop *g = gloop();
  return g ? (int) g->width() : (int) wtags.size();
}

// =========================================================================
//  Struct (single-value category = a set of tag/value pairs)
// =========================================================================
class Struct {
 public:
  Data       *owner = nullptr;
  std::string cat;                  // WITH trailing dot
  bool        write_mode = false;
  std::vector<std::pair<std::string, std::string>> wpairs;
  std::deque<std::string> sret;

  Struct() {}

  pstr GetCategoryName() {
    sret.push_back(detail::strip_dot(cat));
    return (pstr) sret.back().c_str();
  }

  int  GetNofTags();
  pstr GetTag(int tagNo);
  pstr GetField(int tagNo);
  pstr GetString  (cpstr TName, int &RC);
  int  GetReal    (realtype &R, cpstr TName, bool Remove = false);
  int  GetInteger (int &I, cpstr TName, bool Remove = false);

  void PutString  (cpstr S, cpstr TName, bool /*Concatenate*/ = false) {
    write_mode = true; wpairs.emplace_back(TName, S ? S : ".");
  }
  void PutReal    (realtype R, cpstr TName, int prec = 8) {
    char b[64]; std::snprintf(b, sizeof b, "%.*f", prec, (double) R);
    write_mode = true; wpairs.emplace_back(TName, b);
  }
  void PutReal    (realtype R, cpstr TName, cpstr /*format*/) { PutReal(R, TName, 8); }
  void PutInteger (int I, cpstr TName) {
    write_mode = true; wpairs.emplace_back(TName, std::to_string(I));
  }

  std::vector<std::string> collect_tags();   // read: sub-tags present in block
  void flush(gemmi::cif::Block &b);
};

// =========================================================================
//  Data (a data_ block)
// =========================================================================
class Data {
 public:
  gemmi::cif::Document *doc = nullptr;   // resolve block by INDEX (blocks vector reallocs)
  size_t                idx = 0;
  std::unique_ptr<gemmi::cif::Document> owned_doc;  // for a standalone `new Data()`

  std::deque<Loop>     loops;
  std::deque<Struct>   structs;
  std::deque<Category> cats_pool;
  std::unordered_map<std::string, Loop*>   loop_by_cat;
  std::unordered_map<std::string, Struct*> struct_by_cat;
  std::vector<std::string> cat_names;    // WITH dot
  bool cats_built = false;
  std::deque<std::string> sret;

  Data() {}

  gemmi::cif::Block &blk() { return doc->blocks[idx]; }

  // standalone read (Coot: `Data d; d.ReadMMCIFData(fname)`) — own a Document and
  // point at its first block. Used for small-molecule CIFs.
  int SetFlag(int /*flag*/) { return 0; }   // parse flags are gemmi-internal — no-op
  int ReadMMCIFData(cpstr fname) {
    try { owned_doc.reset(new gemmi::cif::Document(gemmi::cif::read_file(fname ? fname : ""))); }
    catch (const std::exception &) { return CIFRC_CantOpenFile; }
    if (owned_doc->blocks.empty()) return CIFRC_NoDataLine;
    doc = owned_doc.get(); idx = 0; cats_built = false;
    return CIFRC_Ok;
  }
  // find the loop containing tags[0] (a null-terminated tag array; core-CIF flat
  // tags). Binds the gemmi loop directly (cat="" so GetString uses full tags).
  // Coot passes both `pstr[]` and `const char*[]`, so accept cpstr.
  PLoop FindLoop(cpstr *tags) {
    if (!tags || !tags[0]) return nullptr;
    gemmi::cif::Loop *gl = blk().find_loop(tags[0]).get_loop();
    if (!gl) return nullptr;
    loops.emplace_back();
    Loop &l = loops.back();
    l.owner = this; l.cat = ""; l.direct = gl;
    return &l;
  }
  PLoop FindLoop(pstr *tags) { return FindLoop((cpstr *) tags); }

  void build_cats() {
    if (cats_built) return;
    cat_names = blk().get_mmcif_category_names();   // returns WITH trailing dot
    cats_built = true;
  }

  pstr GetDataName() {
    sret.push_back(blk().name);
    return (pstr) sret.back().c_str();
  }
  void GetDataName(pstr &dname, bool /*Remove*/ = false) {
    sret.push_back(blk().name);
    dname = (pstr) sret.back().c_str();
  }
  void PutDataName(cpstr dname) { blk().name = dname ? dname : ""; }

  int GetNumberOfCategories() { build_cats(); return (int) cat_names.size(); }

  PCategory GetCategory(int categoryNo) {
    build_cats();
    if (categoryNo < 0 || (size_t) categoryNo >= cat_names.size()) return nullptr;
    cats_pool.emplace_back();
    Category &c = cats_pool.back();
    c.cat = cat_names[categoryNo];
    gemmi::cif::Table t = blk().find_mmcif_category(c.cat);
    c.kind = t.get_loop() ? MMCIF_Loop : MMCIF_Struct;
    return &c;
  }

  PLoop GetLoop(cpstr CName) {
    std::string key = detail::with_dot(CName);
    auto it = loop_by_cat.find(key);
    if (it != loop_by_cat.end()) return it->second;
    if (!blk().find_mmcif_category(key).get_loop()) return nullptr;  // absent or a struct
    loops.emplace_back();
    Loop &l = loops.back();
    l.owner = this; l.cat = key; l.write_mode = false;
    loop_by_cat[key] = &l;
    return &l;
  }

  PStruct GetStructure(cpstr CName) {
    std::string key = detail::with_dot(CName);
    auto it = struct_by_cat.find(key);
    if (it != struct_by_cat.end()) return it->second;
    if (!blk().has_mmcif_category(key)) return nullptr;
    if (blk().find_mmcif_category(key).get_loop()) return nullptr;   // it's a loop
    structs.emplace_back();
    Struct &s = structs.back();
    s.owner = this; s.cat = key; s.write_mode = false;
    struct_by_cat[key] = &s;
    return &s;
  }

  int GetLoopLength(cpstr CName) {
    PLoop l = GetLoop(CName);
    return l ? l->GetLoopLength() : CIFRC_NoCategory;
  }

  // full mmCIF tag from (CName, TName): if CName is empty, TName is already the
  // full tag (small-molecule CIFs pass "" + "_cell_length_a").
  std::string full_tag(cpstr CName, cpstr TName) {
    std::string t = TName ? TName : "";
    return (CName && CName[0]) ? detail::with_dot(CName) + t : t;
  }
  // struct-style direct access (Data::GetString(CName, TName, RC) etc.)
  pstr GetString(cpstr CName, cpstr TName, int &RC) {
    const std::string *v = blk().find_value(full_tag(CName, TName));
    if (!v) { RC = CIFRC_NoTag; return nullptr; }
    RC = CIFRC_Ok;
    if (gemmi::cif::is_null(*v)) return nullptr;
    sret.push_back(gemmi::cif::as_string(*v));
    return (pstr) sret.back().c_str();
  }
  // pstr& form: sets S to the value, returns a CIFRC code (Coot: ierr += ...)
  int GetString(pstr &S, cpstr CName, cpstr TName, bool /*Remove*/ = false) {
    int rc = 0;
    S = GetString(CName, TName, rc);
    return rc;
  }
  int GetReal(realtype &R, cpstr CName, cpstr TName, bool /*Remove*/ = false) {
    R = 0;
    const std::string *v = blk().find_value(full_tag(CName, TName));
    if (!v) return CIFRC_NoTag;
    if (gemmi::cif::is_null(*v)) return CIFRC_NoData;
    try { R = std::stod(gemmi::cif::as_string(*v)); } catch (...) { return CIFRC_WrongFormat; }
    return CIFRC_Ok;
  }
  int GetInteger(int &I, cpstr CName, cpstr TName, bool /*Remove*/ = false) {
    I = 0;
    const std::string *v = blk().find_value(full_tag(CName, TName));
    if (!v) return CIFRC_NoTag;
    if (gemmi::cif::is_null(*v)) return CIFRC_NoData;
    try { I = std::stoi(gemmi::cif::as_string(*v)); } catch (...) { return CIFRC_WrongFormat; }
    return CIFRC_Ok;
  }

  int AddLoop(cpstr CName, PLoop &cifLoop) {
    std::string key = detail::with_dot(CName);
    auto it = loop_by_cat.find(key);
    if (it != loop_by_cat.end()) { cifLoop = it->second; return CIFRC_Ok; }
    loops.emplace_back();
    Loop &l = loops.back();
    l.owner = this; l.cat = key; l.write_mode = true;
    loop_by_cat[key] = &l;
    cifLoop = &l;
    return CIFRC_Created;
  }
  int AddStructure(cpstr CName, PStruct &cifStruct) {
    std::string key = detail::with_dot(CName);
    auto it = struct_by_cat.find(key);
    if (it != struct_by_cat.end()) { cifStruct = it->second; return CIFRC_Ok; }
    structs.emplace_back();
    Struct &s = structs.back();
    s.owner = this; s.cat = key; s.write_mode = true;
    struct_by_cat[key] = &s;
    cifStruct = &s;
    return CIFRC_Created;
  }

  void flush() {
    for (auto &l : loops)   if (l.write_mode) l.flush(blk());
    for (auto &s : structs) if (s.write_mode) s.flush(blk());
  }
};

// ---- Loop methods that need a complete Data -----------------------------
inline gemmi::cif::Loop *Loop::gloop() const {
  if (direct) return direct;                    // FindLoop-bound (core CIF)
  if (write_mode || !owner) return nullptr;
  return owner->blk().find_mmcif_category(cat).get_loop();
}

inline pstr Loop::GetString(cpstr TName, int nrow, int &RC) {
  gemmi::cif::Loop *g = gloop();
  if (!g) { RC = CIFRC_NotALoop; return nullptr; }
  int col = g->find_tag(cat + TName);
  if (col < 0) { RC = CIFRC_NoTag; return nullptr; }
  if (nrow < 0 || (size_t) nrow >= g->length()) { RC = CIFRC_WrongIndex; return nullptr; }
  const std::string &raw = g->val(nrow, col);
  RC = CIFRC_Ok;
  if (gemmi::cif::is_null(raw)) return nullptr;
  sret.push_back(gemmi::cif::as_string(raw));
  return (pstr) sret.back().c_str();
}
inline int Loop::GetReal(realtype &R, cpstr TName, int nrow, bool /*Remove*/) {
  R = 0;
  gemmi::cif::Loop *g = gloop();
  if (!g) return CIFRC_NotALoop;
  int col = g->find_tag(cat + TName);
  if (col < 0) return CIFRC_NoTag;
  if (nrow < 0 || (size_t) nrow >= g->length()) return CIFRC_WrongIndex;
  const std::string &raw = g->val(nrow, col);
  if (gemmi::cif::is_null(raw)) return CIFRC_NoData;
  try { R = std::stod(gemmi::cif::as_string(raw)); } catch (...) { return CIFRC_WrongFormat; }
  return CIFRC_Ok;
}
inline int Loop::GetInteger(int &I, cpstr TName, int nrow, bool /*Remove*/) {
  I = 0;
  gemmi::cif::Loop *g = gloop();
  if (!g) return CIFRC_NotALoop;
  int col = g->find_tag(cat + TName);
  if (col < 0) return CIFRC_NoTag;
  if (nrow < 0 || (size_t) nrow >= g->length()) return CIFRC_WrongIndex;
  const std::string &raw = g->val(nrow, col);
  if (gemmi::cif::is_null(raw)) return CIFRC_NoData;
  try { I = std::stoi(gemmi::cif::as_string(raw)); } catch (...) { return CIFRC_WrongFormat; }
  return CIFRC_Ok;
}
inline pstr Loop::GetTag(int tagNo) {
  gemmi::cif::Loop *g = gloop();
  if (!g) { if (tagNo < 0 || (size_t) tagNo >= wtags.size()) return nullptr;
            sret.push_back(wtags[tagNo]); return (pstr) sret.back().c_str(); }
  if (tagNo < 0 || (size_t) tagNo >= g->tags.size()) return nullptr;
  std::string t = g->tags[tagNo];
  if (t.size() > cat.size() && t.compare(0, cat.size(), cat) == 0) t = t.substr(cat.size());
  sret.push_back(t);
  return (pstr) sret.back().c_str();
}
inline pstr Loop::GetField(int rowNo, int tagNo) {
  gemmi::cif::Loop *g = gloop();
  if (!g) return nullptr;
  if (tagNo < 0 || (size_t) tagNo >= g->width())  return nullptr;
  if (rowNo < 0 || (size_t) rowNo >= g->length()) return nullptr;
  const std::string &raw = g->val(rowNo, tagNo);
  if (gemmi::cif::is_null(raw)) return nullptr;
  sret.push_back(gemmi::cif::as_string(raw));
  return (pstr) sret.back().c_str();
}

// ---- Struct methods that need a complete Data ---------------------------
inline std::vector<std::string> Struct::collect_tags() {
  std::vector<std::string> out;
  if (!owner) return out;
  for (const gemmi::cif::Item &it : owner->blk().items)
    if (it.type == gemmi::cif::ItemType::Pair &&
        it.pair[0].compare(0, cat.size(), cat) == 0)
      out.push_back(it.pair[0].substr(cat.size()));
  return out;
}
inline int Struct::GetNofTags() { return (int) collect_tags().size(); }
inline pstr Struct::GetTag(int tagNo) {
  std::vector<std::string> t = collect_tags();
  if (tagNo < 0 || (size_t) tagNo >= t.size()) return nullptr;
  sret.push_back(t[tagNo]);
  return (pstr) sret.back().c_str();
}
inline pstr Struct::GetField(int tagNo) {
  std::vector<std::string> t = collect_tags();
  if (tagNo < 0 || (size_t) tagNo >= t.size()) return nullptr;
  int rc = 0;
  return GetString(t[tagNo].c_str(), rc);
}
inline pstr Struct::GetString(cpstr TName, int &RC) {
  if (!owner) { RC = CIFRC_NoTag; return nullptr; }
  const std::string *v = owner->blk().find_value(cat + TName);
  if (!v) { RC = CIFRC_NoTag; return nullptr; }
  RC = CIFRC_Ok;
  if (gemmi::cif::is_null(*v)) return nullptr;
  sret.push_back(gemmi::cif::as_string(*v));
  return (pstr) sret.back().c_str();
}
inline int Struct::GetReal(realtype &R, cpstr TName, bool /*Remove*/) {
  R = 0;
  if (!owner) return CIFRC_NoTag;
  const std::string *v = owner->blk().find_value(cat + TName);
  if (!v) return CIFRC_NoTag;
  if (gemmi::cif::is_null(*v)) return CIFRC_NoData;
  try { R = std::stod(gemmi::cif::as_string(*v)); } catch (...) { return CIFRC_WrongFormat; }
  return CIFRC_Ok;
}
inline int Struct::GetInteger(int &I, cpstr TName, bool /*Remove*/) {
  I = 0;
  if (!owner) return CIFRC_NoTag;
  const std::string *v = owner->blk().find_value(cat + TName);
  if (!v) return CIFRC_NoTag;
  if (gemmi::cif::is_null(*v)) return CIFRC_NoData;
  try { I = std::stoi(gemmi::cif::as_string(*v)); } catch (...) { return CIFRC_WrongFormat; }
  return CIFRC_Ok;
}
inline void Struct::flush(gemmi::cif::Block &b) {
  for (auto &tv : wpairs)
    b.set_pair(cat + tv.first, gemmi::cif::quote(tv.second));
}

// =========================================================================
//  File (a whole CIF document)
// =========================================================================
class File {
 public:
  gemmi::cif::Document doc;
  std::deque<Data>     datas;
  std::deque<std::string> sret;

  File() {}

  void rebuild() {
    datas.clear();
    for (size_t i = 0; i < doc.blocks.size(); ++i) {
      datas.emplace_back();
      datas.back().doc = &doc;
      datas.back().idx = i;
    }
  }

  int ReadMMCIFFile(cpstr FName, int /*flags*/ = 0) {
    try { doc = gemmi::cif::read_file(FName ? FName : ""); }
    catch (const std::exception &) { return CIFRC_CantOpenFile; }
    rebuild();
    return CIFRC_Ok;
  }

  int WriteMMCIFFile(cpstr FName, int /*flags*/ = 0) {
    for (auto &d : datas) d.flush();
    std::ofstream os(FName ? FName : "");
    if (!os) return CIFRC_CantOpenFile;
    gemmi::cif::write_cif_to_stream(os, doc, gemmi::cif::WriteOptions());
    return 0;
  }

  int   GetNofData() { return (int) datas.size(); }
  int   GetNumberOfData() { return (int) datas.size(); }

  PData GetCIFData(int dataNo) {
    if (dataNo < 0 || (size_t) dataNo >= datas.size()) return nullptr;
    return &datas[dataNo];
  }
  PData GetCIFData(cpstr name) {
    std::string n(name ? name : "");
    for (auto &d : datas)
      if (d.blk().name == n) return &d;
    return nullptr;
  }

  int AddCIFData(cpstr name) {
    std::string n(name ? name : "");
    if (doc.find_block(n)) return CIFRC_Ok;
    doc.add_new_block(n);
    // index-based Data wrappers survive blocks-vector reallocation, so just
    // append one for the new block (do NOT rebuild — that would drop buffered
    // write state of earlier Data objects).
    datas.emplace_back();
    datas.back().doc = &doc;
    datas.back().idx = doc.blocks.size() - 1;
    return CIFRC_Created;
  }
};

// ---- free helpers --------------------------------------------------------
inline pstr GetCIFMessage(pstr buffer, int rc) {
  const char *m = "unknown mmCIF return code";
  switch (rc) {
    case CIFRC_Ok:             m = "no errors"; break;
    case CIFRC_NoCategory:     m = "category not found"; break;
    case CIFRC_NoTag:          m = "tag not found"; break;
    case CIFRC_NoField:        m = "field not found"; break;
    case CIFRC_WrongFormat:    m = "wrong value format"; break;
    case CIFRC_WrongIndex:     m = "row index out of range"; break;
    case CIFRC_NotALoop:       m = "category is not a loop"; break;
    case CIFRC_NotAStructure:  m = "category is not a structure"; break;
    case CIFRC_NoData:         m = "no data"; break;
    case CIFRC_CantOpenFile:   m = "cannot open file"; break;
    case CIFRC_NoDataLine:     m = "no data_ line"; break;
    case CIFRC_Created:        m = "category created"; break;
    default: break;
  }
  std::strcpy(buffer, m);
  return buffer;
}

}  // namespace mmcif
}  // namespace mmdb

#endif  // COOT_MMDB_SHIM_MMCIF_IMPL_HH
