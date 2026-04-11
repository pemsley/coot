
   TECHNICAL REPORT: MMDB GRAPH MATCHING ANALYSIS FOR GEMMI IMPLEMENTATION

   1. MMDB'S GRAPH MATCHING API

   1.1 Core Classes (mmdb::math namespace)
   - Vertex: Represents graph nodes (atoms)
     * int type: Encoded chemical element (1-based element index or custom encoding)
     * int type_ext: Extended vertex properties (ring info, chirality)
     * int property: User-defined flags (flagwise, for bond counts, chiral marks)
     * std::string name: Atom name
     * Constructors: Vertex(cpstr chem_elem), Vertex(int vtype, cpstr vname)
     * Methods: SetVertex(), SetType(), SetTypeExt(), GetType(), GetTypeExt(),
       GetName(), GetNBonds(), AddBond(), SaveType()/RestoreType()

   - Edge: Represents graph connections (bonds)
     * int v1, v2: 1-indexed vertex numbers
     * int type: Bond type (1=BOND_SINGLE, 2=BOND_DOUBLE, 3=BOND_AROMATIC, 4=BOND_TRIPLE)
     * int property: Additional edge properties
     * Constructors: Edge(int vx1, int vx2, int btype), Edge(int vx1, int vx2, cpstr btype)
     * Methods: SetEdge(), SetType(), GetVertex1(), GetVertex2(), GetType()

   - Graph: Container for vertices and edges representing molecular structure
     * Constructor: Graph(PResidue R, cpstr altLoc=NULL) — directly from mmdb residue
     * Methods:
       - AddVertex(PVertex V), AddEdge(PEdge G): Takes ownership
       - SetVertices(PPVertex V, int vlen), SetEdges(PPEdge G, int glen)
       - MakeGraph(PResidue R, cpstr altLoc): Builds graph from residue atoms and bonds
         Returns: MKGRAPH_Ok, MKGRAPH_NoAtoms, MKGRAPH_ChangedAltLoc, MKGRAPH_MaxOccupancy
       - Build(bool bondOrder): Constructs 2D adjacency matrix (imatrix graph[n][n])
         Sets graph[i][j] = bond_type (or 1 if bondOrder=false)
       - MakeSymmetryRelief(bool noCO2): Assigns symmetry relief modifiers to equivalent vertices
       - IdentifyRings(), IdentifyConnectedComponents()
       - GetVertex(int vertexNo), GetEdge(int edgeNo)
       - GetNofVertices(), GetNofEdges()
       - MakeVertexIDs(): Numbers vertices sequentially

   - GraphMatch: Implements subgraph isomorphism algorithm
     * MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch, bool vertexType=true,
                   VERTEX_EXT_TYPE vertexExt=EXTTYPE_Ignore)
     * Finds maximal common subgraphs of size >= minMatch
     * Returns collection of GMatch objects (one per isomorphism found)
     * Methods:
       - SetFlag(word flag), RemoveFlag(word flag): Controls matching behavior
         GMF_UniqueMatch, GMF_NoCombinations flags available
       - SetMaxNofMatches(int maxNofMatches, bool stopOnMaxN)
       - SetTimeLimit(int maxTimeToRun): Timeout in seconds
       - GetNofMatches(): Returns count of found isomorphisms
       - GetMatch(int MatchNo, ivector &FV1, ivector &FV2, int &nv, realtype &p1, realtype &p2)
         Returns: FV1/FV2 are 1-indexed vertex mappings, nv is match size, p1/p2 are match scores

   - GMatch: Result object for a single match
     * SetMatch(ivector FV1, ivector FV2, int nv, int n, int m)
     * GetMatch(ivector &FV1, ivector &FV2, int &nv, realtype &p1, realtype &p2)
     * GetNofMatches() returns count of internal results
     * isMatch()/isCombination(): Comparison operators

   1.2 Algorithm
   - Uses Ullman's algorithm for subgraph isomorphism detection
   - Citation: Krissinel, E. and Henrick, K. (2004) "Common subgraph isomorphism detection
     by backtracking search" Software - Practice and Experience, 34, 591-607
   - Two implementations available: recursive (Backtrack) and iterative (Ullman)
   - Forward checking during backtracking to prune search space
   - Handles vertex type matching with extended type modifiers via bitwise operations
   - Optional time limiting to prevent long searches

   1.3 Key Properties Encoded in Vertex Type
     Format: 0xCHSSTTTT where:
     - TTTT (bits 0-15): Element type (1-indexed element number or custom)
     - SS (bits 16-23): Symmetry relief modifiers
     - H (bits 24-27): Bond/hydrogen count (stored in upper nibbles)
     - C (bits 28-31): Chirality flags (CHIRAL_LEFT, CHIRAL_RIGHT)
     - TYPE_MASK = 0x00FFFFFF filters to element type only

   2. HOW COOT USES MMDB GRAPH MATCHING

   2.1 Primary Use Case: Dictionary Matching (geometry/dict-utils.cc)
   File: coot/geometry/dict-utils.cc
   Class: coot::dictionary_residue_restraints_t
   Method: match_to_reference()

   Workflow:
     1. Create two graphs from chemical component dictionaries
        - g_1 = working residue dictionary: make_graph(use_hydrogens)
        - g_2 = reference residue dictionary: ref.make_graph(use_hydrogens)

     2. Build internal representation
        - g_1->SetName("working-residue"), g_1->MakeVertexIDs()
        - g_2->SetName("reference-residue"), g_2->MakeVertexIDs()
        - g_1->MakeSymmetryRelief(false)
        - g_2->MakeSymmetryRelief(false)
        - g_1->Build(true), g_2->Build(true)  [bond_order=true]

     3. Configure matching parameters
        - Compute minMatch: ~75% of smaller molecule, clamped to [3, 14]
        - match.SetTimeLimit(2) — limit to 2 seconds
        - bool vertex_type = true

     4. Execute matching
        - match.MatchGraphs(g_1, g_2, minMatch, vertex_type)
        - n_match = match.GetNofMatches()

     5. Select best match and extract atom mappings
        Loop through matches, selecting one with most atoms (nv):
          - match.GetMatch(imatch_best, FV1, FV2, nv, p1, p2)
          - For each matched pair: FV1[ipair] and FV2[ipair] are 1-indexed vertex positions
          - V1 = g_1->GetVertex(FV1[ipair]), V2 = g_2->GetVertex(FV2[ipair])
          - Extract: V1->GetName(), V1->GetType(), V2->GetName()
          - Verify V1->GetType() == V2->GetType() (atoms must be same element)
          - Accumulate atom name mappings for renaming

   2.2 Graph Construction in Coot
   Method: coot::dictionary_residue_restraints_t::make_graph(bool use_hydrogens)

   Implementation (dict-utils.cc:955-1016):
     1. Create empty graph: mmdb::math::Graph *graph = new mmdb::math::Graph
     2. Build name_map: string -> vertex_index (0-indexed) for later edge lookup
     3. Add vertices:
        For each atom in atom_info:
          - Skip hydrogens if use_hydrogens=false
          - Create: new mmdb::math::Vertex(ele.c_str(), name.c_str())
          - graph->AddVertex(v)
          - Track in name_map[atom_id_4c] = i_atom
     4. Add edges from bond restraints:
        For each bond_restraint:
          - Get bond type: mmdb_bond_type = bond_restraint.mmdb_bond_type()
          - Look up both atoms in name_map
          - Skip if hydrogens and use_hydrogens=false
          - Create: new mmdb::math::Edge(it_1->second+1, it_2->second+1, mmdb_bond_type)
            [Note: +1 for 1-indexing]
          - graph->AddEdge(e)
     5. Return graph pointer

   2.3 Other Usage: Connectivity and Ring Detection (ccp4mg-utils/mmut/)
   - mmut_connectivity.cc: Connectivity2::MatchGraphs() for residue matching
     Similar pattern but with HideType/ExcludeType for hydrogen filtering
     Multiple attempts with decreasing minMatch (90%, 80%, ... 10%)

   - mmut_sbase.cc: CMGSBase::MatchGraphs() for SBASE fragment matching
     Creates graphs with MakeGraph(PPAtom atom, int nAtoms) variant
     Iterates through multiple matches, selects by LSQ fit quality

   2.4 Use of Extended Type (type_ext)
   - MakeSymmetryRelief() encodes symmetry information in type_ext
     - Used to distinguish equivalent atoms (e.g., O atoms in CO2)
   - IdentifyRings() marks atoms by ring membership in type_ext
     - ring_mask[n] bits indicate membership in n-membered rings (3-10)

   3. WHAT GEMMI CURRENTLY HAS

   3.1 ChemComp Structure (chemcomp.hpp)
   - BondType enum: Unspec, Single, Double, Triple, Aromatic, Deloc, Metal
   - Restraints::Bond struct with atom1/atom2 (AtomId) and bond_type
   - Restraints::Angle, Torsion, Chirality, Planarity for other restraints
   - ChemComp contains vectors of atoms (Atom) and bonds
   - BUT: No graph representation, no subgraph matching

   3.2 Existing Graph Helper (ace_graph.hpp)
   - AceGraphView: Lightweight adjacency structure
     * std::map<std::string, size_t> atom_index
     * AceBondAdjacency: std::vector<std::vector<AceBondNeighbor>>
     * std::vector<std::vector<int>> neighbors
   - Helper functions for ring detection, aromaticity, CIP ranking
   - NOT a full graph matching implementation
   - Focused on ring/aromaticity analysis for AceDRG, not substructure search

   3.3 Relevant Types
   - Element enum for element types
   - ChemComp::Atom with element, atom_id, and other properties
   - BondType enum readily convertible to/from mmdb bond codes
   - Residue struct in model.hpp with atom array

   4. DESIGN FOR GEMMI SUBGRAPH MATCHING

   4.1 New Core Classes (gemmi/subgraph.hpp or gemmi/graph_match.hpp)

   class Vertex {
     Element element;      // element type
     std::string name;     // atom name
     int type_ext;         // extended properties (rings, chirality)

     // Construction from ChemComp::Atom
     Vertex(const ChemComp::Atom& atom);
     Vertex(Element el, const std::string& name);
   };

   class Edge {
     size_t v1, v2;        // 0-indexed vertex references
     BondType type;        // or int for mmdb compatibility

     Edge(size_t idx1, size_t idx2, BondType btype);
   };

   class Graph {
     std::string name;
     std::vector<Vertex> vertices;    // 0-indexed (differ from mmdb's 1-indexed)
     std::vector<Edge> edges;
     std::vector<std::vector<int>> adjacency;  // adjacency matrix for Build()

     // Construction
     Graph() = default;
     explicit Graph(const ChemComp& cc, bool include_hydrogens = true);

     // Building
     int Build(bool bond_type_matters = true);
     void add_vertex(Element el, const std::string& name);
     void add_edge(size_t idx1, size_t idx2, BondType type);

     // Utilities
     size_t num_vertices() const;
     size_t num_edges() const;
     const Vertex& get_vertex(size_t idx) const;
     const Edge& get_edge(size_t idx) const;
   };

   class SubgraphMatch {
     std::vector<int> F1, F2;    // 0-indexed vertex mappings
     size_t size() const;         // number of matched atoms
     int F1[i], F2[i];           // matched vertex indices
   };

   class SubgraphMatcher {
     std::vector<SubgraphMatch> matches;

     // Main API
     void match(const Graph& g1, const Graph& g2,
                size_t min_match_size,
                bool match_element_type = true,
                bool match_bond_type = true,
                int time_limit_seconds = 0);

     size_t num_matches() const { return matches.size(); }
     const SubgraphMatch& get_match(size_t idx) const;

     // Configuration
     void set_max_matches(size_t max_matches) { ... }
     void set_time_limit(int seconds) { ... }
   };

   4.2 Integration Points with Gemmi
   - Constructor: Graph(const ChemComp& cc, bool include_hydrogens)
     Iterate ChemComp::atoms and ChemComp::bonds
     Use atom_id as vertex name, element for type

   - For Residue matching: Graph(const Residue& res, char alt_loc)
     Similar to mmdb's Graph(PResidue R, cpstr altLoc)
     Handle alternate locations, occupancy filtering

   - BondType conversion: map gemmi::BondType to integer codes
     Single=1, Double=2, Aromatic=3, Triple=4 (matches mmdb)

   4.3 Algorithm Implementation Strategy
   - Adopt Ullman's algorithm (already proven, cited, fast)
   - Use recursive backtracking with forward checking
   - Support optional time limiting via std::chrono
   - Optional: Add iterative version if needed
   - Key optimization: Build adjacency matrix in Build() for O(1) edge lookup

   4.4 API Compatibility Considerations
   - Use 0-indexed internally (more C++-idiomatic)
   - Match results use 0-indexed arrays
   - Document difference from mmdb (1-indexed)
   - Provide helper: graph.get_vertex(match.F1[i])->name equivalent

   4.5 Configuration Options
   - Element type matching (on/off) — default: on
   - Bond order matching (on/off) — default: on
   - Extended type matching (ignore, equal, AND, OR, XOR, etc.) — default: ignore
   - Time limit (seconds) — default: 0 (no limit)
   - Max matches (return first N) — default: unlimited

   5. EXPECTED EFFORT AND COMPLEXITY

   5.1 Core Implementation
   - Vertex, Edge, Graph classes: ~150 lines
   - SubgraphMatch result class: ~30 lines
   - SubgraphMatcher with Ullman's algorithm: ~400-500 lines
   - Build() adjacency matrix: ~50 lines
   - Initialize() forward checking setup: ~100 lines
   - Backtrack() recursive search: ~80 lines
   - Total core: ~800-900 lines C++

   5.2 Integration with Existing Gemmi Code
   - Graph construction from ChemComp: ~50 lines
   - Graph construction from Residue: ~80 lines
   - Bond type conversion table: ~20 lines
   - Helper functions: ~100 lines
   - Total integration: ~250 lines

   5.3 Testing
   - Unit tests for vertex/edge construction
   - Unit tests for Build() adjacency matrix
   - Integration tests: match known similar compounds
   - Regression: match Coot's current results on test cases
   - Performance tests: time vs. graph size

   5.4 Documentation
   - API documentation (doxygen-style)
   - Algorithm citation and explanation
   - Usage examples (dict matching, residue matching)
   - Performance notes and time limit guidance

   6. MIGRATION PATH FOR COOT

   6.1 Minimal Changes Required
   Replace:
     mmdb::math::Graph* g = new mmdb::math::Graph();
   With:
     gemmi::Graph g;

   Replace:
     mmdb::math::GraphMatch match;
     match.MatchGraphs(g1, g2, min_match, vertex_type);
   With:
     gemmi::SubgraphMatcher matcher;
     matcher.match(g1, g2, min_match, vertex_type);

   6.2 Data Type Changes
   - Change mmdb::ivector to std::vector<int>
   - Change mmdb::math::Vertex* to direct Graph::get_vertex(idx) access
   - Update FV1, FV2 array indexing: 1-based → 0-based

   6.3 Compatibility Layer (Optional)
   Could provide wrapper class:
     class GraphMatchAdapter {
       gemmi::SubgraphMatcher impl;
       // Provide mmdb-compatible API
     };

   7. KEY DIFFERENCES FROM MMDB

   Advantage in Gemmi:
   - Modern C++ (STL containers instead of raw pointers)
   - No binary serialization (not needed for in-memory matching)
   - Tighter integration with ChemComp/Residue types
   - Better memory safety (RAII, no manual new/delete)
   - 0-indexed arrays (standard C++)

   Disadvantage vs mmdb:
   - No ChirailtyMask support (could add)
   - No ring identification (separate function)
   - No symmetry relief modifiers (feature for future)

   8. IMPLEMENTATION PRIORITY

   High Priority:
   - Vertex, Edge, Graph classes
   - Build() adjacency matrix
   - SubgraphMatcher with Backtrack algorithm
   - Graph construction from ChemComp
   - Basic element type matching

   Medium Priority:
   - Bond order matching
   - Time limiting
   - Graph construction from Residue (with alt locs)
   - Performance optimization

   Low Priority:
   - Extended type matching (type_ext)
   - Ring identification
   - Symmetry relief modifiers
   - Iterative Ullman() variant
   - Serialization (read/write)

