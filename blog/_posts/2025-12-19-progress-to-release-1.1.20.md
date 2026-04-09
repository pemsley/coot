# What’s New in Coot – October 2025 and Beyond

The [Coot](https://github.com/pemsley/coot) project has seen a steady stream of updates since October 2025, spanning new features, bug fixes, improvements to usability, and under-the-hood tweaks. Below we summarize the highlights and directions from more than 250 changes, as reflected in commit messages and RELEASE note updates. (Tip: [browse all recent commits](https://github.com/pemsley/coot/search?q=committer-date%3A%3E%3D2025-10-01&type=commits).)

---

## 🧪 Integration and Extensions

- **RDKit Experimentation**: Several attempts were made to add new methods for generating or manipulating RDKit molecules directly within Coot, although challenges remain ([details](https://github.com/pemsley/coot/commit/30910e270b0280dce19e3f559637ed012e474eab)).
- **coot-headless-api**: Now links with RDKit directly, improving computational ligand handling away from the GUI ([details](https://github.com/pemsley/coot/commit/7286743a7faf774b1256006d4b58cfc8f3fde913)).

---

## 🧰 GUI, Usability, and Workflow

- **Improved Translation Gizmo**: The translation gizmo's visibility logic was overhauled—aimed at delivering a more robust and predictable user experience ([details](https://github.com/pemsley/coot/commit/a183f4e56f458d5f14034bb928576864a9adcb6f)).
- **Rotamer Dialog**: UI spacing for the close button was improved for more professional dialog layouts.
- **Contour Level Output**: The `contour_level_scroll_scrollable_map()` method now prints additional information, making map interpretation easier ([details](https://github.com/pemsley/coot/commit/91a42bd4363d008145c6f6036507e02f4099c6b2)).
- **Drag and Drop**: "File" support was added for drag-and-drop operations in the UI ([details](https://github.com/pemsley/coot/commit/c4ac45443265b8fc8547465c01b7a56c52267259)).
- **Whitespace Clean-up**: Multiple files received whitespace and style clean-up, reflecting a push for more readable code.

---

## 🧬 Ligand Validation and Chemistry

- **InChIKey Store**: Added an internal store for InChIKeys and a key-binding for "go-to-ligand," directly improving ligand navigation ([details](https://github.com/pemsley/coot/commit/11383bfceeab5a91f30d7318272c54d411c72fe9)).
- **Ligand Validation Updates**: Several internal updates made ligand validation and debugging more transparent for those developing or maintaining chemistry features ([details](https://github.com/pemsley/coot/commit/a035b6b7e0ca354ab863859537a928a26b210c6a)).

---

## 📝 Documentation & Release Notes

- **Doxygen Integration**: Added a warning for a missing Doxygen XML file—helpful for developers generating API docs ([details](https://github.com/pemsley/coot/commit/571757fc36309c8e20872d15b596e607dffe8d79)).
- **RELEASE Notes**: Trivial but frequent updates and typo fixes; the team is keeping record-keeping current ([example](https://github.com/pemsley/coot/commit/c3534e62676ae11b0a4fb54c0425ee663a995b0f)).  
- **New Release (v1.1.19)**: A new tagged release pushed in October ([details](https://github.com/pemsley/coot/commit/de1f0ffea43d2df89e16a6395ae8e6a222f61d8d)).

---

## 🛠️ DevOps and Build

- **Executable Names**: Several executables were renamed in `CMakeLists.txt` for clarity and maintainability.
- **Workflow Actions**: Github workflow actions were temporarily disabled, perhaps to triage stability or build issues.

---

## 🧩 Merge Activity and Collaboration

- Numerous merges of branches and remote-tracking branches signal active collaboration.
- Community developers, such as **Jakub Smulski**, have contributed via merges.

---

## 🧼 Miscellaneous

- **.gitignore**: Removed at one point, possibly as part of project clean-up.
- **GMM Files**: Gaussian mixture model files added to the repository, prefiguring future map/model processing features, although not yet integrated into the build ([details](https://github.com/pemsley/coot/commit/e7d4b6f8946677448caacbfae409cec4b38e75fa)).
- **Parser and GUI**: Added for GPhL BUSTER screen output and additional glycosylation debugging ([details](https://github.com/pemsley/coot/commit/293cc2da917f81ad3b21c87079a9cf603e4ea216)).

---

## 🚧 Note: This is only a partial summary

**GitHub indicates there are 278 commits since October 2025**, but only the first 30 could be programmatically summarized here. **More detailed review and changes can be viewed directly in the [commit history since October 2025](https://github.com/pemsley/coot/search?q=committer-date%3A%3E%3D2025-10-01&type=commits).**

---

_Compiled December 2025_

