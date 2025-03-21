# Changelog
All notable changes to this project will be documented in this file. 

## [1.0.0] - 2025-03-21

### Maintenance
- **Resolved Pytorch version conflict**: Kept version 2.2.0 for consistency across the codebase.
- **Updated the mapper version**: Note that the reproducibility of the paper cannot be guaranteed with the latest version. For reproducibility, please refer to version 0.0.6 and below.
- **Removed `MØD`** from this version for improved functionality.
- **Replaced `synutility` with `synkit`**: This update ensures consistency and better integration across the modules.

### Deprecated
- **Deprecation Warning**: The following features will be removed in version 1.0.1:
  - The `SynUtils`, `SynComp`, and `SynChemistry` modules have been moved to the `synutility` repository.
  - `ITSConstruction` will be removed, as it has been moved to `synutility`.

### Additional Notes
- Users are encouraged to transition to the updated modules and take note of upcoming removals to ensure smooth upgrades in future versions.


 

## [0.0.6] - 2024-12-09

### Added
- **New feature**: Added `ITSArbitrary` to generate all possible ITS Graphs from AAM (Atom-Atom Mapping).  
  - **Note**: This feature is still under development and may experience slow performance and potential memory issues when the number of combinations exceeds 8! (i.e., 40,320).

- **Dependencies split**: Separate lightweight installation options:
  - `pip install syntemp`: Installs the basic version (without Atom-Atom Mapping tools).
  - `pip install syntemp[all]`: Installs the full set of dependencies, including tools for Atom-Atom Mapping (`rxnmapper`, `localmapper`, and `graphormermapper`).

### Changed
- **Enhanced rule clustering**: Improved rule clustering functionality to handle batch processing and mitigate the combinatorial explosion problem.
  - **Note**: This change does not yet integrate with Hierarchical Clustering.
  
- **Improved Isomorphic Filter**: Integrated a new graph signature filter, reducing the ITS clustering time from 78 minutes to just 2 minutes.

- **Hydrogen Inference Refactor**: 
  - Refactored the hydrogen inference function for better readability and efficiency.
  - Integrated graph signature to reduce isomorphism checks, improving processing time by 50%.
  - Removed the `timeout` option (planned for removal in version 0.0.10 due to redundancy).
  - **Issue**: The process remains slow, and there may be potential memory explosion if the number of combinations exceeds 8! (i.e., 40,320).

- **ITSExtraction refactor**: Removed redundant `deepcopy` calls, resulting in a 50% reduction in processing time.

- **Code Cleanup**: 
  - Removed redundant Python functions, which have now been moved to the `synutility` repository.
  - Removed unused variables in the `syntemp` command-line interface (CLI).

### Deprecated
- **Deprecation Warning**: The following features will be removed in version 0.0.10:
  - `SynUtils` and `SynChemistry` modules have been moved to the `synutility` repository.
  - `ITSConstruction` will also be removed as it has been moved to `synutility`.

### Fixed
- **Memory Usage & Performance**: Improved memory management and processing performance in several functions, especially in hydrogen inference and ITS extraction.

### Security
- No security updates in this release.

---
