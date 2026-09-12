# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-09-10

### Added
- Initial release with Clustal Omega support
- Fixed trailing number parsing bug
- 12-point alignment validation system
- Publication-quality entropy plots (600 DPI)
- CSV and Excel output formats
- Position classification (Conserved/Moderate/Highly Variable)
- Comprehensive validation report
- Command-line interface with argparse
- Auto-detection of protein vs. nucleotide sequences

### Fixed
- Clustal Omega trailing position numbers causing validation failures
- `__name__ == "__main__"` syntax error
- Character validation for ambiguous residues
- Gap percentage calculations

### Changed
- Improved error messages with troubleshooting guidance
- Enhanced plot styling for publication use
- Better statistics annotation on figures

### Known Issues
- VPS format support planned for v1.0.1
- Sliding window analysis planned for v1.1.0

---

## [Unreleased]

### Planned
- VPS format file support
- CSV/TSV input support
- Sliding window smoothing
- Interactive Plotly plots
- Batch processing for multiple files
- Web interface option