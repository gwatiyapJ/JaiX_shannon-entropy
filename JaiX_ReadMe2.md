JaiX Shannon Entropy Analysis Tool for Multiple Sequence Alignments
___________________________________________________________________

This ia a Python/Biopython pipeline for calculating Shannon entropy at every position of a multiple sequence alignment (MSA) and identifying conserved and highly variable residues.

Script: shannon_entropy_clustalO_v1.py
Input: ClustalO-formatted multiple sequence alignment
Analysis: Shannon entropy, residue conservation and sequence variability
Output: CSV, Excel, residue lists, entropy plots and publication-quality figures

Note: Although the current script filename contains clustal, their is the FASTA version of the pipeline which can accepts a gapped FASTA multiple sequence alignment produced by Clustal Omega, MAFFT, MUSCLE, ClustalW or another MSA program.

#============================================
## Installation
#============================================
### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

### Install Dependencies

bash
# Create virtual environment (recommended)
python -m venv entropy_env
source entropy_env/bin/activate  # On Windows: entropy_env\Scripts\activate

# Install required packages
pip install biopython pandas matplotlib numpy openpyxl

# Or install from requirements.txt
pip install -r requirements.txt

#================================
## Basic Usage
#================================

# Analyze a Clustal alignment file
python shannon_entropy_clustalO_v1.py -i alignment.aln -o results

# Specify sequence type explicitly
python shannon_entropy_clustalO_v1.py -i alignment.aln -o results --type protein

# For nucleotide alignments
python shannon_entropy_clustalO_v1.py -i dna_alignment.aln -o results --type nucleotide

#================================
##Command-Line Options
#================================
usage: shannon_entropy_clustalO_v1.py [-h] -i INPUT [-o OUTPUT] [--type {auto,protein,nucleotide}]

Calculate Shannon entropy from a Clustal Omega alignment with validation.

optional arguments:
  -h, --help            Show this help message and exit
  -i INPUT, --input INPUT
                        Clustal alignment file (required)
  -o OUTPUT, --output OUTPUT
                        Output directory (default: entropy_results)
  --type {auto,protein,nucleotide}
                        Sequence type. Default: auto (default: auto)


#=======================================
## Input File Format
#=======================================

Supported Formats

1. Clustal/Clustal Omega (.aln, .clustal)

CLUSTAL O(1.2.4) multiple sequence alignment

AVO03647.1      -----------------------------------MGQIITFFQEVPHVIEEVMNIVLIA    25
AYM51835.1      -----------------------------------MGQIITFFQEVPHVIEEVMNIVLIA    25

2. VPS Format (like Lprotein_vps_data.txt)

POSITION CONSENSUS_SEQUENCE VARIABILITY(H)
______________________________________________
1            -        0.254
2            -        0.254
3            M       -0.000

3. CSV/TSV with position and entropy columns

#=====================================
## Validation Checks 
#=====================================

File existence
✅ Clustal header detection
✅ Minimum 2 sequences
✅ Equal sequence lengths
✅ No empty sequences
✅ Unique sequence IDs
✅ Valid characters (amino acids or nucleotides)
✅ No excessive gaps (>50%)
✅ No gap-only columns
✅ Sequence type consistency
✅ Minimum valid residue content (>50%)
✅ Overall validation status

#===================================
## 📤 Output Files
#==================================

File								Description
___________________________________________________________________________
alignment_validation_report.txt		Detailed validation results
shannon_entropy_results.csv			Complete entropy data (all positions)
shannon_entropy_results.xlsx		Excel workbook with multiple sheets
conserved_residues.csv				Positions with H ≤ 0.5
highly_variable_residues.csv		Positions with H ≥ 1.5
shannon_entropy_profile.png			Basic entropy plot (600 DPI)
publication_shannon_entropy.png		Publication-quality figure (600 DPI)
entropy_summary.txt					Analysis summary statistics

Example Output Structure
entropy_results/
├── alignment_validation_report.txt
├── shannon_entropy_results.csv
├── shannon_entropy_results.xlsx
├── conserved_residues.csv
├── highly_variable_residues.csv
├── shannon_entropy_profile.png
├── publication_shannon_entropy.png
└── entropy_summary.txt

📊 Understanding the Results
Shannon Entropy Interpretation

Entropy (H)		Classification				Biological Meaning
___________________________________________________________________________
0.0 - 0.5		**Conserved**				Critical for structure/function
0.5 - 1.5		**Moderately Variable**		Tolerant to substitution
≥ 1.5			**Highly Variable**			Potential antigenic sites, flexible regions


A typical entropy profile can be interpreted as:

Entropy
  ^
4 |                          /\       highly variable
3 |              /\         /  \
2 |        /\   /  \_______/    \
1 | ______/_ \_/_________________\____
0 | ___________________________________> Position
        conserved      variable

Positions close to zero indicate strong sequence conservation.

High entropy indicates greater residue diversity.

However, entropy should not be interpreted as direct evidence of biological function.

A conserved residue may be structurally or functionally important, but this requires independent evidence.

Similarly, a highly variable residue is not necessarily biologically unimportant.


Maximum Entropy Values
Protein: log₂(20) = 4.32 bits (20 amino acids)
Nucleotide: log₂(4) = 2.00 bits (4 nucleotides)

#=========================================
##  Troubleshooting
#=========================================

Common Errors

1. "FAIL: Invalid characters detected"
Cause: Trailing numbers from Clustal Omega or non-standard characters

Solution:
# The fixed script handles this automatically
# If error persists, check for unusual characters:
grep -o '[^ACDEFGHIKLMNPQRSTVWY\-]' alignment.aln | sort | uniq

2. "ERROR: Unable to read alignment"
Cause: Not a valid Clustal format or BioPython not installed

Solution:
# Verify BioPython installation
pip install --upgrade biopython

# Check file header
head -1 alignment.aln
# Should start with "CLUSTAL"

3. "FAIL: Sequences do not have equal lengths"
Cause: Alignment file corrupted or mixed formats

Solution:
# Re-run alignment with Clustal Omega
clustalo -i input.fasta -o alignment.aln --force

4. "WARNING: Sequences with more than 50% gaps"
Cause: Poor alignment quality or too divergent sequences

Solution:
Remove highly divergent sequences
Adjust alignment parameters
Consider using a different alignment algorithm

#=======================================
## Validation Report Example
#=======================================

CLUSTAL ALIGNMENT VALIDATION REPORT
===================================

PASS: Valid Clustal header detected.
PASS: 15 sequences detected.
PASS: Alignment length = 2289 positions.
PASS: Sequence IDs are unique.
PASS: All sequences have equal length.
PASS: No empty sequences detected.
PASS: No unexpected characters detected.
PASS: No sequence contains excessive gaps.
PASS: No gap-only columns detected.
INFO: Sequence type used for analysis: protein
PASS: All sequences contain sufficient valid residues.

OVERALL VALIDATION: PASSED ✓

#====================================
## 🔬 Scientific Background
#====================================

Shannon Entropy Formula

H = -Σ(p_i × log₂(p_i))

Where:
	H = Shannon entropy (bits)
	p_i = Frequency of residue i at a given position
	Sum is over all residue types at that position

Applications:-
	-Epitope prediction: High variability regions may indicate antigenic sites
	-Functional site identification: Conserved regions often critical for function
	-Phylogenetic analysis: Measure sequence divergence
	-Protein engineering: Identify mutable vs. constrained positions
	-Drug target validation: Conserved sites across strains/species


#==========================================
##  Development Setup
#==========================================

# Clone your fork
git clone https://github.com/yourusername/shannon-entropy-clustal.git
cd shannon-entropy-clustal

# Create virtual environment
python -m venv dev_env
source dev_env/bin/activate

# Install in development mode
pip install -e .
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/

# Check code style
flake8 .
black --check .

#========================================
## 📝 Citation
#========================================
If you use this tool in your research, please cite:

@software{JaiX shannon_entropy_clustal,
  author = James SA,
  title = {Shannon Entropy Analysis for Clustal Omega Alignments},
  year = {2026},
  url = {https://github.com/yourusername/shannon-entropy-clustal},
  version = {1.0.0}
}
