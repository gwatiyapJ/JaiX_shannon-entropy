# Contributing to Shannon Entropy Clustal

Thank you for considering contributing to this project! Here's how you can help.

## 🎯 How to Contribute

### Reporting Bugs

Before creating bug reports, please check existing issues. When creating a bug report, include:

- **Clear title and description**
- **Steps to reproduce** the behavior
- **Expected vs. actual behavior**
- **Input file** (if possible, or sample)
- **Error messages** and traceback
- **Python version** and OS information

**Example Bug Report:**
```markdown
**Describe the bug**
Validation fails even though alignment looks correct.

**To Reproduce**
1. Run: `python shannon_entropy_clustal_fixed.py -i my_alignment.aln -o results`
2. See error: "FAIL: Invalid characters detected"

**Expected behavior**
Validation should pass for valid Clustal Omega files.

**Screenshots**
[If applicable]

**Environment:**
- Python: 3.10
- OS: Ubuntu 22.04
- BioPython: 1.79