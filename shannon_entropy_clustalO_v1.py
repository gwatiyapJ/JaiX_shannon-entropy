###JaiX Shannon Entropy Analysis Tool
#!/usr/bin/env python3
"""
SHANNON ENTROPY ANALYSIS PIPELINE - FIXED FOR CLUSTAL OMEGA
Handles Clustal Omega output with trailing position numbers.
"""

import os
import argparse
import re
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO

# CONFIGURATION
# Constants
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
NUCLEOTIDES = set("ACGT")
GAP_CHARACTERS = set("-?.")
NUCLEOTIDE_AMBIGUOUS = set("RYSWKMBDHVN")
PROTEIN_AMBIGUOUS = set("BXZJUO")
CONSERVATION_THRESHOLD = 0.5
VARIABLE_THRESHOLD = 1.5
MAX_GAP_PERCENTAGE = 50.0
MIN_VALID_RESIDUE_PERCENTAGE = 50.0


def read_clustal_alignment(input_file):
    """Read Clustal alignment with proper handling of Omega format."""
    print("\nReading Clustal alignment...")
    print(f"Input file: {input_file}")
    
    try:
        alignment = AlignIO.read(input_file, "clustal")
    except Exception as error:
        print("\nERROR: Unable to read alignment.")
        print("Make sure this is a standard Clustal/Clustal Omega file.")
        print(f"\nBioPython error:\n{error}")
        raise SystemExit(1)
    
    if len(alignment) == 0:
        raise SystemExit("ERROR: No sequences were found in the alignment.")
    if alignment.get_alignment_length() == 0:
        raise SystemExit("ERROR: Alignment contains no residues.")
    
    print(f"Number of sequences : {len(alignment)}")
    print(f"Alignment length    : {alignment.get_alignment_length()}")
    return alignment


def clean_sequence_string(seq_string):
    """
    Remove trailing position numbers from Clustal Omega output.
    Example: '---MGQIITFFQEVPHVIEEVMNIVLIA25' → '---MGQIITFFQEVPHVIEEVMNIVLIA'
    """
    # Remove trailing digits and whitespace
    cleaned = re.sub(r'\d+\s*$', '', seq_string)
    # Remove any remaining whitespace
    cleaned = cleaned.replace(' ', '').replace('\t', '')
    return cleaned


def check_clustal_header(input_file):
    """Check for Clustal header (handles Omega format)."""
    with open(input_file, "r", encoding="utf-8", errors="replace") as file:
        first_line = file.readline().strip()
        # Accept both "CLUSTAL" and "CLUSTAL O" formats
        if first_line.upper().startswith("CLUSTAL"):
            return True, "PASS: Valid Clustal header detected."
        return False, (
            "WARNING: File does not begin with a standard CLUSTAL header."
        )


def detect_sequence_type(alignment):
    """Detect if alignment contains protein or nucleotide sequences."""
    sequence = ""
    for record in alignment:
        seq_str = clean_sequence_string(str(record.seq))
        sequence += seq_str.upper()
    
    for character in GAP_CHARACTERS:
        sequence = sequence.replace(character, "")
    
    if not sequence:
        return "unknown"
    
    nucleotide_like = NUCLEOTIDES | NUCLEOTIDE_AMBIGUOUS
    nucleotide_count = sum(
        character in nucleotide_like for character in sequence
    )
    nucleotide_fraction = nucleotide_count / len(sequence)
    
    if nucleotide_fraction >= 0.95:
        return "nucleotide"
    return "protein"


def validate_alignment(alignment, input_file, sequence_type):
    """Validate alignment with proper sequence cleaning."""
    results = []
    validation_passed = True
    
    # Header check
    header_ok, header_message = check_clustal_header(input_file)
    results.append(header_message)
    
    # Sequence count
    number_sequences = len(alignment)
    if number_sequences < 2:
        results.append("FAIL: Alignment contains fewer than 2 sequences.")
        validation_passed = False
    else:
        results.append(f"PASS: {number_sequences} sequences detected.")
    
    # Alignment length
    alignment_length = alignment.get_alignment_length()
    if alignment_length == 0:
        results.append("FAIL: Alignment contains zero positions.")
        validation_passed = False
    else:
        results.append(f"PASS: Alignment length = {alignment_length} positions.")
    
    # Duplicate IDs
    sequence_ids = [record.id for record in alignment]
    duplicate_ids = [
        sequence_id
        for sequence_id, count in Counter(sequence_ids).items()
        if count > 1
    ]
    if duplicate_ids:
        results.append(
            "FAIL: Duplicate sequence IDs detected: "
            + ", ".join(duplicate_ids)
        )
        validation_passed = False
    else:
        results.append("PASS: Sequence IDs are unique.")
    
    # Sequence lengths
    sequence_lengths = {record.id: len(record.seq) for record in alignment}
    unique_lengths = set(sequence_lengths.values())
    if len(unique_lengths) != 1:
        results.append("FAIL: Sequences do not have equal lengths.")
        for sequence_id, length in sequence_lengths.items():
            results.append(f"    {sequence_id}: {length} residues")
        validation_passed = False
    else:
        results.append("PASS: All sequences have equal length.")
    
    # Empty sequences
    empty_sequences = []
    for record in alignment:
        sequence = clean_sequence_string(str(record.seq)).strip()
        if not sequence:
            empty_sequences.append(record.id)
    if empty_sequences:
        results.append(
            "FAIL: Empty sequences detected: "
            + ", ".join(empty_sequences)
        )
        validation_passed = False
    else:
        results.append("PASS: No empty sequences detected.")
    
    # Character validation (CRITICAL FIX)
    if sequence_type == "protein":
        allowed_characters = AMINO_ACIDS | GAP_CHARACTERS | PROTEIN_AMBIGUOUS
    elif sequence_type == "nucleotide":
        allowed_characters = NUCLEOTIDES | GAP_CHARACTERS | NUCLEOTIDE_AMBIGUOUS
    else:
        allowed_characters = AMINO_ACIDS | NUCLEOTIDES | GAP_CHARACTERS
    
    invalid_characters = {}
    for record in alignment:
        # Clean sequence before validation
        sequence = clean_sequence_string(str(record.seq)).upper()
        invalid = sorted(set(sequence) - allowed_characters)
        if invalid:
            invalid_characters[record.id] = invalid
    
    if invalid_characters:
        results.append("FAIL: Unexpected characters detected.")
        for sequence_id, characters in invalid_characters.items():
            results.append(
                f"    {sequence_id}: {', '.join(characters)}"
            )
        validation_passed = False
    else:
        results.append("PASS: No unexpected characters detected.")
    
    # Gap percentage check
    sequences_with_excessive_gaps = []
    for record in alignment:
        sequence = clean_sequence_string(str(record.seq)).upper()
        gap_count = sum(character in GAP_CHARACTERS for character in sequence)
        if len(sequence) > 0:
            gap_percentage = (gap_count / len(sequence)) * 100
            if gap_percentage > MAX_GAP_PERCENTAGE:
                sequences_with_excessive_gaps.append((record.id, gap_percentage))
    
    if sequences_with_excessive_gaps:
        results.append(
            f"WARNING: Sequences with more than {MAX_GAP_PERCENTAGE}% gaps:"
        )
        for sequence_id, percentage in sequences_with_excessive_gaps:
            results.append(f"    {sequence_id}: {percentage:.2f}% gaps")
    else:
        results.append("PASS: No sequence contains excessive gaps.")
    
    # Gap-only columns
    gap_only_positions = []
    for position in range(alignment_length):
        column = alignment[:, position]
        if all(residue.upper() in GAP_CHARACTERS for residue in column):
            gap_only_positions.append(position + 1)
    
    if gap_only_positions:
        results.append(
            "WARNING: Gap-only alignment positions detected: "
            + ", ".join(map(str, gap_only_positions[:10]))  # Limit display
        )
    else:
        results.append("PASS: No gap-only columns detected.")
    
    # Sequence type info
    results.append(f"INFO: Sequence type used for analysis: {sequence_type}")
    
    # Valid residue content
    low_quality_sequences = []
    valid_characters = AMINO_ACIDS if sequence_type == "protein" else NUCLEOTIDES
    for record in alignment:
        sequence = clean_sequence_string(str(record.seq)).upper()
        valid_count = sum(character in valid_characters for character in sequence)
        if len(sequence) > 0:
            valid_percentage = (valid_count / len(sequence)) * 100
            if valid_percentage < MIN_VALID_RESIDUE_PERCENTAGE:
                low_quality_sequences.append((record.id, valid_percentage))
    
    if low_quality_sequences:
        results.append(
            f"WARNING: Sequences with less than {MIN_VALID_RESIDUE_PERCENTAGE}% valid residues:"
        )
        for sequence_id, percentage in low_quality_sequences:
            results.append(f"    {sequence_id}: {percentage:.2f}% valid")
    else:
        results.append("PASS: All sequences contain sufficient valid residues.")
    
    # Overall result
    if validation_passed:
        results.append("\nOVERALL VALIDATION: PASSED ✓")
    else:
        results.append("\nOVERALL VALIDATION: FAILED ✗")
    
    return results, validation_passed


def save_validation_report(validation_results, output_dir):
    """Save validation report to file."""
    output_file = os.path.join(output_dir, "alignment_validation_report.txt")
    with open(output_file, "w", encoding="utf-8") as file:
        file.write("CLUSTAL ALIGNMENT VALIDATION REPORT\n")
        file.write("===================================\n\n")
        for result in validation_results:
            file.write(result + "\n")
    print(f"Validation report saved: {output_file}")
    return output_file


def calculate_shannon_entropy(column, sequence_type):
    """Calculate Shannon entropy for a single alignment column."""
    valid_residues = []
    for residue in column:
        residue = clean_sequence_string(residue).upper()
        if residue in GAP_CHARACTERS:
            continue
        if sequence_type == "protein":
            if residue in AMINO_ACIDS:
                valid_residues.append(residue)
        else:
            if residue in NUCLEOTIDES:
                valid_residues.append(residue)
    
    if not valid_residues:
        return np.nan, Counter(), {}
    
    counts = Counter(valid_residues)
    total = len(valid_residues)
    frequencies = {residue: count / total for residue, count in counts.items()}
    
    entropy = 0.0
    for probability in frequencies.values():
        entropy -= probability * np.log2(probability)
    
    return entropy, counts, frequencies


def analyze_alignment(alignment, sequence_type):
    """Calculate entropy for all positions in alignment."""
    results = []
    alignment_length = alignment.get_alignment_length()
    total_sequences = len(alignment)
    
    print("\nCalculating Shannon entropy...")
    
    for position in range(alignment_length):
        column = alignment[:, position]
        entropy, counts, frequencies = calculate_shannon_entropy(column, sequence_type)
        
        valid_sequences = sum(counts.values())
        gap_count = total_sequences - valid_sequences
        gap_percentage = (gap_count / total_sequences) * 100 if total_sequences > 0 else 0
        
        if counts:
            most_common_residue, most_common_count = counts.most_common(1)[0]
            conservation = most_common_count / valid_sequences if valid_sequences > 0 else np.nan
        else:
            most_common_residue = "-"
            conservation = np.nan
        
        maximum_entropy = np.log2(20) if sequence_type == "protein" else np.log2(4)
        relative_entropy = entropy / maximum_entropy if not np.isnan(entropy) else np.nan
        
        results.append({
            "Position": position + 1,
            "Entropy": entropy,
            "Relative_Entropy": relative_entropy,
            "Conservation": conservation,
            "Most_Common_Residue": most_common_residue,
            "Valid_Sequences": valid_sequences,
            "Gap_Count": gap_count,
            "Gap_Percentage": gap_percentage
        })
    
    return pd.DataFrame(results)


def classify_positions(df):
    """Classify positions by variability level."""
    classifications = []
    for entropy in df["Entropy"]:
        if pd.isna(entropy):
            classifications.append("No_Data")
        elif entropy <= CONSERVATION_THRESHOLD:
            classifications.append("Conserved")
        elif entropy >= VARIABLE_THRESHOLD:
            classifications.append("Highly_Variable")
        else:
            classifications.append("Moderately_Variable")
    df["Classification"] = classifications
    return df


def save_csv(df, output_dir):
    """Save results to CSV."""
    output_file = os.path.join(output_dir, "shannon_entropy_results.csv")
    df.to_csv(output_file, index=False)
    print(f"CSV saved: {output_file}")
    return output_file


def save_excel(df, output_dir):
    """Save results to Excel with multiple sheets."""
    output_file = os.path.join(output_dir, "shannon_entropy_results.xlsx")
    try:
        with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="Entropy_Analysis", index=False)
            df[df["Classification"] == "Conserved"].to_excel(
                writer, sheet_name="Conserved", index=False
            )
            df[df["Classification"] == "Moderately_Variable"].to_excel(
                writer, sheet_name="Moderately_Variable", index=False
            )
            df[df["Classification"] == "Highly_Variable"].to_excel(
                writer, sheet_name="Highly_Variable", index=False
            )
        print(f"Excel saved: {output_file}")
    except Exception as e:
        print(f"Warning: Excel export failed: {e}")
    return output_file


def save_residue_lists(df, output_dir):
    """Save conserved and variable residue lists."""
    conserved_file = os.path.join(output_dir, "conserved_residues.csv")
    variable_file = os.path.join(output_dir, "highly_variable_residues.csv")
    
    df[df["Classification"] == "Conserved"].to_csv(conserved_file, index=False)
    df[df["Classification"] == "Highly_Variable"].to_csv(variable_file, index=False)
    
    print(f"Conserved residues saved: {conserved_file}")
    print(f"Highly variable residues saved: {variable_file}")


def create_entropy_plot(df, output_dir):
    """Create basic entropy plot."""
    output_file = os.path.join(output_dir, "shannon_entropy_profile.png")
    
    plt.figure(figsize=(14, 6))
    plt.plot(df["Position"], df["Entropy"], linewidth=1.5)
    plt.axhline(CONSERVATION_THRESHOLD, linestyle="--", linewidth=1, color='green', alpha=0.7)
    plt.axhline(VARIABLE_THRESHOLD, linestyle="--", linewidth=1, color='red', alpha=0.7)
    plt.xlabel("Alignment Position", fontsize=12)
    plt.ylabel("Shannon Entropy (bits)", fontsize=12)
    plt.title("Shannon Entropy Profile", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"Entropy plot saved: {output_file}")
    return output_file


def create_publication_figure(df, output_dir):
    """Create publication-quality figure."""
    output_file = os.path.join(output_dir, "publication_shannon_entropy.png")
    
    plt.figure(figsize=(15, 7))
    plt.plot(df["Position"], df["Entropy"], linewidth=1.5, label="Shannon entropy", color='#2E86AB')
    plt.axhline(CONSERVATION_THRESHOLD, linestyle="--", linewidth=1, label="Conserved threshold (0.5)", color='green')
    plt.axhline(VARIABLE_THRESHOLD, linestyle="--", linewidth=1, label="Variable threshold (1.5)", color='red')
    
    variable = df[df["Classification"] == "Highly_Variable"]
    if not variable.empty:
        plt.scatter(variable["Position"], variable["Entropy"], s=20, 
                   label="Highly variable", color='#E94F37', alpha=0.7)
    
    conserved = df[df["Classification"] == "Conserved"]
    if not conserved.empty:
        plt.scatter(conserved["Position"], conserved["Entropy"], s=12, 
                   label="Conserved", color='#06A77D', alpha=0.5)
    
    plt.xlabel("Alignment Position", fontsize=13)
    plt.ylabel("Shannon Entropy (bits)", fontsize=13)
    plt.title("Sequence Variability Based on Shannon Entropy", fontsize=15)
    plt.legend()
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"Publication figure saved: {output_file}")
    return output_file


def generate_summary(df, alignment, sequence_type, output_dir):
    """Generate analysis summary."""
    output_file = os.path.join(output_dir, "entropy_summary.txt")
    
    valid_entropy = df["Entropy"].dropna()
    mean_entropy = valid_entropy.mean() if not valid_entropy.empty else np.nan
    min_entropy = valid_entropy.min() if not valid_entropy.empty else np.nan
    max_entropy = valid_entropy.max() if not valid_entropy.empty else np.nan
    
    with open(output_file, "w", encoding="utf-8") as file:
        file.write("SHANNON ENTROPY ANALYSIS SUMMARY\n")
        file.write("================================\n\n")
        file.write(f"Number of sequences: {len(alignment)}\n")
        file.write(f"Alignment length: {alignment.get_alignment_length()}\n")
        file.write(f"Sequence type: {sequence_type}\n\n")
        file.write(f"Mean entropy: {mean_entropy:.4f}\n")
        file.write(f"Minimum entropy: {min_entropy:.4f}\n")
        file.write(f"Maximum entropy: {max_entropy:.4f}\n\n")
        file.write(f"Conserved positions: {sum(df['Classification'] == 'Conserved')}\n")
        file.write(f"Moderately variable positions: {sum(df['Classification'] == 'Moderately_Variable')}\n")
        file.write(f"Highly variable positions: {sum(df['Classification'] == 'Highly_Variable')}\n\n")
        file.write("Thresholds used:\n")
        file.write(f"Conserved <= {CONSERVATION_THRESHOLD}\n")
        file.write(f"Highly variable >= {VARIABLE_THRESHOLD}\n")
    
    print(f"Summary saved: {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Shannon entropy from a Clustal Omega alignment with validation."
    )
    parser.add_argument("-i", "--input", required=True, help="Clustal alignment file.")
    parser.add_argument("-o", "--output", default="entropy_results", help="Output directory.")
    parser.add_argument("--type", choices=["auto", "protein", "nucleotide"], 
                       default="auto", help="Sequence type. Default: auto.")
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.input):
        raise SystemExit(f"\nERROR: Input file not found:\n{args.input}")
    
    os.makedirs(args.output, exist_ok=True)
    
    alignment = read_clustal_alignment(args.input)
    
    if args.type == "auto":
        sequence_type = detect_sequence_type(alignment)
    else:
        sequence_type = args.type
    
    print(f"\nDetected sequence type: {sequence_type}")
    print("\n" + "=" * 55)
    print("ALIGNMENT VALIDATION")
    print("=" * 55)
    
    validation_results, validation_passed = validate_alignment(
        alignment, args.input, sequence_type
    )
    
    for result in validation_results:
        print(result)
    
    save_validation_report(validation_results, args.output)
    
    if not validation_passed:
        print("\nERROR: Alignment validation failed.")
        print("Entropy analysis was NOT performed.")
        print("Review alignment_validation_report.txt")
        raise SystemExit(1)
    
    results = analyze_alignment(alignment, sequence_type)
    results = classify_positions(results)
    
    save_csv(results, args.output)
    save_excel(results, args.output)
    save_residue_lists(results, args.output)
    create_entropy_plot(results, args.output)
    create_publication_figure(results, args.output)
    generate_summary(results, alignment, sequence_type, args.output)
    
    print("\n" + "=" * 55)
    print("SHANNON ENTROPY ANALYSIS COMPLETED")
    print("=" * 55)
    print(f"\nResults directory:\n{os.path.abspath(args.output)}")
    print("\nGenerated files:")
    print("1. alignment_validation_report.txt")
    print("2. shannon_entropy_results.csv")
    print("3. shannon_entropy_results.xlsx")
    print("4. conserved_residues.csv")
    print("5. highly_variable_residues.csv")
    print("6. shannon_entropy_profile.png")
    print("7. publication_shannon_entropy.png")
    print("8. entropy_summary.txt")


if __name__ == "__main__":
    main()