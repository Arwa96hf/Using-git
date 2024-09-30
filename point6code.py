from argparse import ArgumentParser
from Bio import SeqIO

# Genetic code dictionary to translate codons into amino acids
genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  # Isoleucine and Methionine
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  # Threonine
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  # Asparagine, Lysine
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  # Serine, Arginine
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  # Leucine
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  # Proline
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  # Histidine, Glutamine
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  # Arginine
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  # Valine
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  # Alanine
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  # Aspartic acid, Glutamic acid
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  # Glycine
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  # Serine
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  # Phenylalanine, Leucine
    'TAC':'Y', 'TAT':'Y',                        # Tyrosine
    'TGC':'C', 'TGT':'C',                        # Cysteine (C)
    'TGG':'W',                                   # Tryptophan
    'TAA':'_', 'TAG':'_', 'TGA':'_'              # Stop codons
}

# Function to translate a nucleotide sequence into an amino acid sequence
def translate_dna_to_protein(dna_sequence):
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        amino_acid = genetic_code.get(codon, '?')  # '?' for unknown codons
        if amino_acid == '_':  # Stop codon
            break
        protein_sequence.append(amino_acid)
    return ''.join(protein_sequence)

# Function to parse a FASTA file
def parse_fasta_file(file_path):
    """Read sequences from a FASTA file."""
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Function to find ORFs with RBS and apply ORF length filter
def find_orfs(sequence, min_length, rbs_length):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    rbs_sequences = ["AGGAGG", "GGAGG", "AGG", "GAG"]  # Common RBS sequences
    orfs = []

    # Check three different reading frames
    for frame in range(3):
        start_positions = []  # To store start codon positions
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i + 3]
            if codon == start_codon:
                start_positions.append(i)  # Found a start codon
            elif codon in stop_codons and start_positions:
                # Found a stop codon; process ORF
                while start_positions:
                    start_pos = start_positions.pop(0)
                    orf = sequence[start_pos:i + 3]
                    orf_length = (i + 3 - start_pos) // 3  # ORF length in codons

                    # Check for RBS upstream of the start codon
                    rbs_start = max(0, start_pos - rbs_length)
                    rbs_region = sequence[rbs_start:start_pos]

                    # If the ORF length is sufficient and an RBS is found
                    if orf_length >= min_length and any(rbs in rbs_region for rbs in rbs_sequences):
                        protein_sequence = translate_dna_to_protein(orf)  # Translate ORF to protein
                        orfs.append((frame + 1, start_pos + 1, i + 3, orf, orf_length, protein_sequence))
    
    return orfs

# Main function
def main_task(filename, min_length, rbs_length):
    sequences = parse_fasta_file(filename)
    
    for seq_id, sequence in sequences.items():
        print(f"Sequence ID: {seq_id}")
        orfs = find_orfs(sequence, min_length, rbs_length)
        
        if orfs:
            for orf in orfs:
                print(f"Frame {orf[0]}: Start at {orf[1]}, Stop at {orf[2]}, Length: {orf[4]} codons")
                print(f"DNA Sequence: {orf[3]}")
                print(f"Protein Sequence: {orf[5]}")
                print("-" * 60)
        else:
            print("No valid ORFs found.")

if __name__ == "__main__":
    parser = ArgumentParser(description="Find ORFs in a FASTA file with RBS filtering.")
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    parser.add_argument("-l", "--min-length", type=int, default=100, help="Minimum ORF length in codons (default: 100)")
    parser.add_argument("-r", "--rbs-length", type=int, default=20, help="Length to search for RBS upstream (default: 20)")
    args = parser.parse_args()

    main_task(args.file, args.min_length, args.rbs_length)
