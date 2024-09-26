import sys

# Define stop codons
STOP_CODONS = ["TAA", "TAG", "TGA"]

def find_genes(sequence):
    """Find genes between start (ATG) and stop codons in three reading frames."""
    sequence = sequence.upper()  # Convert the sequence to uppercase
    for frame in range(3):
        print(f"\nReading frame {frame + 1}:")
        start = -1
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i + 3]
            if codon == "ATG":  # Start codon
                start = i
            elif codon in STOP_CODONS and start != -1:  # Stop codon
                gene = sequence[start:i + 3]
                print(f"Gene found from position {start} to {i + 3}: {gene}")
                start = -1  # Reset to find the next gene

def read_fasta(filename):
    """Read a FASTA (or .fna) file and return the sequence."""
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(">"):  # Ignore header lines starting with '>'
                sequence += line.strip()
    return sequence

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    # Check if the file has a proper extension
    if not fasta_file.endswith(('.fna', '.fasta', '.fa')):
        print("Error: Please provide a valid FASTA file (.fna, .fasta, .fa)")
        sys.exit(1)
    
    # Read the sequence and find genes
    genome_sequence = read_fasta(fasta_file)
    find_genes(genome_sequence)

if __name__ == "__main__":
    main()


