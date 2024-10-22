import random

def generate_random_dna_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of a given length.

    Args:
        length (int): Length of the DNA sequence.
    
    Returns:
        str: Random DNA sequence.
    """
    return ''.join(random.choices('ATCG', k=length))

def create_fasta_file(num_sequences: int, sequence_length: int, output_file: str) -> None:
    """
    Create a FASTA file with a given number of DNA sequences.

    Args:
        num_sequences (int): Number of sequences to generate.
        sequence_length (int): Length of each DNA sequence.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for i in range(1, num_sequences + 1):
            seq_id = f">sequence_{i}"
            dna_sequence = generate_random_dna_sequence(sequence_length)
            f.write(f"{seq_id}\n")
            f.write(f"{dna_sequence}\n")

# Parameters for the FASTA file
num_sequences = 1000
sequence_length = 10000 * 3  # 1000 codons = 3000 nucleotides
output_fasta = "generated_sequences.fasta"

# Create the FASTA file
create_fasta_file(num_sequences, sequence_length, output_fasta)

output_fasta
