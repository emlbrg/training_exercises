import random

def generate_dna_sequence(length, repeat_unit='AT', min_repeats=6, max_repeats=30):
    """Generate a random DNA sequence with specified tandem repeats."""
    sequence = []
    while len(sequence) < length:
        # Randomly decide to insert a tandem repeat or a random base
        if random.random() < 0.3:  # 30% chance to insert a tandem repeat
            repeats = random.randint(min_repeats, max_repeats)  # Random repeat count
            sequence.append(repeat_unit * repeats)
        else:
            bases = ['A', 'T', 'C', 'G']
            sequence.append(''.join(random.choices(bases, k=100)))  # Add 100 random bases

    return ''.join(sequence)[:length]  # Ensure the length is exactly the specified

fasta_file_path = "test.fasta"
with open(fasta_file_path, 'w') as fasta_file:
    for i in range(1, 1111):
        seq_id = f">seq{i}"
        sequence = generate_dna_sequence(1000)
        fasta_file.write(f"{seq_id}\n{sequence}\n")

fasta_file_path
