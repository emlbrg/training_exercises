import hashlib
import json
from Bio import SeqIO
from parse import parse_fasta
from kmer import naive_kmer_count, hashmap_kmer_count, optimized_counter_kmer_count

def calculate_md5_checksum(kmer_dict: dict) -> str:
    """
    Calculate the MD5 checksum of the k-mer counts dictionary.
    
    Args:
        kmer_dict (dict): Dictionary of k-mer counts.
    
    Returns:
        str: The MD5 checksum as a hexadecimal string.
    """
    # Convert the dictionary to a JSON string and sort it to ensure consistent ordering
    kmer_str = json.dumps(kmer_dict, sort_keys=True)
    checksum = hashlib.md5(kmer_str.encode('utf-8')).hexdigest()
    
    return checksum

sequences = parse_fasta("test.fasta")
k = 3
all_results_same = True
md5_checksums = [] 

for seq_id, sequence in sequences.items():
    print(f"Processing sequence {seq_id}")
    naive_output = naive_kmer_count(sequence, k)
    hashmap_output = hashmap_kmer_count(sequence, k)
    optimized_output = optimized_counter_kmer_count(sequence, k)

    # checksums
    naive_checksum = calculate_md5_checksum(naive_output)
    hashmap_checksum = calculate_md5_checksum(hashmap_output)
    optimized_checksum = calculate_md5_checksum(optimized_output)

    md5_checksums.append((seq_id, naive_checksum, hashmap_checksum, optimized_checksum))
    if not (naive_checksum == hashmap_checksum == optimized_checksum):
            all_results_same = False


for seq_id, naive_checksum, hashmap_checksum, optimized_checksum in md5_checksums:
    print(f"\nSequence ID: {seq_id}")
    print(f"Naive MD5 Checksum: {naive_checksum}")
    print(f"Hashmap MD5 Checksum: {hashmap_checksum}")
    print(f"Optimized MD5 Checksum: {optimized_checksum}")

if all_results_same:
    print("All implementations produce identical results for all sequences.")
else:
    print("Results are different.")


# Compare to gold standard 
def parse_fasta_rhetorical(fasta_file: str) -> dict:
    """
    Parse a FASTA file and create a dictionary where sequences are the keys
    and sequence IDs (seq_id) are the values.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with sequences as keys and seq_id as values.
    """
    seq_dict = {}
    with open(fasta_file, 'r') as file:
        seq_id = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    seq_dict[sequence] = seq_id
                seq_id = line[1:]  # Remove the '>' and take the sequence ID
                sequence = ""  # Reset sequence
            else:
                sequence += line
        if sequence:
            seq_dict[sequence] = seq_id

    return seq_dict

fasta_file = "jellyfish_kmers.fasta"
jellyfish_output = parse_fasta_rhetorical(fasta_file)
print(jellyfish_output)
jellyfish_checksum = calculate_md5_checksum(jellyfish_output)

# Compare the Jellyfish results with Python k-mer counting implementations
for seq_id, naive_checksum, hashmap_checksum, optimized_checksum in md5_checksums:
    print(f"\nSequence ID: {seq_id}")
    print(f"Naive MD5 Checksum: {naive_checksum}")
    print(f"Hashmap MD5 Checksum: {hashmap_checksum}")
    print(f"Optimized MD5 Checksum: {optimized_checksum}")
    print(f"Jellyfish MD5 Checksum: {jellyfish_checksum}")

    # Compare against Jellyfish
    if naive_checksum == hashmap_checksum == optimized_checksum == jellyfish_checksum:
        print(f"All implementations and Jellyfish produce identical results for sequence {seq_id}!")
    else:
        print(f"Mismatch in results for sequence {seq_id} compared to Jellyfish.")
