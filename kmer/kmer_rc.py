from collections import defaultdict, Counter
from typing import Dict
from parse import parse_fasta
from decorators import time_it

# Helper function to calculate reverse complement
def reverse_complement(sequence: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

@time_it
def naive_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Naive sliding window approach to count k-mers, including reverse complement.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts (forward and reverse complement).
    """
    kmer_counts = {}
    reverse_seq = reverse_complement(sequence)
    
    # Count k-mers in forward sequence
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    
    # Count k-mers in reverse complement sequence
    for i in range(len(reverse_seq) - k + 1):
        kmer = reverse_seq[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
            
    return kmer_counts

@time_it
def hashmap_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Hash map-based k-mer counting, including reverse complement.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts (forward and reverse complement).
    """
    kmer_counts = defaultdict(int)
    reverse_seq = reverse_complement(sequence)
    
    # Count k-mers in forward sequence
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] += 1

    # Count k-mers in reverse complement sequence
    for i in range(len(reverse_seq) - k + 1):
        kmer = reverse_seq[i:i+k]
        kmer_counts[kmer] += 1
    
    return kmer_counts

@time_it
def optimized_counter_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Optimized k-mer counting using Counter, including reverse complement.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts (forward and reverse complement).
    """
    kmers = (sequence[i:i+k] for i in range(len(sequence) - k + 1))
    reverse_seq = reverse_complement(sequence)
    reverse_kmers = (reverse_seq[i:i+k] for i in range(len(reverse_seq) - k + 1))
    
    return Counter(kmers) + Counter(reverse_kmers)

# Example usage
sequences = parse_fasta("test.fasta")
k = 3
for seq_id, sequence in sequences.items():
    print(f"Processing sequence {seq_id}")
    naive_output = naive_kmer_count(sequence, k)
    hashmap_output = hashmap_kmer_count(sequence, k)
    optimized_output = optimized_counter_kmer_count(sequence, k)

    print(optimized_output)
