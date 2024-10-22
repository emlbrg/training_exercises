from collections import defaultdict
from collections import Counter
from typing import Dict
from parse import parse_fasta
from decorators import time_it

@time_it
def naive_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Naive sliding window approach to count k-mers.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts.
    """
    kmer_counts = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    return kmer_counts

@time_it
def hashmap_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Hash map-based k-mer counting.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts.
    """
    kmer_counts = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

@time_it
def optimized_counter_kmer_count(sequence: str, k: int) -> Dict[str, int]:
    """
    Optimized k-mer counting using Counter.
    
    Args:
        sequence (str): The DNA sequence.
        k (int): Length of the k-mers.
    
    Returns:
        Dict[str, int]: Dictionary of k-mer counts.
    """
    kmers = (sequence[i:i+k] for i in range(len(sequence) - k + 1))
    return Counter(kmers)


sequences = parse_fasta("test.fasta")
k = 3
for seq_id, sequence in sequences.items():
    print(f"Processing sequence {seq_id}")
    naive_output = naive_kmer_count(sequence, k)
    hashmap_output = hashmap_kmer_count(sequence, k)
    optimized_output = optimized_counter_kmer_count(sequence, k)

    print(optimized_output)