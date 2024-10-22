import hashlib
from parse import parse_fasta
from palindromic import find_palindromes_bruteforce, find_sliding_window_palindromes, find_palindromes_dict


def compute_checksum(palindromes) -> str:
    """
    Compute the MD5 checksum for a list or dictionary of palindromes.

    Args:
        palindromes (list or dict): Palindromes found by the function.
    
    Returns:
        str: MD5 checksum of the palindromes.
    """
    if isinstance(palindromes, dict):
        combined_palindromes = ''.join(palindromes.values())  # If its a dic
    else:
        combined_palindromes = ''.join(palindromes)

    return hashlib.md5(combined_palindromes.encode("utf-8")).hexdigest()

def benchmark_palindrome_finding_with_checksum(seq_dict: dict, length: int) -> None:
    """
    Benchmark the three palindrome-finding functions and compute their checksums.

    Args:
        seq_dict (dict): Dictionary of sequence IDs and sequences.
        length (int): Length of palindromes to find.
    """
    print("Computing checksum for Bruteforce method...")
    bruteforce_palindromes = []
    for seq_id, seq in seq_dict.items():
        bruteforce_palindromes.extend(find_palindromes_bruteforce(seq, length))
    bruteforce_checksum = compute_checksum(bruteforce_palindromes)
    print(f"Bruteforce checksum: {bruteforce_checksum}")

    print("Computing checksum for Sliding window method...")
    sliding_window_palindromes = []
    for seq_id, seq in seq_dict.items():
        sliding_window_palindromes.extend(find_sliding_window_palindromes(seq, length))
    sliding_window_checksum = compute_checksum(sliding_window_palindromes)
    print(f"Sliding window checksum: {sliding_window_checksum}")

    print("Computing checksum for Dictionary method...")
    dict_palindromes = {}
    for seq_id, seq in seq_dict.items():
        dict_palindromes.update(find_palindromes_dict(seq, length))
    dict_checksum = compute_checksum(dict_palindromes)
    print(f"Dictionary checksum: {dict_checksum}")

    print("\nComparing checksums:")
    if bruteforce_checksum == sliding_window_checksum == dict_checksum:
        print("All methods produced the same results.")
    else:
        print("Methods produced different results.")

if __name__ == "__main__":
    fasta_file = "test.fasta"
    seq_dict = parse_fasta(fasta_file)
    
    palindrome_length = 6
    
    benchmark_palindrome_finding_with_checksum(seq_dict, palindrome_length)




# import time
# from palindromic import find_palindromes_bruteforce, find_sliding_window_palindromes, find_palindromes_dict
# from parse import parse_fasta

# def benchmark_palindrome_search(dna_sequences: dict, length: int):
#     for seq_id, dna_sequence in dna_sequences.items():
#         print(f"\nBenchmarking for sequence ID: {seq_id}")
        
#         # Brute-force
#         start = time.time()
#         result_bruteforce = find_palindromes_bruteforce(dna_sequence, length)
#         time_bruteforce = time.time() - start

#         # Sliding window
#         start = time.time()
#         result_sliding = find_sliding_window_palindromes(dna_sequence, length)
#         time_sliding = time.time() - start

#         # Dictionary-based
#         start = time.time()
#         result_dict = find_palindromes_dict(dna_sequence, length)
#         time_dict = time.time() - start

#         # benchmark results for each sequence
#         print(f"Brute-force time: {time_bruteforce:.4f}s, Palindromes found: {len(result_bruteforce)}, Palindromes: {result_bruteforce}")
#         print(f"Sliding window time: {time_sliding:.4f}s, Palindromes found: {len(result_sliding)}, Palindromes: {result_sliding}")
#         print(f"Dictionary-based time: {time_dict:.4f}s, Palindromes found: {len(result_dict)}, Palindromes: {result_dict}")


# dna_seq = parse_fasta('test.fasta')
# benchmark_palindrome_search(dna_seq, 6)