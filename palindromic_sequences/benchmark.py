import time
from parse import parse_fasta
from palindromic import find_palindromes_bruteforce, find_sliding_window_palindromes, find_palindromes_dict

def benchmark_palindrome_finding(seq_dict: dict, length: int) -> None:
    """
    Benchmark three palindrome-finding functions on the provided sequence dictionary.
    
    Args:
        seq_dict (dict): Dictionary of sequence IDs and sequences (pardsed FASTA).
        length (int): Length of palindromes to find.
    """
    # Bruteforce
    start_time = time.time()
    for seq_id, seq in seq_dict.items():
        find_palindromes_bruteforce(seq, length)
    bruteforce_duration = time.time() - start_time
    print(f"Bruteforce method took: {bruteforce_duration:.4f} seconds")

    # Sliding window
    start_time = time.time()
    for seq_id, seq in seq_dict.items():
        find_sliding_window_palindromes(seq, length)
    sliding_window_duration = time.time() - start_time
    print(f"Sliding window method took: {sliding_window_duration:.4f} seconds")

    # Dictionary
    start_time = time.time()
    for seq_id, seq in seq_dict.items():
        find_palindromes_dict(seq, length)
    dict_method_duration = time.time() - start_time
    print(f"Dictionary method took: {dict_method_duration:.4f} seconds")

if __name__ == "__main__":
    fasta_file = "generated_sequences.fasta"
    seq_dict = parse_fasta(fasta_file)
    
    benchmark_palindrome_finding(seq_dict, 6)


















### PREVIOUS VERSIONS ###
# import time
# from palindromic import find_palindromes_bruteforce, sliding_window_palindromes, find_palindromes_dict
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
#         result_sliding = sliding_window_palindromes(dna_sequence, length)
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



# import time

# def benchmark_palindrome_methods(fasta_file: str, length: int):
#     seq_dict = parse_fasta(fasta_file)

#     # Benchmark each method
#     start_time = time.time()
#     for seq in seq_dict.values():
#         find_palindromes_bruteforce(seq, length)
#     bruteforce_time = time.time() - start_time

#     start_time = time.time()
#     for seq in seq_dict.values():
#         sliding_window_palindromes(seq, length)
#     sliding_window_time = time.time() - start_time

#     start_time = time.time()
#     for seq in seq_dict.values():
#         find_palindromes_dict(seq, length)
#     dict_time = time.time() - start_time

#     # Print results
#     print(f"Brute force method took: {bruteforce_time:.4f} seconds")
#     print(f"Sliding window method took: {sliding_window_time:.4f} seconds")
#     print(f"Dictionary method took: {dict_time:.4f} seconds")

# # Example usage:
# benchmark_palindrome_methods("test.fasta", length=6)  # Adjust length as needed
