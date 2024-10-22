from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import numpy as np
from decorators import time_it
import time
import re

### TASK 1 ###

@time_it
def parse_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Sequence header
                current_sequence = line[1:]  # Remove '>'
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line  # Append the sequence
    return sequences

@time_it
def parse_fasta_biopython(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences[record.id] = str(record.seq)
    return sequences


fasta_file = 'antibody_sequences_updated.fasta'

print('#### Approach 1: Using Basic String Manipulation ####')
sequences = parse_fasta(fasta_file)
for seq_id, sequence in sequences.items():
    print(f'{seq_id}: {sequence}')

print('#### Approach 2: Using Biopython ####')
sequences_biopython = parse_fasta_biopython(fasta_file)
for seq_id, sequence in sequences_biopython.items():
    print(f'{seq_id}: {sequence}')

### TASK 2 ###

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) * 100

def sequence_statistics(sequences):
    for seq_id, sequence in sequences.items():
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)
        print(f'Sequence ID: {seq_id}')
        print(f'  Length: {length}')
        print(f'  GC Content: {gc_content:.2f}%')
        print('----')


def calculate_gc_numpy(sequence):
    seq_array = np.array(list(sequence))
    gc_count = np.sum((seq_array == 'G') | (seq_array == 'C'))
    return gc_count / len(seq_array) * 100

def sequence_statistics_numpy(sequences):
    stats = {}
    for seq_id, sequence in sequences.items():
        length = len(sequence)
        gc_content = calculate_gc_numpy(sequence)
        stats[seq_id] = {'length': length, 'gc_content': gc_content}
    return stats

print('#### Approach 1: Using string methods ####')
sequence_statistics(sequences)
print('#### Approach 2: Using Numpy ####')
stats_numpy = sequence_statistics_numpy(sequences)
print(stats_numpy)




def find_orfs_biopython(sequences_dict, min_len=3):
    """
    Find all ORFs in both forward and reverse strands using manual codon checking,
    similar to the 'find_orfs_basic' approach, while leveraging Biopython for sequence handling.
    
    Args:
        sequences_dict (dict): Dictionary with record IDs as keys and sequences as values.
        min_len (int): Minimum length of the ORF to be considered valid (default is 3).
    
    Returns:
        list: List of ORFs, with details on record ID, strand, frame, start, end, and sequence.
    """
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    all_orfs = []
    
    for record_id, sequence in sequences_dict.items():
        seq = Seq(sequence)
        
        # Forward strand ORFs (+1)
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                codon = str(seq[i:i+3])
                if codon == start_codon:
                    for j in range(i + 3, len(seq), 3):
                        stop_codon = str(seq[j:j+3])
                        if stop_codon in stop_codons:
                            orf = seq[i:j+3]
                            if len(orf) >= min_len:
                                all_orfs.append((record_id, +1, frame, i, j + 3, str(orf)))
                            break
        
        # Reverse strand ORFs (-1)
        reverse_seq = seq.reverse_complement()
        for frame in range(3):
            for i in range(frame, len(reverse_seq), 3):
                codon = str(reverse_seq[i:i+3])
                if codon == start_codon:
                    for j in range(i + 3, len(reverse_seq), 3):
                        stop_codon = str(reverse_seq[j:j+3])
                        if stop_codon in stop_codons:
                            orf = reverse_seq[i:j+3]
                            if len(orf) >= min_len:
                                all_orfs.append((record_id, -1, frame, i, j + 3, str(orf)))
                            break
    
    return all_orfs


def find_orfs_basic(sequences_biopython, min_len=3):
    """
    Find all ORFs in both forward and reverse strands by manually checking codons.
    
    Args:
        sequences_biopython (dict): Dictionary with record IDs as keys and sequences as values.
        min_len (int): Minimum length of the ORF to be considered valid (default is 3).
    
    Returns:
        list: List of ORFs with details on record ID, strand, frame, start, end, and sequence.
    """
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    all_orfs = []
    seen_orfs = set()  # Track unique ORFs based on their positions

    for record_id, sequence in sequences_biopython.items():
        # Forward strand ORFs (+1)
        for frame in range(3):
            for i in range(frame, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon == start_codon:
                    for j in range(i + 3, len(sequence), 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf = sequence[i:j+3]
                            orf_tuple = (record_id, +1, frame, i, j + 3)  # Use just the positions for uniqueness
                            if len(orf) >= min_len and orf_tuple not in seen_orfs:
                                all_orfs.append((record_id, +1, frame, i, j + 3, orf))
                                seen_orfs.add(orf_tuple)  # Mark this ORF as seen
                            break

        # Reverse strand ORFs (-1)
        reverse_seq = str(Seq(sequence).reverse_complement())
        for frame in range(3):
            for i in range(frame, len(reverse_seq), 3):
                codon = reverse_seq[i:i+3]
                if codon == start_codon:
                    for j in range(i + 3, len(reverse_seq), 3):
                        stop_codon = reverse_seq[j:j+3]
                        if stop_codon in stop_codons:
                            orf = reverse_seq[i:j+3]
                            orf_tuple = (record_id, -1, frame, i, j + 3)  # Use just the positions for uniqueness
                            if len(orf) >= min_len and orf_tuple not in seen_orfs:
                                all_orfs.append((record_id, -1, frame, i, j + 3, orf))
                                seen_orfs.add(orf_tuple)  # Mark this ORF as seen
                            break

    return all_orfs


def print_longest_orf(orfs):  # Function to find and print the longest ORF for each sequence
    for record_id in set(orf[0] for orf in orfs):
        record_orfs = [orf for orf in orfs if orf[0] == record_id]
        if record_orfs:
            # Ensure we are accessing the correct part of the tuple (ORF sequence at index 5)
            longest_orf = max(record_orfs, key=lambda x: len(x[5]) if isinstance(x[5], str) else 0)
            translated = Seq(longest_orf[5]).translate()
            print(f"Record ID: {record_id}, Longest ORF: {longest_orf}, Translated: {translated} ")


### TESTING TASK 3 ###

# Run all three ORF-finding functions
fasta_file = fasta_file = 'mrgd_qf_fasta-2/test.fasta'
sequences_biopython = parse_fasta_biopython(fasta_file)

orfs_basic = find_orfs_basic(sequences_biopython)
orfs_biopython = find_orfs_biopython(sequences_biopython)
# orfs_regex = find_all_orfs_with_regex(sequences_biopython)


# Run some more detailed comparison if the methods find different ORFs
def compare_orfs(orfs_basic, orfs_biopython, orfs_regex):
    # Convert ORF lists to sets for easy comparison
    orfs_basic_set = set(orfs_basic)
    orfs_biopython_set = set(orfs_biopython)
    # orfs_regex_set = set(orfs_regex)

    # Find ORFs that are in basic but not in biopython or regex
    basic_not_in_biopython = orfs_basic_set - orfs_biopython_set
    # basic_not_in_regex = orfs_basic_set - orfs_regex_set

    # Find ORFs that are in biopython but not in basic or regex
    biopython_not_in_basic = orfs_biopython_set - orfs_basic_set
    # biopython_not_in_regex = orfs_biopython_set - orfs_regex_set

    # Find ORFs that are in regex but not in basic or biopython
    # regex_not_in_basic = orfs_regex_set - orfs_basic_set
    # regex_not_in_biopython = orfs_regex_set - orfs_biopython_set

    # Print differences
    print("ORFs in find_orfs_basic but not in find_orfs_biopython:")
    for orf in basic_not_in_biopython:
        print(orf)

    # print("\nORFs in find_orfs_basic but not in find_all_orfs_with_regex:")
    # for orf in basic_not_in_regex:
    #     print(orf)

    print("\nORFs in find_orfs_biopython but not in find_orfs_basic:")
    for orf in biopython_not_in_basic:
        print(orf)

    # print("\nORFs in find_orfs_biopython but not in find_all_orfs_with_regex:")
    # for orf in biopython_not_in_regex:
    #     print(orf)

    # print("\nORFs in find_all_orfs_with_regex but not in find_orfs_basic:")
    # for orf in regex_not_in_basic:
    #     print(orf)

    # print("\nORFs in find_all_orfs_with_regex but not in find_orfs_biopython:")
    # for orf in regex_not_in_biopython:
    #     print(orf)

# Compare and print the differences between methods
compare_orfs(orfs_basic, orfs_biopython)

# Print the longest ORF found by each method
print("Longest ORF from find_orfs_basic:")
print_longest_orf(orfs_basic)

print("Longest ORF from find_orfs_biopython:")
print_longest_orf(orfs_biopython)

# print("Longest ORF from find_all_orfs_with_regex:")
# print_longest_orf(orfs_regex)

# Check for differences
if orfs_basic == orfs_biopython:
    print("All methods found the same ORFs!")
else:
    print("Methods found different ORFs.")