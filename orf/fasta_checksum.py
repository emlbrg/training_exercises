import hashlib
from Bio import SeqIO
from fasta_task import parse_fasta, parse_fasta_biopython, find_orfs_basic, find_orfs_biopython, find_all_orfs_with_regex
from Bio.Seq import Seq

# Function to calculate the checksum of parsed sequences
def calculate_sequences_checksum(sequences, algorithm='md5'):
    """
    Calculate a checksum for a dictionary of sequences.

    Args:
        sequences (dict): Dictionary with sequence record IDs as keys and sequences as values.
        algorithm (str): Hash algorithm to use (default is 'md5').

    Returns:
        str: The computed checksum.
    """
    hash_func = hashlib.new(algorithm)
    for seq in sequences.values():
        hash_func.update(seq.encode())
    return hash_func.hexdigest()

# Function to calculate the checksum for ORFs
def calculate_orfs_checksum(orfs, algorithm='md5'):
    """
    Calculate a checksum for a list of ORFs.

    Args:
        orfs (list): List of tuples representing ORFs.
        algorithm (str): Hash algorithm to use (default is 'md5').

    Returns:
        str: The computed checksum.
    """
    hash_func = hashlib.new(algorithm)
    for orf in orfs:
        # Convert the ORF tuple into a string and update the hash
        orf_str = ''.join(map(str, orf))
        hash_func.update(orf_str.encode())
    return hash_func.hexdigest()

# Existing code for parsing sequences
fasta_file = 'antibody_sequences_updated.fasta'

# Parse using the first function
sequences_1 = parse_fasta(fasta_file)

# Parse using the second function (Biopython-based)
sequences_2 = parse_fasta_biopython(fasta_file)

# Calculate checksums for both parsing functions
checksum_1 = calculate_sequences_checksum(sequences_1, algorithm='md5')
checksum_2 = calculate_sequences_checksum(sequences_2, algorithm='md5')

print(f"Checksum for parse_fasta: {checksum_1}")
print(f"Checksum for parse_fasta_biopython: {checksum_2}")

# Compare checksums of parsing results
if checksum_1 == checksum_2:
    print("The outputs of both parsing functions are identical!")
else:
    print("The outputs of the parsing functions differ.")

# Now find ORFs for both methods and calculate their checksums
orfs_basic = find_orfs_basic(sequences_1)  # Using find_orfs_basic on the first parsed sequences
orfs_biopython = find_orfs_biopython(sequences_2)  # Using find_orfs_biopython on the second parsed sequences

# Calculate checksums for ORFs from both methods
checksum_orfs_basic = calculate_orfs_checksum(orfs_basic, algorithm='md5')
checksum_orfs_biopython = calculate_orfs_checksum(orfs_biopython, algorithm='md5')

print(f"Checksum for find_orfs_basic: {checksum_orfs_basic}")
print(f"Checksum for find_orfs_biopython: {checksum_orfs_biopython}")

# After calculating checksums, add this code to compare the outputs
if checksum_orfs_basic == checksum_orfs_biopython:
    print("The outputs of both ORF finding functions are identical!")
else:
    print("The outputs of the ORF finding functions differ.")

    # Convert ORF lists to dictionaries for easier comparison
    orfs_dict_basic = {(orf[0], orf[2], orf[3]): orf for orf in orfs_basic}  # (Record ID, Frame, Start) as key
    orfs_dict_biopython = {(orf[0], orf[2], orf[3]): orf for orf in orfs_biopython}  # (Record ID, Frame, Start) as key

    # Find unique ORFs in both dictionaries
    unique_to_basic = {key: orfs_dict_basic[key] for key in orfs_dict_basic if key not in orfs_dict_biopython}
    unique_to_biopython = {key: orfs_dict_biopython[key] for key in orfs_dict_biopython if key not in orfs_dict_basic}

    # Print unique ORFs found in each method with Record ID
    if unique_to_basic:
        print("Unique ORFs found in find_orfs_basic:")
        for key, orf in unique_to_basic.items():
            record_id, strand, frame, start, end, orf_seq = orf
            print(f"Record ID: {record_id}, Strand: {strand}, Frame: {frame}, Start: {start}, End: {end}, ORF: {orf_seq}")

    if unique_to_biopython:
        print("Unique ORFs found in find_orfs_biopython:")
        for key, orf in unique_to_biopython.items():
            record_id, strand, frame, start, end, orf_seq = orf
            print(f"Record ID: {record_id}, Strand: {strand}, Frame: {frame}, Start: {start}, End: {end}, ORF: {orf_seq}")

    # # Optional: Print the entire lists for manual comparison
    # print("\nORFs from find_orfs_basic:")
    # for orf in orfs_basic:
    #     print(orf)

    # print("\nORFs from find_orfs_biopython:")
    # for orf in orfs_biopython:
    #     print(orf)

