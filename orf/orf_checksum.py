import hashlib
from orffinderA import parse_fasta, find_orfs, find_orf_translate, find_longest_orf

def calculate_checksum(proteins):
    """
    Calculate the checksum for a list of proteins.

    Args:
        proteins (list): List of protein details (ID, strand, frame, start, end, sequence).

    Returns:
        str: Hexadecimal checksum of the concatenated protein sequences.
    """
    concatenated_proteins = ''.join(sorted(str(pro[-1]) for pro in proteins))  # Convert Seq to str
    # print(len(concatenated_proteins))
    return hashlib.md5(concatenated_proteins.encode()).hexdigest()


def compare_results(all_orfs, all_proteins):
    """
    Compare the results of find_orfs and find_orf_translate by checksum.

    Args:
        all_orfs (list): Results from find_orfs.
        all_proteins (list): Results from find_orf_translate.

    Returns:
        bool: True if checksums match, False otherwise.
    """
    # Calculate checksums
    checksum_orfs = calculate_checksum(all_orfs)
    checksum_proteins = calculate_checksum(all_proteins)

    # Print checksums
    print(f"Checksum for find_orfs: {checksum_orfs}")
    print(f"Checksum for find_orf_translate: {checksum_proteins}")

    # Compare checksums
    return checksum_orfs == checksum_proteins


if __name__ == "__main__":
    fasta_file = 'mrgd_qf_fasta-2/test.fasta'
    
    sequences_biopython = parse_fasta(fasta_file)

    # Get results from both functions
    all_orfs = find_orfs(sequences_biopython, translate=True, trim_trailing=False)
    all_prots = find_orf_translate(sequences_biopython)
    longest_orf = find_longest_orf(all_orfs)
    longest_prot = find_longest_orf(all_prots)
    print(longest_orf)
    print(longest_prot)

    # Compare the results
    if compare_results(longest_orf, longest_prot):
        print("The results are identical.")
    else:
        print("The results differ.")

