from Bio.Seq import Seq
from decorators import time_it

def parse_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences



#### VERSION A ####
### This is the second fastest ###
@time_it
def find_orf_translate(sequences, table=1, min_pro_len=30):
    """
    Find long proteins from a dictionary of sequences and return detailed information.

    Args:
        sequences (dict): Dictionary of sequences where keys are IDs and values are sequences.
        table (int): Translation table to use (default is 1).
        min_pro_len (int): Minimum protein length to consider (default is 100).
    
    Returns:
        list: List of proteins with details on record ID, strand, frame, length, and sequence.
    """
    all_proteins = []

    for seq_id, seq in sequences.items():
        nuc = Seq(seq)  # Essentially, it's a string?

        for strand, nuc in [(+1, nuc), (-1, nuc.reverse_complement())]:  # Forward + reverse
            for frame in range(3):
                protein_sequence = nuc[frame:].translate(table)
                for pro in protein_sequence.split("*"):
                    start_index = pro.find("M")
                    if start_index != -1:  # There is at least one 'M'
                        trimmed_pro = pro[start_index:]
                        start = frame + start_index * 3
                        end = start + len(trimmed_pro) * 3
                        end += 3 # Stop codon 
                        # Check if the trimmed protein meets the min len
                        if len(trimmed_pro) >= min_pro_len:
                            all_proteins.append((seq_id, strand, frame, start, end, trimmed_pro))

    return all_proteins


#### VERSION B ####
### This is the slowest ###
@time_it
def find_orfs(seq_dic, translate=False, trim_trailing=True):
    """
    Find all ORFs in the given sequences and return DNA or translated proteins.

    Args:
        seq_dic (dict): Dictionary containing sequence IDs as keys and sequences as values.
        translate (bool): If True, return translated proteins, else return DNA ORFs.
        trim_trailing (bool): If True, exclude ORFs that don't end with a stop codon. 
                              If False, include ORFs that don't end with a stop codon.

    Returns:
        list: A list of tuples with ORF information (seq_id, strand, frame, start, end, sequence).
    """
    all_orfs = []

    for seq_id, sequence in seq_dic.items():
        seq_obj = Seq(sequence)
        for strand, seq in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
            for frame in range(3):
                for i in range(frame, len(seq), 3):
                    codon = seq[i:i + 3]
                    if codon == 'ATG':
                        orf_found = False
                        for j in range(i, len(seq), 3):
                            stop_codon = seq[j:j + 3]
                            if stop_codon in ['TAA', 'TAG', 'TGA']:
                                orf_seq = seq[i:j + 3]
                                if translate:
                                    translated_protein = orf_seq.translate(table=1, to_stop=True)
                                    all_orfs.append((seq_id, strand, frame, i, j + 3, str(translated_protein)))
                                else:
                                    all_orfs.append((seq_id, strand, frame, i, j + 3, str(orf_seq)))
                                orf_found = True
                                break
                        # Handle ORFs without stop codons if trim_trailing is False
                        if not orf_found and not trim_trailing:
                            orf_seq = seq[i:]
                            if translate:
                                translated_protein = orf_seq.translate(table=1, to_stop=False)  # Translate up to the end
                                all_orfs.append((seq_id, strand, frame, i, j + 3, str(translated_protein)))  # Coordinates are wrong?
                            else:
                                all_orfs.append((seq_id, strand, frame, i, j + 3, str(orf_seq)))

    return all_orfs


def find_longest_orf(all_orfs):
    """
    Find the longest ORF for each sequence ID from a list of ORFs.

    Args:
        all_orfs (list): List of ORF details (ID, strand, frame, start, end, sequence).

    Returns:
        list: A list of the longest ORFs for each sequence ID.
    """
    longest_orfs = {}
    
    # Group ORFs by sequence ID and find the longest for each ID
    for orf in all_orfs:
        seq_id = orf[0]
        if seq_id not in longest_orfs or len(orf[-1]) > len(longest_orfs[seq_id][-1]):
            longest_orfs[seq_id] = orf  # Update with the longer ORF

    return list(longest_orfs.values())


@time_it
def save_orf(all_orfs, out_file_path):
    """
    Print all proteins from a list of long proteins.

    Args:
        all_orfs (list): List of protein details (ID, strand, frame, start, end, sequence).
        out_file_path (str): Path to the output file.
    """
    if not all_orfs:
        return print("No proteins found.")

    with open(out_file_path, 'w+') as f:
        for protein in all_orfs:
            seq_id, strand, frame, start, end, pro = protein
            f.write(f"{seq_id}: {pro} - length {len(pro)}, strand {strand}, frame {frame}, start {start}, end {end}\n")
    
    print(f"Protein details written to {out_file_path}")


#### EXAMPLE ####

fasta_file = 'mrgd_qf_fasta-2/test.fasta'
out_file_path = 'mrgd_qf_fasta-2/test_frame.txt'

sequences_biopython = parse_fasta(fasta_file)
all_orfs = find_orfs(sequences_biopython, translate=True, trim_trailing=False)
all_prots = find_orf_translate(sequences_biopython)
# print(all_orfs)
longest_orf = find_longest_orf(all_orfs)
longest_prot = find_longest_orf(all_prots)
print(longest_orf)
print(longest_prot)
save_orf(longest_orf, out_file_path)



