from parse import parse_fasta
from typing import List, Dict
import re

def find_str_sliding_window(seq_dict: Dict, min_repeat_count: int=14) -> List:
    """Finds simple tandem repeats using a sliding window approach with fixed repeat unit length.

    Args:
        seq_dict (dict): Dictionary containing sequence IDs as keys and DNA sequences as values.
        min_repeat_count (int): Minimum number of repeats to consider an STR.

    Returns:
        list: List of STRs (sequence ID, start position, end position, repeat unit, repeat count).
    """
    str_list = []

    for seq_id, sequence in seq_dict.items():
        seq_length = len(sequence)
        k = 2
        found_positions: list = []

        for i in range(seq_length - k + 1):
            repeat_unit = sequence[i:i + k]
            count = 0
            j = i

            while j < seq_length and sequence[j:j + k] == repeat_unit:  # Count consecutive?
                count += 1
                j += k

            if count >= min_repeat_count:
                if not any(start < j and end > i for start, end in found_positions):  # Longest only
                    str_list.append((seq_id, i + 1, j, repeat_unit, count))
                    found_positions.append((i + 1, j))

    return str_list

def find_str_regex(seq_dict: Dict, min_repeat_count: int=14) -> List:
    """Finds simple tandem repeats using regular expressions with fixed repeat unit length.

    Args:
        seq_dict (dict): Dictionary containing sequence IDs as keys and DNA sequences as values.
        min_repeat_count (int): Minimum number of repeats to consider an STR.

    Returns:
        list: List of STRs (sequence ID, start position, end position, repeat unit, repeat count).
    """
    str_list = []

    for seq_id, sequence in seq_dict.items():
        k = 2  # Repeat length (AT for example)
        pattern = f"(.{{{k}}})\\1{{{min_repeat_count - 1},}}"  # e.g., (AT){14,} for 14 or more repeats
        matches = re.finditer(pattern, sequence)

        for match in matches:
            repeat_unit = match.group(1)
            count = match.group(0).count(repeat_unit)
            start = match.start() + 1
            end = match.end()  # End position is inclusive!
            str_list.append((seq_id, start, end, repeat_unit, count))

    return str_list

def find_str_hashmap(seq_dict: Dict, min_repeat_count: int=14, max_repeat_length: int=30) -> List:
    """Finds simple tandem repeats using a hashmap for counting. It stores positions of k-mers
    and counts consecutive repeats.

    Args:
        seq_dict (dict): Dictionary containing sequence IDs as keys and DNA sequences as values.
        min_repeat_count (int): Minimum number of repeats to consider an STR.
        max_repeat_length (int): Maximum length of the repeat unit.

    Returns:
        list: List of STRs (sequence ID, start position, end position, repeat unit, repeat count).
    """
    str_list = []

    for seq_id, sequence in seq_dict.items():
        seq_length = len(sequence)
        found_positions = []  # start-end of previously found repeats

        for k in range(2, max_repeat_length + 1):
            kmer_map = {}
            for i in range(seq_length - k + 1):
                kmer = sequence[i:i + k]
                if kmer in kmer_map:
                    kmer_map[kmer].append(i)
                else:
                    kmer_map[kmer] = [i]

            for kmer, positions in kmer_map.items():  # check for STRs in the hashmap
                count = 1
                prev_pos = positions[0]

                for pos in positions[1:]:
                    if pos == prev_pos + k:  # consecutive k-mer?
                        count += 1
                    else:
                        if count >= min_repeat_count:
                            start_pos = prev_pos - (count - 1) * k + 1
                            end_pos = prev_pos + k
                            if not any(start < end_pos and end > start_pos for start, end in found_positions):
                                str_list.append((seq_id, start_pos, end_pos, kmer, count))
                                found_positions.append((start_pos, end_pos))

                        count = 1

                    prev_pos = pos

                if count >= min_repeat_count:
                    start_pos = prev_pos - (count - 1) * k + 1
                    end_pos = prev_pos + k
                    if not any(start < end_pos and end > start_pos for start, end in found_positions):
                        str_list.append((seq_id, start_pos, end_pos, kmer, count))
                        found_positions.append((start_pos, end_pos))

    return str_list


seq_dict = parse_fasta("test.fasta")
sliding_window = find_str_sliding_window(seq_dict)
regex = find_str_regex(seq_dict)
hashmap = find_str_hashmap(seq_dict)
print(f"Sliding window results: {sliding_window}")
print(f"Regex method results: {regex}")
print(f"Hashmap method results: {hashmap}")
