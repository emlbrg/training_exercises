### Challenge: Finding Simple Tandem Repeats (STRs)

Simple tandem repeats (STRs) are short sequences of DNA (typically 2-6 nucleotides) that repeat consecutively. These sequences are important in areas like population genetics, forensic biology, and personal genome analysis. The task here is to identify STRs in a DNA sequence.

### Problem Specifications:
- Input: A dictionary of sequences from a FASTA file (`seq_dict = parse_fasta(fasta_file)`).
- Output: A list of STRs (sequence IDs, start and end positions, repeat unit, and the number of repeats).

### Implementations:

#### 1. **Sliding Window with Direct Comparison**
This method uses a sliding window approach to scan for consecutive repeats of a motif. You define the repeat unit size (e.g., 2-6 nucleotides) and check whether adjacent windows contain the same subsequence.

**Steps**:
- For each sequence, slide a window of size `k` (where `k` is the size of the repeat unit).
- Check whether the subsequence repeats consecutively.
- Record regions where the subsequence repeats a specified number of times (e.g., 3 or more).

```python 
def find_str_sliding_window(seq_dict: Dict, min_repeat_count: int=14) -> List:
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
```
**Pros**:
- Easy to implement and straightforward to understand.
- Works for small sequences or specific motifs.

**Cons**:
- Inefficient for large sequences due to repeated comparisons.
- Hard to handle larger or complex STR patterns.

#### 2. **Regular Expression Matching**
Use regular expressions to search for STRs by defining a pattern that captures consecutive repeats of a short motif.

**Steps**:
- For each sequence, construct a regular expression pattern for possible repeat units (e.g., `(A[TG]{2,6})`).
- Use Python’s `re` module to find all matches of the repeat pattern in the sequence.
- Record the position and number of repeats for each match.

```python
def find_str_regex(seq_dict: Dict, min_repeat_count: int=14) -> List:
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
```
**Pros**:
- Powerful and flexible for capturing various repeat patterns.
- Efficient for finding STRs across large sequences.

**Cons**:
- Complex to construct regular expressions for different repeat units.
- May miss complex STR patterns not easily captured by a single regular expression.

#### 3. **Hashmap-Based Counting**
Use a hashmap (dictionary) to store k-mers and their counts in the sequence. Check for consecutive occurrences of the same k-mer and identify STRs based on this count.

**Steps**:
- For each sequence, extract all k-mers and store their positions in a hashmap.
- Check for consecutive k-mers at adjacent positions to identify STRs.
- Record regions with repeated k-mers as STRs.

```python
def find_str_hashmap(seq_dict, min_repeat_count=16, max_repeat_length=30):
    str_list = []

    for seq_id, sequence in seq_dict.items():
        seq_length = len(sequence)
        kmer_map = {}

        # Populate the hashmap with k-mers
        for k in range(2, max_repeat_length + 1):  # Repeat unit length from 2 to max_repeat_length
            for i in range(seq_length - k + 1):
                kmer = sequence[i:i + k]
                if kmer in kmer_map:
                    kmer_map[kmer].append(i)
                else:
                    kmer_map[kmer] = [i]

        # Check for STRs in the hashmap
        for kmer, positions in kmer_map.items():
            count = 0
            prev_pos = -1

            for pos in positions:
                if prev_pos == -1 or pos == prev_pos + len(kmer):
                    count += 1
                else:
                    if count >= min_repeat_count:
                        str_list.append((seq_id, prev_pos - len(kmer) + 1, prev_pos, kmer, count))
                    count = 1
                prev_pos = pos
            
            # Check for last count
            if count >= min_repeat_count:
                str_list.append((seq_id, prev_pos - len(kmer) + 1, prev_pos, kmer, count))

    return str_list
```
**Pros**:
- Efficient for counting k-mers and identifying repeats.
- Can be extended to handle more complex repeat patterns.

**Cons**:
- Memory-intensive for large sequences and k-mer sizes.
- Does not easily capture mismatches or imperfect repeats.

### Gold Standard Approach:
Compare your results with a tool like **Tandem Repeats Finder (TRF)**, which is widely used for detecting STRs in DNA sequences.

### Checksum Benchmarking:
To compare the STRs found by your implementations with the gold standard, you can compute a checksum.

```python
import hashlib

def compute_str_checksum(strs: list) -> str:
    """Compute checksum for a list of STRs.

    Args:
        strs (list): List of STRs (ID, start, end, repeat unit, repeat count).

    Returns:
        str: The checksum of the concatenated STR details.
    """
    str_data = ''.join([f"{id}:{start}-{end}:{unit}:{count}" for id, start, end, unit, count in strs])
    return hashlib.sha256(str_data.encode()).hexdigest()

# Example usage
strs = [('seq1', 0, 10, 'AT', 5), ('seq2', 50, 70, 'GCG', 7)]
checksum = compute_str_checksum(strs)
print(f"Checksum: {checksum}")
```

### Benchmarking and Performance Metrics:
1. **Accuracy**: Compare the STRs detected by your methods against the output from Tandem Repeats Finder.
2. **Time Efficiency**: Measure the time taken by each method to find STRs.
3. **Memory Usage**: Monitor the memory consumption of each approach using `memory_profiler`.

### Summary:

- **Sliding Window**: Simple but inefficient for large sequences.
- **Regular Expression**: Powerful and flexible but may require complex patterns.
- **Hashmap-Based**: Efficient for counting repeats but memory-intensive.

This problem is straightforward but offers plenty of room for comparing performance and accuracy across different methods.