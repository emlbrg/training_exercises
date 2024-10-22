String manipulation and parsing are critical in bioinformatics, especially when dealing with sequence data, alignments, or formats like FASTA, SAM, etc etc.
As an example, I want to parse and analyze DNA sequences from a **FASTA** file using different methodologies.

### Problem Statement
You are given a FASTA file containing DNA sequences, and you need to:
1. Parse the sequences.
2. Calculate basic statistics such as sequence length and GC content.
3. Identify any sequences containing specific motifs (e.g., a promoter region) to find ORF.
4. Optimization for larger datasets

Example FASTA format:
```
>Sequence1
ATGCGTAGCTAGTGCATGCTAGTCGTAGCTAGCTAGCGTACG
>Sequence2
GCTAGCTAGCTAGGCGATCGTAGCATCGATGCATGCATGCTAGC
>Sequence3
ATGCATGCATCGATCGTAGCTAGCGTAGCTAGCTAGC
```

### Task 1: Parsing the FASTA File
#### Approach 1: Using Basic String Manipulation
You can read the file line by line and process the sequences. This is the most intuitive method IMHO.

```python
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

# Example usage
fasta_file = 'example.fasta'
sequences = parse_fasta(fasta_file)
for seq_id, sequence in sequences.items():
    print(f'{seq_id}: {sequence}')
```

**Pros**:
- Simple and direct approach.
- Efficient for small to medium-sized FASTA files.
- No dependencies or external libraries needed.

**Cons**:
- Not suitable for extremely large datasets where memory management becomes an issue.
- You need to ensure sequences are properly formatted (e.g., no multiline headers).
- Loads the entire thing into memory.

#### Approach 2: Using Biopython
```python
from Bio import SeqIO

def parse_fasta_biopython(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences[record.id] = str(record.seq)
    return sequences

fasta_file = 'example.fasta'
sequences = parse_fasta(fasta_file)
for seq_id, sequence in sequences.items():
    print(f'{seq_id}: {sequence}')
```

**Pros:**
- Efficient and scalable: Handles large FASTA files better than the basic approach.
- Provides additional functionality (e.g., handling multiple sequence formats, validation).
- Handles errors and inconsistencies more gracefully.

**Cons:**
- Requires installing and importing the biopython library.
- May be overkill for simple tasks.

### Task 2: Calculating Sequence Statistics (Length and GC Content)
Once the sequences are parsed, you can easily calculate statistics such as the length of each sequence and its GC content.

#### Approach 1: Using string methods

```python
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

# Example usage
sequence_statistics(sequences)
```

**Pros**:
- Straightforward to calculate various statistics like length, GC content, or other sequence properties.
- Easily extendable to compute other properties (e.g., AT content, dinucleotide frequencies).

**Cons**:
- No parallelization; could be slow for very large datasets.
- Assumes sequences are already parsed into memory.
- Slower for larger datasets due to manual looping
- Does not handle ambiguous bases (like N or X)

#### Approach 2: Using Numpy
Using `numpy` arrays for handling sequences can significantly speed up operations on large datasets.

```python
import numpy as np

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

sequence_statistics_numpy(sequences)
```

**Pros:**
- Fast and efficient: Handles larger datasets quickly by leveraging vectorized operations in numpy.
- Scalable for large sequence datasets.

**Cons:**
- Requires numpy, an external library.
- Not as intuitive for someone unfamiliar with numpy arrays.

<!-- ### Task 3: Finding Specific Motifs in Sequences
Now, let’s identify sequences containing a specific motif (e.g., a promoter region like "TATA").

#### Approach 1: Using Basic String Search
```python
def find_motifs(sequences, motif):
    motif_results = {}
    for seq_id, sequence in sequences.items():
        if motif in sequence:
            motif_results[seq_id] = sequence
    return motif_results

# Example usage: Find sequences containing the "TATA" motif
motif = 'TATA'
motif_sequences = find_motifs(sequences, motif)
for seq_id, sequence in motif_sequences.items():
    print(f"Motif '{motif}' found in {seq_id}: {sequence}")
```

**Pros**:
- Quick and easy to implement.
- Works well for simple motifs or patterns.

**Cons**:
- Not suitable for finding more complex patterns (e.g., degenerate motifs, which allow for mismatches).
- Requires exact matches.

#### Approach 2: Using Regular Expressions for Complex Patterns
Bioinformatics often deals with motifs that are not exact, such as consensus sequences with ambiguous nucleotides. Regular expressions can be used to find more flexible patterns.

```python
import re

def find_motifs_regex(sequences, motif_regex):
    motif_results = {}
    pattern = re.compile(motif_regex)
    for seq_id, sequence in sequences.items():
        if pattern.search(sequence):
            motif_results[seq_id] = sequence
    return motif_results

# Example usage: Find sequences containing a flexible "TATA" box (e.g., "TATA[AT]A")
motif_regex = r'TATA[AT]A'
motif_sequences = find_motifs_regex(sequences, motif_regex)
for seq_id, sequence in motif_sequences.items():
    print(f"Motif '{motif_regex}' found in {seq_id}: {sequence}")
```

**Pros**:
- Flexible for finding motifs with ambiguous nucleotides or allowing for mismatches.
- Useful for searching for consensus sequences.

**Cons**:
- Regular expressions can become complex.
- Slower for very large datasets due to regex matching overhead. -->

### **Task 3: Identifying ORFs (Open Reading Frames)**

Identifying open reading frames (ORFs) in DNA sequences involves finding stretches of nucleotides that start with a start codon (e.g., `ATG`) and end with a stop codon (e.g., `TAA`, `TAG`, or `TGA`). Generally, the correct ORF is the longest one that starts with `ATG` and ends with a stop codon.
If you have expected protein lengths (based on known antibodies), you can use that to filter your candidates.

**IMPORTANT: Validate that the extracted ORFs are composed of triplet codons and that they do not contain frameshifts or premature stop codons.**

#### **Approach 1: Basic String Search**

This method uses a basic string search to scan through the sequence in all three reading frames and identify ORFs.

```python
def find_orfs_basic(sequence, min_len=3):  # 6 frames
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    
    # Forward strand (3 frames)
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                # Find the nearest stop codon
                for j in range(i + 3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = sequence[i:j+3]
                        if len(orf) >= min_len:
                            orfs.append((1, frame, i, j + 3, orf))  # Strand +1 for forward
                        break
    
    # Reverse complement strand (3 frames)
    reverse_seq = str(Seq(sequence).reverse_complement())  # Get reverse complement
    for frame in range(3):
        for i in range(frame, len(reverse_seq), 3):
            codon = reverse_seq[i:i+3]
            if codon == start_codon:
                # Find the nearest stop codon
                for j in range(i + 3, len(reverse_seq), 3):
                    stop_codon = reverse_seq[j:j+3]
                    if stop_codon in stop_codons:
                        orf = reverse_seq[i:j+3]
                        if len(orf) >= min_len:
                            orfs.append((-1, frame, i, j + 3, orf))  # Strand -1 for reverse
                        break
    
    return orfs
```

**Pros**:
- Simple and easy to understand.
- No external libraries needed.

**Cons**:
- Can be slow for large sequences because of nested loops.
- Reverse complement strand requires additional nested loops.
- Not optimized for large datasets.

#### **Approach 2: Using Biopython’s CodonTable and Seq Object**

Biopython provides useful tools for working with genetic sequences, including codon tables and `Seq` objects, which can help identify ORFs more efficiently and handle both forward and reverse strands.

```python
from Bio.Seq import Seq
from Bio.Data import CodonTable
def find_orf_translate(sequence, table=1, min_pro_len=100):
    all_proteins = []

    for seq_id, seq in sequence.items():
        nuc = Seq(seq)  # Essentially, it's a string

        for strand, nuc in [(+1, nuc), (-1, nuc.reverse_complement())]:  # Forward + reverse
            for frame in range(3):
                protein_sequence = nuc[frame:].translate(table)
                # print(protein_sequence)
                for pro in protein_sequence.split("*"):
                    # print(pro)
                    if len(pro) >= min_pro_len:
                        start = frame
                        end = start + len(pro) * 3
                        all_proteins.append((seq_id, strand, frame, start, end, pro))
                        # print(all_proteins)

    return all_proteins
```

**Pros**:
- **Efficient**: Uses Biopython's internal translation and reverse complement methods.
- **Comprehensive**: Handles both forward and reverse complement strands.
- **Codon Table Integration**: Uses Biopython’s codon tables, allowing for easy changes to genetic code.

**Cons**:
- Requires Biopython library.
- Slightly more complex to understand for those unfamiliar with Biopython’s `Seq` object and codon tables.

### Task 4: Optimizing for Large Datasets (Streaming Parsing)
For very large datasets, we might want to avoid loading the entire FASTA file into memory. Instead, we can use a streaming approach.

#### Approach 1: Streaming with Python Generators
Using Python’s generator functions (via `yield`) allows you to parse sequences one at a time, avoiding the need to load everything into memory.

```python
def parse_fasta_stream(fasta_file):
    current_sequence = None
    sequence = []
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Sequence identifier
                if current_sequence:
                    yield current_sequence, ''.join(sequence)
                current_sequence = line[1:]  # Remove '>'
                sequence = []
            else:
                sequence.append(line)  # Append the sequence
        if current_sequence:
            yield current_sequence, ''.join(sequence)

# Example usage: Streaming parsing and calculating sequence lengths
for seq_id, sequence in parse_fasta_stream(fasta_file):
    print(f'Sequence ID: {seq_id}, Length: {len(sequence)}')
```

**Pros**:
- Memory-efficient, suitable for large datasets and doesn't load it all into memory.
- Only processes one sequence at a time.

**Cons**:
- More complex to implement and use compared to in-memory parsing.
- Harder to perform operations requiring multiple sequences simultaneously. For example, cannot calculate sequence-wide statistics in one pass unless all sequences are parsed (or handled in chunks).

#### Approach 2: Using Biopython with Iterators
Biopython provides a built-in method for streaming sequences using iterators, which achieves the same memory efficiency as the generator-based approach.

```python
from Bio import SeqIO

def parse_fasta_biopython_stream(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        yield record.id, str(record.seq)
```

**Pros:**
- Memory-efficient and simple to implement with Biopython.
- Error handling and other Biopython features are still available.

**Cons:**
- Requires Biopython to be installed.
- Slightly more overhead than the generator method for simple use cases.


### **Pros and Cons Overview (Including Task 3)**

| **Task**                | **Approach 1 (Basic)**                         | **Approach 2 (Biopython/Numpy)**                       |
|-------------------------|------------------------------------------------|--------------------------------------------------------|
| **Parsing (Task 1)**     | - Simple and easy to implement<br> - No dependencies | - Faster for large datasets<br> - More robust error handling |
| **GC Content (Task 2)**  | - Easy to understand<br> - No dependencies | - Optimized for performance using vectorized operations |
| **ORF Finding (Task 3)** | - Simple string search<br> - No dependencies | - Handles both strands<br> - Efficient codon handling |
| **Streaming Parsing (Task 4)** | - Memory-efficient<br> - Customizable generator | - Memory-efficient<br> - Biopython features and error handling |

### Conclusion
This example illustrates how string manipulation and parsing are common tasks in bioinformatics, with typical operations including parsing sequence data, calculating statistics, and identifying motifs. Different approaches (basic string methods, regex, streaming) offer varying levels of performance and flexibility depending on the size and complexity of the data.
The choice depends on the size of the data, performance requirements, and whether you want additional features like format handling, validation, and error checking. While the basic approaches are intuitive and easy to implement, Biopython and Numpy provide more scalability and efficiency, especially for larger datasets.
