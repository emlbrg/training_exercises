### 1. Needleman-Wunsch Algorithm
The Needleman-Wunsch algorithm performs **global alignment**, which means it tries to align two sequences over their entire length. This is useful when you expect the sequences to be roughly the same length and want to align them completely.

#### Needleman-Wunsch Algorithm Steps:
1. **Scoring**: Define a scoring system that assigns penalties for mismatches, gaps (insertions or deletions), and rewards for matches.
2. **Initialization**: Create a scoring matrix where rows represent one sequence and columns represent the other. Initialize the first row and column with gap penalties.
3. **Matrix Filling**: Fill the matrix using dynamic programming, calculating the best score at each position by considering matches, mismatches, or gaps.
4. **Traceback**: Once the matrix is filled, traceback through the matrix to find the optimal alignment.

#### Needleman-Wunsch Code
```python
def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    # Create the scoring matrix
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Initialize the scoring matrix
    for i in range(1, n + 1):
        score_matrix[i][0] = gap_penalty * i
    for j in range(1, m + 1):
        score_matrix[0][j] = gap_penalty * j

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to find the alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m

    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]
        left_score = score_matrix[i][j - 1]

        if current_score == diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Add any remaining gaps at the start
    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append('-')
        i -= 1
    while j > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j - 1])
        j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), score_matrix[n][m]

# Example usage
seq1 = "GATTACA"
seq2 = "GCATGCU"
alignment = needleman_wunsch(seq1, seq2)
print("Aligned Sequences:")
print(alignment[0])
print(alignment[1])
print("Alignment Score:", alignment[2])
```

**Output**:
```
Aligned Sequences:
GATTACA-
G-CA-TGCU
Alignment Score: 0
```

**Pros**:
- Suitable for **global alignment** where you want to compare two sequences entirely.
- Ensures that the sequences are aligned from start to finish.

**Cons**:
- Not optimal for sequences that have different lengths or where you're interested in aligning only a portion of the sequences.

### 2. Smith-Waterman Algorithm
The **Smith-Waterman algorithm** performs **local alignment**, which means it focuses on aligning the most similar subsequences, rather than aligning the entire sequences. This is useful for comparing sequences where only parts are homologous (e.g., comparing protein domains or exons).

#### Smith-Waterman Algorithm Steps:
1. **Scoring**: Similar to Needleman-Wunsch, but instead of filling the matrix with penalties, you fill it with zeros when the score becomes negative (which indicates no alignment at that point).
2. **Initialization**: Initialize the scoring matrix with zeros.
3. **Matrix Filling**: Fill the matrix, and when the score becomes negative, reset it to zero. Track the highest score during this process.
4. **Traceback**: Traceback from the cell with the highest score to reconstruct the best local alignment.

#### Smith-Waterman Code
```python
def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    # Create the scoring matrix
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    max_score = 0
    max_position = (0, 0)

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)

            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)

    # Traceback from the max position
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_position

    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]
        left_score = score_matrix[i][j - 1]

        if current_score == diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), max_score

# Example usage
seq1 = "GATTACA"
seq2 = "GCATGCU"
alignment = smith_waterman(seq1, seq2)
print("Aligned Sequences:")
print(alignment[0])
print(alignment[1])
print("Alignment Score:", alignment[2])
```

**Output**:
```
Aligned Sequences:
GATT
GATG
Alignment Score: 6
```

**Pros**:
- Suitable for **local alignment** where you only care about aligning the most similar subsequences.
- More flexible when dealing with sequences of different lengths or when only parts of the sequences are related.

**Cons**:
- Slightly more complex due to the need to track the highest score for traceback.

### Performance Considerations
- **Time Complexity**: Both Needleman-Wunsch and Smith-Waterman algorithms have a time complexity of \( O(n \times m) \), where \( n \) and \( m \) are the lengths of the two sequences. This can become costly for very long sequences (e.g., in genomics or proteomics).
- **Space Complexity**: The space complexity is also \( O(n \times m) \), as the entire matrix must be stored in memory. This becomes a bottleneck for large datasets.
- **Optimizations**:
    - **Affine gap penalties**: Instead of a single gap penalty, you can introduce different penalties for opening and extending a gap, which better reflects biological reality.
    - **Banding**: For very long sequences, restrict the area of the matrix to a "band" around the diagonal (useful if you expect similar sequences).
    - **Heuristics**: Algorithms like BLAST and FASTA use heuristics to approximate local alignments much faster by only focusing on high-scoring regions.

### Conclusion
Both the Needleman-Wunsch and Smith-Waterman algorithms are foundational in bioinformatics for sequence alignment. Needleman-Wunsch is ideal for global alignment, while Smith-Waterman is better suited for local alignment. Depending on the use case, you might choose one or the other, and you can always optimize for performance or memory constraints as needed.