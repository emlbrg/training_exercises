Compare the results of your Python k-mer counting implementations against **Jellyfish** as a "gold standard" tool.

### Step 1: Run Jellyfish for k-mer Counting
Once installed, you can use Jellyfish to count k-mers on your FASTA file. Here’s the basic command to run Jellyfish:

```bash
jellyfish count -m 3 -s 100M -o jellyfish_output.jf test.fasta
```

- `-m 3`: Specifies the k-mer size (in this case, k=3).
- `-s 100M`: Sets the hash size, which can be adjusted based on the file size.
- `-o jellyfish_output.jf`: The output file for the k-mer counts.
- `test.fasta`: The input FASTA file.

After counting, you’ll need to dump the results into a readable format (like a text file):

```bash
jellyfish dump jellyfish_output.jf > jellyfish_kmers.fasta
```

This will output a file `jellyfish_kmers.fasta` with k-mer counts in the following format:
```
>45
AAA
>32
AAG
>27
AAC
...
```

### Step 3: Parse Jellyfish Output in Python
Load Jellyfish’s results into your Python script to compare them with the outputs from your Python k-mer counting implementations.

Here’s a helper function to parse the Jellyfish output:

```python
def parse_jellyfish_output(file_path: str) -> dict:
    """
    Parse the Jellyfish output file into a dictionary.

    Args:
        file_path (str): Path to the Jellyfish k-mer output file.

    Returns:
        dict: Dictionary where keys are k-mers and values are counts.
    """
    kmer_counts = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                kmer, count = parts
                kmer_counts[kmer] = int(count)
    return kmer_counts
```

### Step 4: Compare Python Outputs with Jellyfish
Now, you can compare the k-mer counts from your Python functions (`naive_kmer_count`, `hashmap_kmer_count`, `optimized_counter_kmer_count`) with the results from Jellyfish.

Here’s how you can integrate it into your existing script:

```python
# Parse the Jellyfish k-mer counts
jellyfish_output = parse_jellyfish_output('jellyfish_kmers.txt')

# Calculate the MD5 checksum of the Jellyfish output
jellyfish_checksum = calculate_md5_checksum(jellyfish_output)

# Output the Jellyfish checksum
print(f"Jellyfish MD5 Checksum: {jellyfish_checksum}")

# Compare the Jellyfish results with Python k-mer counting implementations
for seq_id, naive_checksum, hashmap_checksum, optimized_checksum in md5_checksums:
    print(f"\nSequence ID: {seq_id}")
    print(f"Naive MD5 Checksum: {naive_checksum}")
    print(f"Hashmap MD5 Checksum: {hashmap_checksum}")
    print(f"Optimized MD5 Checksum: {optimized_checksum}")

    # Compare against Jellyfish
    if naive_checksum == hashmap_checksum == optimized_checksum == jellyfish_checksum:
        print(f"All implementations and Jellyfish produce identical results for sequence {seq_id}!")
    else:
        print(f"Mismatch in results for sequence {seq_id} compared to Jellyfish.")
```

### Step 5: Validate Results
1. **Run your Python script** with the k-mer counting functions and compare the MD5 checksums against the Jellyfish output.
2. **Analyze any differences**: If there are mismatches, you can further investigate the Python implementations or fine-tune Jellyfish’s settings.

### Summary:
- **Install Jellyfish**.
- **Run Jellyfish** to generate k-mer counts.
- **Parse Jellyfish’s output** into a dictionary.
- **Compare** the MD5 checksums of your Python results with Jellyfish’s results.
