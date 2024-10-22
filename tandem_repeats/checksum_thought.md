The differences in performance between the **sliding window**, **regular expression**, and **hashmap-based** methods are due to how each approach handles sequence processing and memory usage. Here's a detailed explanation of why you're seeing these results:

### 1. **Sliding Window:**
- **Why is it surprisingly fast?**
    - **Direct Iteration**: The sliding window method directly scans the DNA sequence character by character. It uses a fixed-length window (`k`), sliding across the sequence, and checks for repeating units. This avoids complex backtracking or pattern matching like regex.
    - **Simple Looping Logic**: Its logic is straightforward and doesn't involve hashing, pattern matching, or additional structures, keeping the runtime per character fairly low.
    - **Efficiency in Code Structure**: Even though the sliding window method loops through each position and counts consecutive repeats, the simplicity of iteration ensures it runs efficiently.
  
- **Why is it more memory intensive?**
    - **Tracking Repeat Positions**: The sliding window approach needs to keep track of positions where STRs are found to avoid overlaps. This list grows as more STRs are found, consuming more memory.
    - **Full Sequence Processing**: Since it processes the entire sequence without shortcuts or optimizations like hashing, the memory usage accumulates more quickly, especially for large sequences or multiple repeat lengths.

### 2. **Regular Expression:**
- **Why is it less memory intensive?**
    - **Compact Representation**: Regex engines are highly optimized for matching patterns. When using regular expressions, the engine constructs an efficient internal representation of the pattern, often using less memory than explicitly iterating through every possible position like in the sliding window method.
    - **Built-in Optimization**: Regular expressions have built-in optimizations that automatically discard intermediate states that aren’t relevant, reducing memory usage compared to manual iteration approaches.
  
- **Why can it be challenging to write?**
    - **Complex Patterns**: Crafting the correct regex pattern for more complex tandem repeats can be tricky. For example, adjusting patterns for variable-length repeat units or handling edge cases (like partial overlaps) can lead to very complicated expressions.
    - **Pattern Matching vs Simple Counting**: Regex engines excel at pattern matching but require fine-tuning to ensure correct counting of repeats, especially when working with DNA sequences.

### 3. **Hashmap-Based:**
- **Why does it consume similar memory to sliding window?**
    - **K-mer Hashmap**: The hashmap-based approach constructs a dictionary to store the positions of each k-mer (repeat unit). While this helps speed up lookups, it increases memory usage because each k-mer is stored as a key with a list of positions. This can use more memory for sequences with many k-mers or repeats.
    - **Storing Multiple Repeats**: The method stores every possible repeat length (from 2 up to `max_repeat_length`), so for each sequence, the hashmap keeps track of all positions for every repeat unit. This results in increased memory usage for larger repeat lengths and longer sequences.

- **Why is it much slower?**
    - **Repeated Hash Lookups**: The hashmap-based approach has an additional overhead from hashing every k-mer and looking it up in the hashmap. While hashing can be efficient, the overhead accumulates, especially for large sequences with many potential k-mers.
    - **Sequential Position Check**: After finding a k-mer in the hashmap, the method still needs to verify if the k-mers are consecutive by iterating through the positions, which can become computationally expensive.
    - **Multiple Repeat Lengths**: Since it calculates and checks k-mers of varying lengths (`k` from 2 to `max_repeat_length`), the method has to iterate multiple times over the sequence, further slowing down execution.

### Summary of Performance:
- **Sliding Window**: It’s faster because it directly scans the sequence without additional structure (like regex engines or hashmaps), but higher memory usage comes from keeping track of repeat positions.
- **Regular Expression**: More optimized in memory due to efficient internal handling by the regex engine. However, developing the right regex for different cases can be tricky and may not always be the most intuitive solution.
- **Hashmap-Based**: While using hashmaps can reduce repeated work in some cases, it introduces overhead with hash lookups and managing position lists, making it slower than direct iteration or regex. Memory usage remains high due to the storage of k-mers and their positions.

### Choosing the Best Approach:
- If your primary concern is speed and your sequences are large, **sliding window** might be the best approach.
- For memory efficiency, especially with simpler patterns, **regular expressions** could be a good choice.
- If you need to handle more complex repeat unit lengths and have sufficient memory, **hashmap-based** methods can offer flexibility, albeit with slower runtime.