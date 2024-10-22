def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))


#### OPTION A ####
### Brure force ###
def find_palindromes_bruteforce(dna_sequence: str, length: int) -> list:
    palindromes = []
    for i in range(len(dna_sequence) - length + 1):
        substring = dna_sequence[i:i+length]
        if substring == reverse_complement(substring):
            palindromes.append(substring)
    return palindromes


#### OPTION B ####
### Sliding window ###
def find_sliding_window_palindromes(dna_sequence: str, length: int) -> list:
    palindromes = []
    reverse_complements = {i: reverse_complement(dna_sequence[i:i+length])
                           for i in range(len(dna_sequence) - length + 1)}
    
    for i in range(len(dna_sequence) - length + 1):
        substring = dna_sequence[i:i+length]
        if substring == reverse_complements[i]:
            palindromes.append(substring)
    return palindromes

#### OPTION C ####
### Dicitonary ###
def find_palindromes_dict(dna_sequence: str, length: int) -> dict:
    palindromes = {}
    
    for i in range(len(dna_sequence) - length + 1):
        substring = dna_sequence[i:i+length]
        if substring == reverse_complement(substring):
            palindromes[i] = substring
    
    return palindromes
# def find_palindromes_dict(dna_sequence: str, length: int) -> dict:
#     reverse_complements = {}
#     palindromes = {}

#     for i in range(len(dna_sequence) - length + 1):
#         substring = dna_sequence[i:i+length]
#         # rev_comp = reverse_complement(substring)

#         # only count palindromes, not just reverse complements
#         if substring == reverse_complement(substring) and substring not in reverse_complements:
#             palindromes[i] = substring
        
#         reverse_complements[substring] = i  # Track the reverse complement seen
#     return palindromes


