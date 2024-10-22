from parse import parse_fasta
import hashlib
import time
from memory_profiler import memory_usage
from typing import List, Dict
import re
from tandem_repeats import find_str_sliding_window, find_str_regex, find_str_hashmap


def calculate_md5(data) -> str:
    """Calculate MD5 checksum for a list of STRs."""
    md5 = hashlib.md5()
    for entry in data:
        md5.update(str(entry).encode('utf-8'))
    return md5.hexdigest()

def profile_method(method, seq_dict, min_repeat_count=14, max_repeat_length=30):
    """Profiles time and memory usage for the given method."""
    start_time = time.time()

    if method == find_str_hashmap:
        mem_usage = memory_usage((method, (seq_dict,), {'min_repeat_count': min_repeat_count, 'max_repeat_length': max_repeat_length}))
        result = method(seq_dict, min_repeat_count=min_repeat_count, max_repeat_length=max_repeat_length)
    else:
        mem_usage = memory_usage((method, (seq_dict,), {'min_repeat_count': min_repeat_count}))
        result = method(seq_dict, min_repeat_count=min_repeat_count)
    
    end_time = time.time()
    checksum = calculate_md5(result)
    time_taken = end_time - start_time
    peak_memory = max(mem_usage)
    
    return {
        'checksum': checksum,
        'time': time_taken,
        'peak_memory': peak_memory,
        'result': result
    }

def compare_results(sliding_window_results, regex_results, hashmap_results):
    """Compares the results from three different methods and prints if they are the same or differ."""
    
    sliding_window_set = set(sliding_window_results)
    regex_set = set(regex_results)
    hashmap_set = set(hashmap_results)

    if sliding_window_set == regex_set:
        print("Sliding window and regex results are the same.")
    else:
        print("Sliding window and regex results differ.")
    
    if sliding_window_set == hashmap_set:
        print("Sliding window and hashmap results are the same.")
    else:
        print("Sliding window and hashmap results differ.")
    
    if regex_set == hashmap_set:
        print("Regex and hashmap results are the same.")
    else:
        print("Regex and hashmap results differ.")

if __name__ == "__main__":
    seq_dict = parse_fasta("test.fasta")
    
    print("Profiling sliding window...")
    sliding_window_profile = profile_method(find_str_sliding_window, seq_dict)
    print(f"Sliding Window - MD5: {sliding_window_profile['checksum']}, Time: {sliding_window_profile['time']}s, Peak Memory: {sliding_window_profile['peak_memory']} MB")

    print("Profiling regex...")
    regex_profile = profile_method(find_str_regex, seq_dict)
    print(f"Regex - MD5: {regex_profile['checksum']}, Time: {regex_profile['time']}s, Peak Memory: {regex_profile['peak_memory']} MB")
    
    print("Profiling hashmap...")
    hashmap_profile = profile_method(find_str_hashmap, seq_dict)
    print(f"Hashmap - MD5: {hashmap_profile['checksum']}, Time: {hashmap_profile['time']}s, Peak Memory: {hashmap_profile['peak_memory']} MB")
    
    compare_results(sliding_window_profile['result'], regex_profile['result'], hashmap_profile['result'])