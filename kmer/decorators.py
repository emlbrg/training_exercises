"""
author: @emlbrg
"""

import time
from functools import wraps
from typing import Callable, Any

def time_it(func: Callable) -> Callable:
    """
    Decorator to measure the execution time of a function.

    Args:
        func (function): The function to be timed.

    Returns:
        function: Wrapped function with timing logic.
    """
    @wraps(func)
    def timer(*args: Any, **kwargs: Any) -> Any:
        print("--Starting timer--")
        start_time = time.time()
        results = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        print(f'Total execution time: {int(hours)} hours, {int(minutes)} minutes, {int(seconds)} seconds')
        return results
    return timer

# usage:
# @time_it
# def your_function():