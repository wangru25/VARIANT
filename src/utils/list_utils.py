# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-20 11:13:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/utils/list_utils.py
Description: List processing utilities and helper functions.
'''

from typing import Any, List, Tuple


def delete_element_from_list(input_list: List[Any], index: int) -> List[Any]:
    '''
    Delete an element from a list at a given index.

    Args:
        input_list: List to modify
        index: Index of element to delete

    Returns:
        Modified list with element removed

    Note:
        This function modifies the input list in-place and returns it.
        If the index is out of bounds, the list is returned unchanged.
    '''
    if 0 <= index < len(input_list):
        del input_list[index]
    return input_list


def delete_char_from_string(sequence: str, index: int) -> str:
    '''
    Delete a character from a string at a given index.

    Args:
        sequence: String to modify
        index: Index of character to delete

    Returns:
        String with character removed

    Raises:
        IndexError: If index is out of bounds
    '''
    if index < 0 or index >= len(sequence):
        raise IndexError(
            f"Index {index} out of bounds for string of length {len(sequence)}"
        )

    return sequence[:index] + sequence[index + 1 :]


def safe_delete_element(input_list: List[Any], index: int) -> Tuple[List[Any], bool]:
    '''
    Safely delete an element from a list without raising exceptions.

    Args:
        input_list: List to modify
        index: Index of element to delete

    Returns:
        Tuple of (modified_list, success_flag)
    '''
    if 0 <= index < len(input_list):
        result_list = input_list.copy()
        del result_list[index]
        return result_list, True
    else:
        return input_list.copy(), False


def find_element_indices(input_list: List[Any], target: Any) -> List[int]:
    '''
    Find all indices of a target element in a list.

    Args:
        input_list: List to search in
        target: Element to find

    Returns:
        List of indices where the target element appears
    '''
    return [i for i, element in enumerate(input_list) if element == target]


def remove_duplicates_preserve_order(input_list: List[Any]) -> List[Any]:
    '''
    Remove duplicate elements from a list while preserving order.

    Args:
        input_list: List to deduplicate

    Returns:
        List with duplicates removed
    '''
    seen = set()
    result = []

    for element in input_list:
        if element not in seen:
            seen.add(element)
            result.append(element)

    return result


def split_list_by_condition(
    input_list: List[Any], condition_func: callable
) -> Tuple[List[Any], List[Any]]:
    '''
    Split a list into two parts based on a condition function.

    Args:
        input_list: List to split
        condition_func: Function that returns True/False for each element

    Returns:
        Tuple of (elements_where_true, elements_where_false)
    '''
    true_elements = []
    false_elements = []

    for element in input_list:
        if condition_func(element):
            true_elements.append(element)
        else:
            false_elements.append(element)

    return true_elements, false_elements


def chunk_list(input_list: List[Any], chunk_size: int) -> List[List[Any]]:
    '''
    Split a list into chunks of specified size.

    Args:
        input_list: List to chunk
        chunk_size: Size of each chunk

    Returns:
        List of chunks

    Raises:
        ValueError: If chunk_size is less than 1
    '''
    if chunk_size < 1:
        raise ValueError("Chunk size must be at least 1")

    return [
        input_list[i : i + chunk_size] for i in range(0, len(input_list), chunk_size)
    ]


def flatten_nested_list(nested_list: List[List[Any]]) -> List[Any]:
    '''
    Flatten a nested list structure.

    Args:
        nested_list: List of lists to flatten

    Returns:
        Flattened list
    '''
    return [element for sublist in nested_list for element in sublist]


def count_occurrences(input_list: List[Any]) -> dict:
    '''
    Count occurrences of each element in a list.

    Args:
        input_list: List to count elements in

    Returns:
        Dictionary mapping elements to their counts
    '''
    counts = {}
    for element in input_list:
        counts[element] = counts.get(element, 0) + 1
    return counts


def find_most_common(input_list: List[Any], n: int = 1) -> List[Tuple[Any, int]]:
    '''
    Find the n most common elements in a list.

    Args:
        input_list: List to analyze
        n: Number of most common elements to return

    Returns:
        List of tuples (element, count) sorted by frequency
    '''
    counts = count_occurrences(input_list)
    sorted_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    return sorted_items[:n]


def rotate_list(input_list: List[Any], positions: int) -> List[Any]:
    '''
    Rotate a list by a specified number of positions.

    Args:
        input_list: List to rotate
        positions: Number of positions to rotate (positive = right, negative = left)

    Returns:
        Rotated list
    '''
    if not input_list:
        return input_list

    n = len(input_list)
    positions = positions % n

    return input_list[-positions:] + input_list[:-positions]


def contains_element(test_string: str, test_list: List[str]) -> bool:
    '''
    Check if a string contains any element from a list.

    Args:
        test_string: String to search in
        test_list: List of strings to search for

    Returns:
        True if any element from test_list is found in test_string
    '''
    return any(element in test_string for element in test_list)


def is_keys_subset(sub_list: List[Any], test_list: List[Any]) -> Tuple[bool, List[Any]]:
    '''
    Check if the first list is a proper subset of the second list.

    Args:
        sub_list: List to check if it's a subset
        test_list: List to check against

    Returns:
        Tuple of (is_subset, difference_list)

    Example:
        >>> is_keys_subset([1, 2], [1, 2, 3, 4])
        (True, [3, 4])
    '''
    is_subset_result = set(sub_list).issubset(set(test_list)) and len(test_list) > len(
        sub_list
    )

    if is_subset_result:
        difference = [x for x in test_list if x not in sub_list]
    else:
        difference = []

    return is_subset_result, difference
