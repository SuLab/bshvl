"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Refactored by: Tong Shu Li (tongli@scripps.edu)
Last updated: 2016-02-25

This file parses the dependency tree and returns the path to the common root.
"""
from itertools import islice

def get_path(node, deptree):
    """Given a starting node, walk the dependency tree and return the path walked."""
    path = []
    while node in deptree:
        path.append({
            "current": node,
            "parent": deptree[node]["parent"],
            "label": deptree[node]["label"]
        })

        node = deptree[node]["parent"]

    path.append({
        "current": node,
        "parent": -1,
        "label": "ROOT"
    })

    return path

def get_common_root(pathA, pathB):
    """This function __should__ compare two lists of potentially unequal length
    and look for the longest suffix list that is the same in the two lists.
    It should return the first element in that suffix list.

    Example: if given [1, 2, 3, 4, 4, 5] and [1, 1, 1, 3, 4, 5], the function
    should return the element 4, because [4, 5] is the longest suffix list.

    WARNING:

    However, in the original code for this function, Emily Mallory made an
    off-by-one error. The original code does not check the first element of
    the second list, and therefore the results are one smaller than what they
    should be.

    The original code was:

    commonroot = None
    for i in range(0, len(path1)):
        j = len(path1) - 1 - i
        if -i-1 <= -len(path2) or path1[j]["current"] != path2[-i-1]["current"]:
            break

        commonroot = path1[j]["current"]

    The first part of the if statement is used to prevent accessing path2
    out of bounds when the list path1 is longer than path2. However, when
    -i-1 == len(path2), this is still within bounds of path2. As in,
    list[0] == list[-len(list)]. However, the original code breaks at this
    point instead of checking the last element due to the <= operator. The
    correct operator should use < and not <=.

    Since the refactored code is trying to match the original code exactly,
    this off-by-one bug has been intentionally left the same.

    To remove this off-by-one error, simply change "len(pathB) - 1" to
    "len(pathB)" in the while loop below.

    Tong Shu Li
    Last updated: 2016-02-24
    """
    pos = 1
    while (pos <= len(pathA) and pos <= len(pathB) - 1
        and pathA[-pos]["current"] == pathB[-pos]["current"]):
        pos += 1

    return pathA[-pos+1]["current"] if pos > 1 else None

def walk_path(path, direction, common_root, lemma):
    assert direction in ["left", "right"], "Wrong walk path direction."

    route = []
    arrow = "->" if direction == "right" else "<-"

    if path[0]["current"] != common_root:
        route += ["--", path[0]["label"], arrow]

        for node in islice(path, 1, None):
            if node["current"] == common_root:
                break

            route += ["--", node["label"], arrow]
            if node["parent"] != common_root and node["parent"] != -1:
                route.append(lemma[node["parent"]].lower())

        route.append("|")

    return "".join(route if direction == "right" else route[::-1])

def dep_path(deptree, lemma, start1, start2):
    """
    Name: dep_path
    Input: dependency tree, lemma, and start word postions
    Return: Dependency path between two words

    Simplified version of Sentence class dependency path code
    """
    # length of the dependency tree should be the length of the sentence, which
    # should always be greater than zero and <= 50
    assert 0 < len(deptree) <= 50, "Deptree length out of bounds."

    path1 = get_path(start1, deptree)
    path2 = get_path(start2, deptree)

    common_root = get_common_root(path1, path2)

    left_path = walk_path(path1, "right", common_root, lemma)
    right_path = walk_path(path2, "left", common_root, lemma)

    if common_root is None:
        return left_path + "NONEROOT" + right_path

    if common_root == start1 or common_root == start2:
        return left_path + "SAMEPATH" + right_path

    return left_path + lemma[common_root].lower() + right_path
