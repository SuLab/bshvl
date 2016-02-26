# Tong Shu Li
# Created on: 2016-02-24
# Last updated: 2016-02-25
"""
Simple utilities for helping the gene-gene extractor.
"""
import re
import string

from itertools import combinations

def is_verb(s):
    assert isinstance(s, str)
    return re.search(r'VB\w*', s) is not None

def no_comma(s):
    return "," not in s

def remove_underscores(s):
    assert isinstance(s, str)
    return s.replace("_", "")

def is_neg_word(s):
    assert isinstance(s, str)
    return s in set(["no", "not", "neither", "nor"])

def is_plural_noun(s):
    assert isinstance(s, str)
    return s in ["NNS", "NNPS"]

#-------------------------------------------------------------------------------

def is_gene_list(kind, ws, dict_gene_symbols_all):
    assert kind in ["A", "B"], "Invalid gene list type."

    if kind == "A":
        vals = ["_", ",", "_"]

        return all(
            vals[i%4] == v if i%4 < 3 else v in dict_gene_symbols_all
            for i, v in enumerate(ws)
        )

    # gene list type B
    return all(
        v == "," if i%2 == 0 else v in dict_gene_symbols_all
        for i, v in enumerate(ws)
    )

#-------------------------------------------------------------------------------

def get_short_sentences(document):
    MAX_WORDS_IN_SENTENCE = 50
    for sentence in document.sents:
        if len(sentence.words) <= MAX_WORDS_IN_SENTENCE:
            yield sentence

def get_gene_pairs(genes):
    for geneA, geneB in combinations(genes, 2):
        if remove_underscores(geneA.word) != remove_underscores(geneB.word):
            yield (geneA, geneB)

def get_dependency_tree(sentence):
    lemma = [word.lemma for word in sentence.words]
    deptree = {
        word.insent_id: {"label": word.dep_label, "parent": word.dep_par}
        for word in sentence.words
    }
    return (lemma, deptree)

#-------------------------------------------------------------------------------

def start_cap_letter(s):
    return len(s) > 0 and s[0] in string.ascii_uppercase

def appear_in_same_doc(geneA, geneB, dict_gene_pmid):
    if start_cap_letter(geneA) and start_cap_letter(geneB) and geneA in dict_gene_pmid:
        for pmid in dict_gene_pmid[geneA]:
            if geneB in dict_pmid_gene[pmid]:
                return True

    return False

def slice(a_list, N):
    for i in range(len(a_list) - N + 1):
        yield a_list[i:i+N]

def no_interact_phrase(ws):
    # check if not interact/bind phrase is in ws and not just the words
    """
    WARNING:

    This code is left intentionally wrong. The original code by Mallory fails
    to check the last possible pair of words in the ws list due to the incorrect
    out-of-bounds checking condition "j + 1 < len(ws) - 1".

    To check all of the two word pairs in the ws list, she should have written
    "j + 1 < len(ws)" instead.

    To replicate this bug, we only slice ws[:-1] instead of the correct ws.

    Tong Shu Li
    2016-02-25
    """
    return any(
        i == "not" and j in set(["bind", "interact", "interacts"])
        for i, j in slice(ws[:-1], 2)
    )

def no_interact_words(ws):
    KEYWORDS = set([
        "bind", "binds", "bound",
        "interact", "interacts", "interacted",
        "associates", "associated",
        "complex"
    ])

    return set(ws).isdisjoint(KEYWORDS)
