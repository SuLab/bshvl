# Tong Shu Li
# Created on: 2016-02-24
# Last updated: 2016-02-25
"""
Simple utilities for helping the gene-gene extractor.
"""

import re

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
