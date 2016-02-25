# Tong Shu Li
# Created on: 2016-02-24
# Last updated: 2016-02-24
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
