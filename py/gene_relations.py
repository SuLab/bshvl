#!/usr/bin/python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Extractor for gene-gene relations, compiled in a single script.

Refactored by: Tong Shu Li (tongli@scripps.edu)
Last updated: 2016-02-25
"""
import csv
import random
import re
import sys

from itertools import islice

from helper.easierlife import *

from util import is_verb
from util import no_comma
from util import is_neg_word
from util import is_plural_noun

from util import is_gene_list
from util import appear_in_same_doc
from util import no_interact_phrase

from util import get_short_sentences
from util import get_gene_pairs
from util import get_dependency_tree

from dependency_tree import dep_path

#extractor dictionaries
dict_gene_symbols_all = {}
dict_gene_pruned = {}
dict_interact = {}
dict_no_interact = {}
dict_pmid_gene = {}
dict_gene_pmid = {}
dict_drug_names = {}
dict_compound_bio_roles = set()
dict_geneid2name = {}
dict_geneid2official = {}
dict_name2geneid = {}
dict_snowball = {}
dict_exclude_dist_sup = {}
dict_english = {}
dict_abbv = {}
dict_domains = {}
dict_y2h = {}
dict_pmid2plos = {}
dict_gs_docids = set()

def load_dict():

    """
    Name: load_dict
    Input: None
    Return: None

    Store all relevant dictionaries for the gene-gene extractor
    """

    csv.field_size_limit(sys.maxsize)
    GENE_DICT = "/dicts/genes_pruned.tsv"
    GENE_DICT_ALL = "/dicts/genes.tsv"
    NEG_INT_DICT = "/dicts/negatome_combined_stringent_names.txt"
    DRUG_DICT = "/dicts/drugs.tsv"
    DICT_DIALECT = "excel-tab"
    SNOWBALL_DICT = "/dicts/genegene_snowball.txt"
    SUPERVSION_EXCLUDE_DICT = "/dicts/genegene_exclusion_distant_supervision.txt"
    TF_DICT = "/dicts/chea-background.csv"
    PLOS2PMID_BIOGRID_DICT = "/dicts/plos_journals_BioGRID_pmids.txt"

    #GS filter
    INPUT_FILE_GS_SKIP = "/data/plos_journals_dip_mint_pmids.txt"
    INPUT_FILE_10K_SKIP = "/data/plos_docids_sample_10000.txt"

    #gold standard docids
    for x in [x.strip().split("\t")[0] for x in open(BASE_FOLDER + INPUT_FILE_GS_SKIP).readlines()]:
        dict_gs_docids.add(x.rstrip().split(".pdf")[0].lower())

    for x in [x.strip() for x in open(BASE_FOLDER + INPUT_FILE_10K_SKIP).readlines()]:
        dict_gs_docids.add(x.rstrip().split(".pdf")[0].lower())

    #dictionary with all gene symbols, not pruned
    with open(BASE_FOLDER + GENE_DICT_ALL) as tsv:
        r = csv.reader(tsv, dialect=DICT_DIALECT)
        headers = r.next()
        for line in r:

            if line[1] not in dict_geneid2name:
                dict_geneid2name[line[1]] = {}

            if len(line[4]) > 2 and line[4] != "":
                if line[4] not in dict_name2geneid:
                    dict_name2geneid[line[4]] = {}

                dict_gene_symbols_all[line[4].rstrip()] = "symbol"#line[3]
                dict_geneid2name[line[1]][line[4].rstrip()] = 1
                dict_geneid2official[line[1]] = line[4].rstrip()
                dict_name2geneid[line[4].rstrip()][line[1]] = 1
            alt_names = line[6].rstrip()
            alt_names = re.sub("\"", "", alt_names)
            alt_names = re.sub(",$", "", alt_names)
            alt_names = alt_names.split(",")
            for alt in alt_names:
                if len(alt) > 2:
                    if alt not in dict_name2geneid:
                        dict_name2geneid[alt] = {}

                    dict_gene_symbols_all[alt] = "alt name"
                    dict_geneid2name[line[1]][alt] = 1
                    dict_name2geneid[alt][line[1]] = 1

    with open(BASE_FOLDER + GENE_DICT) as tsv:
        r = csv.reader(tsv, dialect=DICT_DIALECT)
        headers = r.next()
        for line in r:

            if len(line[0]) > 2 and line[0] != "" and not re.search("^\d+\.?\d+$",line[0]):
                dict_gene_pruned[line[0]] = "symbol"#line[3]
            alt_names = line[1].rstrip()
            #alt_names = re.sub("\"", "", alt_names)
            #alt_names = re.sub(",$", "", alt_names)
            alt_names = alt_names.split(",")
            for alt in alt_names:
                if len(alt) > 2 and not re.search("^\d+\.?\d+$",alt):
                    dict_gene_pruned[alt] = line[1]



    for l in open(BASE_FOLDER + NEG_INT_DICT):
        ss = l.rstrip().split('\t')
        w1 = ss[3]
        w2 = ss[4]

        if w1 == "NULL" or w2 == "NULL" : continue

        if w1 not in dict_no_interact:
            dict_no_interact[w1] = {}

        dict_no_interact[w1][w2] = 1

    for l in open(BASE_FOLDER + SNOWBALL_DICT):
        ss = l.rstrip().split('\t')
        w1 = ss[0]
        w2 = ss[1]

        if w1 not in dict_interact:
            dict_interact[w1] = {}

        dict_interact[w1][w2] = 1

        if w2 not in dict_interact:
            dict_interact[w2] = {}

        dict_interact[w2][w1] = 1

    for l in open(BASE_FOLDER + SUPERVSION_EXCLUDE_DICT):
        ss = l.rstrip().split("\t")

        if ss[0] not in dict_exclude_dist_sup:
            dict_exclude_dist_sup[ss[0]] = {}

        dict_exclude_dist_sup[ss[0]][ss[1]] = 1

        if ss[1] not in dict_exclude_dist_sup:
            dict_exclude_dist_sup[ss[1]] = {}

        dict_exclude_dist_sup[ss[1]][ss[0]] = 1

    for l in open(BASE_FOLDER + PLOS2PMID_BIOGRID_DICT):
        ss = l.split("\t")
        plos = ss[0]
        pmid = ss[1].rstrip()
        dict_pmid2plos[pmid] = plos

    for l in open(BASE_FOLDER + "/dicts/BIOGRID-ALL-3.2.112.tab.txt"):
        ss = l.split('\t')
        pmids = ss[8]

        if re.search(r";", pmids): print "WARNING"
        if re.search(r";", pmids): log('&&&&&'+"WARNING")
        for pmid in pmids.split(';'):
            if pmid in dict_pmid2plos:
                plos_pmid = dict_pmid2plos[pmid]
                if plos_pmid not in dict_pmid_gene:
                    dict_pmid_gene[plos_pmid] = {}

        y2h = 0
        if ss[6] == "Two-hybrid":
            y2h = 1
        elif ss[6] != "Co-crystal Structure" and ss[6] != "Reconstituted Complex" and ss[6] != "Co-purification":
            continue

        skip_genes = ["p38", "PI3K"]

        if y2h == 0:
            g1 = ss[2]
            g2 = ss[3]
            if g1 not in skip_genes and g2 not in skip_genes:
                if g1 not in dict_interact:
                    dict_interact[g1] = {}
                dict_interact[g1][g2] = 1
                if g2 not in dict_interact:
                    dict_interact[g2] = {}
                dict_interact[g2][g1] = 1

            #alias
            alt_list_1 = ss[4].split("|")
            alt_list_2 = ss[5].split("|")

            for a1 in alt_list_1:
                for a2 in alt_list_2:

                    if a1 not in skip_genes and a2 not in skip_genes:
                        if a1 not in dict_interact:
                            dict_interact[a1] = {}
                        if a2 not in dict_interact:
                            dict_interact[a2] = {}
                        dict_interact[a1][a2] = 1
                        dict_interact[a2][a1] = 1

                        if pmids.rstrip() in dict_pmid2plos:
                            plos_pmid = dict_pmid2plos[pmids.rstrip()]
                            if plos_pmid in dict_pmid_gene:
                                dict_pmid_gene[plos_pmid][a1] = 1
                                dict_pmid_gene[plos_pmid][a2] = 1

            if pmids.rstrip() in dict_pmid2plos:
                plos_pmid = dict_pmid2plos[pmids.rstrip()]
                if plos_pmid in dict_pmid_gene:
                    dict_pmid_gene[plos_pmid][g1] = 1
                    dict_pmid_gene[plos_pmid][g2] = 1
        else:
            g1 = ss[2]
            g2 = ss[3]
            if g1 not in dict_y2h:
                dict_y2h[g1] = {}
            dict_y2h[g1][g2] = 1
            if g2 not in dict_y2h:
                dict_y2h[g2] = {}
            dict_y2h[g2][g1] = 1

            #alias
            alt_list_1 = ss[4].split("|")
            alt_list_2 = ss[5].split("|")

            for a1 in alt_list_1:
                for a2 in alt_list_2:

                    if a1 not in dict_y2h:
                        dict_y2h[a1] = {}
                    if a2 not in dict_y2h:
                        dict_y2h[a2] = {}
                    dict_y2h[a1][a2] = 1
                    dict_y2h[a2][a1] = 1


    with open(BASE_FOLDER + TF_DICT, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[7] == "human":
                if row[1] not in dict_interact:
                    dict_interact[row[1]] = {}
                dict_interact[row[1]][row[3]] = 1

                if row[3] not in dict_interact:
                    dict_interact[row[3]] = {}
                dict_interact[row[3]][row[1]] = 1

    for pmid in dict_pmid_gene:
        for gene in dict_pmid_gene[pmid]:
            if gene not in dict_gene_pmid:
                dict_gene_pmid[gene] = {}
            dict_gene_pmid[gene][pmid] = 1

    for l in open(BASE_FOLDER + "dicts/compounds_bio_roles.txt"):
        for w in l.split(";"):
            dict_compound_bio_roles.add(w.rstrip().lower())


    with open(BASE_FOLDER + DRUG_DICT) as tsv:
        r = csv.reader(tsv, dialect=DICT_DIALECT)
        headers = r.next()
        for line in r:
            if line[5] == 'Drug/Small Molecule':
                if line[1].lower() not in dict_compound_bio_roles:
                    if line[1].lower != "nitric oxide":
                        dict_drug_names[line[1].lower()] = line[1]

    for l in open(BASE_FOLDER + "/dicts/med_acronyms_pruned.txt"):
        word = l.rstrip().split("\t")[0]
        dict_abbv[word] = 1

    for l in open(BASE_FOLDER + "/dicts/words"):
        dict_english[l.rstrip().lower()] = 1

    with open(BASE_FOLDER + "/dicts/smart_domain_list.txt") as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            dict_domains[row[0].rstrip()] = 1


def normalize(word):
    """
    Name: normalize
    Input: word
    Return: normalized word

    Normalization for word characters
    """
    return word.encode("ascii", "ignore").replace("'", '_').replace('{', '-_-').replace('}','-__-').replace('"', '-___-').replace(', ,', ',__')


def normalize_utf(word):
    """
    Name: normalize_utf
    Input: text word
    Return: normalized word

    Replaces common UTF codes with appropriate characters
    """

    word = re.sub('\\xe2\\x80\\x94', '-', word)
    word = re.sub('\\xef\\xac\\x81', 'fi', word)
    word = re.sub('\\xc2\\xb0', "DEGREE", word)
    word = re.sub('\\xe2\\x80\\x99', "\'", word)
    word = re.sub('\\xef\\xac\\x82', "fl", word)
    word = re.sub('\\xc2\\xa3', 'POUND', word)
    word = re.sub('\\xe2\\x80\\x98', "\'", word)
    word = re.sub('\\xe2\\x80\\x9c', "\"", word)
    word = re.sub('\\xe2\\x80\\x9d', "\"", word)
    word = re.sub('\\xe2\\x80\\x93', "-", word)
    word = re.sub('\\xe2\\x80\\x94', "--", word)
    word = re.sub('\\xe2\\x80\\xa6', "...", word)
    word = re.sub('\\xc2\\x82', ',', word)       # High code comma
    word = re.sub('\\xc2\\x84', ',,', word)       # High code double comma
    word = re.sub('\\xc2\\x85', '...', word)     # Tripple dot
    word = re.sub('\\xc2\\x88', '^', word)       # High carat
    word = re.sub('\\xc2\\x91', '\'', word)    # Forward single quote
    word = re.sub('\\xc2\\x92', '\'', word)    # Reverse single quote
    word = re.sub('\\xc2\\x93', '\"', word)    # Forward double quote
    word = re.sub('\\xc2\\x94', '\"', word)    # Reverse double quote
    word = re.sub('\\xc2\\x95', '_', word)
    word = re.sub('\\xc2\\x96', '-', word)       # High hyphen
    word = re.sub('\\xc2\\x97', '--', word)      # Double hyphen
    word = re.sub('\\xc2\\x99', '_', word)
    word = re.sub('\\xc2\\xa0', '_', word)
    word = re.sub('\\xc2\\xa6', '|', word)       # Split vertical bar
    word = re.sub('\\xc2\\xab', '<<', word)      # Double less than
    word = re.sub('\\xc2\\xbb', '>>', word)      # Double greater than
    word = re.sub('\\xc2\\xbc', '1/4', word)     # one quarter
    word = re.sub('\\xc2\\xbd', '1/2', word)     # one half
    word = re.sub('\\xc2\\xbe', '3/4', word)     # three quarters
    word = re.sub('\\xca\\xbf', '\'', word)    # c-single quote
    word = re.sub('\\xcc\\xa8', '', word)        # modifier - under curve
    word = re.sub('\\xcc\\xb1', '', word)         # modifier - under line
    word = re.sub('\\xc2\\xa7', 'CODE', word)

    return word

def extract(doc):
    """
    Name: extract
    Input: Document object
    Return: relation and features for genegene database table

    Extractor code to generate truth label and features
    """

    def get_genes(sentence):
        # list of ambiguous gene symbols to exclude from detected genes
        EXCLUDED_GENES = set([
            "ECM", "EMT", "AML", "CLL", "ALL", "spatial", "PDF", "ANOVA", "MED",
            "gamma", "San", "RSS", "2F1", "ROS", "zeta", "ADP", "ALS", "GEF", "GAP"
        ])

        genes = []
        for word in sentence.words:
            if word.word in dict_gene_symbols_all and word.word not in EXCLUDED_GENES:
                genes.append(word)

        return genes

    def not_brace(s):
        assert isinstance(s, str), "Can't check if {} has a brace!".format(s)
        return s not in "{}"

    def output(label):
        """Remember, the inner function can access the variables in the outer function!"""
        vals = [
            doc.docid, mid1, mid2, w1.word, w2.word,
            label,
            feature, sent_text, "\\N"
        ]

        print "\t".join(vals)
        log("&&&&&{}".format("\t".join(vals)))

#-------------------------------------------------------------------------------

    for sent in get_short_sentences(doc):
        genes = get_genes(sent)
        lemma, deptree = get_dependency_tree(sent)

        words = sent.words

        for w1, w2 in get_gene_pairs(genes):
            assert w1.insent_id < w2.insent_id

            minindex = w1.insent_id
            maxindex = w2.insent_id

            # whitespace?
            ws = [lemma[i] for i in range(minindex+1, maxindex) if no_comma(lemma[i])]
            ## Do not include as candidates ##
            if not set(ws).isdisjoint(set(["while", "whereas", "but", "where", "however"])):
                continue

            features = []

            ############## FEATURE EXTRACTION ##################################

            ####### FEATURE: WORD SEQUENCE BETWEEN MENTIONS AND VERB PATHS #####
            verbs_between = []
            minl_w1 = 100
            minp_w1 = None
            minw_w1 = None
            mini_w1 = None
            minl_w2 = 100
            minp_w2 = None
            minw_w2 = None
            mini_w2 = None
            neg_found = False

            # for each word between the current gene pair
            for i in range(minindex+1, maxindex):
                if is_verb(words[i].pos) and not_brace(words[i].word) and no_comma(words[i].word):
                    p_w1 = sent.get_word_dep_path(minindex, words[i].insent_id)
                    p_w2 = sent.get_word_dep_path(words[i].insent_id, maxindex)

                    if len(p_w1) < minl_w1:
                        minl_w1 = len(p_w1)
                        minp_w1 = p_w1
                        minw_w1 = lemma[i]
                        mini_w1 = words[i].insent_id

                    if len(p_w2) < minl_w2:
                        minl_w2 = len(p_w2)
                        minp_w2 = p_w2
                        minw_w2 = lemma[i]
                        mini_w2 = words[i].insent_id

                    if is_neg_word(lemma[i-1]):
                        if i < maxindex - 2:
                            neg_found = True
                            features.append("NEG_VERB_BETWEEN_with[%s]" % words[i-1].word + "-" + lemma[i])
                    else:
                        verbs_between.append(sent.words[i].lemma)


            ##### FEATURE: HIGH QUALITY PREP INTERACTION PATTERNS ##############
            high_quality_verb = False
            if len(verbs_between) == 1 and not neg_found:
                features.append("SINGLE_VERB_BETWEEN_with[%s]" % verbs_between[0])

                if verbs_between[0] in set(["interact", "associate", "bind", "regulate", "phosporylate", "phosphorylated"]):
                    high_quality_verb = True
            else:
                for verb in verbs_between:
                    features.append("VERB_BETWEEN_with[%s]" % verb)


            if len(ws) == 1 and ws[0] == "and" and minindex > 1:
                if lemma[minindex - 1] in ["of", "between"] and words[minindex - 2].word in ["interaction", "binding"]:
                    if minindex == 2 or minindex > 2 and not is_neg_word(lemma[minindex - 3]):
                        high_quality_verb = True


            if len(ws) == 1 and maxindex < len(words) - 1:
                if ((ws[0] == "-" and lemma[maxindex + 1] == "complex")
                    or (ws[0] == "and" and words[maxindex + 1].word in ["interaction", "interactions"])):
                    high_quality_verb = True


            ##### FEATURE: WORDS BETWEEN MENTIONS ##############################
            if 0 < len(ws) < 7 and set(ws).isdisjoint(set(["{", "}", "\"", "/", "\\", ","])):
                 if " ".join(ws) not in set(["_ and _", "and", "or",  "_ or _"]):
                     features.append("WORDS_BETWEEN_with[%s]" % " ".join(ws))


            ##### FEATURE: 3-GRAM WORD SEQUENCE ################################
            """
            WARNING:

            Mallory messed up this block of code. When checking the 3-gram, she
            accidentally wrote no_comma(ws) instead of no_comma(ws[i+2]).
            Therefore this does not correctly check the 3-gram. To fix this
            bug, change no_comma(ws) to no_comma(ws[i+2]).

            Tong Shu Li
            2016-02-24
            """
            BAD_CHARS = set(["\'", "\"", "{", "}", "[", "]", "-", ","])
            if 4 < len(ws) < 15:
                for i in range(len(ws) - 3 + 1):
                    if BAD_CHARS.isdisjoint(set(ws[i:i+3])):
                        if no_comma(ws[i]) and no_comma(ws[i+1]) and no_comma(ws):
                            features.append("WS_3_GRAM_with[{}]".format("-".join(ws[i:i+3])))


            ##### FEATURE: PREPOSITIONAL PATTERNS #####
            if minindex > 1:
                if sent.words[minindex - 2].lemma.lower() in ["association", "interaction", "complex", "activation", "bind", "binding"]:
                    if sent.words[minindex - 1].word.lower() in ["of", "between"] and ("with" in ws or "and" in ws or "to" in ws) and len(ws) ==1:
                        features.append("PREP_PATTERN[{0}_{1}_{2}]".format(sent.words[minindex-2].lemma.lower(), sent.words[minindex-1].word.lower(), sent.words[minindex+1].word.lower()))
                        high_quality_verb = True


            ##### FEATURE: NEGATED GENES #####
            if sent.words[maxindex-1].lemma in ["no", "not", "neither", "nor"]:
                features.append("NEG_SECOND_GENE[%s]" % sent.words[maxindex - 1].lemma)

            if minindex > 0:
                if sent.words[minindex-1].lemma in ["no", "not", "neither", "nor"]:
                    features.append("NEG_FIRST_GENE[%s]" % sent.words[minindex - 1].lemma)


            if mini_w2 == mini_w1 and mini_w1 != None and len(minp_w1) < 100: # and "," not in minw_w1:
                feature2 = 'DEP_PAR_VERB_CONNECT_with[' + minw_w1 + ']'
                features.append(feature2)

            ##### FEATURE: DEPENDENCY PATH #####
            p = dep_path(deptree, lemma, w1.insent_id, w2.insent_id)
            word1_parent_idx = w1.dep_par
            word1_parent_path = w1.dep_label

            if len(p) < 100:
                try:
                    word1_parent_path = normalize_utf(word1_parent_path)
                    p = normalize_utf(p)
                    p.decode('ascii')
                    norm_p = re.sub(",", "_", normalize(p))

                    if word1_parent_idx == -1:
                        features.append("ROOT_'" + norm_p + "'")

                    else:
                        par_word = sent.words[word1_parent_idx]
                        par_word_lemma = normalize_utf(par_word.lemma)

                        if "," not in par_word_lemma:
                            feature = 'DEP_PAR[' + normalize(par_word_lemma) + '--' + word1_parent_path + '--' + norm_p + ']'
                            features.append(feature)
                except UnicodeDecodeError:
                    pass

            ##### FEATURE: WINDOW FEATURES #####
            bad_char = ["\'", "}", "{", "\"", "-", ",", "[", "]"]
            flag_family = 0

            if minindex > 0:
                if sent.words[minindex - 1].lemma not in bad_char and "," not in sent.words[minindex - 1].lemma:
                    if sent.words[minindex-1].lemma in dict_gene_pruned:
                        features.append('WINDOW_LEFT_M1_1_with[GENE]')
                    else:
                        features.append('WINDOW_LEFT_M1_1_with[%s]' % sent.words[minindex-1].lemma)

            if minindex > 1:
                if sent.words[minindex-2].word in dict_gene_pruned:
                    left_phrase = "GENE"+"-"+sent.words[minindex-1].lemma
                else:
                    left_phrase = sent.words[minindex-2].lemma+"-"+sent.words[minindex-1].lemma

                if sent.words[minindex - 2].lemma not in bad_char and sent.words[minindex - 1].lemma not in bad_char and "," not in left_phrase:
                    features.append('WINDOW_LEFT_M1_PHRASE_with[%s]' % left_phrase)
                elif sent.words[minindex - 2].lemma not in bad_char and "," not in sent.words[minindex - 2].lemma:
                    if sent.words[minindex-2].word in dict_gene_pruned:
                        features.append('WINDOW_LEFT_M1_2_with[GENE]')
                    else:
                        features.append('WINDOW_LEFT_M1_2_with[%s]' % sent.words[minindex-2].lemma)

            if maxindex < len(sent.words) - 1:
                if sent.words[maxindex + 1].lemma not in bad_char and "," not in sent.words[maxindex + 1].lemma:
                    if sent.words[maxindex+1].word in dict_gene_pruned:
                        features.append('WINDOW_RIGHT_M2_1_with[GENE]')
                    else:
                        if sent.words[maxindex+1].lemma in ["family", "superfamily"]: flag_family = 1
                        features.append('WINDOW_RIGHT_M2_1_with[%s]' % sent.words[maxindex+1].lemma)

            if maxindex < len(sent.words) - 2:
                if sent.words[maxindex + 2].word in dict_gene_pruned:
                    right_phrase = "GENE"+"-"+sent.words[maxindex+1].lemma
                else:
                    right_phrase = sent.words[maxindex+2].lemma+"-"+sent.words[maxindex+1].lemma

                if sent.words[maxindex + 2].lemma not in bad_char and sent.words[maxindex + 1].lemma not in bad_char and "," not in right_phrase:
                    features.append('WINDOW_RIGHT_M2_PHRASE_with[%s]' % right_phrase)
                elif sent.words[maxindex + 2].lemma not in bad_char and "," not in sent.words[maxindex + 2].lemma:
                    if sent.words[maxindex + 2].word in dict_gene_pruned:
                        features.append('WINDOW_RIGHT_M2_2_with[GENE]')
                    else:
                        features.append('WINDOW_RIGHT_M2_2_with[%s]' % sent.words[maxindex+2].lemma)

            if len(ws) > 4:
                if sent.words[minindex + 1].lemma not in bad_char and "," not in sent.words[minindex + 1].lemma:
                    if sent.words[minindex + 1].word in dict_gene_pruned:
                        features.append('WINDOW_RIGHT_M1_1_with[GENE]')
                    else:
                        if sent.words[minindex+1].lemma in ["family", "superfamily"]: flag_family = 1
                        features.append('WINDOW_RIGHT_M1_1_with[%s]' % sent.words[minindex+1].lemma)

                if sent.words[minindex+2].word in dict_gene_pruned:
                    m1_right_phrase = "GENE"+"-"+sent.words[minindex+1].lemma
                else:
                    m1_right_phrase = sent.words[minindex+2].lemma+"-"+sent.words[minindex+1].lemma

                if sent.words[minindex + 2].lemma not in bad_char and sent.words[minindex + 1].lemma not in bad_char and "," not in m1_right_phrase:
                    features.append('WINDOW_RIGHT_M1_PHRASE_with[%s]' % m1_right_phrase)
                elif sent.words[minindex + 2].lemma not in bad_char and "," not in sent.words[minindex + 2].lemma:
                    if sent.words[minindex+2].word in dict_gene_pruned:
                        features.append('WINDOW_RIGHT_M1_2_with[GENE]')
                    else:
                        features.append('WINDOW_RIGHT_M1_2_with[%s]' % sent.words[minindex+2].lemma)

                if sent.words[maxindex - 1].lemma not in bad_char and "," not in sent.words[maxindex - 1].lemma:
                    if sent.words[maxindex - 1].word in dict_gene_pruned:
                        features.append('WINDOW_LEFT_M2_1_with[GENE]')
                    else:
                        features.append('WINDOW_LEFT_M2_1_with[%s]' % sent.words[maxindex-1].lemma)

                if sent.words[maxindex-2].word in dict_gene_pruned:
                    m2_left_phrase = "GENE"+"-"+sent.words[maxindex-1].lemma
                else:
                    m2_left_phrase = sent.words[maxindex-2].lemma+"-"+sent.words[maxindex-1].lemma

                if sent.words[maxindex - 2].lemma not in bad_char and sent.words[maxindex - 1].lemma not in bad_char and "," not in m2_left_phrase:
                    features.append('WINDOW_LEFT_M2_PHRASE_with[%s]' % m2_left_phrase)
                elif sent.words[maxindex - 2].lemma not in bad_char and "," not in sent.words[maxindex - 2].lemma:
                    if sent.words[maxindex-2].word in dict_gene_pruned:
                        features.append('WINDOW_LEFT_M2_2_with[GENE]')
                    else:
                        features.append('WINDOW_LEFT_M2_2_with[%s]' % sent.words[maxindex-2].lemma)


            ##### FEATURE: DOMAIN ##############################################
            DOMAIN_WORDS = set([
                "domain", "domains", "motif", "motifs", "site", "sites",
                "region", "regions", "sequence", "sequences", "elements"
            ])

            found_domain = False
            if ((minindex > 0 and words[minindex + 1].word in DOMAIN_WORDS)
                or (maxindex < len(words) - 1 and words[maxindex + 1].word in DOMAIN_WORDS)):
                features.append("GENE_FOLLOWED_BY_DOMAIN_WORD")
                found_domain = True


            ##### FEATURE: PLURAL GENES ########################################
            """
            WARNING:

            There is a typo for the second feature name. There should be no
            right parenthesis after "NOUN". This typo was kept for consistency.
            """
            if minindex > 0 and is_plural_noun(words[minindex + 1].pos):
                features.append('GENE_M1_FOLLOWED_BY_PLURAL_NOUN_with[%s]' % words[minindex + 1].word)
            elif maxindex < len(words) - 1 and is_plural_noun(words[maxindex + 1].pos):
                features.append('GENE_M2_FOLLOWED_BY_PLURAL_NOUN)_with[%s]' % words[maxindex + 1].word)


            ##### FEATURE: GENE LISTING ########################################
            if len(ws) > 0:
                # why do we need to add the same feature twice?
                # is this intentional?
                if ws[-1] == "_" and is_gene_list("A", ws, dict_gene_symbols_all):
                    features.append("GENE_LISTING")

                if ws[-1] == "," and is_gene_list("B", ws, dict_gene_symbols_all):
                    features.append("GENE_LISTING")


            #### GENERATE FEATURE ARRAY FOR POSTGRESQL #####
            feature = "{" + ','.join(features) + '}'

            mid1 = "{}_{}_{}".format(doc.docid, sent.sentid, w1.insent_id)
            mid2 = "{}_{}_{}".format(doc.docid, sent.sentid, w2.insent_id)

            ############## DISTANT SUPERVISION #################################

            sent_text = sent.__repr__()
            if sent_text.endswith("\\"):
                sent_text = sent_text[:-1]

            if w1.word in dict_exclude_dist_sup and w2.word in dict_exclude_dist_sup[w1.word]:
                output("\\N")
            elif words[0].word == "Abbreviations" and words[1].word == "used":
                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                    output("false")

                output("\\N")
            elif w1.word not in dict_abbv and w1.word not in dict_english and w2.word not in dict_english and w2.word not in dict_abbv and w1.word not in dict_domains and w2.word not in dict_domains:
                if w1.word in dict_interact and w2.word in dict_interact[w1.word] and "mutation" not in sent_text and "mutations" not in sent_text and "variant" not in sent_text and "variants" not in sent_text and "polymorphism" not in sent_text and "polymorphisms" not in sent_text:

                    if not found_domain and flag_family == 0:
                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            output("true")

                        output("\\N")
                    else:
                        output("\\N")
                else:
                    if w1.word in dict_no_interact and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws):
                        if w2.word in dict_no_interact[w1.word] and not high_quality_verb: 
                            if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                output("false")

                            output("\\N")
                        else:
                            output("\\N")
                    elif w2.word in dict_no_interact and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws):
                        if w1.word in dict_no_interact[w2.word] and not high_quality_verb:
                            if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                output("false")

                            output("\\N")
                        else:
                            output("\\N")
                    elif appear_in_same_doc(w1.word, w2.word, dict_pmid_gene) and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws ) and not high_quality_verb:
                        # Negative Example: Mention appear in KB in same doc, but no interaction extracted in KB

                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            output("false")

                        output("\\N")
                    elif no_interact_phrase(ws) and not high_quality_verb:
                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            output("false")

                        output("\\N")
                    elif w1.ner == "Person" or w2.ner == "Person":
                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            output("false")

                        output("\\N")
                    elif random.random() < .08 and not high_quality_verb:
                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            output("false")

                        output("\\N")
                    else:
                        output("\\N")
            else:
                output("\\N")

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    log("START of gene_relations.py!")

    random.seed(12345) # to keep logs consistent

    load_dict()

    for row in fileinput.input():
        doc = deserialize(row.rstrip('\n'))
        try:
            extract(doc)
        except:
            pass
