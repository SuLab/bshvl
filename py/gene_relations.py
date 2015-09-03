############################################
# Author: Emily K Mallory and Ce Zhang
# Date: 08/24/15
# Contact: emily.mallory@stanford.edu
#
# Extractor for gene-gene relations, compiled
#   in a single script
# 
############################################

#!/usr/bin/python
import sys
from helper.easierlife import *
import csv
import re
import random

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

######
# Name: dep_path
# Input: dependency tree, sentence, lemma, and start and end word postions
# Return: Dependency path between two words
#
# Simplified version of Sentence class dependency path code
######

def dep_path(deptree, sent, lemma, start1, end1, start2, end2):
    if len(deptree) > 0:
        path1 = []
        end = end1 - 1
        ct = 0
        while True:
            ct = ct + 1
            if ct > 100:
                break
            if end not in deptree:
                path1.append({"current":end, "parent": -1, "label":"ROOT"})
                break
            path1.append({"current":end, "parent": deptree[end]["parent"], "label":deptree[end]["label"]})
            end = deptree[end]["parent"]

        path2 = []
        end = end2 - 1
        ct = 0
        while True:
            ct = ct + 1
            if ct > 100:
                break
            if end not in deptree:
                path2.append({"current":end, "parent": -1, "label":"ROOT"})
                break
            path2.append({"current":end, "parent": deptree[end]["parent"], "label":deptree[end]["label"]})
            end = deptree[end]["parent"]

        commonroot = None
        for i in range(0, len(path1)):
            j = len(path1) - 1 - i
            #plpy.notice(path1[j])  
            #plpy.notice(path2[-i-1])  
            if -i-1 <= -len(path2) or path1[j]["current"] != path2[-i-1]["current"]:
                break
            commonroot = path1[j]["current"]

        left_path = ""
        lct = 0
        for i in range(0, len(path1)):
            lct = lct + 1
            if path1[i]["current"] == commonroot:
                break
            if path1[i]["parent"] == commonroot or path1[i]["parent"]==-1:
                left_path = left_path + ("--" + path1[i]["label"] + "->" + '|')
            else:
                w = lemma[path1[i]["parent"]].lower()
                if i == 0: 
                    w = ""
                left_path = left_path + ("--" + path1[i]["label"] + "->" + w)

        right_path = ""
        rct = 0
        for i in range(0, len(path2)):
            rct = rct + 1
            if path2[i]["current"] == commonroot:
                break
            if path2[i]["parent"] == commonroot or path2[i]["parent"]==-1:
                right_path = ('|' + "<-" + path2[i]["label"] + "--") + right_path
            else:
                w = lemma[path2[i]["parent"]].lower()
                if i == 0:
                    w = ""
                right_path = (w + "<-" + path2[i]["label"] + "--" ) + right_path

        path = ""
        if commonroot == end1-1 or commonroot == end2-1:
            path = left_path + "SAMEPATH" + right_path
        else:
            if commonroot != None:
                path = left_path + lemma[commonroot].lower() + right_path
            else:
                path = left_path + "NONEROOT" + right_path
        if path != "":
            return path
        else:
            return None

######
# Name: load_dict
# Input: None
# Return: None
#
# Store all relevant dictionaries for the gene-gene extractor
######

def load_dict():

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

    bowdict = {'-agent--interfere': 2, 'lasted': 6, 'gatekeeper': 2, 'lrp-binding': 2, 'assembles': 4, 'abundances': 2, 'ad--dobj-': 2, '-dobj--jnk': 6, 'controversial': 2, 'for--maintain': 2, '-agent--lp': 2, 'assembled': 4, 'astressin': 4, 'replication--nsubj-': 2, 'at--inhibitor': 2, 'ventral': 2, 'c-flip': 6, 'c-jun-n-terminal': 6, 'polyubiquitin': 6, 'of--hiv-': 2, 'recycle': 2, 'myocardin-related': 2, 'deoxycytidine': 4, '-agent--arabidopsis': 2, 'immuno': 4, '-advcl--establish': 3, 'delc': 4, 'tat-activated': 2, 'delg': 4, 'histone-binding': 2, 'regional': 2, 'connects': 6, 'delt': 25, '-advcl--generate': 2, 'co-precipitation--nsubjpass-': 2, 'bringing': 10, 'monkeys': 2, 'bms-': 4, 'to--member': 2, 'ki-ras': 4, 'mean--ccomp-': 2, 'profile--ccomp-': 3, 'of--dsb': 4, 'without--have': 4, 'method--dobj-': 2, 'disc': 4, '-advcl--recruit': 2, 'perivitelline': 2, 'glued': 2, 'run-on': 10, '-dep--sensitive': 6, 'disclose': 2, 'senescence': 4, 'poorly-differentiated': 2, '-nsubj--fibroblast': 2, 'in--support': 5, 'naip': 2, 'direct': 9, 'accumulation--agent-': 2, '-dobj--ser': 2, 'erk-mediated': 4, '-ccomp--mechanism': 2, 'homeologous': 6, 'ghrs': 2, 'sterile': 4, 'incompatible': 2, 'tethers': 2, 'regulator--rcmod-': 2, '-ccomp--independent': 3, 'and--over-expression': 3, 'lie--rcmod-': 3, 'fingers': 2, 'liberated': 4, '-phosphorylation': 2, 'mviic': 2, 'cancer--conj': 2, '-rcmod--lead': 8, '-dep--require': 3, '-nsubj--rsk': 2, 'med': 8, 'mec': 2, 'menage': 2, 'mbl-associated': 2, 'mek': 15, 'with--d-site': 2, 'alpha-globin': 2, 'by--sek': 2, 'two--nsubj-': 5, 'hyperoxia': 2, 'xfs': 6, 'of--trail': 5, 'smrad': 4, 'fos-b': 2, 'amino-terminal': 24, 'polyphenol': 2, 'drs': 4, 'early-onset': 2, 'dmap': 4, 'site--rcmod-': 2, 'specific--dep-': 7, 'growth-related': 18, 'to-ala': 6, 'sequestering': 2, 'j-ghrhr': 2, 'adolescents': 2, 'neu-transfected': 2, 'c-myc-induced': 4, '-amod--determine': 2, 'rsai': 2, 'det': 2, '-ser-deficient': 2, 'carm': 2, '-tmod--alone': 2, 'hypophysiotropic': 2, '-advmod--remarkably': 3, 'decline--nsubj-': 2, 'chaperones': 2, '-amod--analogous': 2, 'circumvent': 2, '-dobj--atm': 2, 'convertase': 3, 'conformation--dobj-': 2, 'nls-deficient': 4, 'beta-pix': 2, '-dobj--entry': 2, 'substitutions': 2, 'er-': 2, 'parkin': 2, 'of--copy': 2, 'while--dissociation': 3, 'fanconi': 6, 'autophosphorylation--nsubjpass-': 3, 'successful': 3, '-xcomp--bring': 5, 'classify--parataxis-': 2, 'subepithelial': 2, 'localisation--dobj-': 2, 'compensate--dep-': 2, '-amod--signal-regulated': 2, 'including--homodimerization': 2, 'confer--rcmod-': 2, 'podosome': 4, 'multiplexes': 2, 'xpa--nsubj-': 2, 'kit-mediated': 2, 'diphosphorylated': 2, 'account--ccomp-': 2, 'erk--prep': 2, '-rcmod--unable': 2, 'translocation--nsubj-': 5, 'and--potentiate': 2, 'halphacgrp': 2, 'phospho-cdk': 4, 'family--xcomp-': 4, 'alpha-adrenergic-stimulated': 4, 'open--xcomp-': 2, 'bisindolylmaleimide': 4, 'lipopeptides': 2, 'and--impair': 3, 'elimination--nsubj-': 2, 'degenerate': 2, '-xcomp--prevent': 2, 'base-line': 2, 'pro-apoptotic': 4, 'hog': 5, 'endometrial': 6, 'cilia': 2, 'regrowth': 2, 'taf--nsubj-': 2, 'tpa-stimulated': 2, 'dsb--prep': 4, 'and--recognize': 2, 'downregulate--parataxis-': 2, 'coactivating': 2, '-tmod--inhibitor': 4, 'phosphoserine-': 4, 'ids': 2, 'minority': 2, 'transcription--dep-': 4, 'mesenteries': 2, 'withdrawal--nsubj-': 8, 'acetylation-induced': 2, 'phasic': 2, 'mechanism--ccomp-': 2, 'mesf': 2, 'of--gp': 3, 'attempt': 2, 'acid-treated': 2, 'that--dobj-': 3, '-nsubjpass--autophosphorylation': 3, 'alteration--nsubjpass-': 3, 'ttf-': 2, 'delinst': 4, 'with--sb': 2, 'decoy': 4, 'neuropathological': 2, 'rsv-tat': 8, 'tgfalpha': 12, '-appos--snp': 3, '-day-old': 2, 'lysines': 2, '-agent--fact': 2, 'in--woman': 4, 'array--prep': 4, 'depletion--nsubj-': 3, 'amphipathic': 2, 'chromosome-linked': 2, 'apoptosis-inducing': 4, 'innovative': 8, 'april---': 3, 'fiv': 4, 'p-her-': 2, 'than--complex': 4, '-nsubj--profile': 3, 'as-odns': 2, '-ccomp--compensate': 4, 'v-snare': 2, 'modify--conj': 2, 'reprogramming': 4, 'compelling': 2, 'with--up-regulate': 2, 'hhv-': 2, 'maximal': 2, 'inter-sh': 2, '-sek': 2, '-appos--overexpression': 2, 'neurological': 2, 'oxidant': 4, 'wehi-': 3, 'decrease--partmod-': 2, 'convertase--prep': 2, '-acomp--mdm': 2, 'pre-incubated': 2, 'amidating': 2, '-conj--model': 2, 'd-site': 4, 'bcl-xs': 2, 'jnk-specific': 6, 'four-fold': 12, 'c-jun-dependent': 2, 'array--dobj-': 2, 'gene-targeted': 4, 'signal--dep-': 5, 'subset--dobj-': 5, 'male-specific': 2, 'mdm': 69, 'molecule-associated': 2, 'versus--': 2, '-agent--sb': 2, 'kinase-sensitive': 4, 'but--phosphorylate': 2, 'gating': 2, 'form--dep-': 2, 'protein-truncating': 2, 'rt': 11, 'and--related': 2, 'clade': 2, 'kinase-specific': 4, 'ras-induced': 4, 'smn': 42, 'independent--nsubj-': 3, 'hec-': 2, 'adenylate': 18, 'sma': 14, 'trif-related': 2, 'implicate--parataxis-': 3, 'trans-interaction': 5, 'ntr': 4, 'homologue--dobj-': 2, 'ascending': 2, 'aporeceptors': 2, 'by--accumulation': 2, '-dibutyrate': 8, 'angiomyolipoma--prep': 2, 'of--mdm': 3, 'kinase--agent-': 3, 'fret': 2, 'depolarizing': 2, 'baff': 9, 'anti-t': 2, 'virus--nsubjpass-': 2, 'fancc--prep': 7, '-csubj--contain': 2, 'replication-associated': 2, 'harvey-ras': 2, '-acomp--pol': 3, 'fas-dd': 2, 'separable': 6, 'demonstrate--conj': 8, 'punctate': 6, 'gamma-radiation': 2, 'trka': 6, 'nmdar': 4, 'vittoria': 2, 'non-informative': 8, 'on--surface': 2, 'mineralocorticoids': 2, '-nsubj--step': 2, '-dobj--hyper-proliferation': 2, 'necl-': 6, 'importance--appos-': 2, 'b-inducing': 4, 'tos': 2, 'co-regulated': 2, 'wrn--prep': 2, 'phospho-c-jun': 4, 'needed': 8, 'o-linked': 2, 'too': 4, 'of--mutant': 3, 'cc-type': 2, 'including--p': 2, 'dilated': 2, 'pif-pocket-interacting': 2, 'molecular-weight': 2, 'estimates': 2, 'ischemia-induced': 4, 'tata-binding': 10, 'rsv-tat-induced': 2, 'reinforced': 2, 'mutation-associated': 2, '-partmod--inactivate': 5, 'phosphorylation-deficient': 16, '-dobj--egf': 2, 'over-express': 2, 'on--presence': 2, '-rcmod--propose': 5, 'peculiar': 2, 'coassemble': 3, 'metaphase': 8, 'facilitate--rcmod-': 2, 'entry--dobj-': 2, 'syncytium-inducing': 4, '-nsubj--importance': 6, 'modulators': 8, 'smn--tmod-': 2, 'epithelial-to-mesenchymal': 2, 'paved': 2, 'jnk-derived': 2, 'and--bind': 6, 'dna-protein': 2, 'nbn': 2, 'bridge': 12, 'rad': 66, 'by--overexpression': 2, 'core-binding': 2, 'uv--prep': 2, 'plzf': 4, '-appos--melanoma': 2, 'with--pd': 3, 'limitation--dobj-': 2, 'superinduction': 2, 'raw': 2, 'trans-interacts': 2, '-nsubjpass--co-precipitation': 2, 'and--destroy': 2, 'proliferation-inducing': 2, 'and--akt': 2, 'initiate--xcomp-': 2, 'eq': 2, 'of--permeability': 3, 'egfr-mrna': 4, 'predominantly': 6, 'envelope-expressing': 2, 'from--lack': 2, 'multifactor-dimensionality': 4, 'sapk--prep': 2, 'snp--prep': 3, 'wild-type-conformation': 2, 'chemotropic': 6, 'membrane-localized': 4, '-to-v': 2, 'predisposing': 6, 'fact--agent-': 2, 'betacatenin-tcf': 2, 'complexed': 12, 'tyk': 2, 'model--conj-': 2, 'bso': 4, 'cycloheximide': 10, 'cool-': 2, '-dobj--form': 6, 'laagq': 2, 'block--rcmod-': 3, 'hyperproliferative': 18, 'gtpase--prep': 2, 'flox': 4, 'observation': 7, 'dominant-interfering': 10, 'doc': 2, 'bnip': 2, 'cysteine-independent': 2, 'introns': 2, 'epistasis--prep': 2, 'as--adapter': 4, 'itself--ccomp-': 3, 'derive--csubj-': 2, 'cv-n': 4, 'note--parataxis-': 2, 'sphingomyelinase': 4, 'of--pifithrin-alpha': 2, 'pkr': 2, 'd-': 8, '-agent--level': 2, 'manifestation': 2, 'pathway--tmod-': 2, 'serum-stimulated': 2, 'delactt': 4, 'subtle': 4, 'cyclase-activating': 30, 'prevent--xcomp-': 2, 'to--ability': 2, 'after--inhibition': 2, 'palindrome': 2, 'guidance': 2, 'kinase': 458, 'hrs': 2, 'mpirh': 2, 'dnmt': 21, 'unleashed': 4, 'kainate-induced': 8, 'by--phosphorylate': 5, 'judged': 8, 'allo-proliferation': 2, 'tropic': 2, 'with--change': 3, 'sauvagine': 2, 't-til': 2, 'transactivator': 3, 'subcloned': 2, '-mmp': 6, 'manner--nsubj-': 2, '-dep--mapk': 2, 'altogether': 2, 'co-activator': 9, 'with--association': 3, 'dr': 15, '-nsubj--pd': 4, 'block--dep-': 2, 'infrequent': 4, 'beads': 4, 'mnli': 2, 'domain-mediated': 2, 'accordance': 10, 'bax': 27, 'peroxide': 2, 'emt': 2, 'bad': 7, '-rcmod--detect': 2, 'versus-': 2, 'fertility': 10, 'adapters': 8, 'hrs-binding': 2, 'y-box': 2, 'antibody--nsubj-': 2, 'capan-': 4, '-nsubj--co-association': 2, 'sra-': 2, 'nut': 4, '-ccomp--consist': 3, 'synergize--ccomp-': 3, 'cytokine-dependent': 2, 'rock-independent': 2, 'of--mek': 3, 'nub': 6, 'hydrocarbon': 6, 'nua': 2, 'lps-binding': 4, 'mut': 6, 'interactor': 5, '-dep--recruitment': 2, 'ad-nbk-induced': 2, 'feline': 4, 'methyltransferase--tmod-': 2, 'kinase-deficient': 2, 'cooperation': 21, 'hmdm': 4, 'c-jun': 111, 'deltaegf-transfected': 2, 'wts': 2, '-nsubj--oncoprotein': 2, 'of--drug': 2, '-jnk': 8, 'herg': 4, 'tbp-': 4, 'biphasically': 4, 'phosphorylate--ccomp-': 2, 'to--transcription': 2, 'with--another': 4, 'fragment--dobj-': 2, 'generates': 2, '-xcomp--capable': 6, 'phospho': 4, 'degree--prep': 7, 'country': 4, 'extracts': 16, 'adhesion--prep': 3, 'distinction': 4, 'cell-cell': 8, 'rptk': 2, 'tone': 20, 'ethanol': 2, 'detectable--ccomp-': 3, 'and--consequently': 2, 'without--that': 3, 'compose--rcmod-': 2, 'kappab-dependent': 6, 'inactivate': 3, 'dna-bound': 2, '-rcmod--abrogate': 2, 'tif': 3, 'binge': 2, 'nf-interleukin': 2, 'subfamily--rcmod-': 2, '-rcmod--confer': 2, 'sl-': 2, 'vinblastine-induced': 4, 'incompetent': 2, 'quinoxalinecarboxamide': 2, 'stresses': 2, 'placentas': 2, 'mutant--dep-': 4, 'c-jun--pobj-': 2, 'ore-like': 2, 'erse': 2, '-dep--only': 2, 'tumor-suppressor': 4, 'chk-': 6, 'enhance--parataxis-': 3, 'at--location': 2, 'lifr': 6, 'harboured': 4, 'scntfr': 2, 'to--binding': 4, 'fviii-vwf': 2, 'due--conj-': 2, 'comprise--rcmod-': 2, 'ikk-phosphorylated': 2, 'tand': 2, 'two-fold': 2, 'alone--tmod-': 2, 'slk': 2, 'aural': 2, 'helodermin': 4, 'three--dep-': 2, 'with--tbk': 2, 'ig-like': 8, 'raralpha': 9, 'cand': 8, 'on--dep-': 3, 'with--level': 2, 'conditioned': 8, 'skipping': 4, 'r-deficient': 2, '-rcmod--titrate': 2, 'tumorigenesis': 10, 'erbitux': 2, 'ie': 12, 'polypeptide--prep': 2, '-appos--regulation': 4, 'bhlh-pas': 4, 'lef-': 2, '-xcomp--arrest': 3, 'no--prep': 2, 'cumulative': 2, '-npadvmod--activity': 3, 'c-jun--agent-': 2, 'd-phe': 2, 'accumulating': 4, 'lithium': 4, 'modulator': 23, 'hepatogenesis': 2, 'evolve': 2, 'chromosome-positive': 2, 'nic': 4, 'macrophage-colony': 2, 'nix': 2, 'and--ki-ra': 2, 'nip': 2, 'whereby': 10, 'nf-at-ap-': 4, 'rfxank--dobj-': 2, 'drip': 4, 'to--enhancement': 3, 'by--irradiation': 2, 'client': 2, '-ccomp--profile': 3, 'thf': 8, 'cytosine': 4, 'nf-kappa': 14, '-ccomp--link': 2, 'and--reveal': 2, 'trif--dobj-': 2, 'charcot-marie-tooth': 2, 'scoparone': 2, '-appos--prb': 12, 'bronchospasm': 2, 'wip--tmod-': 2, '-nsubj--sb': 7, '-mdm': 6, 'regarded': 2, 'coactivate--dep-': 2, 'for--growth': 2, 'neisserial': 8, 'from--those': 5, 'stoichiometry': 3, 'intraluminal': 6, 'nocodazole-induced': 2, 'hat': 6, 'mncs': 2, 'of--presence': 2, 'splice-sites': 2, 'rfxb': 4, 'lcp': 2, 'in--enhancement': 2, 'antiestrogens': 2, 'present--prepc': 2, 'activity--npadvmod-': 3, 'gabarap': 2, 'arac': 2, 'kcl': 2, 'rnr': 4, '-pobj--modulation': 2, 'mapk-specific': 6, 'give--partmod-': 3, 'corepress': 2, 'apotag': 2, 'temporally': 2, '-nsubj--five': 3, '-kcne': 2, 't-tropic': 4, 'msi-positive': 2, 'arbitrary': 2, '-wt': 2, 'trail--nsubj-': 2, '-pobj--whi-p': 2, 'fsh-r-mediated': 2, 'for--enhancement': 2, '-dep--specific': 7, '-nsubj--disruption': 4, 'begin': 2, 'elderly': 2, 'unnecessary': 4, 'polycystin-': 4, 'necessary--amod-': 6, 'permeabilization': 4, 'desensitize': 2, 'plscr': 3, 'hdaci': 2, 'rasgap': 6, 'neomycin': 2, '-dep--rpsk': 2, 'fancf': 36, 'monoubiquitination--prep': 3, 'only--': 2, 'lxxll': 2, 'to--kinase': 2, 'fog': 4, 'myotilin': 2, 'translocator': 20, 'jnk--prep': 10, 'support--conj': 3, 'hsc': 2, 'epigenetic': 2, 'two-component': 2, 'hgcn': 2, 'proline-arginine': 2, 'd-tyr': 4, 'shifting': 2, 'filaments': 2, 'replace--rcmod-': 5, 'including--activity': 5, '-parataxis--note': 2, 'determination--amod-': 2, '-advcl--grow': 3, 'mek-erk': 2, 'ar-treated': 2, 'gtp-ase': 2, 'domain--conj': 4, 'terminal--amod-': 7, 'and--between': 10, 'ensure': 2, 'sept': 2, 'conformations': 2, 'decrease--xcomp-': 2, 'cc-rcc': 8, 'down-regulation--dobj-': 4, 'soh': 10, 'sparing': 2, 'nonpenetrance': 4, 'uncontrolled': 2, 'membrane-bound': 2, 'tfiih': 2, 'versions': 2, 'ppx': 4, 'rounds': 6, 'rely': 3, '-dep--three': 2, 'and--support': 5, 'connote': 4, 'cul-': 2, 'waf': 61, 'diferuloylmethane': 2, 'substantiate': 2, 'map-kinases': 2, 'accordingly': 6, 'triacyl': 2, 'for--sik': 2, 'mt-': 2, '-dep--c-fo': 7, 'rela': 45, 'hmc-': 2, 'atf-': 22, 'term--csubj-': 3, 'attempted': 2, 'becoming': 4, 'coprecipitated': 4, 'with--view': 2, 'noncycling': 2, 'only--dep-': 2, 'tam-': 2, 'of--gtpase': 2, '-regulated': 12, 'sm-d': 4, 'cancers': 2, 'non-lysine': 4, 'terminally': 2, 'rtt-like': 2, 'renaturation': 6, 'adopt': 2, '-nsubj--stimulus': 7, 'trid': 2, '-dobj--ad': 2, 'analyse--ccomp-': 4, 'epiregulin': 2, 'trail-sensitive': 2, '-ccomp--c': 2, '-nsubj--depletion': 3, 'co-association--nsubj-': 2, 'nh': 175, 'ring-containing': 2, 'aipl': 2, 'nn': 2, 'proenzyme': 2, 'nb': 2, 'camkiialpha': 14, 'monoubiquitinylate': 2, '-amod--mitogen-activated': 3, 'irradiation--prep': 2, 'consist--xcomp-': 6, 'brings': 4, 'tid': 4, '-agent--show': 4, 'by--analyze': 5, 'egf-dependent': 4, 'innervates': 2, 'ccas': 2, 'convincing': 2, '-nsubjpass--mapk': 4, '-dobj--bax': 2, 'sumo-': 2, '-dep--on': 3, 'action--agent-': 3, 'ikappabs': 4, 're-localization': 2, 'tcl-tropic': 2, 'mitogen-activated': 72, 'assume': 2, 'tar-cat': 4, 'modules': 2, 'of--subset': 2, 'fss': 4, 'hdelag': 2, 'excitotoxic': 8, 'fsx': 2, 'isolated-enzyme': 4, 'immunocomplex--prep': 2, 'skip': 2, 'relieve': 4, 'rog': 2, 'including--nf-kappab': 4, 'tripartite': 4, 'of--lung': 6, 'paf-stimulated': 2, 'functioned': 6, '-ccomp--means': 2, 'clarification': 8, 'cell-surface': 2, 'novo--prep': 2, 'phalloidin': 2, 'modulation--pobj-': 2, 'parenchyma': 2, '-phosphoinositide-dependent': 2, 'render--dep-': 2, 'sequestered': 6, 'aspp': 8, 'stabilizing': 4, 'analogous': 6, 'with--potential': 2, 'pcna-dependent': 2, 'ubiquitive': 2, 'benzodiazepines': 4, 'bcgf-ii': 2, 'oligodeoxynucleotides': 4, 'and--enhancement': 3, 'erk--conj': 2, '-dep--role': 4, 'including--cc-rcc': 4, 'as--mechanism': 4, 'responders': 4, 'considering-': 3, 'determine--amod-': 2, 'scf-type': 4, 'synovial': 2, 'fgfs': 2, 'abrogate--rcmod-': 2, 'skewed': 2, 'opposite--parataxis-': 2, 'comparable--advcl-': 3, 'fftd': 2, 'glb': 2, 'bhlh': 4, '-nsubjpass--alteration': 3, 'ceh--prep': 2, 'mn': 2, 'hvem': 2, 'scchn': 4, 'analogous--amod-': 2, 'trailr': 4, 'hindered': 2, 'pacap-initiated': 2, 'heterodimerizes': 2, '-d-pgj': 2, 'eating': 2, 'tuberous': 2, 'delcct': 4, '-nsubj--evaluation': 2, '-dobj--c-flip': 2, 'domain-': 2, 'scaffold': 22, 'necrotic': 4, 'f-mediated': 2, '-nsubj--withdrawal': 8, 'bax--dobj-': 2, 'transactivation-mediated': 2, '-partmod--detect': 7, 'efficacious': 2, 'to--fibroblast': 4, 'ligand-bound': 2, 'unleash--partmod-': 2, 'ht-': 4, 'cytokine-induced': 2, '-dep--hdac': 2, 'proasthmatic': 2, 'understanding': 5, '-mediated--acomp-': 2, '-agent--kinase': 3, 'activator--nsubj-': 2, 'overlap--dep-': 2, 'accomplishes': 2, 'support--prep': 6, 'heteromeric': 4, 'muscle-like': 2, 'ph-induced': 2, 'cross-react': 2, 'accomplished': 4, 'renatured': 8, 'caspase-resistant': 2, 'hbz': 2, 'anoxia': 4, 'xpf-binding': 2, 'alpha-granules': 2, 'heterodimerized': 2, 'tasks': 2, 'prommp': 4, 'subunit--nsubj-': 9, 'telomerase': 6, 'actin-microspike': 2, 'anisomycin': 4, 'src-': 13, 'msci': 4, 'of--mkk': 3, 'd-box': 2, 'transplant': 6, 'oxide-cgmp-induced': 2, 'terminal': 98, 'characterization--nsubj-': 2, '-dep--overlap': 2, 'cooperates': 6, 'estrogen-treated': 4, 'and--tlr': 4, 'camkii': 10, 'msl': 4, '-appos--gln': 2, 'upregulate--partmod-': 3, 'before--transmit': 2, 'cooperated': 4, 'ubiquitinate': 3, 'lung--prep': 5, 'for--mechanism': 3, 'heterodimers': 30, 'macrophage-tropic': 6, 'rfx-type': 4, 'absence--nsubj-': 2, 'hbds': 2, 'copy--prep': 2, 'hyper-proliferation': 4, 'cytoplasmically': 2, 'besides--maintain': 2, 'whatever': 2, 'phospholipases': 2, 'complex--appos-': 2, 'translocator--prep': 5, 'map-kinase': 2, 'ai-mediated': 2, 'rcc-associated': 2, 'two--npadvmod-': 2, 'thr-': 10, 'ubiquitin-b': 2, 'facilitating': 2, 'fas-r': 4, '-partmod--unleash': 2, 'cak': 10, 'ulk': 2, 'fold-induction': 8, 'ligand-activated': 4, 'modifier': 2, 'of--role': 3, 'flg': 2, 'social': 8, 'and--maintain': 2, 'disruption--nsubj-': 4, 'isoproterenol-stimulated': 2, 'jun-family': 2, 'to--uv': 2, '-dep--form': 2, 'phosphorylate': 69, 'ectopically': 2, 'spc': 2, 'snap-': 2, 'c-fo--dep-': 7, 'sblp': 4, 'neddylation': 2, 'd-dependent': 6, 'rarbeta': 2, 'septin': 4, 'snp--appos-': 3, 'hbrm': 2, 'short-patch': 2, '-rcmod--refer': 2, 'oropharyngeal': 4, '-nsubj--activator': 2, '-nsubj--model': 2, '-nsubj--two': 5, 'ebna': 8, '-ccomp--mean': 2, 'inputs': 2, 'mar': 2, 'and--suppose': 2, 'through--formation': 4, 'and--cancer': 2, 'cardiak': 2, 'applications': 2, 'syndrome-associated': 2, 'refolded': 2, 'methylate': 3, 'improving': 2, 'srp': 11, 'antagonise': 2, 'dynamics--dobj-': 2, 'maybe': 4, 'ataxia-telangiectasia-mutated': 2, 'sv': 20, '-ucla': 2, 'phosphotyrosine-independent': 2, 'antiparallel': 2, '-dep--transcription': 4, 'co-expressions': 2, 'pkd-s': 18, 'erk--agent-': 4, 'tagged': 4, 'bbr': 4, 'r--amod-': 2, 'regulatory-factors': 4, 'chromaffin': 4, 'holocytochrome': 2, 'dissociation--prep': 3, 'c-flice': 4, 'and--extend': 2, '-oh': 2, 'sirna-treated': 2, 'and--allow': 2, 'beta-catenin-dependent': 2, 'tta': 2, '-dep--thr-': 3, 'coreceptor': 39, 'implicating': 14, 'conjunction': 3, 'gtf': 2, '-trk': 7, 'immuno-precipitation': 4, 'for--ccr': 2, 'jip': 8, 'lrp--dep-': 2, 'bezafibrate': 2, '-nsubjpass--action': 2, 'creb-binding': 6, 'paralog--prep': 2, '-xcomp--phosphorylate': 4, 'lpa-mediated': 4, 'monensin': 4, 'qualified': 2, 'immunoblot--nsubj-': 3, 'proteinaceous': 4, '-dep--inhibition': 5, '-rcmod--regulator': 2, 'tanycytes': 4, 'sub-pool': 2, '-nsubj--correlation': 3, 'between--g': 2, 'play--advcl-': 3, 'of--convertase': 2, 'revert': 2, 'apoptogenic': 2, 'egf-il': 2, 'synchronize': 2, 'between--p': 2, 'physiologic': 4, 'bzlf': 2, 'for--absence': 2, '-dobj--rcc': 6, 'coinduce': 2, '-dobj--dsptp': 2, 'hypoxic-ischemic': 2, 'in--tumour': 2, 'alpha-amidating': 2, 'mhv-infected': 2, 'entering': 2, 'see--ccomp-': 8, 'in--mutant': 2, 'ma-': 2, '-dobj--down-regulation': 4, '-transformed': 8, 'ep': 2, 'beta-defensins': 2, 'hcr': 2, 'bly': 2, 'hct': 12, 'bubr': 2, 'rounding': 2, 'sensitive--dep-': 6, 'dihydropyrimidine': 2, 'kcnq': 5, 'hcc': 13, 'inactivate--partmod-': 5, 'osmr': 6, 'kcne': 9, 'hcn': 27, '-nsubj--mutlalpha': 6, 'mitochondria-dependent': 4, 'cooperatively': 2, 'pathogenic': 2, 'paralog--appos-': 2, '-driven': 4, 'intensely': 4, 'independent--rcmod-': 6, 'ore': 16, 'remainder': 2, 'baf': 4, 'gain-of-function': 2, 'pvwf-dependent': 4, 'spatial--dobj-': 3, 'trans-interact': 2, 'turn': 10, 'but--independent': 4, 'nodal': 2, 'c-jun--nsubj-': 7, 'signature-like': 2, 'n-sh': 2, 'antisenses': 2, 'nepp': 4, 'surviving': 2, '-nsubj--mkk': 4, 'xy': 2, 'n-smase': 6, 'skrp': 15, 'array': 6, 'six-hour': 2, 'antiserum--nsubj-': 2, 'btp': 2, 'necessarily': 2, 'ionomycin': 6, 'participates': 8, 'gid': 2, 'ngf-stimulated': 2, 'distinct--amod-': 3, 'jnk-': 2, 'ikki': 4, '-ccomp--exhibit': 4, 'betaark-ct': 2, 'fos-jun': 2, 'mtbp': 2, 'hub': 2, 'as-jnk': 2, 'recruitment--dep-': 2, 'accurately': 14, 'of--cytotoxicity': 2, 'crevice': 2, 'ush': 2, 'titrate--rcmod-': 2, 'relb': 18, 'ebv-infected': 2, 'cre-binding': 4, '-nsubj--manner': 2, 'infiltrate': 2, 'lucent': 2, '-partmod--decrease': 2, '-ccomp--regulate': 3, 'jnks': 25, 'extend--ccomp-': 2, 'of--mapk': 8, 'sknsh': 4, 'with--sp': 3, 'oligonucleotides': 4, 'co-transfection': 5, 'tfiie': 4, 'acetyllysine': 2, 'dual-specific': 2, 'hbeta': 2, 'upstream--dobj-': 2, 'jnkk': 10, 'on--subset': 2, 'er-beta': 14, 'accelerative': 2, 'sam': 2, 'prb': 20, 'dynamically': 2, 'turned': 2, 'beta-tc': 2, 'through--stabilization': 2, 'uv-induced': 4, 'factor-mediated': 8, 'distorted': 4, '-pobj--antagonist': 4, '-nsubj--mapk': 2, 'l-jnki': 2, 'with--expression': 8, 'lead--rcmod-': 8, 'c-jun--nsubjpass-': 14, 'lps-lbp': 2, 'erse-ii': 2, '-dobj--translocator': 2, 'efficient': 15, 'glucocorticoid-receptor-interacting': 4, 'overcoming': 2, 'converges': 4, 'inhibition--dep-': 5, '-advcl--abolish': 3, 'rcc--dobj-': 6, 'hyperacetylation': 2, '-nsubj--met': 3, 'jund': 4, 'ultimate': 4, 'fcr': 2, 'trichomes': 4, '-ccomp--polarize': 3, 'paralogs': 46, 'probed': 4, 'cotransfection--nsubj-': 2, 'yotiao': 2, 'autophosphorylation': 2, 'link--ccomp-': 2, '-agent--number': 2, '-xcomp--consist': 6, 'c-jun--prep': 76, 'refeeding': 2, 'stimulus--nsubj-': 7, 'deltg': 4, 'later': 2, 'adrenals': 2, 'zo-': 7, 'in--snp': 4, 'in--trafficking': 3, 'regulating-molecules': 2, 'jurkat': 6, 'tnf-rp': 2, 'mrg': 6, 'c-flip--dobj-': 3, 'precise': 4, 'drastically': 8, 'ccrl': 2, 'fas-disc': 4, 'implicate--partmod-': 2, 'nle': 2, 'compose--partmod-': 2, '-dobj--status': 2, 'feedback': 8, 'kinesin-i': 2, '-nsubj--uptake': 2, '-dobj--death': 2, 'homologs': 4, 'mutation-harboring': 2, 'liberation': 2, '-dep--facilitate': 2, 'lipids': 2, 'ddldpy': 2, 'activator--rcmod-': 2, 'slam-associated': 2, 'although--': 3, 'mutlalpha': 12, 'epo-r': 2, 'inhibit--dobj-': 4, 'assistant': 4, 'in--pathway': 5, '-dobj--kda': 2, '-pobj--activator': 2, '-dep--jnk': 4, 'in--present': 2, 'truncate--partmod-': 2, 'mechanism--advcl-': 3, 'over-expression--conj': 3, 'bifurcation': 4, 'multiprotein': 4, 'mapk--nsubjpass-': 4, 'summary': 2, 'cell-matrix': 2, 'capable--xcomp-': 6, 'signal-induced': 2, 'mapk--dobj-': 2, 'discussing': 2, 'giantin': 2, '-dep--block': 2, 'irradiation': 17, 'bisphosphatase': 2, 'concert': 2, 'rankl--prep': 2, 'physically': 8, 'of--c-jun': 7, 'tpbeta-mediated': 2, '-dep--immunize': 2, '-sensitive': 2, 'co-exists': 2, 'arabidopsis--agent-': 2, 'rmc': 2, 'intratracheally': 2, 'basis--dobj-': 3, 'manganese-induced': 4, 'combined-antisense': 2, 'hrgr': 6, '-agent--reduce': 2, 'surfaces': 2, 'beta-sheet': 2, 'tnf-r': 18, 'n-terminus': 14, 'or--stabilize': 2, 'infrequent--ccomp-': 2, 'alpha-naphthoflavone': 2, 'acylation': 4, 'cytoprotective': 2, 'nuclear-speckles': 2, 'import': 7, '-ccomp--account': 2, 'oeas': 6, 'therapies': 16, 'hypophosphorylated': 4, 'since--call': 10, 'tandemly': 2, 'jp': 6, 'dome': 2, 'cys-free': 6, 'killing': 2, 'phosphoglycerate': 2, 'and--demonstrate': 9, 'stages': 2, 'lp--agent-': 2, 'soleus': 2, 'propagation': 4, 'to--d': 3, '-appos--c-jun': 7, 'inward': 2, 'in--stabilize': 2, 'siii': 2, 'means--ccomp-': 2, 'melanocortin-': 2, 'sccs': 2, 'destabilizes': 2, 'suppressor-like': 2, 'n-terminal-iapp': 2, 'facilitate--prepc': 2, 'than--rcmod-': 2, 'uv-driven': 6, 'analogue--prep': 2, 'pkc--prep': 4, '-dobj--jip-': 2, 'a-related': 8, 'gynecomastia': 4, '-nsubj--therapy': 2, 'immunodeficiency': 2, '-ccomp--infrequent': 2, 'caspase-activated': 2, 'igm-april': 6, 'deposited': 4, '-advmod--not': 4, 'at--ser': 4, 'arrest--xcomp-': 3, '-ccomp--define': 3, 'sblp--tmod-': 2, 'transactivations': 10, 'with--variety': 5, '-parataxis--classify': 2, 'action--nsubjpass-': 2, 'egfr-devoid': 4, 'gist': 4, 'disrupt--conj': 3, 'receptor-ras-erk--prep': 2, 'kda--dobj-': 2, 'facilitator': 2, '-agent--p': 2, 'through--modulation': 3, 'dicer--nsubj-': 2, 'subconfluent': 2, 'folding': 12, 'c-fes': 2, 'proconvertases': 2, 'cell-death': 2, 'adrenoleukodystrophy': 2, 'understand--partmod-': 4, 'since--p': 6, 'eif-': 2, 'duplexes': 2, 'atm--dobj-': 2, 'hdm': 3, 'level--agent-': 2, '-nsubjpass--protein-': 2, 'basis--nsubj-': 2, 'pneumococcal': 2, 'pelle-like': 2, 'egf-induced': 2, 'topoisomerase': 2, 'upperstream': 4, 'angiotensinogen': 4, '-ccomp--phosphorylate': 2, 'ontogeny': 2, 'spleens': 2, 'illustrates': 4, 'complexation--nsubj-': 4, 'inducibility--dobj-': 2, 'egf--dobj-': 2, 'erlotinib': 2, '-ccomp--colocalize': 2, 'phosphorylating': 22, 'dsptp--dobj-': 2, 'beta-adaptin': 6, 'distributions': 4, 'in--immunocomplex': 2, 'effectively--advmod-': 2, 'ifn-gamma-mediated': 2, 'envelope-mediated': 2, '-slm': 2, 'subunit': 54, 'nadph-producing': 2, 'nab': 4, 'adaptin': 4, 'kinase-inhibitors': 2, 'assist': 6, '-dobj--mapk': 2, 'chain-restricted': 2, 'alpha--dobj-': 2, 'without--activate': 2, 'tumor-necrosis': 2, 'ceh': 4, 'crhr': 2, 'not--conj': 2, 'tracking': 2, 'roughly': 2, 'restrain': 2, 'pathway--dep-': 2, 'trimerization': 3, 'downmodulation': 2, 'ccg': 4, '-nsubj--fd': 2, '-dobj--calcitonin': 2, 'consisted': 10, 'redirects': 2, 'crk-associated': 2, 'anti-rat': 2, 'cocapping': 9, 'transmitting': 2, 'spliced': 14, 'sp-induced': 4, 'virus-specific': 2, 'tif-ia': 6, 'from--macrophage': 2, 'heterophilic': 4, 'downregulate--ccomp-': 2, 'mimicking': 4, 'blos': 2, 'of--plasmid': 2, 'heparin-binding': 10, 'atra': 2, '-nsubj--antibody': 2, 'retinoid-mediated': 2, 'to--release': 4, 'n-wasp-dependent': 2, 'decomposed': 2, 'independent--ccomp-': 3, 'amd': 4, 'c-raf-': 4, 'two--conj': 2, 'gemin': 2, 'extracellular-related': 2, 'mof': 4, 'clan': 2, 'ki-ra--conj': 2, 'detect--partmod-': 7, 'hyper-proliferation--dobj-': 2, 'let-': 2, 'hccsmc--prep': 3, '-parataxis--similar': 2, 'cfos': 2, 'flow-induced': 2, 'calcium-activated': 2, 'pdgf-rbeta': 2, 'activity--nsubj-': 2, 'a-mediated': 2, 'gnrh': 2, '-dep--accomplish': 3, 'lps-dependent': 4, 'cyclase--prep': 2, '-nsubj--independent': 3, 'e-bp': 2, 'dsred': 2, 't-snare': 4, 'x-box-binding': 2, 'gvl': 2, 'cavity': 2, 'nherf-': 2, '-tmod--cell': 3, 'ad-dn-c-jun': 4, 'in--distinction': 2, 'ddb': 10, 'using--': 5, 'vcp-interacting': 2, 'synergize': 4, 'duration--dobj-': 5, 'r-fc': 2, 'mutant--pobj-': 7, 'preincubate': 2, 'derlin-': 3, 'continental': 4, 'coexpressing': 2, 'tdu': 2, 'tdt': 2, 'srb': 8, 'om': 2, 'hampers': 2, 'oi': 2, '-xcomp--antagonise': 2, 'action--prep': 3, 'acm': 2, 'c-loop': 2, 'au-rich': 2, 'discrimination': 2, 'mixed': 2, 'hyperactivation': 4, 'sepcr': 2, '-rcmod--potentiate': 2, 'don-exposed': 4, 'filamin': 2, 'inclusion': 2, 'communication': 2, '-nsubj--pretreatment': 5, 'mutant--prep': 6, 'enhancement--prep': 6, 'accounts': 4, 'anti-phospho-mitogen-activated': 2, 'hla-wtg': 2, 'dna-damage-inducible': 2, 'cyclins': 4, 'such--dep-': 2, 'rb-cdk': 8, 'colocalizes': 2, '-dobj--method': 2, 'utilize--rcmod-': 3, '-appos--process': 2, 'have--prepc': 2, 'dnase': 4, 'independent--conj': 2, 'c-rel': 88, 'transition--prep': 2, 'recognize--conj': 2, 'to--carcinogenesis': 2, 'osteochondroma': 4, 'in--ap-': 2, 'ptpn': 2, 'compensate': 2, 'granules': 4, 'telomere-driven': 4, 'atrophy': 6, 'baffr': 4, '-nsubjpass--region': 9, 'phospho-pecam-': 2, 'without--cc-rcc': 2, 'survived': 4, '-nsubjpass--phosphatase': 2, 'bring--xcomp-': 5, 'transgenes': 2, 'mapk--prep': 9, 'applying': 2, 'over-expression': 4, 'and--interfere': 4, 'reasoned': 2, 'vhr': 6, 'fsh': 2, 'hstaf': 3, 'in--degradation': 4, 'peaking': 2, 'role--dep-': 4, 'phosphatase--nsubjpass-': 2, '-nsubj--subunit': 9, 'betatrcp': 4, 'augmentative': 2, 'uva-mediated': 4, '-conj--due': 2, 'immunoprecipitation': 6, 'helicase': 2, '-dep--pathway': 2, 'conversely': 10, 'with--oligonucleotide': 5, 'again': 2, 'cancer-associated': 2, 'suppress--rcmod-': 5, 'as--family': 2, 'cp-': 2, '-parataxis--downregulate': 2, 'xpb--dobj-': 2, 'actin-disrupting': 2, '-advcl--tend': 2, 'homodimerization': 10, 'carbohydrate-rich': 2, 'spatial': 6, '-npadvmod--function': 2, 'telangiectasium': 4, 'translocator--dobj-': 2, 'alcl': 2, 'haptens': 4, 'adeoad': 4, '-nsubj--dependence': 2, 'mekk': 40, 'degradation--prep': 2, '-partmod--have': 4, 'sultibi': 2, 'by--improve': 2, 'profile--nsubj-': 3, 'stimulate--partmod-': 5, 'mask': 2, '-kinase-independent': 2, 'amino-terminus': 2, 'ensuring': 2, 'anisomycin-mediated': 2, 'trithorax': 2, 'resembles': 4, 'nutrient': 2, 'mitogen-associated': 8, 'hgf-stimulated': 6, '-advcl--high': 4, 'sci': 2, 'flice-like': 4, 'ptd': 2, 'papilloma': 2, 'differentiation-associated': 4, 'staurosporine-activated': 2, '-dobj--inducibility': 2, 'c-fo': 3, '-dobj--modulation': 5, 'ptt': 2, 'oncoprotein--nsubj-': 2, 'prolongs': 2, 'fa-associated': 4, '-ccomp--downregulate': 2, 'gelatinase': 4, 'destroy--conj': 2, 'neuroregeneration': 2, '-dep--capable': 2, 'gp': 16, 'permeability': 18, 'aogen': 4, 'on--activity': 5, 'nmbr': 2, 'deltant': 2, 'crest': 2, 'axonal': 10, 'elk-': 14, 'ftd': 2, 'dupat': 4, 'returning': 2, '-agent--erk': 2, 'subunits': 66, 'maintain--prepc': 2, 'of--antibody': 5, 'jnk--appos-': 2, 'heterodimerize': 5, 'participate--advcl-': 8, 'icad': 8, 'faf': 4, 'of--growth': 4, 'dualtropic': 2, 'but--not': 3, 'corpuscles': 6, 'ageing': 2, '-partmod--precede': 2, 'phkb': 2, 'agonistic': 6, 'gos': 2, 'wild-type--prep': 4, 'another--prep': 4, 'ubiquitination--dep-': 2, 'paralog--conj': 4, 'kinase-': 5, '-dep--sdf-': 3, 'ten': 2, '-nsubj--antiserum': 2, 'supershifted': 2, 'finger-containing': 2, 'tef': 2, 'nbk': 4, 'gammagcs': 2, 'and--implicate': 3, 'displays': 4, 'hormone-binding': 2, 'interfere--conj': 4, '-ccomp--compatible': 2, 'mdm-': 7, 'adapter--prep': 4, 'modulation--dobj-': 5, 'c-jun--tmod-': 2, 'tfeb': 4, 'c-c': 2, 'in--translocator': 5, 'suv': 8, 'to--monoubiquitination': 2, 'sur': 4, 'activations': 4, '-xcomp--initiate': 2, 'calcitonin--dobj-': 2, 'oligonucleotide--prep': 4, 'cohorts': 2, '-rcmod--than': 2, 'trail-r': 23, 'cell-recruiting': 2, 'rejector': 2, 'behaviour': 12, 'wild': 2, '-appos--homolog': 5, 'iras': 14, 'complexation': 17, 'insertion-deletion': 2, 'striatum': 2, 'release--parataxis-': 4, '-amod--determination': 2, 'multiplex': 2, 'pics': 2, 'ski-interacting': 4, 'prickle': 2, 'fsh-induced': 2, 'of--dr': 3, 'osmolarity': 4, 'remarkably--advmod-': 3, 'pathological': 8, 'actin-cytoskeleton': 2, 'serine-': 2, 'octapeptide': 2, 'mediator': 2, 'as--polypeptide': 2, 'lipopolysaccharide-binding': 2, '-ccomp--appear': 2, '-appos--sb': 2, 'analyzing': 8, 'pre-initiation': 4, 'dc-signr': 2, '-advcl--mechanism': 3, 'exostosin-': 2, 'bundle': 2, 'culture--dep-': 2, 'etc': 2, '-rcmod--facilitate': 2, 'of--protein-': 2, 'n-cor': 16, 'tondu': 4, 'and--inactivate': 5, 'with--tetraspanin': 2, '-lz': 2, '-dobj--fragment': 2, 'wxxxf': 2, '-dep--coactivate': 2, 'bax--prep': 2, 'consist--ccomp-': 3, 'baff-r': 4, 'infertile--dep-': 2, 'compose--ccomp-': 3, 'osteoprotegerin': 10, 'for--three': 2, 'unnecessary--ccomp-': 2, 'initiates': 6, 'transform--dep-': 5, 'varying': 2, 'a-raf': 2, 'fragilis': 2, 'serines': 2, '-nsubj--tumor': 3, 'lrp-': 2, 'raf-': 4, 'fas-binding': 2, 'as-odn--prep': 2, 'baff--prep': 2, 'beta-pleated': 2, '-nsubjpass--mt': 2, 'endocytosis': 6, 'syntaxin': 4, 'del': 6, '-dobj--conformation': 2, '-using': 4, 'specialized': 2, 'scramble': 4, 'implicate--conj': 3, 'delt--conj': 3, 'cad--pobj-': 2, 'circulate--csubj-': 2, 'll-z': 2, 'suppression--prep': 5, '-dobj--upstream': 2, '-dobj--subset': 5, 'minimum': 4, 'pten-null': 6, 'loop-deleted': 2, 'concordant': 6, 'scar': 4, 'f-': 4, '-rcmod--record': 2, 'similar--parataxis-': 2, 'needs': 2, 'of--modulation': 2, 'rather': 10, 'mapk': 289, 'actr': 4, 'to--nucleus': 2, '-ccomp--synergize': 3, 'by--binding': 2, 'actn': 6, '-smad': 2, 'fancf--prep': 10, 'preference--nsubjpass-': 2, 'ubiquitin': 8, 'ascertain--rcmod-': 2, 'denaturing': 2, 'fr': 2, 'b-chain-c-peptide': 2, '-appos--polymerase': 2, 'agonist-induced': 2, '-dobj--cooperation': 4, 'questioning': 2, 'function--npadvmod-': 2, '-year-old': 2, 'fd': 6, 'of--no': 2, 'anneal--rcmod-': 3, 'fo': 2, 'trois': 2, 'chiefly': 2, '-dep--tyr-': 3, 'mapkkk': 16, '-advcl--consist': 2, 'dispensable': 12, 'trafficking--prep': 3, 'with--therapy': 2, 'phosphorylate--csubj-': 2, '-rcmod--lie': 3, 'on--kinase': 3, 'ishikawa': 2, 'with--phosphorylation': 3, 'with--prominent': 2, 'on--phosphorylation': 4, 'paralogues': 4, 'cause--prepc': 2, 'ffh': 8, '-rcmod--replace': 5, 'uvc': 14, 'ires-dependent': 2, 'orexigenic': 2, '-xcomp--open': 2, 'non-responders': 2, 'but--diminish': 2, 'polymerases': 2, 'nonrecurrent': 4, 'actually': 2, 'jhco': 3, 'on--jnk': 2, 'pcl': 2, 'mefs': 2, 'uv-inducible': 2, 'osmosensitive': 2, '-xcomp--desensitize': 2, 'nephrin': 4, 'propose': 12, 'permeability--prep': 3, 'direct--ccomp-': 2, 'tgf-beta-induced': 4, 'mage-a': 2, 'due--conj': 3, 'refer--rcmod-': 2, 'bcl-xl-transfected': 2, 'with--paralog': 2, 'tyr-phosphorylation': 6, 'impair--ccomp-': 3, 'non-peptidic': 2, 'c-jun--dobj-': 20, 'racp': 2, 'monoubiquitinylated': 2, 'ikappabepsilon': 2, 'fance': 58, 'fancd': 22, 'arf-mdm': 2, 'phosphoinositide-dependent': 2, 'inherent': 2, 'ifh': 4, '-ccomp--dispensable': 2, 'investigators': 2, '-agent--accumulation': 2, 'czeta': 6, 'consistent--nsubj-': 2, 'loxp-mediated': 4, 'level-egfr': 8, 'at--degree': 7, 'pol--acomp-': 3, 'oncoproteins': 4, 'loops': 8, 'melanocytes': 2, 'with--plzf': 3, 'tbp-tafi': 2, 'sb-': 7, 'inactivates': 10, 'remodelling': 2, 'tfiid--nsubj-': 2, 'on--trial': 2, 'interfere--agent-': 2, 'inactivated': 18, 'tfiih---': 4, 'nf-y-independent': 4, 'idurd': 4, '-rcmod--comprise': 2, 'dysregulation': 3, '-subunit': 2, 'in--apc': 2, 'asv': 2, 'f-transfected': 2, 'antisera': 2, '-partmod--upregulate': 3, '-partmod--give': 3, '-appos--importance': 2, 'abnormally': 2, 'tetrameric': 2, 'asm': 2, 'pprb': 2, 'pairs': 4, 'retained': 10, 'mdl-': 4, 'and--two': 2, 'titrates': 4, 'tcdd': 4, 'not--advmod-': 4, 'cssna': 2, 'by--subunit': 2, '-dobj--versatility': 2, 'atc': 4, 'assemblage': 2, 'tend--advcl-': 2, 'co-operates': 2, 'as--kinase': 4, 'utilizes': 4, 'aso-mediated': 2, 'radioresistant': 2, 'jeg-': 2, 'fancg': 34, 'hyaluronan': 2, 'sequester': 18, 'immunoreactivities': 2, 'fancc': 58, 'fanca': 50, 'anti-phosphotyrosine': 2, 'define--ccomp-': 3, 'fancl': 6, 'splice-site': 2, 'vinorelbine': 2, 'homologues': 2, 'impairing': 2, '-nsubj--complexation': 4, 'antisense-treated': 2, 'gpd': 2, 'pms': 40, '-macroglobin': 2, 'kinesin-i-driven': 2, 'hspc': 2, 'dff': 2, '-nsubj--isoform': 2, '-nsubj--cyclin': 2, 'kinase--appos-': 4, '-rcmod--activator': 2, 't-all': 4, 'co-culture': 4, 'bowel': 2, 'with--chromogranin': 3, 'corticotroph': 2, 'concordance': 3, 'ets-domain': 4, 'focusing': 2, 'ebp-atf': 2, '-dep--waf': 4, 'immunoprecipitations': 4, 'polarized': 2, 'co-receptor': 9, 'gnrh-induced': 4, '-dep--culture': 2, 'aak': 4, 'brct-motif-containing': 4, 'pairwise': 2, 'phosphoserine': 4, 'lysate': 2, 'sb': 85, '-nsubj--immunoblot': 3, 'transcriptase-polymerase': 2, 'dc-sign--prep': 2, '-deltac': 2, 'intraocular': 2, 'allele--dep-': 2, 'wafi': 4, '-nsubjpass--virus': 2, 'teleost': 2, 'kinase-inactive': 4, 'by--compare': 3, 'antagonizes': 2, 'a-associated': 4, '-dobj--growth': 3, 'destabilizing': 2, 'but--affect': 3, '-dep--adaptor': 2, 'entirety': 6, 'src-mitogen-activated': 4, 'lip': 6, '-nsubj--cotransfection': 2, 'selp': 2, 'centrosomes': 2, 'mpl-induced': 4, 'hypophosphorylation': 5, '-ccomp--form': 8, 'i-vegf': 2, 'dna-pk': 7, 'trbp': 2, 'plap': 2, 'as--subunit': 2, 'trail--prep': 7, 'cep-': 3, 'phosphorylate--conj': 3, 'nfkb': 6, 'patterned': 4, 'chondrosarcomas': 4, 'estimate--prep': 2, 'tfiih--nsubj-': 4, 'kd': 4, 'of--array': 4, 'ext': 3, '-advcl--resist': 2, 'kv': 2, 'mineralocorticoid': 4, '-dep--such': 2, 'tc-ptp': 4, 'heterodimerization': 24, 'three--nsubj-': 4, 'in--stabilization': 3, 'hcap-c': 2, 'conjugation--nsubj-': 2, 'hcap-e': 2, 'hcap-g': 2, 'scattered': 2, 'pbsf': 2, 'rtt': 2, 'menage--nsubj-': 3, 'label--dep-': 2, 'denatured': 8, 'utilization': 2, 'carefully': 2, 'of--suppression': 3, 'pka-anchoring': 2, 'neuroendocrine': 2, 'versatility--dobj-': 2, 'tvl-': 6, 'detect--rcmod-': 2, 'bfa': 4, 'prb--appos-': 12, 'nfy-b': 2, 'translocates': 2, '-appos--paralog': 2, '-ccomp--compose': 3, 'destabilize': 2, 'experiment': 12, 'dysfunctional': 6, 'through--receptor': 5, 'jjhan': 2, 'importance--nsubj-': 6, '-dobj--pak': 2, 'begins': 4, 'inducibly': 2, '-tmod--sblp': 2, 'suppose--conj': 2, 'adcc': 2, 'hiv': 22, 'meanwhile': 4, 'component--prep': 2, 'attenuation--prep': 4, 'recruit--advcl-': 2, '-pobj--pd': 2, 'cumulate': 2, 'excitotoxin-induced': 2, '-gondii-mediated': 4, 'bifurcation--dobj-': 2, 'bona': 12, 'co-operate': 3, 'guinea-pigs': 2, 'rac-jnk': 2, 'death--dobj-': 2, 'noncanonical': 4, 'aryl': 12, 'withdrawal': 12, 'overexpression--agent-': 9, 'nonspecific': 6, 'hinder': 8, 'receptor-blocking': 2, 'of--wild-type': 3, 'and--event': 3, 'unneddylated': 4, 'of--induction': 10, 'homophilically': 2, 'undergoes': 2, 'heterotopias': 2, 'antiproliferative': 2, 'therapy--nsubj-': 2, 'n-aryl': 4, 'affinity--prep': 2, 'arf': 34, 'obtain--conj': 2, 'phosphogluconate': 2, 'ark': 2, 'positive--ccomp-': 3, 'retinopathy': 2, 'hassall': 6, 'and--dependent': 5, '-coupled': 2, 'wip-wasp': 2, 'wox': 2, 'distinctive': 4, 'currently': 2, 'quasi-native': 4, 'chlorpyrifos-induced': 2, 'between--conj': 8, 'woodchuck': 2, 'consist--advcl-': 2, 'beta--conj': 2, 'arf-bp': 2, 'from--death': 4, 'europe': 4, 'contain--csubj-': 2, '-dep--signal': 5, 'signal-regulated--amod-': 2, 'campestri': 6, 'complementary': 2, 'inhibitor--tmod-': 4, 'interfaces': 6, 'ulcerative': 2, 'succeed': 2, 'pcat': 4, '-dobj--basis': 3, 'jnk-binding': 8, 'heat-induced': 4, 'event--conj': 3, '-tmod--pathway': 2, 'chromogranin--prep': 3, 'context': 2, 'retarded': 4, 'emsa': 8, 'stress-activated': 38, 'experimental': 4, 'rag': 5, 'lauryl': 2, 'harbour': 4, 'with--deacetylase': 2, 'with--hsp': 2, 'vessel': 2, 'cumulated': 2, '-creb-dependent': 6, 'platinum': 2, 'rb-related': 2, 'delga': 6, 'ubiquination': 2, '-advcl--comparable': 3, 'eda--dep-': 2, '-pobj--mapkkk': 2, 'dcis': 4, 'lymphopenic': 2, '-appos--abnormality': 2, 'coimmunization': 2, 'smn--prep': 3, 'squirrel': 2, 'heregulin': 4, 'silenced': 2, 'from--stimulate': 2, '-advcl--link': 2, 'uva-induced': 10, 'resist--advcl-': 2, 'pacemaker': 4, 'and--pd': 4, 'fp-mediated': 2, 'shps-': 4, 'mlk': 2, 'in--affinity': 2, 'at--site': 2, 'marbp': 4, 'ras-mutant': 6, 'transgenic--prepc': 2, '-dep--mlh': 2, 'c-fo--prep': 2, 'utp-induced': 2, 'abnormality--appos-': 2, 'relum': 2, 'assay--rcmod-': 2, 'rin': 2, 'destabilization': 2, 'rim': 4, 'sgii': 4, 'terminator': 2, 'simian': 2, 'smcd': 2, 'chaffeensis': 4, '-nsubj--brca': 2, 'podocin': 2, 'jips': 4, 'hb-egf': 2, 'maximising': 2, 'trans-activation': 2, 'predominated': 6, 'edg-': 5, 'but--attenuate': 3, 'composition': 2, 'prerequisite': 2, 'gmmsd': 2, '-nsubj--relocalization': 2, 'attenuation': 6, 'forced': 20, 'restraint': 2, 'for--deletion': 2, 'restrains': 2, 'dca-induced': 2, 'targetting': 2, 'psd': 14, 'neuraminidase': 4, 'syt': 4, 'tata-box': 4, '-terminal': 74, '-parataxis--display': 3, '-dep--infertile': 2, 'zinc-superoxide': 2, 'interneurones': 4, 'pstat': 3, 'machinery': 6, 'jip-': 20, 'entry': 14, 'analyze--prepc': 5, 'doxorubicin-induced': 4, 'c-fos': 28, 'with--promoter': 2, 'precede--partmod-': 2, 'within--promoter': 2, 'nk-kappab': 2, 'unable--rcmod-': 2, '-rcmod--suppress': 5, 'irritable': 2, 'with--mutant': 3, 'colocalize--ccomp-': 2, 'rgr': 4, 'sensitive--advcl-': 2, 'melanoma--appos-': 2, 'transcriptase': 6, 'heterophilically': 2, 'tgf-alpha': 4, 'omega-conotoxin': 2, 'cue': 4, '-amod--terminal': 7, 'immunohistology-localized': 4, 'lack--parataxis-': 6, 'snap': 2, 'dysregulate': 2, 'transformation': 3, 'hne-mediated': 2, '-rcmod--ascertain': 2, 'and--disrupt': 2, 'big': 2, 'abolish--advcl-': 3, 'generate--advcl-': 2, 'emigration': 2, 'and--paralog': 4, 'in--gp': 2, 'k-c': 2, 'collagen--nsubj-': 3, '--lack': 2, '-rcmod--assay': 2, 'heterodimer--dep-': 2, '-dobj--array': 2, 'sequesters': 4, 'interfered': 6, 'isopeptide': 2, 'deleterious': 2, 'reveal--rcmod-': 5, 'with--c-jun': 4, '-nsubj--translocation': 5, 'cntfr': 2, '-dobj--alpha': 4, 'pretreating': 4, 'anti-ccr': 2, 'absolutely': 2, 'translational': 10, 'transactivators': 2, 'element-binding': 2, 'understood': 4, 'protrh': 2, 'with--as-odn': 2, '-dobj--brca': 2, 'since--': 2, 'gradually': 2, 'contacts': 4, '-ccomp--impair': 3, 'pretreatment--nsubj-': 5, 'pelp': 3, 'conditionally': 2, 'inoculate': 4, 'calmodulin-dependent': 4, 'in--transgenic': 2, 'ra-effect': 4, '-ccomp--see': 8, 'mitotic': 2, 'compatible--ccomp-': 2, '-csubj--derive': 2, 'histone-methylating': 2, '-nsubj--irradiation': 3, 'avp-stimulated': 2, 'through--regulate': 2, 'mof--prep': 2, 'cgc': 2, 'cgb': 4, 'perpendicular': 2, 'equipotent': 2, 'assists': 2, 'co-chaperone': 4, 'jnk--dobj-': 17, '-tmod--methyltransferase': 2, 'pifithrin-alpha--prep': 2, 'transactivator--nsubj-': 2, 'extravescicular': 2, 'potentiate--rcmod-': 2, 'ndel': 2, 'phosphatase-treated': 6, 'nfkappa-b': 6, '-partmod--compose': 2, 'dsb-related': 4, '-advcl--ensure': 2, 'cajal': 2, 'connection--prep': 2, 'phosphoserine-phosphorylated': 2, '-nsubj--conjugation': 2, '-nsubj--tfiih': 2, 'rpsk': 4, 'because--majority': 2, 'integral': 2, 'cdk-': 2, 'msr': 14, '-nsubj--elimination': 2, 'require--dep-': 3, 'dependence--nsubj-': 2, 'for--ubiquitination': 2, '-rcmod--site': 2, 'tall-': 10, 'arylhydrocarbon': 4, 'link--advcl-': 2, 'soce': 8, '-nsubj--absence': 2, 'ptcl': 4, 'mkk': 137, 'perp': 4, '-ccomp--key': 3, 'mkp': 10, 'cell--tmod-': 3, 'wildtype': 12, 'ra-induced': 2, 'c-mediated': 4, 'pkg-i': 2, '-rcmod--aid': 2, 'c-jun---': 5, 'chx': 2, 'unstimulated': 6, 'resistant-nature': 2, '-dobj--interaction': 8, 'pirh': 2, 'msc': 6, 'transformation-suppressor': 2, 'cin': 2, 'with--chromosome': 2, 'chemosensitivity': 2, 'uo': 6, 'to--estimate': 2, 'yes-associated': 2, 'highlight': 19, '-rcmod--compose': 2, '-dep--heterodimer': 2, 'region--nsubjpass-': 9, 'spectrum--nsubj-': 3, 'l-ohp': 4, 'scavenging': 4, 'hepatomas': 6, '-nsubj--replication': 2, 'monoubiquitination': 12, 'uptake--nsubj-': 2, '-dep--allele': 2, 'naes': 8, 'prevail': 2, 'aid--rcmod-': 2, 'acid-induced': 2, 'judge--dep-': 4, 'to--inhibit': 2, 'grx': 2, 'endosomes': 8, 'gre': 2, 'homodimerization--prep': 2, 'extinction': 2, 'hat-deficient': 2, 'form--rcmod-': 2, 'blys': 10, 'nonactivated': 12, 'in--adhesion': 3, 'envelope': 13, 'treatment--conj': 2, 'to--heterodimer': 2, 'of--phosphorylation': 6, 'suppression': 21, 'fishes': 2, '-xcomp--compensate': 2, 'drai': 6, '-npadvmod--two': 2, 'stabilization--prep': 5, '-nsubj--spectrum': 3, '-rcmod--subfamily': 2, '-hiv-': 2, 'and--sp': 2, '-partmod--hypophosphorylate': 2, 'enrolled': 2, 'fresh': 2, 'sftd': 2, 'having': 2, 'eda': 4, 'edg': 2, 'because--': 4, 'form--ccomp-': 8, 'gbetagamma-induced': 2, 'signatures': 4, 'broader': 4, 'tumours': 2, 'discordant': 4, '-nsubj--menage': 3, 'gk': 2, 'model--nsubj-': 2, 'migg': 4, 'in--cause': 2, 'hiv-': 12, 'supershift': 4, 'acetylates': 2, 'heterotrimer': 2, 'tert-butyl-hydroperoxide': 2, 'process--appos-': 2, 'concomitantly': 18, 'of--pkc': 4, 'as--activity': 2, 'npm-alk-mediated': 4, 'clathrin-dependent': 2, 'confirmation': 6, 'taxol-induced': 2, 'itc-induced': 4, 'tel-aml': 2, 'including--receptor-ras-erk': 2, 'continues': 2, '-rcmod--form': 2, 'engineered': 2, 'wave': 4, '-specific': 2, 'tet-on': 2, '-dobj--erk': 7, 'continued': 2, 'signal-regulating': 4, 'in--attenuation': 5, 'connect--advcl-': 2, '-ic': 4, '-advcl--assemble': 3, 'to--screen': 4, 'extend--conj': 2, 'phe--dep-': 3, 'dimer': 2, '-dobj--homologue': 2, 'pma-treated': 2, 'rrm': 2, 'oddball': 2, 'odns': 4, 'ddavp-induced': 2, 'in--loss': 2, 'rdna': 10, 'regulation--appos-': 4, 'abin-': 10, 'for--peptide': 4, 'noncarriers': 2, 'hungary': 2, 'polyubiquitination': 4, 'supressing': 2, 'sillence': 2, 'transplanted': 4, 'sapk': 45, 'calcineurin-regulated': 2, 'subtype-': 2, '-csubj--term': 3, 'access': 12, 'exogenously-expressed': 2, 'bdp': 2, 's-phase': 2, 'ct-': 6, 'consequently--conj': 2, 'desensitize--xcomp-': 2, '-partmod--understand': 4, 'stabilizer': 2, '-regulated--amod-': 2, 'appear--ccomp-': 2, '-dobj--spatial': 3, '-dep--ubiquitination': 2, 'mapk--dep-': 2, 'ha-ras': 4, 'nacht': 2, 'sik--prep': 2, 'intermediate-intensity': 2, '-dobj--dynamics': 2, 'engage': 5, 'dna-repair': 2, '-ccomp--clear': 2, 'monosomy': 6, 'telangiectasia': 2, 'caf': 10, 'rpoa': 4, 'vbp': 6, 'co-immunoprecipitation': 6, 'calcitonin-receptor-like': 2, 'prominent--prepc': 2, '-appos--complex': 2, 'chang': 2, 'makes': 8, 'delays': 2, 'composed': 14, '-aza-': 2, 'stabilization': 36, 'ago': 4, 'win': 2, 'of--analogue': 2, 'clara': 2, 'sgt': 4, 'proteasome-dependent': 2, 'lt-e': 2, 'pxc': 2, '-nsubjpass--preference': 2, 'apf': 8, 'bard': 9, 'rpsd': 4, '-prep--on': 3, 'telangiectasia-rad': 2, 'chm': 4, 'next': 14, 'zif': 6, 'duplicate': 2, 'and--v': 5, 'of--adenovirus': 3, 'invariably': 2, '-ccomp--direct': 2, 'msi-h': 4, 'and--essential': 2, 'nimodipine': 2, 'neuroendocrine-specific': 2, '-dobj--localisation': 2, 'damage-induced': 8, 'lactogen': 2, 'cooperate': 4, 'trail': 33, 'ra-inducible': 2, 'b-cll': 4, 'vip-': 2, 'detects': 2, 'b-infected': 2, 'unwinding': 2, 'fos-related': 4, 'kinase-alpha': 2, 'origins': 6, 'augment': 2, 'signal-related': 10, 'sre': 4, 'akt--conj': 2, 'other--prep': 2, 'stabilizes': 4, 'ubiquitin-mediated': 2, 'regulation--dep-': 2, '-amod--necessary': 6, 'rpsk--dep-': 2, 'vgl-': 2, 'enlargement': 2, 'weaker': 16, 'overexpression--appos-': 2, 'dispensable--ccomp-': 2, 'gamma-tubulin': 2, 'ra-mediated': 4, 'fra-': 18, 'taf': 6, 'adce': 2, '-generated': 2, 'tab': 10, 'related--conj': 2, 'co-purifies': 2, 'provoke': 2, 'docking': 22, 'tgfbeta-mediated': 2, 'defective': 17, 'delaa': 2, 'distinction--prep': 2, 'and--due': 2, 'delag': 14, 'epistasis': 4, 'gcp': 2, 'varies': 4, 'animal': 4, 'sik': 4, 'filament': 2, 'sweat': 4, 'pro-caspase-': 4, 'co-translational': 2, '-advcl--sensitive': 2, 'tp-mediated': 2, 'masks': 4, 'taci': 2, 'd-site--prep': 2, 'inst': 6, 'met-transfected': 2, 'co-transfect': 3, 'fd--nsubj-': 2, 'condensin': 4, 'ebp-alpha': 2, 'p-erk': 2, '-amod--distinct': 3, 'tethering': 6, 'tfiiagamma': 2, 'cbp--amod-': 2, 'hla-mg': 4, 'redundant': 3, 'beta-adrenoreceptor': 2, 'and--obtain': 2, 'substrate--conj': 3, '-partmod--stimulate': 5, 'saturated': 2, 'hydropathy': 4, 'hypophosphorylate--partmod-': 2, 'same--prep': 4, 'crosstalk': 4, 'kappabalpha': 8, 'mutl': 2, 'gcl': 2, 'sperm': 2, 'muts': 6, 'line-adapted': 2, 'fas-death-inducing': 4, 'jnk': 522, 'facilitates': 16, 'falsely': 2, 'meis': 7, 'mafg': 8, 'equimolar': 4, 'considering--fact': 2, 'integrators': 2, 'mapkk': 26, 'temporal': 2, 'hypothermic': 2, 'ikk-gamma': 2, 'prohb-egf': 2, 'degrees': 14, 'reduce--agent-': 2, 'evaluation--nsubj-': 2, 'lz': 2, 'induction--conj': 2, 'microfilament-disrupting': 4, 'using-': 4, 'perikarya': 2, 'repletion': 2, 'ensure--advcl-': 2, '-rcmod--independent': 6, 'hemochromatosis': 2, '-hsp': 14, '-parataxis--reduce': 3, 'lzk': 2, 'mitigated': 2, 'accumulates': 2, 'satellites': 2, 'grpr': 2, 'ea-induced': 2, 'd-kityf': 2, 'compensate--xcomp-': 2, 'caspase-independent': 4, 'propose--rcmod-': 5, 'mere': 4, 'form--dobj-': 6, 'heteromer': 2, 'concave': 4, '-dep--phe': 3, 'c-jun--appos-': 19, 'clinicopathologic': 2, 'nft': 4, 'topbp': 2, 'nor--modify': 2, 'identifies': 8, 'in--event': 2, 'presenilin--prep': 10, '-ccomp--analyse': 4, 'insa': 2, 'fusion--nsubjpass-': 2, '-mkk': 4, 'cck-induced': 2, 'to--cyclase': 2, 'salicylates': 4, 'pir-b': 2, 'tnf-r-associated': 4, 'cryopyrin': 2, 'and--stimulate': 9, 'investigate--partmod-': 4, 'rplq': 4, 'follicle': 4, 'gene-related': 2, 'rel': 26, 'of--fancf': 8, '-dobj--rad': 3, 'of--fancc': 6, 'dose-dependently': 4, 'ecd': 2, '-pobj--ag': 2, 'chaperone-like': 4, 'n-cad': 2, 'regulate--ccomp-': 3, 'have--partmod-': 4, 'rhok': 2, 'approaches': 2, 'pkc-zeta': 2, 'segregated': 2, 'five--nsubj-': 3, 'cyclin-dependant': 2, 'with--member': 2, 'timer': 2, 'ltc': 2, 'capable--dep-': 2, 'and--domain': 4, '-partmod--implicate': 2, 'predominate': 4, 'hivenv': 6, 'sek': 23, '-dep--mutant': 4, 'to--c-jun': 4, '-dep--label': 2, 'establish--advcl-': 3, 'stimulus-induced': 4, 'rectifying': 2, '-dep--phosphatase': 3, '-ccomp--itself': 3, 'myelogenous': 6, 'sterically': 4, 'hef': 2, 'fibroblast--nsubj-': 2, 'rfxank': 2, 'vitronectin-exposed': 2, 'immunize--dep-': 2, 'laddering': 2, '-dep--compensate': 2, 'neurotransmitter': 14, 'fosb': 2, 'adopts': 2, 'to--cul': 3, 'cyt': 10, '-rcmod--map': 2, 'heterodimeric': 2, 'beta-arrestins': 2, 'manipulation': 3, 'ogd': 2, '-nsubj--nine': 2, 'erk-': 2, 'endothelial-dependent': 20, 'sphingosine-': 2, 'display--parataxis-': 3, 'n-ras': 2, '-csubj--circulate': 2, 'bed': 2, 'pander-treated': 2, 'assumes': 4, 'alpha-helical-crf': 2, 'sr': 8, 'ovalbumin': 2, 'tfiid': 4, 'n-terminal': 540, 'hpr-induced': 2, 'prohormone-convertase': 2, 'mapk-phosphatase-': 2, 'bind--conj': 5, 'trpc': 94, 'mediators': 6, 'poisons': 2, 'iron-induced': 2, 'k--nsubjpass-': 5, '-rcmod--act': 17, 'methyl-cpg': 2, 'status--dobj-': 2, 'unlikely': 2, 'brg-': 3, 'need': 5, '-mediated': 4, 'to--n-terminus': 2, 'relb-containing': 2, 'viewed': 4, 'ubiquitination--prep': 2, 'minimally': 4, '-nsubj--basis': 2, 'and--treatment': 2, 'exhibit--ccomp-': 4, 'akt-dependent': 4, '-promoter': 2, 'nbk-mediated': 4, 'incorporating': 2, '-rcmod--reveal': 5, 'provoked': 6, 'mutant--xcomp-': 2, 'facilitate--dep-': 2, 'competence': 2, 'isoform--nsubj-': 2, 'nf-yc': 2, 'homomultimer': 2, 'level--conj': 2, '-agent--overexpression': 9, 'and--overcome': 2, 'cyclin--nsubj-': 2, 'dependent--conj': 3, 'agglutinins': 2, 'telomere-binding': 2, 'with--the': 2, '-nsubj--tip': 2, 'centriolar': 2, 'so-called': 6, '-dep--regulation': 2, 'antagonise--xcomp-': 2, '-nsubj--three': 4, 'ikk-beta': 8, 'nonfunctional': 4, '-dobj--inhibit': 4, 'looked': 6, 'cooperation--dobj-': 4, '-nsubj--taf': 2, 'fact': 4, 'p-jnk': 4, 'restore--conj': 3, 'activator--pobj-': 2, 'circulates': 2, 'terminate': 2, '-dobj--tnf-r': 2, 'convergence': 2, '-me-induced': 2, 'focus-forming': 2, 'phosphatase--dep-': 3, 'with--fragment': 2, 'polarize--ccomp-': 3, '-ccomp--detectable': 3, 'nitroprusside': 2, '-appos--kinase': 4, 'producer': 2, 'qt': 2, 'lipofectin': 4, 'derepress': 2, '-advmod--effectively': 2, 'prolong': 4, 'co-association': 4, 'egfr--dobj-': 2, 'equivalents': 4, 'in--component': 2, 'mapk--nsubj-': 2, 'tapl': 6, 'bcma': 4, 'of--sequence': 3, 'antagonism': 8, 'mgl': 2, 'vim': 2, 'tbk': 3, 'for--epistasis': 2, 'matrigel': 4, 'pecam-': 4, 'dmbta': 2, 'essential--conj': 2, '-dep--render': 2, 'tunel': 4, 'translocalization': 2, 'thalamic': 4, 'in--context': 3, 'pronase': 2, '-ccomp--extend': 2, 'view': 6, 'strengthened': 2, 'record--rcmod-': 2, 'rscd': 4, 'tumor--nsubj-': 3, 'coreceptors': 24, 'coprecipitation': 4, 'correlations': 12, 'cc-rcc--prep': 4, 'operates': 2, 'hut-': 2, 'homolog--appos-': 5, 'accomplish--dep-': 3, 'ptracer-cmv': 2, 'stimulate--conj': 6, 'brs-': 2, '-nsubj--mekk': 2, 'dlg': 8, 'pedigree': 4, 'dle': 2, 'cullin': 5, 'tnf-alpha-related': 4, 'tend': 4, 'key--ccomp-': 3, '-parataxis--lack': 6, 'about--same': 3, 'miep': 4, 'signalosome': 2, 'show--agent-': 4, 'rephosphorylated': 4, 'interfere': 2, 'hypothesis': 12, 'cross-regulation': 2, 'context--prep': 3, 'grow--advcl-': 3, 'kt-': 2, 'comparably': 2, 'to--action': 4, 'homozygosity': 4, 'pka-induced': 2, 'buster': 4, '-mitogen': 2, 'apo-': 3, 'sham-irradiated': 8, 'april': 6, 'to--novo': 2, '-partmod--truncate': 2, 'eth': 2, 'pkd--prep': 3, 'jnk--dep-': 11, 'masses': 2, 'tyr-': 9, 'harmonin': 2, 'rpt': 12, 'co-precipitation': 7, 'nine--nsubj-': 2, 'store-operated': 10, 'in--chromophore': 2, 'influencing': 2, 'reconstitution': 2, '-csubj--phosphorylate': 2, 'aggregation-prone': 2, 'creb-mediated': 2, '-nsubj--decline': 2, 'ucn': 2, 'connotes': 4, 'rpm': 8, 'rpn': 16, 'cbp--dep-': 3, 'dimerized': 2, 'nik--nsubj-': 2, 'ikappabbeta': 4, 'high--advcl-': 4, 'cell-tropic': 2, 'rb-': 2, 'corepressors': 4, 'certain': 8, 'eightfold': 2, 'number--agent-': 2, '-dobj--phosphatidylinositol-': 2, 'nnos': 12, 'protein-protein': 6, 'b--nsubj-': 3, 'recq-like': 2, '-ccomp--ask': 3, 'controller': 2, 'ccaat-enhancer-binding': 2, 'of--erk': 5, 'fa-a': 4, 'hmg': 2, 'fa-d': 4, 'hmb': 4, 'snbs': 2, 'hepatocellular': 8, '-dep--judge': 4, 'fsh-r': 2, '-ccomp--unnecessary': 2, 'fructose-': 2, 'sirna': 19, 'pheromone': 12, 'biphasic': 4, 'va': 2, 'bal-': 8, 'anneals': 4, '-parataxis--enhance': 3, 'vi': 6, '-parataxis--implicate': 3, '-nsubj--activity': 2, 'disengagement': 2, 'b-lymphocyte': 4, 'insulin-stimulated': 8, 'of--transfect': 2, 'adaptor--dep-': 2, 'alpha-actin-positive': 2, '-parataxis--opposite': 2, 'fkn': 2, 'intranasally': 2, 'immunocomplexes': 8, 'appreciable': 6, 'map--rcmod-': 2, 'tccs': 46, '-activating': 2, '-xcomp--family': 4, 'protein-interacting': 4, 'fide': 12, 'c-rel--dep-': 2, 'plzf--prep': 6, 'with--tail': 2, 'unconventional': 2, '-agent--action': 3, '-advcl--connect': 2, 'h-ras-induced': 2, 'phd': 6, 'on--prep-': 3, 'threonine-': 2, 'ser--dobj-': 2, '-dep--transform': 5, 'stnf-r': 28, '-rcmod--utilize': 3, 'preimmune': 12, 'retard': 2, 'nevertheless': 2, 'fmlp': 2, '-nsubj--consistent': 2, 'pdb': 2, 'glycerol-': 2, 'mutlalpha--nsubj-': 6, 'comparisons': 2, 'sensing': 4, 'reintroduce': 2, 'subgroup': 2, 'deltap': 4, 'melanocyte': 6, 'onda': 2, 'mfd': 2, 'compensate--ccomp-': 4, 'reduce--parataxis-': 3, 'sptb': 2, 'on--polymorphism': 3, '-partmod--investigate': 4, 'deltac': 4, 'growth--dobj-': 3, 'divergently': 2, 'bmp-': 4, 'tcc': 26, 'spaciotemporal': 2, 'neurodegeneration': 2, 'extralesional': 2, 'ectopic': 4, 'tumour--prep': 2, 'arresting': 2, 'extends': 6, 'attractive': 4, 'antagonist--pobj-': 4, 'usage': 7, 'nata': 6, 'relocalization--nsubj-': 2, 'unprecedented': 4, 'rfxap': 4, 'dimers': 6, 'mdf-': 2, 'repair-associated': 2, 'clear--ccomp-': 2, 'cpan': 2, 'quercetin-induced': 2, 'sequencing': 2, 'immediately': 2, 'angiotensin-i': 4, 'cadmium-induced': 4, 'egf--nsubjpass-': 2, 'no-induced': 4, 'round--prep': 3, 'with--component': 5, 'alternatively': 14, 'between--level': 4, 'susceptibilities': 4, 'migration--prep': 3, 'displaces': 2, '-advcl--participate': 8, 'pdgf-beta': 4, 'works': 2, 'meleu': 2, 'box-containing': 2, 'a-infected': 2, 'shuttling': 2, 'phenomena': 2, '-dobj--limitation': 2, '-ccomp--positive': 3, 'polymerase--appos-': 2, 'immunodepletion': 2, 'mitomycin': 4, 'oligodeoxynucleotide': 4, 'proper': 2, 'mapkkk--pobj-': 2, '-xcomp--decrease': 2, 'kdr': 8, 'xiap-binding': 2, 'dd-dependent': 2, 'lean': 2, 'receptor-ras-erk': 2, 't-cell-tropic': 4, 'of--hccsmc': 3, 'rstn': 2, '-dobj--that': 3, 'assemble--advcl-': 3, 'limitations': 4, '-rcmod--anneal': 3, 'cas-deficient': 2, 'ensures': 2, 'metalloproteinase-': 2, '-dobj--duration': 5, 'ngfr': 2, 'act--rcmod-': 17, 'approximating': 2, 'dcr': 10, 'chemoattracting': 2, 'jun-d': 2, 'of--apoptosis': 5, '-amino-quinazoline': 2, '-parataxis--release': 4, 'jun-b': 2, 'mspi': 2, 'ndc': 4, '-nsubj--characterization': 2, 'organelle': 2, 'podocyte': 10, 'cochaperone': 4, '-advcl--play': 3, 'rrrcatgyyy': 2, '-nsubj--nr': 2, 'own': 4, 'methoxyverapamil': 2, 'cbfns': 2, 'or--treat': 3, 'pathway-regulating': 10, 'immunostained': 2, 'wortmannin-sensitive': 2, 'cancerous': 4, 'binding-deficient': 4, 'prkdc': 2, '-nsubj--erk': 2, 'c-iap': 15, 'actin-associated': 4, 'reactivate': 6, 'neuronally': 4, 'transcription-independent': 2, 'vac': 2, 'mptp-induced': 8, 'cycle-regulatory': 4, 'for--round': 3, 'delag--conj': 2, 'impair--conj': 3, '-rcmod--block': 3, 'tetraspanin--prep': 2, '-nsubjpass--fusion': 2, '-null': 4, 'eae': 2, 'hu': 2, 'trimester': 2, 'c-related': 2, '-xcomp--mutant': 2, 'upre': 2, '-nsubj--transactivator': 2, 'truncation': 6, 'dma': 2, 'phosphorylate--xcomp-': 4, 'glutamate-mediated': 2, '-dobj--bifurcation': 2, 'from--angiomyolipoma': 2, 'and--induction': 4, 'chaperone': 2, '-amod--r': 2, 'decays': 2, 'glioblastomas': 2, 'a-containing': 10, 'under--condition': 2, 'exon-by-exon': 8, '-pobj--mutant': 7, 'betapix': 2, 'but--restore': 2, 'inactivate--conj': 4, 'hgfr': 2, 'by--kinase': 3, 'galaag': 8, 'pic': 8, 'interaction--dobj-': 8, 'ahrr': 6, 'and--substrate': 4, 'c-sh': 2, 'paralog': 20, 'by--isoform': 3, 'ink': 14, 'to--accumulation': 2, '-dobj--mkk': 3, 'pir': 2, 'ing': 4, 'identities': 2, 'inc': 4, 'since-': 13, 'serine--prep': 2, '-nsubj--collagen': 3, 'trap-positive': 2, 'for--presenilin': 8, 'l-sign': 8, 'hiv-gp': 4, 'c-jun--conj': 16, 'phosphorylation--prep': 3, 'call--prep': 10, 'epo--prep-': 2, 'insults': 2, 'overcome--conj': 2, 'stag': 2, 'ubf': 6, 'of--dag': 2, 'stoichiometric': 2, 'additionally': 4, 'mapkkks': 2, 'correlation--nsubj-': 3, 'pharmacologic': 10, 'anti-tapasin': 2, 'hne': 2, 'auxiliary': 2, 'with--migration': 2, 'irradiation--nsubj-': 3, 'deregulated': 2, 'mitogen-activated--amod-': 3, 'in--polypeptide': 2}

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

######
# Name: normalize
# Input: word
# Return: normalized word
#
# Normalization for word characters
######

def normalize(word):
    return word.encode("ascii", "ignore").replace("'", '_').replace('{', '-_-').replace('}','-__-').replace('"', '-___-').replace(', ,', ',__')

######
# Name: normalize_utf
# Input: text word
# Return: normalized word
#
# Replaces common UTF codes with appropriate characters
######

def normalize_utf(word):
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

######
# Name: extract
# Input: Document object
# Return: relation and features for genegene database table
#
# Extractor code to generate truth label and features
######

def extract(doc):

    for sent in doc.sents:
        genes = []
        lemma = []
        deptree = {}

        if len(sent.words) > 150: continue

        for word in sent.words:
            if word.word in dict_gene_symbols_all:
                genes.append(word)
            deptree[word.insent_id] = {"label":word.dep_label, "parent":word.dep_par}
            lemma.append(word.lemma)

        seen_pair_index = {}

        if len(sent.words) <= 50:
            for w1 in genes:
                for w2 in genes:
                    
                    #list of ambiguous gene symbols to exclude from detected genes
                    gene_exclusion = ["ECM", "EMT", "AML", "CLL", "ALL", "spatial", "PDF", "ANOVA", "MED", "gamma", "San", "RSS", "2F1", "ROS", "zeta", "ADP", "ALS", "GEF", "GAP"]

                    if w1 == w2: continue
                    if w1.word == w2.word: continue
                    if w1.word in gene_exclusion or w2.word in gene_exclusion: continue

                    minindex = min(w1.insent_id, w2.insent_id)
                    maxindex = max(w1.insent_id, w2.insent_id)

                    #don't store pairs twice
                    if minindex in seen_pair_index:
                        if maxindex in seen_pair_index[minindex]:
                            continue
                        else:
                            seen_pair_index[minindex][maxindex] = 1
                    else:
                        seen_pair_index[minindex] = {}
                        seen_pair_index[minindex][maxindex] = 1


                    features = []

                    if w1 == w2:
                        features.append('SAMEOBJ')

                    ## Don't include same genes as mentions ##
                    word1_tmp = re.sub("_", "", w1.word)
                    word2_tmp = re.sub("_", "", w2.word)
                    if w1.word == w2.word or word1_tmp == word2_tmp:
                        mention = None
                        continue

                    ## Don't include same genes as mentions ##
                    if w1.word in dict_name2geneid and w2.word in dict_name2geneid:
                        geneids1 = dict_name2geneid[w1.word]
                        geneids2 = dict_name2geneid[w2.word]

                        for geneid1 in geneids1:
                            if geneid1 in geneids2:
                                mention = None
                                continue

                    ############## FEATURE EXTRACTION ####################################
                    
                    # ##### FEATURE: WORD SEQUENCE BETWEEN MENTIONS AND VERB PATHS #####
                    ws = []
                    verbs_between = []
                    minl_w1 = 100
                    minp_w1 = None
                    minw_w1 = None
                    mini_w1 = None
                    minl_w2 = 100
                    minp_w2 = None
                    minw_w2 = None
                    mini_w2 = None
                    neg_found = 0

                    high_quality_verb = 0
                    for i in range(minindex+1, maxindex):
                        if "," not in sent.words[i].lemma:
                            ws.append(sent.words[i].lemma)
                        if re.search('VB\w*', sent.words[i].pos): # and sent.words[i].lemma != "be":
                            if sent.words[i].word != "{" and sent.words[i].word != "}" and "," not in sent.words[i].word:
                                p_w1 = sent.get_word_dep_path(minindex, sent.words[i].insent_id)
                                p_w2 = sent.get_word_dep_path(sent.words[i].insent_id, maxindex)

                                if len(p_w1) < minl_w1:
                                    minl_w1 = len(p_w1)
                                    minp_w1 = p_w1
                                    minw_w1 = sent.words[i].lemma
                                    mini_w1 = sent.words[i].insent_id

                                if len(p_w2) < minl_w2:
                                    minl_w2 = len(p_w2)
                                    minp_w2 = p_w2
                                    minw_w2 = sent.words[i].lemma
                                    mini_w2 = sent.words[i].insent_id
                                
                                if i > 0:
                                    if sent.words[i-1].lemma in ["no", "not", "neither", "nor"]:
                                        if i < maxindex - 2:
                                            neg_found = 1
                                             features.append("NEG_VERB_BETWEEN_with[%s]" % sent.words[i-1].word + "-" + sent.words[i].lemma)
                                    else:
                                        if sent.words[i] != "{" and sent.words[i] != "}":
                                            verbs_between.append(sent.words[i].lemma)

                    ## Do not include as candidates ##
                    if "while" in ws or "whereas" in ws or "but" in ws or "where" in ws or "however" in ws:
                        mention = None
                        continue

                    ##### FEATURE: HIGH QUALITY PREP INTERACTION PATTERNS #####
                    high_quality_verb = 0
                    if len(verbs_between) == 1 and neg_found == 0:
                        features.append("SINGLE_VERB_BETWEEN_with[%s]" % verbs_between[0])
                        if verbs_between[0] in ["interact", "associate", "bind", "regulate", "phosporylate", "phosphorylated"]:
                            high_quality_verb = 1
                    else:
                        for verb in verbs_between:
                            features.append("VERB_BETWEEN_with[%s]" % verb)

                    if len(ws) == 1 and ws[0] == "and" and minindex > 1:
                        if minindex > 2:
                            if sent.words[minindex - 3].lemma not in ["no", "not", "neither", "nor"] and \
                            sent.words[minindex - 1].lemma in ["of", "between"] and sent.words[minindex - 2].word in ["interaction", "binding"]:
                                high_quality_verb = 1
                        elif sent.words[minindex - 1].lemma in ["of", "between"] and sent.words[minindex - 2].word in ["interaction", "binding"]:
                            high_quality_verb = 1

                    if len(ws) == 1 and ws[0] == "-" and maxindex < len(sent.words) - 1:
                        if sent.words[maxindex + 1].lemma == "complex":
                            high_quality_verb = 1

                    if len(ws) == 1 and ws[0] == "and" and maxindex < len(sent.words) - 1:
                        if sent.words[maxindex + 1].word in ["interaction", "interactions"]:
                            high_quality_verb = 1
                    

                    ##### FEATURE: WORDS BETWEEN MENTIONS #####
                    if len(ws) < 7 and len(ws) > 0 and "{" not in ws and "}" not in ws and "\"" not in ws and "/" not in ws and "\\" not in ws and "," not in ws:
                         if " ".join(ws) not in ["_ and _", "and", "or",  "_ or _"]:
                             features.append("WORDS_BETWEEN_with[%s]" % " ".join(ws))

                    ##### FEATURE: 3-GRAM WORD SEQUENCE #####
                    bad_char = ["\'", "}", "{", "\"", "-", ",", "[", "]"] #think about adding parens
                    if len(ws) > 4 and len(ws) < 15:
                        for i in range(2,len(ws)):
                            if ws[i-2] not in bad_char and ws[i - 1] not in bad_char and ws[i] not in bad_char:
                                if "," not in ws[i-2] and "," not in ws[i-1] and "," not in ws:
                                    features.append("WS_3_GRAM_with[" + ws[i - 2] + "-" + ws[i - 1] + "-" + ws[i]+"]")

                    ##### FEATURE: PREPOSITIONAL PATTERNS #####
                    if minindex > 1:
                        if sent.words[minindex - 2].lemma.lower() in ["association", "interaction", "complex", "activation", "bind", "binding"]:
                            if sent.words[minindex - 1].word.lower() in ["of", "between"] and ("with" in ws or "and" in ws or "to" in ws) and len(ws) ==1:
                                features.append("PREP_PATTERN[{0}_{1}_{2}]".format(sent.words[minindex-2].lemma.lower(), sent.words[minindex-1].word.lower(), sent.words[minindex+1].word.lower()))
                                high_quality_verb = 1


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
                    p = dep_path(deptree, sent, lemma, w1.insent_id, w1.insent_id+1,w2.insent_id,w2.insent_id+1)
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

                    ##### FEATURE: DOMAIN #####
                    domain_words = ["domains", "motif", "motifs", "domain", "site", "sites", "region", "regions", "sequence", "sequences", "elements"]
                    found_domain = 0
                    if minindex > 0:
                        if sent.words[minindex + 1].word in domain_words: 
                            found_domain = 1
                            features.append('GENE_FOLLOWED_BY_DOMAIN_WORD')
                    
                    if maxindex < len(sent.words) - 1 and found_domain == 0:
                        if sent.words[maxindex + 1].word in domain_words: 
                            features.append('GENE_FOLLOWED_BY_DOMAIN_WORD')
                            found_domain = 1

                    
                    ##### FEATURE: PLURAL GENES #####
                    found_plural = 0
                    if minindex > 0:
                        if sent.words[minindex + 1].pos == "NNS" or sent.words[minindex + 1].pos == "NNPS": 
                            found_plural = 1
                            features.append('GENE_M1_FOLLOWED_BY_PLURAL_NOUN_with[%s]' % sent.words[minindex + 1].word)
                    
                    if maxindex < len(sent.words) - 1 and found_plural == 0:
                        if sent.words[maxindex + 1].pos == "NNS" or sent.words[maxindex + 1].pos == "NNPS": 
                            found_plural = 1
                            features.append('GENE_M2_FOLLOWED_BY_PLURAL_NOUN)_with[%s]' % sent.words[maxindex + 1].word)

                    ##### FEATURE: GENE LISTING #####
                    if len(ws) > 0:
                        count = 0
                        flag_not_list = 0
                        for w in ws:

                            #should be comma
                            if count % 4 == 0:
                                if w != '_':
                                    flag_not_list = 1
                            elif count % 4 == 1:
                                if w != ",":
                                    flag_not_list = 1
                            elif count %4 == 2:
                                if w != "_":
                                    flag_not_list = 1
                            #should be a gene
                            else:
                                if w not in dict_gene_symbols_all:
                                    flag_not_list = 1
                            count = count + 1

                        if ws[-1] != '_':
                            flag_not_list = 1

                        if flag_not_list == 0:
                            features.append('GENE_LISTING')

                    if len(ws) > 0:
                        count = 0
                        flag_not_list = 0
                        for w in ws:

                            #should be comma
                            if count % 2 == 0:
                                if w != ',':
                                    flag_not_list = 1

                            #should be a gene
                            else:
                                if w not in dict_gene_symbols_all:
                                    flag_not_list = 1
                            count = count + 1

                        if ws[-1] != ',':
                            flag_not_list = 1

                        if flag_not_list == 0:
                            features.append('GENE_LISTING')

                    #### GENERATE FEATURE ARRAY FOR POSTGRESQL #####
                    feature = "{" + ','.join(features) + '}'

                    mid1 = doc.docid + '_' + '%d' % sent.sentid + '_' + '%d' % w1.insent_id
                    mid2 = doc.docid + '_' + '%d' % sent.sentid + '_' + '%d' % w2.insent_id

                    ############## DISTANT SUPERVISION ###################
                    
                    sent_text = sent.__repr__()
                    if sent_text.endswith("\\"):
                        sent_text = sent_text[0:len(sent_text) - 1]
                    
                    if w1.word in dict_exclude_dist_sup and w2.word in dict_exclude_dist_sup[w1.word]:
                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])

                    elif sent.words[0].word == "Abbreviations" and sent.words[1].word == "used":
                        if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                            print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                            print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                        else:
                            print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])

                    elif w1.word not in dict_abbv and w1.word not in dict_english and w2.word not in dict_english and w2.word not in dict_abbv and w1.word not in dict_domains and w2.word not in dict_domains:
                        if w1.word in dict_interact and w2.word in dict_interact[w1.word] and "mutation" not in sent_text and "mutations" not in sent_text and "variant" not in sent_text and "variants" not in sent_text and "polymorphism" not in sent_text and "polymorphisms" not in sent_text:

                            if found_domain == 0 and flag_family == 0:
                                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "true", feature, sent_text, "\\N"])
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])

                            else:
                                print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                        
                        else:
                            # Negative Example: Mention appear in KB in same doc, but no interaction extracted in KB
                            appear_in_same_doc = False
                            if re.search('^[A-Z]', w1.word) and re.search('^[A-Z]', w2.word):
                                if w1.word in dict_gene_pmid:
                                    for pmid in dict_gene_pmid[w1.word]:
                                        if w2.word in dict_pmid_gene[pmid]:
                                            appear_in_same_doc = True

                            #check if not interact/bind phrase is in ws and not just the words
                            no_interact_phrase = False
                            for j, var in enumerate(ws):
                                if var == 'not' and j + 1 < len(ws) - 1:
                                    if ws[j+1] == "interacts" or ws[j+1] == "interact" or ws[j+1] == "bind":
                                        no_interact_phrase = True

                            if w1.word in dict_no_interact and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws):
                                if w2.word in dict_no_interact[w1.word] and high_quality_verb == 0: 
                                    if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                    else:
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            elif w2.word in dict_no_interact and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws):
                                if w1.word in dict_no_interact[w2.word] and high_quality_verb == 0:
                                    if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                    else:
                                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            elif appear_in_same_doc == True and ("binds" not in ws and "interacts" not in ws and "interacted" not in ws and "bound" not in ws and "complex" not in ws and "associates" not in ws and "associated" not in ws and "bind" not in ws and "interact" not in ws ) and high_quality_verb == 0:
                                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            elif no_interact_phrase == True and high_quality_verb == 0: #("binds" in ws or "interacts" in ws or "bind" in ws or "interact" in ws) and "not" in ws:
                                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            elif w1.ner == "Person" or w2.ner == "Person":
                                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])    
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            elif random.random() < .08 and high_quality_verb == 0:
                                if doc.docid.split(".pdf")[0] not in dict_gs_docids:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "false", feature, sent_text, "\\N"])
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])                        
                                else:
                                    print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                            else:
                                print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])
                    else:
                        print '\t'.join([doc.docid, mid1, mid2, w1.word, w2.word, "\\N", feature, sent_text, "\\N"])


if __name__ == '__main__':
    log("START!")

    load_dict()

    for row in fileinput.input():
        doc = deserialize(row.rstrip('\n'))
        try:
            extract(doc)
        except:
            wrong = True
