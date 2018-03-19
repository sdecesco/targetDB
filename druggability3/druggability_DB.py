#!/usr/bin/env python
import argparse
import time

try:
    import urllib.request as urllib
    from urllib.error import *
except ImportError:
    import urllib
    from urllib2 import *
# To remove danger of using input in python2 environments.
try:
    input = raw_input
except NameError:
    pass

import sys, os, mygene, requests, xlsxwriter, subprocess, re
from Bio import ExPASy
from Bio import SwissProt
from intermine.webservice import Service
from intermine.errors import ServiceError, WebserviceError
from druggability3 import pocket_finder as pocket
from druggability3 import db_connection as db
from druggability3 import pdb_parser
from protein_atlas_api import proteinatlas as patlas
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from operator import itemgetter
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
from druggability3 import cns_mpo as mpo
from druggability3 import drugg_errors

# ===================# SETTING UP PATHS #============================#

dbase_file_path = '/data/sdecesco/databases/druggability/'
output_lists_path = '/data/sdecesco/databases/druggability/outputs/lists/'
output_single_path = '/data/sdecesco/databases/druggability/outputs/single_targets/'

# dbase_file_path = '/Users/stephanedecesco/PycharmProjects/Druggability_3/druggability3/temp/druggability/'
# output_lists_path = '/Users/stephanedecesco/PycharmProjects/Druggability_3/druggability3/temp/druggability/outputs/lists/'
# output_single_path = '/Users/stephanedecesco/PycharmProjects/Druggability_3/druggability3/temp/druggability/outputs/single_targets/'

if not os.path.exists(dbase_file_path):
    os.makedirs(dbase_file_path)
if not os.path.exists(output_lists_path):
    os.makedirs(output_lists_path)
if not os.path.exists(output_single_path):
    os.makedirs(output_single_path)

# ===================# ChembL DATABASE VERSION IN USE #============================#

chembL_db_version = 'chembl_23'
druggability_db = 'druggability'


# ===================# RECOVER THE LIST OF ENTRIES IN THE DBASE #===================#

def get_list_entries():
    dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
    res = dbase.get("SELECT Target_id,Gene_name,Synonyms FROM Targets")
    list_of_entries = {i['Target_id']: i['Gene_name'] for i in res}
    gene_in_db = [
        {'name': i['Gene_name'], 'syn': list(map(str.strip, i['Synonyms'].upper().split(','))), 'ID': i['Target_id']}
        for i
        in res]
    dbase.close()
    return list_of_entries, gene_in_db


# ===================# INITIATE HUMANMINE SERVICE #============================#


service = Service("http://www.humanmine.org/humanmine/service")


class NetworkError(RuntimeError):
    pass


def retryer(max_retries=10, timeout=5):
    def wraps(func):
        request_exceptions = (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.HTTPError, urllib.URLError, urllib.HTTPError, ServiceError, WebserviceError
        )

        def inner(*args, **kwargs):
            for i in range(max_retries):
                try:
                    result = func(*args, **kwargs)
                except request_exceptions:
                    time.sleep(timeout * i)
                    continue
                else:
                    return result
            else:
                raise NetworkError

        return inner

    return wraps


@retryer(max_retries=10, timeout=10)
def gene_to_uniprotid(list_of_gene_name):
    mg = mygene.MyGeneInfo()
    gene_id = {}
    request = mg.querymany(list_of_gene_name, scopes="symbol", fields=['uniprot'], species=9606, as_dataframe=True)

    try:
        uniprot_dict = request['uniprot'].to_dict()
    except KeyError:
        print(
            'None of the genes in the list has a uniprot accession number, make sure you entered the gene name properly')
        sys.exit()
    for gene in list_of_gene_name:
        gene = gene.rstrip('\n')
        gene = gene.rstrip('\r')
        try:
            if type(uniprot_dict[gene]['Swiss-Prot']) is list:
                counter = 0
                for entry in uniprot_dict[gene]['Swiss-Prot']:
                    name = gene + '_' + str(counter)
                    gene_id[name] = entry
                    counter += 1
            else:
                gene_id[gene] = uniprot_dict[gene]['Swiss-Prot']
        except (TypeError, KeyError):
            gene_id[gene] = "No Match"
    return gene_id


@retryer(max_retries=10, timeout=10)
def get_uniprot(gene_id):
    try:
        handle = ExPASy.get_sprot_raw(gene_id)
    except:
        # print(gene_id)
        try:
            handle = ExPASy.get_sprot_raw(gene_id)
        except:
            return None
    try:
        record = SwissProt.read(handle)
    except ValueError:
        record = None
        return record
    return record


def get_crossref_pdb_go_chembl(record):
    # ------------Cross references----------# cross_references contains all the ID number of multiple database such
    # as PDB, CHEMBL, GO, etc... It data structure is composed of a list of tuples, each tuple starting with the
    # database prefix We want to store all the references, sometime one ref (e.g. PDB) contains multiple entries,
    # those multiple entries will be stored as either list or dict in a dict of database example : Dict{'PDB':[
    # 'CODE1','CODE2',...],'GO':{'GO:0006954':'P:inflammatory response','GO:XXXXXX':'C:LocationBlabla'}}
    CrossRef = {}
    PDB_list = {}
    chembl_id = None
    GO_process = []
    GO_Function = []
    GO_Component = []
    for ref in record.cross_references:
        if ref[0] == 'PDB':
            chain_lengths = str(ref[4]).split(',')
            if len(chain_lengths) > 1:
                chain = []
                for i in chain_lengths:
                    if '/' in i:
                        c_list = i.replace(' ', '').split('=')[0].split('/')
                        for c in c_list:
                            chain.append(c)
                    else:
                        chain.append(i.split('=')[0])
            else:
                chain = str(ref[4]).split(',')[0].split('=')[0].split('/')
            for i in range(len(chain)):
                chain[i] = str(chain[i]).strip(' ')
            chain = list(set(chain))
            PDB_list[ref[1]] = {'PDB_code': ref[1], 'Technique': ref[2], 'Resolution': ref[3],
                                'Chain': chain, 'Domain': []}
        elif ref[0] == 'GO':
            GO = ref[2].split(':')
            if GO[0] == 'P':
                GO_process.append(GO[1])
            if GO[0] == 'F':
                GO_Function.append(GO[1])
            if GO[0] == 'C':
                GO_Component.append(GO[1])
        elif ref[0] == 'Ensembl':
            if ref[0] in CrossRef.keys():
                continue
            CrossRef[ref[0]] = ref[3].split('.')[0]
        elif ref[0] == 'ChEMBL':
            chembl_id = ref[1]
        else:
            CrossRef[ref[0]] = ref[1]
    # pdb = PDB_list
    # CrossRef['GO'] = {'P': GO_process, 'F': GO_Function, 'C': GO_Component}
    return CrossRef, PDB_list, {'P': GO_process, 'F': GO_Function, 'C': GO_Component}, chembl_id


@retryer(max_retries=10, timeout=10)
def get_humanmine_data(gene):
    disease = []
    phenotypes = []
    differential_exp_tissues = []
    differential_exp_diseases = []
    gwas = []
    pathways = []

    query = service.new_query("Gene")
    query.add_view("symbol", "id")
    query.add_constraint("Gene", "LOOKUP", gene, code="A")
    query.add_constraint("organism.species", "=", "sapiens", code="B")

    try:
        primary_id = list(query.rows())[0]["id"]
    except IndexError:
        return disease, phenotypes, differential_exp_diseases, differential_exp_tissues, gwas, pathways

    query = service.new_query("Gene")
    query.add_view("diseases.identifier", "diseases.name")
    query.add_constraint("id", "=", primary_id, code="A")

    for row in query.rows():
        disease.append({'disease': row["diseases.name"], 'disease_id': row["diseases.identifier"]})

    query = service.new_query("Gene")
    query.add_view(
        "symbol", "homologues.homologue.organism.shortName",
        "homologues.homologue.alleles.symbol",
        "homologues.homologue.alleles.primaryIdentifier",
        "homologues.homologue.alleles.genotypes.phenotypeTerms.name",
        "homologues.homologue.alleles.genotypes.phenotypeTerms.description",
        "homologues.homologue.alleles.genotypes.name",
        "homologues.homologue.alleles.genotypes.zygosity",
        "homologues.homologue.alleles.type"
    )
    query.add_constraint("id", "=", primary_id, code="A")
    for row in query.rows():
        phenotypes.append({
            'gene': row["symbol"], 'organism': row["homologues.homologue.organism.shortName"],
            'Allele_id': row["homologues.homologue.alleles.primaryIdentifier"],
            'Phenotype': row["homologues.homologue.alleles.genotypes.phenotypeTerms.name"],
            'Phenotype_desc': row["homologues.homologue.alleles.genotypes.phenotypeTerms.description"],
            'genotype': row["homologues.homologue.alleles.genotypes.name"],
            'zygosity': row["homologues.homologue.alleles.genotypes.zygosity"],
            'Allele_symbol': row["homologues.homologue.alleles.symbol"],
            'Allele_type': row["homologues.homologue.alleles.type"]})

    query = service.new_query("Gene")
    query.add_view(
        "atlasExpression.tStatistic", "atlasExpression.condition",
        "atlasExpression.pValue", "atlasExpression.type",
        "atlasExpression.expression"
    )
    query.add_constraint("id", "=", primary_id, code="A")
    query.add_constraint("atlasExpression.type", "!=", "FPKM value", code="B")
    query.add_constraint("atlasExpression.expression", "!=", "NONDE", code="C")
    for row in query.rows():
        if row["atlasExpression.type"] == 'disease_state':
            differential_exp_diseases.append({
                'T_statistic': row["atlasExpression.tStatistic"], 'Condition': row["atlasExpression.condition"],
                'p_value': row["atlasExpression.pValue"],
                'expression_status': row["atlasExpression.expression"]})
        elif row["atlasExpression.type"] == 'organism_part':
            differential_exp_tissues.append({
                'T_statistic': row["atlasExpression.tStatistic"], 'Tissue': row["atlasExpression.condition"],
                'p_value': row["atlasExpression.pValue"],
                'expression_status': row["atlasExpression.expression"]})
        else:
            continue

    query = service.new_query("GWAS")
    query.add_view(
        "results.pValue", "results.phenotype", "firstAuthor",
        "publication.pubMedId", "results.associatedGenes.organism.shortName",
        "publication.year", "publication.doi"
    )
    query.add_sort_order("GWAS.publication.year", "DESC")
    query.add_constraint("results.associatedGenes.id", "=", primary_id, code="A")
    for row in query.rows():
        gwas.append({'p_value': row["results.pValue"], 'phenotype': row["results.phenotype"],
                     'first_author': row["firstAuthor"], 'pubmed_id': row["publication.pubMedId"],
                     'organism': row["results.associatedGenes.organism.shortName"],
                     'publication_year': row["publication.year"], 'doi': row["publication.doi"]})

    query = service.new_query("Gene")
    query.add_view("pathways.name", "pathways.dataSets.name")
    query.add_constraint("id", "=", primary_id, code="A")
    query.add_constraint("organism.name", "=", "Homo sapiens", code="B")

    for row in query.rows():
        pathways.append({'pathway_name': row["pathways.name"], 'pathway_dataset': row["pathways.dataSets.name"]})

    return disease, phenotypes, differential_exp_diseases, differential_exp_tissues, gwas, pathways


def get_domains(record=None, gene_id=None, chembl_id=None):
    # ------------ domain and binding sites ----------#
    domain = []
    if not record:
        if gene_id:
            record = get_uniprot(gene_id)
            if record is None:
                return None
        else:
            print("[ArgumentError]: Combination of arguments is invalid")
            raise drugg_errors.ArgumentError

    for feature in record.features:
        if feature[0] == 'DOMAIN':
            start = feature[1]
            finish = feature[2]
            while not str(start)[0].isdigit():
                start = str(start)[1:]
            while not str(finish)[0].isdigit():
                finish = str(finish)[1:]
            domain_name = str(feature[3]).split('.')[0]
            domain_id = str(record.accessions[0]) + str(start) + str(finish) + '_uniprot'
            domain.append(
                {'Start': int(start), 'Stop': int(finish), 'name': domain_name, 'length': int(finish) - int(start),
                 'domain_id': domain_id, 'source_name': 'Uniprot', 'Source_id': 'n.a.'})
    if chembl_id:
        dbase = db.open_db(chembL_db_version, pwd=args.db_password, user=args.db_username)
        query = "SELECT CD.start_position,CD.end_position,DOM.domain_name,DOM.source_domain_id,DOM.domain_type FROM target_dictionary " \
                "TD,target_components TC,domains DOM,component_domains CD WHERE TD.chembl_id='%s' AND TD.tid=TC.tid AND " \
                "TC.component_id=CD.component_id AND DOM.domain_id=CD.domain_id GROUP BY TD.chembl_id,CD.domain_id" % \
                chembl_id
        domains = dbase.get(query)
        dbase.close()
        for i in domains:
            domain_id = str(record.accessions[0]) + str(i['start_position']) + str(i['end_position']) + '_' + i[
                'source_domain_id']
            domain.append({'Start': int(i['start_position']), 'Stop': int(i['end_position']), 'name': i['domain_name'],
                           'length': int(i['end_position']) - int(i['start_position']),
                           'domain_id': domain_id, 'source_name': i['domain_type'], 'Source_id': i['source_domain_id']})
        items_to_delete = []
        for i in range(len(domain)):
            for j in range(len(domain)):
                if i == j:
                    continue
                elif j in items_to_delete:
                    continue
                elif domain[i]['domain_id'].split('_')[0] == domain[j]['domain_id'].split('_')[0]:
                    items_to_delete.append(i)
        for i in reversed(items_to_delete):
            domain.pop(i)
    return domain


def get_chembl_info(chembl_id):
    dbase = db.open_db(chembL_db_version, pwd=args.db_password, user=args.db_username)
    query = "SELECT PC.pref_name ,PC.short_name ,PC.protein_class_desc ,group_concat(DISTINCT(" \
            "CSYN.component_synonym)) AS Synonym FROM target_dictionary TD, target_components TC, " \
            "component_sequences CS, component_class CC, protein_classification PC, " \
            "component_synonyms CSYN WHERE TD.chembl_id='%s' AND TD.tid=TC.tid AND " \
            "TC.component_id=CS.component_id AND TC.component_id=CC.component_id AND " \
            "TC.component_id=CSYN.component_id AND CC.protein_class_id=PC.protein_class_id GROUP BY " \
            "TD.chembl_id" % chembl_id
    entry_info = dbase.get(query)
    dbase.close()
    return entry_info


def synonyms(record):
    try:
        synonym = record.gene_name.split(';')[1].split('=')[1]
    except IndexError:
        synonym = 'Not found'
    return synonym


def get_comments(record):
    comments = {}
    for i in record.comments:
        comment = i.split(':', 1)
        comments[comment[0]] = comment[1]
        # for key,value in self.comments.items():
        # print key+": "+value
    return comments


def get_variants(record, domains):
    modifications = {}
    for feat in record.features:
        previous = ''
        new = ''
        action = ''
        comment = ''
        if feat[0] in ['VAR_SEQ', 'VARIANT']:
            if str(feat[3]).startswith('Missing'):
                action = 'remove'
                previous = ''
                new = ''
                try:
                    comment = '(' + str(feat[3]).split('(')[1]
                except IndexError:
                    comment = '(' + str(feat[3]).split('(')[0]
            elif '->' in str(feat[3]):
                mod = str(feat[3]).split('(')[0].replace(' ', '').split('->')
                previous = mod[0]
                try:
                    new = mod[1].split('.')[0]
                except IndexError:
                    new = mod[1]
                action = 'replace'
                try:
                    comment = '(' + str(feat[3]).split('(')[1]
                except IndexError:
                    comment = '(' + str(feat[3]).split('(')[0]
            modifications[feat[4]] = {'start': feat[1], 'stop': feat[2], 'previous': previous,
                                      'new': new, 'action': action, 'comment': comment, 'domains': []}
        if feat[0] == 'MUTAGEN':
            if str(feat[3]).startswith('Missing'):
                action = 'remove'
                previous = ''
                new = ''
                try:
                    comment = str(feat[3]).split(':', 1)[1]
                except IndexError:
                    comment = str(feat[3]).split(':', 1)[0]
            elif '->' in str(feat[3]):
                mod = str(feat[3]).split(':')[0].replace(' ', '').split('->')
                previous = mod[0]
                new = mod[1]
                try:
                    comment = str(feat[3]).split(':', 1)[1]
                except IndexError:
                    comment = str(feat[3]).split(':', 1)[0]
                action = 'replace'
            modifications['MUTAGEN_' + str(feat[1])] = {'start': feat[1], 'stop': feat[2], 'previous': previous,
                                                        'new': new, 'action': action, 'comment': comment, 'domains': []}
    for d in domains:
        for m in modifications.keys():
            if modifications[m]['start'] >= d['Start'] and modifications[m]['stop'] <= d['Stop']:
                modifications[m]['domains'].append(d['name'])

    variants = []
    for item in record.comments:
        if str(item).startswith('ALTERNATIVE'):
            splitted = str(item).split('; Name=')[1:]
            for i in splitted:
                isoform = i.replace(' ', '').split(';')
                name = isoform[0].split('{')[0]
                isoid = ''
                seq_mod = ''
                for j in isoform:
                    if j.startswith('IsoId='):
                        isoid = j.split('=')[1]
                    if j.startswith('Sequence'):
                        seq_mod = j.split('=')[1].split(',')
                variants.append({'name': name, 'isoid': isoid, 'seq_mod': seq_mod})
    if args.verbose:
        print("[ISOFORMS SEQUENCE ALIGNMENT IN PROGRESS]")
    for isoform in variants:
        temp_seq = list(record.sequence)
        if 'Displayed' in isoform['seq_mod']:
            isoform['sequence'] = record.sequence
            isoform['seq_list'] = list(record.sequence)
        elif 'External' in isoform['seq_mod']:
            isoform['sequence'] = get_uniprot(isoform['isoid']).sequence
            isoform['seq_list'] = list(isoform['sequence'])
        else:
            for changes in isoform['seq_mod']:
                if changes == 'Notdescribed':
                    temp_seq = ['no_sequence']
                    break
                if modifications[changes]['action'] == 'remove':
                    start = modifications[changes]['start'] - 1
                    stop = modifications[changes]['stop']
                    for i in range(start, stop):
                        temp_seq[i] = ''
                if modifications[changes]['action'] == 'replace':
                    start = modifications[changes]['start'] - 1
                    stop = modifications[changes]['stop']
                    temp_seq[start] = modifications[changes]['new']
                    for i in range(start + 1, stop):
                        temp_seq[i] = ''

            isoform['sequence'] = ''.join(temp_seq)
            isoform['seq_list'] = temp_seq
        seq_object = Seq(isoform['sequence'])
        isoform['alignment'] = align(record.sequence, seq_object)
    if args.verbose:
        print("[ISOFORMS SEQUENCE ALIGNMENT DONE]")
    return variants, modifications


def align(sequence1, sequence2, end_gaps=True, print_align=False):
    if sequence1 == 'no_sequence' or sequence2 == 'no_sequence' or sequence1 == '' or sequence2 == '':
        return None

    seq1 = sequence1.replace('O', 'X').replace('U', 'X')
    seq2 = str(sequence2).replace('O', 'X').replace('U', 'X')

    alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5, one_alignment_only=True,
                                          penalize_end_gaps=(end_gaps, end_gaps))
    if not alignments:
        alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5, one_alignment_only=True,
                                              penalize_end_gaps=(False, False))
    if not alignments:
        return None
    try:
        total_length = len(alignments[0][1])

        match = 0
        gap = 0
        for i in range(len(alignments[0][0])):
            if alignments[0][0][i] == alignments[0][1][i]:
                match += 1
            elif alignments[0][0][i] == '-' or alignments[0][1][i] == '-':
                gap += 1
        identity = round((match / total_length) * 100, 2)
        gaps = round((gap / total_length) * 100, 2)
        score = alignments[0][2]
        if print_align:
            print(pairwise2.format_alignment(*alignments[0]))
        results = {'score': score, 'identity': identity, 'gaps': gaps}
    except IndexError:
        return None
    return results


@retryer(max_retries=10, timeout=10)
def get_pdb(list_of_pdb, path):
    file_path = path + 'PDB/'
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    for k in list_of_pdb.keys():
        saved_pdb = file_path + str(k) + '.pdb'
        # print('[PDB DOWNLOAD]: ' + str(list_of_pdb[k]['PDB_code']) + '.pdb')
        if os.path.isfile(saved_pdb):
            if args.verbose:
                print('[FILE ALREADY THERE]: ' + str(k) + '.pdb')
            list_of_pdb[k]['path'] = saved_pdb
        else:
            try:
                urllib.urlretrieve('http://files.rcsb.org/download/' + str(k) + '.pdb', saved_pdb)
                list_of_pdb[k]['path'] = saved_pdb
                if args.verbose:
                    print('[PDB DOWNLOAD]: ' + str(k) + '.pdb')
            except HTTPError:
                if args.verbose:
                    print('[ERROR]: ' + str(k) + '.pdb file not found (do not exist or is '
                                                 'too big)')
    if args.verbose:
        print('[PDB DOWNLOAD DONE]')
    return file_path


def get_pdb_seq_info(pdb_list, domain, isoforms):
    if args.verbose:
        print("[PDB INFO]: Extracting PDB information")
    for keys in pdb_list.keys():
        if 'path' not in pdb_list[keys].keys():
            pdb_list[keys]['length'] = 0
            pdb_list[keys]['Domain'] = 'no_domain'
            pdb_list[keys]['sequences'] = {}
            continue
        seq = pdb_parser.get_sequence(keys, pdb_list[keys]['path'],
                                      pdb_list[keys]['Chain'], domain)
        chain_done = []
        for chain in seq.keys():
            if seq[chain]['chain_name'] not in pdb_list[keys]['Chain']:
                continue
            for i in seq[chain]['equal']:
                if i in chain_done:
                    seq[chain]['aligned'] = seq[i]['aligned']
                    chain_done.append(chain)
                    break
            if seq[chain]['aligned']:
                continue
            if isoforms:
                max_score = 0
                isoform_match = []
                list_of_align = []
                for variants in isoforms:
                    partial_seq = ''.join(variants['seq_list'][seq[chain]['start'] - 1:seq[chain]['stop'] - 1])
                    seq_align = align(partial_seq, seq[chain]['sequence'], end_gaps=False)
                    if seq_align is None:
                        continue
                    if seq_align['score'] > max_score:
                        max_score = seq_align['score']
                    list_of_align.append((seq_align, variants['isoid']))
                for item in list_of_align:
                    if item[0]['score'] == max_score:
                        isoform_match.append({'isoform': item[1], 'identity': item[0]['identity']})
                seq[chain]['aligned'] = isoform_match
                chain_done.append(chain)
        pdb_list[keys]['sequences'] = seq
        length = 0
        for chain in seq.keys():
            if seq[chain]['length'] > length:
                length = seq[chain]['length']
            for c in seq[chain]['domain']:
                pdb_list[keys]['Domain'].append(c)

        pdb_list[keys]['Domain'] = list(set(pdb_list[keys]['Domain']))
        pdb_list[keys]['length'] = length
    if args.verbose:
        print("[PDB INFO]: Done")


def get_ligands_to_do(chembl_code):
    # ====================# GETTING THE LIST OF LIGAND WITH BIOACTIVITIES IN THE DB #=========================#
    lig_to_do = []
    if args.verbose:
        print("[LIGANDS]: Extracting ligand informations")
    # =====================# GETTING THE LIST OF LIGAND ASSOCIATED TO  THE TARGET #==========================#
    dbase = db.open_db(chembL_db_version, pwd=args.db_password, user=args.db_username)
    query_lig_target = "SELECT TD.chembl_id AS target_id,MOL.chembl_id AS lig_id,version.name as " \
                       "chembl_version FROM target_dictionary TD,activities BIO,assays AC," \
                       "molecule_dictionary MOL,version WHERE TD.chembl_id='%s' AND TD.tid=AC.tid AND " \
                       "AC.assay_id=BIO.assay_id AND BIO.published_value is not null AND " \
                       "BIO.molregno=MOL.molregno GROUP BY BIO.molregno" % chembl_code

    res_lig_target = dbase.get(query_lig_target)
    dbase.close()
    if not res_lig_target:
        return lig_to_do

    count = 0
    limit = 500
    if len(res_lig_target) <= limit:
        limit = len(res_lig_target)
    n_loop = 0
    excess = len(res_lig_target) % limit
    dbase2 = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
    lig_ids = '('
    query_lig_bioact_db_base = """SELECT DISTINCT lig_id,chembl_version FROM bioactivities WHERE lig_id IN """
    res_lig_in_db = ()
    for i in range(len(res_lig_target)):
        count += 1
        lig_ids += "'" + res_lig_target[i]['lig_id'] + "',"
        if count % limit == 0:
            lig_ids = lig_ids.rstrip(',') + ')'
            full_query = query_lig_bioact_db_base + lig_ids
            tmp_res = dbase2.get(full_query)
            lig_ids = '('
            count = 0
            n_loop += 1
            res_lig_in_db = res_lig_in_db + tmp_res
    if excess != 0:
        lig_ids = lig_ids.rstrip(',') + ')'
        full_query = query_lig_bioact_db_base + lig_ids
        tmp_res = dbase2.get(full_query)
        res_lig_in_db = res_lig_in_db + tmp_res
    dbase2.close()

    # lig_ids = '('
    # for i in res_lig_target:
    #     lig_ids += "'" + i['lig_id'] + "',"
    # lig_ids = lig_ids.rstrip(',') + ')'
    # dbase2 = db.open_db(druggability_db,pwd=args.db_password,user=args.db_username)
    # query_lig_bioact_db = """SELECT DISTINCT lig_id,chembl_version FROM bioactivities WHERE lig_id IN """ + lig_ids
    # res_lig_in_db = dbase2.get(query_lig_bioact_db)
    dict_lig_in_db = {i['lig_id']: i['chembl_version'] for i in res_lig_in_db}
    # dbase2.close()

    # =========# ADDING IN THE TO-DO LIST ONLY LIGAND WITH NO BIOACTIVITY IN THE DB #========================#

    for i in res_lig_target:
        if i['lig_id'] in dict_lig_in_db.keys():
            if int(str(i['chembl_version']).split('_')[1]) <= int(str(dict_lig_in_db[i['lig_id']]).split('_')[1]):
                pass
            else:
                lig_to_do.append(i['lig_id'])
        else:
            lig_to_do.append(i['lig_id'])
    return lig_to_do


def get_assays(target_chembl):
    if args.verbose:
        print("[ASSAYS]: Extracting assays information")
    dbase = db.open_db(chembL_db_version, pwd=args.db_password, user=args.db_username)
    query = "SELECT TD.chembl_id AS target_id, AC.chembl_id AS assay_id FROM target_dictionary TD,assays AC WHERE " \
            "TD.chembl_id = '%s' AND TD.tid = AC.tid GROUP BY AC.chembl_id" % target_chembl

    res = dbase.get(query)
    dbase.close()

    return res


def get_bioactivity(lig_to_do):
    # =========# CREATING THE STRING OF LIGANDS ('CHEMBLXXXX','CHEMBLXXXX',....) #========================#
    if args.verbose:
        print("[BIOACTIVITIES]: Extracting bioactivities")

    dbase = db.open_db(chembL_db_version, pwd=args.db_password, user=args.db_username)

    # ========# FETCHING ALL BIOACTIVITIES (FOR THE LIST OF LIGANDS: ALL TARGETS) #=======================#
    res_bioactivity = ()
    lig_str = ''
    counter = 0
    for i in lig_to_do:
        counter += 1
        lig_str += "'" + str(i) + "',"
        if counter == 1000:
            lig_str = lig_str.rstrip(',')
            query_bioactivity = "SELECT EXP.chembl_id AS assay_id,TD.chembl_id AS Target_id,TD.target_type AS target_type," \
                                "TD.pref_name AS target_name,TD.organism AS target_organism, MOL.chembl_id AS " \
                                "lig_id,ACT.standard_units AS units,ACT.standard_relation AS operator," \
                                "ACT.standard_value AS value_num,DOC.doi AS doi,CONCAT((CASE WHEN " \
                                "ACT.standard_relation!='=' THEN ACT.standard_relation ELSE '' END)," \
                                "ROUND(ACT.standard_value,2),' ',CASE WHEN ACT.standard_units is null THEN '' " \
                                "ELSE ACT.standard_units END) AS value_text,VER.name AS chembl_version," \
                                "ACT.standard_type AS standard_type,ACT.activity_comment," \
                                "ACT.data_validity_comment,ACT.pchembl_value FROM molecule_dictionary MOL," \
                                "activities ACT,assays EXP,target_dictionary TD,docs DOC,version VER WHERE " \
                                "MOL.chembl_id in (%s) AND MOL.molregno = ACT.molregno AND " \
                                "ACT.assay_id=EXP.assay_id AND ACT.published_value is not NULL AND EXP.tid=TD.tid " \
                                "AND ACT.doc_id=DOC.doc_id" % lig_str
            res_bioactivity += dbase.get(query_bioactivity)
            counter = 0
            lig_str = ''

    if counter < 1000:
        lig_str = lig_str.rstrip(',')
        query_bioactivity = "SELECT EXP.chembl_id AS assay_id,TD.chembl_id AS Target_id,TD.target_type AS target_type," \
                            "TD.pref_name AS target_name,TD.organism AS target_organism, MOL.chembl_id AS " \
                            "lig_id,ACT.standard_units AS units,ACT.standard_relation AS operator," \
                            "ACT.standard_value AS value_num,DOC.doi AS doi,CONCAT((CASE WHEN " \
                            "ACT.standard_relation!='=' THEN ACT.standard_relation ELSE '' END)," \
                            "ROUND(ACT.standard_value,2),' ',CASE WHEN ACT.standard_units is null THEN '' " \
                            "ELSE ACT.standard_units END) AS value_text,VER.name AS chembl_version," \
                            "ACT.standard_type AS standard_type,ACT.activity_comment," \
                            "ACT.data_validity_comment,ACT.pchembl_value FROM molecule_dictionary MOL," \
                            "activities ACT,assays EXP,target_dictionary TD,docs DOC,version VER WHERE " \
                            "MOL.chembl_id in (%s) AND MOL.molregno = ACT.molregno AND " \
                            "ACT.assay_id=EXP.assay_id AND ACT.published_value is not NULL AND EXP.tid=TD.tid " \
                            "AND ACT.doc_id=DOC.doc_id" % lig_str
        res_bioactivity += dbase.get(query_bioactivity)
    return res_bioactivity


def blast_launcher(sequence, seq_file, db, output_name, num_core=8):
    if os.path.isfile(seq_file):
        subprocess.check_output(
            ['/ssddata/sdecesco/blast/bin/blastp', '-db', str(db), '-query', str(seq_file), '-out', str(output_name),
             '-num_threads', str(num_core), '-outfmt', str(5), '-max_target_seqs', str(100)],
            env={'BLASTDB': '/data/sdecesco/databases/blastdb'})
    else:
        save_file = open(seq_file, 'w')
        save_file.write(str(sequence))
        save_file.close()
        subprocess.check_output(
            ['/ssddata/sdecesco/blast/bin/blastp', '-db', str(db), '-query', str(seq_file), '-out', str(output_name),
             '-num_threads', str(num_core), '-outfmt', str(5), '-max_target_seqs', str(100)],
            env={'BLASTDB': '/data/sdecesco/databases/blastdb'})


def pdb_blast(sequence, path, gene_id, gene='', pdb_list=[]):
    file_path = path + 'PDB_BLAST/'
    blast_file = file_path + gene + '_' + gene_id + '.xml'
    seq_file = file_path + gene + '_' + gene_id + '.seq'
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    if args.verbose:
        print('[3D BLAST]:' + gene + '(' + gene_id + ')')
    alternate_pdb = {}
    if os.path.isfile(blast_file):
        if ((((time.time() - os.path.getmtime(blast_file)) / 60) / 60) / 24) <= 15:
            if args.verbose:
                print('[3D BLAST FILE FOUND]:' + gene)
            result_handle = open(blast_file)
        else:
            if args.verbose:
                print('[3D BLAST FILE FOUND]: File older than 2 weeks, a new blast will be performed (' + gene + ')')
            os.remove(blast_file)
            blast_launcher(sequence, seq_file, 'pdbaa', blast_file, num_core=args.num_core)
            if os.path.isfile(blast_file):
                result_handle = open(blast_file)
            else:
                print("[3D BLAST][ERROR]: Something went wrong, no blast result generated")
                return alternate_pdb
    else:
        blast_launcher(sequence, seq_file, 'pdbaa', blast_file, num_core=args.num_core)
        if os.path.isfile(blast_file):
            result_handle = open(blast_file)
        else:
            print("[3D BLAST][ERROR]: Something went wrong, no blast result generated")
            return alternate_pdb

    if args.verbose:
        print('[3D BLAST DONE]: Now parsing the data - ' + gene + '(' + gene_id + ')')
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    e_value_treshold = 0.0001
    query_length = float(blast_record.query_length)
    for alignment in blast_record.alignments:
        pdb_code = alignment.title.split('|')[3]
        if pdb_code in pdb_list:
            continue
        chain_id = alignment.accession
        for hsp in alignment.hsps:
            length = float(hsp.align_length)
            percent_seq = round((length / query_length) * 100, 1)
            similarity = round((float(hsp.positives) / length) * 100, 1)
            if hsp.score > 200 and hsp.expect < e_value_treshold and similarity > 60:
                alternate_pdb[pdb_code] = {'pdb_code': pdb_code, 'seq_percent': percent_seq, 'similarity': similarity,
                                           'Chain_id': chain_id, 'chain_letter': chain_id.split('_')[1]}
    get_pdb(alternate_pdb, path)
    for pdb in alternate_pdb.values():
        if 'path' in pdb.keys():
            result = pdb_to_uniprot(pdb['chain_letter'], pdb['path'])
            if result:
                pdb['gene'] = result['gene']
                pdb['organism'] = result['organism']
                pdb['uniprot_id'] = result['uniprot_id']
            else:
                pdb['gene'] = 'NA'
                pdb['organism'] = 'NA'
                pdb['uniprot_id'] = 'NA'
        else:
            pdb['gene'] = 'NA'
            pdb['organism'] = 'NA'
            pdb['uniprot_id'] = 'NA'
            continue

    return alternate_pdb


def pdb_to_uniprot(chain, pdb_file):
    result = {}
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith('DBREF'):
                current_line = list(filter(None, line.strip('\n').replace(' ', '_').split('_')))
                if current_line[2] != chain:
                    continue
                if current_line[0] == 'DBREF1':
                    if current_line[5] == 'UNP':
                        try:
                            result['gene'] = current_line[6]
                            result['uniprot_id'] = current_line[6]
                            result['organism'] = current_line[7]
                        except IndexError:
                            result['gene'] = 'NA'
                            result['organism'] = 'NA'
                            result['uniprot_id'] = 'NA'
                else:
                    if current_line[5] == 'UNP':
                        try:
                            result['gene'] = current_line[7]
                            result['uniprot_id'] = current_line[6]
                            result['organism'] = current_line[8]
                        except IndexError:
                            result['gene'] = 'NA'
                            result['organism'] = 'NA'
                            result['uniprot_id'] = 'NA'
    return result


def proteins_blast(sequence, gene_id, gene, path):
    file_path = path + 'PROTEIN_BLAST/'
    blast_file = file_path + gene + '_' + gene_id + '.xml'
    seq_file = file_path + gene + '_' + gene_id + '.seq'
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    if args.verbose:
        print('[PROTEIN BLAST] ' + gene + '(' + gene_id + ')')
    if os.path.isfile(blast_file):
        if ((((time.time() - os.path.getmtime(blast_file)) / 60) / 60) / 24) <= 15:
            if args.verbose:
                print('[PROTEIN BLAST FILE FOUND]:' + gene)
            result_handle = open(blast_file)
        else:
            if args.verbose:
                print(
                    '[PROTEIN BLAST FILE FOUND]: File older than 2 weeks, a new blast will be performed (' + gene + ')')
            os.remove(blast_file)
            blast_launcher(sequence, seq_file, 'swissprot', blast_file, num_core=args.num_core)
            if os.path.isfile(blast_file):
                result_handle = open(blast_file)
            else:
                print("[PROTEIN BLAST][ERROR]: Something went wrong, no blast result generated")
                return []
    else:
        blast_launcher(sequence, seq_file, 'swissprot', blast_file, num_core=args.num_core)
        if os.path.isfile(blast_file):
            result_handle = open(blast_file)
        else:
            print("[PROTEIN BLAST][ERROR]: Something went wrong, no blast result generated")
            return []
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    if args.verbose:
        print('[PROTEIN BLAST DONE]: Now parsing the data - ' + gene + '(' + gene_id + ')')
    e_value_treshold = 0.0001
    query_length = float(blast_record.query_length)
    list_of_neighbours = []
    list_of_accession_ID = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            length = float(hsp.align_length)
            neighbour_accession_code = alignment.accession
            neighbour_gene_name = alignment.hit_id.split('|')[-1].split('_')[0]
            neighbour_gene_species = alignment.hit_id.split('|')[-1].split('_')[1]
            if neighbour_accession_code in list_of_accession_ID:
                continue
            else:
                list_of_accession_ID.append(neighbour_accession_code)

            percent_seq = round((length / query_length) * 100, 1)
            similarity = round((float(hsp.positives) / length) * 100, 1)
            if hsp.score > 200 and hsp.expect < e_value_treshold and (
                    length / query_length) > 0.5 and similarity >= 40 and neighbour_accession_code != gene_id:
                list_of_neighbours.append(
                    (neighbour_accession_code, percent_seq, similarity, neighbour_gene_name, neighbour_gene_species))
    return list_of_neighbours


def write_to_db(target, db_name):
    if target.record is None:
        return None

    # ========# OPEN DATABASE #========#
    if args.verbose:
        print("[DATABASE]: Start to write info into the database")
    dbase = db.open_db(db_name, pwd=args.db_password, user=args.db_username)

    # ========# FILLING THE TARGETS TABLE #=========#
    if not target.target_in:
        location = ""
        for item in target.go['C']:
            location += str(item) + " / "
        location = location.rstrip('/')
        function = ""
        for item in target.go['F']:
            function += str(item) + '/'
        function = function.rstrip('/')
        process = ""
        for item in target.go['P']:
            process += str(item) + " / "
        process = process.rstrip('/')

        values = [(
            target.swissprotID, target.gene, target.record.organism, target.record.taxonomy_id[0], target.sequence,
            location, process,
            function, str(len(target.variants)), target.synonyms, target.ProteinClass['Class_name'],
            target.ProteinClass['Class_description'], target.ProteinClass['Class_short'], target.chembl_id)]
        dbase.multi_write("""INSERT INTO Targets(Target_id, Gene_name, Species, species_id, Sequence, Cell_location,
        Process, Function, Number_isoforms, Synonyms, Protein_class, Protein_class_desc, Protein_class_short,chembl_id)""",
                          values, update="ON DUPLICATE KEY UPDATE Gene_name=VALUES(Gene_name),Species=VALUES(Species),"
                                         "species_id=VALUES(species_id),Sequence=VALUES(Sequence),Cell_location=VALUES("
                                         "Cell_location),Process=VALUES(Process),Function=VALUES("
                                         "Function),Number_isoforms=VALUES(Number_isoforms),Synonyms=VALUES("
                                         "Synonyms),Protein_class=VALUES(Protein_class),Protein_class_desc=VALUES("
                                         "Protein_class_desc), "
                                         "Protein_class_short=VALUES(Protein_class_short)")
        if args.verbose:
            print("[DATABASE]: Target table populated/updated")

    # ========# FILLING THE DOMAIN TABLE #=========#

    if target.domain:
        dom_values = []
        for dom in target.domain:
            dom_values.append((target.swissprotID, dom['Start'], dom['Stop'], dom['name'], dom['domain_id'],
                               dom['source_name'], dom['Source_id'], dom['length']))
        dbase.multi_write(
            "INSERT INTO Domain_targets(Target_id, Domain_start, Domain_stop, Domain_name, Domain_id,source_name,"
            "Source_id,length)", dom_values,
            update='ON DUPLICATE KEY UPDATE source_name=VALUES(source_name),Source_id=VALUES(Source_id)')
        if args.verbose:
            print("[DATABASE]: Domain table populated/updated")

    # ========# FILLING THE 3D BLAST TABLE #=========#

    if target.alternate_pdb:
        values = []
        for key in target.alternate_pdb.keys():
            values.append((target.swissprotID, target.alternate_pdb[key]['uniprot_id'], key,
                           target.alternate_pdb[key]['Chain_id'], target.alternate_pdb[key]['similarity'],
                           target.alternate_pdb[key]['organism'], target.alternate_pdb[key]['gene'],
                           target.alternate_pdb[key]['chain_letter'],
                           target.swissprotID + '_' + target.alternate_pdb[key]['uniprot_id'] + '_' +
                           target.alternate_pdb[key]['Chain_id']))

        dbase.multi_write(
            "INSERT INTO 3D_Blast(Query_target_id, Hit_gene_id, Hit_PDB_code, Hit_Chain_id, similarity,Hit_gene_species,Hit_gene_name,Chain_Letter,Unique_ID)",
            values, update="ON DUPLICATE KEY UPDATE Hit_gene_species=VALUES(Hit_gene_species),similarity=VALUES("
                           "similarity),Hit_gene_name=VALUES(Hit_gene_name)")
        if args.verbose:
            print("[DATABASE]: PDB Blast table populated/updated")

    # ========# FILLING THE PROTEIN BLAST TABLE #=========#

    if target.neighbours:
        values = []
        for item in target.neighbours:
            values.append(
                (target.swissprotID, item[0], item[1], item[2], item[3], item[4], target.swissprotID + '_' + item[0]))
        dbase.multi_write(
            "INSERT INTO protein_blast(Query_target_id, Hit_gene_id, sequence_in_common, similarity,Hit_gene_name,Hit_gene_species,Unique_ID)",
            values, update="ON DUPLICATE KEY UPDATE sequence_in_common=VALUES(sequence_in_common),"
                           "similarity=VALUES(similarity),Hit_gene_name=VALUES(Hit_gene_name),"
                           "Hit_gene_species=VALUES(Hit_gene_species)")
        if args.verbose:
            print("[DATABASE]: Blast table populated/updated")

    # ========# FILLING THE MODIFICATIONS TABLE #=========#

    if target.modifications:
        mod_values = []
        for mod in target.modifications.keys():
            mod_type = mod.split('_')[0]
            mod_domains = ''
            for i in target.modifications[mod]['domains']:
                mod_domains += str(i) + ','
            mod_domains = mod_domains.rstrip(',')
            mod_values.append((mod, mod_type, target.modifications[mod]['action'], target.modifications[mod]['comment'],
                               target.modifications[mod]['previous'], target.modifications[mod]['new'],
                               target.modifications[mod]['start'], target.modifications[mod]['stop'],
                               str(target.swissprotID) + '_' + mod, str(target.swissprotID), mod_domains))

        dbase.multi_write("INSERT INTO modifications(mod_id,mod_type, action, comment, previous, new, start, stop,"
                          "Unique_modID,Target_id,domains)", mod_values)

    # ========# FILLING THE ISOFORMS + ISOFORMS MODIFICATIONS TABLES #=========#

    if target.variants:
        var_values = []
        var_mod_values = []
        for var in target.variants:
            canonical = 0
            if 'Displayed' in var['seq_mod']:
                canonical = 1
            elif 'External' in var['seq_mod']:
                pass
            else:
                for mod in var['seq_mod']:
                    if mod == 'Notdescribed':
                        continue
                    var_mod_values.append((var['isoid'], target.swissprotID + '_' + mod))
            if var['alignment'] is None:
                var_values.append((target.swissprotID, var['name'], var['isoid'], var['sequence'], len(var['sequence']),
                                   canonical, '', '', ''))
            else:
                var_values.append((target.swissprotID, var['name'], var['isoid'], var['sequence'], len(var['sequence']),
                                   canonical, var['alignment']['identity'], var['alignment']['gaps'],
                                   var['alignment']['score']))

                # ====# ISOFORMS #=====#

        dbase.multi_write(
            "INSERT INTO Isoforms(Target_id, Isoform_name, Isoform_id, Sequence, n_residues, Canonical, Identity, "
            "Gaps,Score)", var_values, update="ON DUPLICATE KEY UPDATE Date=CURRENT_TIMESTAMP")

        # ====# ISOFORMS MOD #=====#

        dbase.multi_write("INSERT INTO isoform_modifications(isoform_id, mod_id)", var_mod_values,
                          update="ON DUPLICATE KEY UPDATE Date=CURRENT_TIMESTAMP")
        if args.verbose:
            print("[DATABASE]: Isoforms tables populated/updated")

    # ========# FILLING THE PDB RELATED TABLES #=========#

    if target.pdb:
        values_PDB = []
        values_PDBChains = []
        values_PDBChains_iso = []
        values_PDBChains_equal = []
        values_PDBChains_dom = []
        for key in target.pdb.keys():
            chain_done = []
            values_PDB.append((target.pdb[key]['PDB_code'], target.pdb[key]['Technique'],
                               target.pdb[key]['Resolution']))
            for chain in target.pdb[key]['sequences'].keys():
                start_stop = ''
                for i in target.pdb[key]['sequences'][chain]['start_stop_pairs']:
                    start_stop += str(i[0]) + '-' + str(i[1]) + ' | '
                start_stop = start_stop.rstrip(' | ')
                values_PDBChains.append((str(key) + '_' + str(chain), key, chain,
                                         target.pdb[key]['sequences'][chain]['sequence'],
                                         target.swissprotID,
                                         target.pdb[key]['sequences'][chain]['start'],
                                         target.pdb[key]['sequences'][chain]['stop'],
                                         target.pdb[key]['sequences'][chain]['length'], start_stop))
                if target.pdb[key]['sequences'][chain]['equal']:
                    for eq in target.pdb[key]['sequences'][chain]['equal']:
                        if eq not in chain_done:
                            values_PDBChains_equal.append((str(key) + '_' + str(chain), str(key) + '_' + str(eq)))
                    chain_done.append(chain)
                if target.pdb[key]['sequences'][chain]['aligned']:
                    for iso in target.pdb[key]['sequences'][chain]['aligned']:
                        values_PDBChains_iso.append((str(key) + '_' + str(chain), iso['isoform'], iso['identity']))
                if 'no_domain' not in target.pdb[key]['sequences'][chain]['domain']:
                    for dom_id in target.pdb[key]['sequences'][chain]['domain_id']:
                        values_PDBChains_dom.append((str(key) + '_' + str(chain), dom_id))

        dbase.multi_write("INSERT INTO PDB(PDB_code, Technique, Resolution)", values_PDB,
                          update="ON DUPLICATE KEY UPDATE Technique=VALUES(Technique),Resolution=VALUES(Resolution),Date=CURRENT_TIMESTAMP")
        dbase.multi_write(
            "INSERT INTO PDB_Chains(Chain_id, PDB_code, Chain, Sequence, Target_id, Start, Stop, n_residues,start_stop)",
            values_PDBChains, update="ON DUPLICATE KEY UPDATE Date=CURRENT_TIMESTAMP")
        dbase.multi_write("INSERT INTO Chain_equality(Chain_id_A, Chain_id_B)", values_PDBChains_equal)
        dbase.multi_write("INSERT INTO PDBChain_Isoforms(Chain_id, Isoform_id, similarity)", values_PDBChains_iso)
        dbase.multi_write("INSERT INTO PDBChain_Domain(Chain_id, Domain_id)", values_PDBChains_dom)
        if args.verbose:
            print("[DATABASE]: PDB tables populated/updated")

    # ========# FILLING THE fPOCKET RELATED TABLES #=========#

    if target.pockets:
        values_pockets = []
        values_pockets_Domain = []
        values_pockets_Chain = []
        for key in target.pockets.keys():
            for p in target.pockets[key].keys():
                pocket_id = target.swissprotID + '_' + key + '_' + target.pockets[key][p].pocket_number
                values_pockets.append((key, target.swissprotID, target.pockets[key][p].pocket_number,
                                       pocket_id,
                                       target.pockets[key][p].score,
                                       target.pockets[key][p].druggability_score,
                                       target.pockets[key][p].apolar_sasa, target.pockets[key][p].polar_sasa,
                                       target.pockets[key][p].total_sasa, target.pockets[key][p].volume,
                                       target.pockets[key][p].druggable))
                if target.pockets[key][p].chain_coverage:
                    for c in target.pockets[key][p].chain_coverage.keys():
                        values_pockets_Chain.append(
                            (pocket_id, key + '_' + c, str(target.pockets[key][p].chain_coverage[c])))
                if target.pockets[key][p].part_of_domain:
                    for d in target.pockets[key][p].part_of_domain:
                        values_pockets_Domain.append((pocket_id, d['domain_id'], d['coverage']))

        dbase.multi_write(
            "INSERT INTO fPockets(PDB_code, Target_id, Pocket_number, Pocket_id, Score, DrugScore, apolar_sasa, "
            "polar_sasa, total_sasa, volume,druggable)",
            values_pockets, update="ON DUPLICATE KEY UPDATE last_updated=CURRENT_TIMESTAMP,Score=VALUES(Score),"
                                   "DrugScore=VALUES(DrugScore),apolar_sasa=VALUES(apolar_sasa),polar_sasa=VALUES("
                                   "polar_sasa),total_sasa=VALUES(total_sasa),volume=VALUES(volume),druggable=VALUES(druggable)")
        dbase.multi_write("INSERT INTO fPockets_Chain(Pocket_id, Chain_id, List_of_contacts)", values_pockets_Chain)
        dbase.multi_write("INSERT INTO fPockets_Domain(Pocket_id, Domain_id, Coverage)", values_pockets_Domain)
        if args.verbose:
            print("[DATABASE]: fPockets tables populated/updated")

    if target.alternate_pockets:
        values_altpock = []
        for key in target.alternate_pockets.keys():
            for p in target.alternate_pockets[key].keys():
                values_altpock.append((key, target.swissprotID, target.alternate_pockets[key][p].pocket_number,
                                       target.swissprotID + '_' + key + '_' + target.alternate_pockets[key][
                                           p].pocket_number,
                                       target.alternate_pockets[key][p].score,
                                       target.alternate_pockets[key][p].druggability_score,
                                       target.alternate_pockets[key][p].apolar_sasa,
                                       target.alternate_pockets[key][p].polar_sasa,
                                       target.alternate_pockets[key][p].total_sasa,
                                       target.alternate_pockets[key][p].volume, 'TRUE',
                                       target.alternate_pockets[key][p].druggable))
        dbase.multi_write(
            "INSERT INTO fPockets(PDB_code, Target_id, Pocket_number, Pocket_id, Score, DrugScore, apolar_sasa,"
            " polar_sasa, total_sasa, volume,blast,druggable)", values_altpock, update="ON DUPLICATE KEY UPDATE "
                                                                                       "last_updated=CURRENT_TIMESTAMP,"
                                                                                       "Score=VALUES(Score),DrugScore=VALUES("
                                                                                       "DrugScore),apolar_sasa=VALUES("
                                                                                       "apolar_sasa),polar_sasa=VALUES("
                                                                                       "polar_sasa),total_sasa=VALUES("
                                                                                       "total_sasa),volume=VALUES(volume),druggable=VALUES(druggable)")
        if args.verbose:
            print("[DATABASE]: alternate fPockets tables populated/updated")

    # ========# FILLING THE ASSAY_TARGET TABLE #=========#

    if target.assay:
        assay_target_values = []

        for assay in target.assay:
            assay_target_values.append((assay['assay_id'], target.swissprotID))

        dbase.multi_write("INSERT INTO assay_target(assay_id, target_id)", assay_target_values)
        if args.verbose:
            print("[DATABASE]: Assay table populated/updated")

    # ========# FILLING THE BIOACTIVITIES TABLE #=========#

    if target.bioactivities:
        bioactivity_values = []
        for entry in target.bioactivities:
            bioactivity_values.append((entry['lig_id'], entry['Target_id'], entry['assay_id'], entry['units'],
                                       entry['operator'], entry['value_num'], entry['doi'], entry['value_text'],
                                       entry['chembl_version'], entry['standard_type'], entry['activity_comment'],
                                       entry['data_validity_comment'], entry['pchembl_value'], entry['target_type'],
                                       entry['target_name'], entry['target_organism']))

        dbase.multi_write(
            "INSERT INTO bioactivities(lig_id, Target_id, assay_id, units, operator, value_num, doi, value_txt, "
            "chembl_version, standard_type, activity_comment, data_validity_comment,pchembl_value,target_type,target_name,target_organism)",
            bioactivity_values,
            update='ON DUPLICATE KEY UPDATE date=CURRENT_TIMESTAMP,chembl_version=VALUES(chembl_version)')
        if args.verbose:
            print("[DATABASE]: Bioactivities table populated/updated")

    # ========# FILLING THE PROTEIN EXPRESSION TABLE #=========#

    if target.protein_expression:
        protein_level_values = []
        protein_selectivity_values = [(target.swissprotID, float(target.protein_expression.selective_entropy),
                                       target.protein_expression.max_organ)]
        for i in target.protein_expression.protein_lvl['Cells']:
            protein_level_values.append((target.swissprotID, i['organ'], i['tissue'], i['cell_type'], i['level'],
                                         target.swissprotID + '_' + i['organ'] + '_' + i['tissue'] + '_' + i[
                                             'cell_type']))

        dbase.multi_write("INSERT INTO protein_expression_levels(Target_id, organ, tissue, cell, value,Entry_id)",
                          protein_level_values,
                          update='ON DUPLICATE KEY UPDATE date=CURRENT_TIMESTAMP,value=VALUES(value)')
        dbase.multi_write("INSERT INTO protein_expression_selectivity(Target_id, Selectivity_entropy, max_organ)",
                          protein_selectivity_values,
                          update='ON DUPLICATE KEY UPDATE date=CURRENT_TIMESTAMP,Selectivity_entropy=VALUES(Selectivity_entropy),max_organ=VALUES(max_organ)')
        if args.verbose:
            print("[DATABASE]: Protein expression levels tables populated/updated")

    # ========# FILLING THE HUMAN MINE DATA TABLES #=========#

    if target.disease:
        value = []
        for d in target.disease:
            value.append(
                (target.swissprotID, d['disease'], d['disease_id'], target.swissprotID + '_' + d['disease_id']))
        dbase.multi_write("INSERT INTO disease(Target_id, disease_name, disease_id, Unique_id)", value,
                          update="ON DUPLICATE KEY UPDATE disease_name=VALUES(disease_name)")
        if args.verbose:
            print("[DATABASE]: Disease table populated/updated")

    if target.differential_exp_tissues:
        value = []
        for exp in target.differential_exp_tissues:
            value.append(
                (target.swissprotID, exp['T_statistic'], exp['Tissue'], exp['expression_status'], exp['p_value']))

        dbase.multi_write("INSERT INTO diff_exp_tissue(Target_id, t_stat, Tissue, expression_status, p_value)", value)
        if args.verbose:
            print("[DATABASE]: Differential expression (tissues) table populated/updated")

    if target.differential_exp_disease:
        value = []
        for exp in target.differential_exp_disease:
            value.append(
                (target.swissprotID, exp['Condition'], exp['T_statistic'], exp['expression_status'], exp['p_value']))

        dbase.multi_write("INSERT INTO diff_exp_disease(Target_id, disease, t_stat, expression_status, p_value)  ",
                          value)
        if args.verbose:
            print("[DATABASE]: Differential expression (disease) table populated/updated")

    if target.gwas:
        value = []
        for g in target.gwas:
            value.append((target.swissprotID, g['doi'], g['first_author'], g['organism'], g['p_value'], g['phenotype'],
                          g['publication_year'], g['pubmed_id']))
        dbase.multi_write(
            "INSERT INTO gwas(Target_id, doi, first_author, organism, p_value, phenotype, publication_year, pubmed_id)",
            value)
        if args.verbose:
            print("[DATABASE]: GWAS table populated/updated")

    if target.phenotypes:
        value = []
        for p in target.phenotypes:
            value.append((target.swissprotID, p['Allele_id'], p['Phenotype'], p['Phenotype_desc'], p['genotype'],
                          p['organism'], p['zygosity'], p['Allele_symbol'], p['Allele_type']))

        dbase.multi_write(
            "INSERT INTO phenotype(Target_id, Allele_id, Phenotype, Phenotype_desc, genotype, organism, zygosity,Allele_symbol,Allele_type)",
            value)
        if args.verbose:
            print("[DATABASE]: phenotype table populated/updated")

    if target.pathways:
        value = []
        for p in target.pathways:
            value.append((target.swissprotID, p['pathway_dataset'], p['pathway_name']))
        dbase.multi_write("INSERT INTO pathways(Target_id, pathway_dataset, pathway_name) ", value)
        if args.verbose:
            print("[DATABASE]: Pathway table populated/updated")

    # ========# FILLING THE CROSSREF TABLE #=========#
    # Only ChemblID at the moment

    if target.chembl_id:
        value = []
        value.append(
            (target.swissprotID, target.chembl_id, target.swissprotID + '_' + target.chembl_id))
        dbase.multi_write("INSERT INTO Crossref(target_id,Chembl_id,unique_id)", value,
                          update="ON DUPLICATE KEY UPDATE date=CURRENT_TIMESTAMP")
        if args.verbose:
            print("[DATABASE]: Cross-references table populated/updated")

    dbase.close()


def get_single_excel(target_id):
    if target_id in list_of_entries.values():
        writer = pd.ExcelWriter(output_single_path + target_id + '_' + list_of_entries[target_id] + '.xlsx',
                                engine='xlsxwriter')

        workbook = writer.book

        # ============================ STYLES ===============================#

        bold_center = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center'})
        red = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'red'})
        green = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'green'})
        col_header = workbook.add_format({'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'vcenter'})

        col_header_greenbg = workbook.add_format(
            {'bold': True, 'bg_color': '#98F5A4', 'align': 'center', 'valign': 'vcenter'})
        col_header_orangebg = workbook.add_format(
            {'bold': True, 'bg_color': '#F5D287', 'align': 'center', 'valign': 'vcenter'})
        col_header_redbg = workbook.add_format(
            {'bold': True, 'bg_color': '#FF9E9E', 'align': 'center', 'valign': 'vcenter'})
        wrap = workbook.add_format({'text_wrap': True, 'valign': 'vcenter'})
        v_center = workbook.add_format({'valign': 'vcenter'})
        link = workbook.add_format(
            {'bold': True, 'valign': 'vcenter', 'align': 'center', 'color': 'blue', 'underline': True})

        wb_general_info = workbook.add_worksheet('General info')
        writer.sheets['General info'] = wb_general_info
        wb_disease = workbook.add_worksheet('diseases')
        writer.sheets['diseases'] = wb_disease
        wb_expression = workbook.add_worksheet('expression')
        writer.sheets['expression'] = wb_expression
        wb_genotypes = workbook.add_worksheet('genotypes')
        writer.sheets['genotypes'] = wb_genotypes
        wb_isoforms = workbook.add_worksheet('isoforms')
        writer.sheets['isoforms'] = wb_isoforms
        wb_var_mut = workbook.add_worksheet('variants_mutants')
        writer.sheets['variants_mutants'] = wb_var_mut
        wb_struct = workbook.add_worksheet('Structure')
        writer.sheets['Structure'] = wb_struct

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query = "SELECT * FROM Targets WHERE Target_id='" + target_id + "'"
        query_disease = "SELECT disease_name,disease_id FROM disease WHERE Target_id='" + target_id + "'"
        query_reactome = "SELECT pathway_name FROM pathways WHERE pathway_dataset='Reactome pathways data set' AND Target_id='" + target_id + "'"
        query_kegg = "SELECT pathway_name FROM pathways WHERE pathway_dataset='KEGG pathways data set' AND Target_id='" + target_id + "'"
        res_gen_info = dbase.get(query)
        res_disease = dbase.get(query_disease)
        res_reactome = dbase.get(query_reactome)
        res_kegg = dbase.get(query_kegg)
        dbase.close()

        # GENERAL INFO HEADER
        header_index = {'Gene_name': (0, 0), 'Synonyms': (1, 0), 'Target_id': (2, 0), 'Protein_class': (3, 0),
                        'Protein_class_desc': (4, 0), 'Species': (5, 0), 'Number_isoforms': (6, 0),
                        'Sequence': (0, 3), 'Cell_location': (0, 4), 'DISEASE': (8, 0, 8, 1), 'disease_id': (9, 0),
                        'disease_name': (9, 1), 'PATHWAYS': (8, 3, 8, 4), 'Reactome': (9, 3), 'KEGG': (9, 4)}
        for head in header_index.keys():
            if len(header_index[head]) == 2:
                row, col = header_index[head]
                wb_general_info.write(row, col, head, col_header)
            elif len(header_index[head]) == 4:
                row, col, last_row, last_col = header_index[head]
                wb_general_info.merge_range(row, col, last_row, last_col, head, col_header)
        sequence = ''
        if res_gen_info:
            for k, v in res_gen_info[0].items():
                if k in header_index:
                    row, col = header_index[k]
                    if k == 'Sequence' or k == 'Cell_location':
                        if k == 'Sequence':
                            sequence = v
                        row = row + 1
                        last_row = row + 5
                        wb_general_info.merge_range(row, col, last_row, col, v, wrap)
                    elif k == 'Number_isoforms':
                        col = col + 1
                        wb_general_info.write_url(row, col, "internal:isoforms", link)
                        wb_general_info.write(row, col, v, link)
                    else:
                        col = col + 1
                        wb_general_info.write(row, col, v, wrap)
        if res_disease:
            for i in range(len(res_disease)):
                for k, v in res_disease[i].items():
                    row, col = header_index[k]
                    row = row + i + 1
                    if k == 'disease_id':
                        wb_general_info.write_url(row, col, 'https://omim.org/entry/' + v.split(':')[1], link)
                        wb_general_info.write(row, col, v, link)
                    else:
                        wb_general_info.write(row, col, v, v_center)
        if res_reactome:
            row, col = header_index['Reactome']
            row = row + 1
            reactome = [i['pathway_name'] for i in res_reactome]
            wb_general_info.write_column(row, col, reactome, v_center)
        if res_kegg:
            row, col = header_index['KEGG']
            row = row + 1
            kegg = [i['pathway_name'] for i in res_kegg]
            wb_general_info.write_column(row, col, kegg, v_center)

        # DISEASE HEADER
        dis_header_index = {'DISEASE REGULATION': (0, 0, 0, 4), 'GWAS': (0, 6, 0, 11), 'disease': (1, 0),
                            't_stat': (1, 1),
                            'std_dev_t': (1, 2), 'n': (1, 3), 'direction': (1, 4), 'phenotype': (1, 6),
                            'organism': (1, 7), 'author': (1, 8), 'year': (1, 9), 'p_value': (1, 10),
                            'pubmed_id': (1, 11)}

        for head in dis_header_index.keys():
            if len(dis_header_index[head]) == 2:
                row, col = dis_header_index[head]
                wb_disease.write(row, col, head, col_header)
            elif len(dis_header_index[head]) == 4:
                row, col, last_row, last_col = dis_header_index[head]
                wb_disease.merge_range(row, col, last_row, last_col, head, col_header)
        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_disease = """SELECT
  disease,
  round(avg(t_stat),1) as t_stat,
  round(stddev(t_stat),1) as std_dev_t,
  count(t_stat) as n,
  max(expression_status) as direction
  FROM diff_exp_disease
  WHERE Target_id='%s'
  GROUP BY Target_id,disease
  ORDER BY t_stat DESC""" % target_id
        query_gwas = """SELECT
  phenotype,
  organism,
  p_value,
  first_author as author,
  publication_year as 'year',
  pubmed_id
FROM gwas
WHERE Target_id='%s'
ORDER BY phenotype""" % target_id
        res_disease_exp = dbase.get(query_disease)
        res_gwas = dbase.get(query_gwas)
        dbase.close()
        if res_disease_exp:
            for i in range(len(res_disease_exp)):
                for k, v in res_disease_exp[i].items():
                    row, col = dis_header_index[k]
                    row = row + i + 1
                    wb_disease.write(row, col, v)
        wb_disease.conditional_format(1, 1, row, 1, {'type': 'data_bar'})
        wb_disease.conditional_format(1, 2, row, 2,
                                      {'type': 'icon_set', 'reverse_icons': True, 'icon_style': '3_traffic_lights'})
        if res_gwas:
            for i in range(len(res_gwas)):
                for k, v in res_gwas[i].items():
                    row, col = dis_header_index[k]
                    row = row + i + 1
                    if k == 'pubmed_id':
                        wb_disease.write_url(row, col, 'https://www.ncbi.nlm.nih.gov/pubmed/' + v, link)
                        wb_disease.write(row, col, v, link)
                    else:
                        wb_disease.write(row, col, v)

        # EXPRESSION HEADER

        expression_header_index = {'Tissue Expression': (0, 0, 0, 3), 'Tissue': (1, 0), 't_stat': (1, 1),
                                   'std_dev_t': (1, 2), 'n': (1, 3), 'Selectivity': (0, 5, 0, 6),
                                   'ORGANS': (1, 5, 1, 8), 'organ_name': (2, 5), 'Total_value': (2, 6),
                                   'n_tissues': (2, 7), 'avg_value': (2, 8)}
        for head in expression_header_index.keys():
            if len(expression_header_index[head]) == 2:
                row, col = expression_header_index[head]
                wb_expression.write(row, col, head, col_header)
            elif len(expression_header_index[head]) == 4:
                row, col, last_row, last_col = expression_header_index[head]
                wb_expression.merge_range(row, col, last_row, last_col, head, col_header)

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_tissue = """SELECT
  Tissue,
  round(avg(t_stat),1) as t_stat,
  round(stddev(t_stat),1) as std_dev_t,
  count(t_stat) as n
  FROM diff_exp_tissue
  WHERE Target_id='%s'
  GROUP BY Tissue
  ORDER BY t_stat DESC""" % target_id
        query_selectivity = """SELECT
Selectivity_entropy
FROM protein_expression_selectivity
WHERE Target_id='%s'""" % target_id
        query_organ_expression = """SELECT
  organ as organ_name,
  sum(value) as Total_value,
  count(value)as n_tissues,
  avg(value) as avg_value
  FROM protein_expression_levels
  WHERE Target_id='%s'
  GROUP BY organ
  ORDER BY avg_value DESC""" % target_id
        query_tissue_expression = """SELECT
  organ,
  tissue,
  cell,
  value
  FROM protein_expression_levels
  WHERE Target_id='%s'""" % target_id
        res_tissue = dbase.get(query_tissue)
        res_selectivity = dbase.get(query_selectivity)
        res_organ_exp = dbase.get(query_organ_expression)
        res_tissue_exp = dbase.get(query_tissue_expression)
        dbase.close()

        if res_tissue:
            for i in range(len(res_tissue)):
                for k, v in res_tissue[i].items():
                    row, col = expression_header_index[k]
                    row = row + i + 1
                    wb_expression.write(row, col, v)
            wb_expression.conditional_format(1, 1, row, 1, {'type': 'data_bar'})
            wb_expression.conditional_format(1, 2, row, 2, {'type': 'icon_set', 'reverse_icons': True,
                                                            'icon_style': '3_traffic_lights'})
        if res_selectivity:
            wb_expression.merge_range(0, 7, 0, 8, res_selectivity[0]['Selectivity_entropy'], col_header)
        if res_organ_exp:
            for i in range(len(res_organ_exp)):
                for k, v in res_organ_exp[i].items():
                    row, col = expression_header_index[k]
                    row = row + i + 1
                    wb_expression.write(row, col, v)
            wb_expression.conditional_format(3, 8, row, 8, {'type': 'data_bar'})

            organ_chart = workbook.add_chart({'type': 'bar'})
            organ_chart.add_series({'values': '=expression!$I$4:$I$16',
                                    'categories': '=expression!$F$4:$F$16',
                                    'name': 'Organ Expression'})
            organ_chart.set_legend({'none': True})
            organ_chart.set_x_axis({'min': 0, 'max': 3, 'major_unit': 1, 'minor_unit_type': 'level',
                                    'major_gridlines': {'visible': True, 'line': {'width': 1.25, 'dash_type': 'dash'}}})
            wb_general_info.insert_chart('G1', organ_chart)

        if res_tissue_exp:
            previous_organ = ''
            row = 0
            col = 10
            organ_count = 0
            for i in range(len(res_tissue_exp)):
                if res_tissue_exp[i]['organ'] != previous_organ:
                    if row >= 65:
                        col += 5
                        row = 0
                    if row == 0:
                        pass
                    else:
                        row += 2
                    wb_expression.merge_range(row, col, row, col + 3, res_tissue_exp[i]['organ'].upper(), col_header)
                    row += 1
                    wb_expression.write(row, col, 'tissue name', col_header)
                    wb_expression.merge_range(row, col + 1, row, col + 2, 'Cell type', col_header)
                    wb_expression.write(row, col + 3, 'Value', col_header)
                    previous_organ = res_tissue_exp[i]['organ']
                    organ_count += 1
                row += 1
                wb_expression.write(row, col, res_tissue_exp[i]['tissue'])
                wb_expression.merge_range(row, col + 1, row, col + 2, res_tissue_exp[i]['cell'])
                wb_expression.write(row, col + 3, res_tissue_exp[i]['value'])
            brain_chart = workbook.add_chart({'type': 'bar'})
            brain_chart.add_series({'values': '=expression!$N$3:$N$13',
                                    'categories': '=expression!$K$3:$L$13',
                                    'name': 'Brain Expression'})
            brain_chart.set_legend({'none': True})
            brain_chart.set_x_axis({'min': 0, 'max': 3, 'major_unit': 1, 'minor_unit_type': 'level',
                                    'major_gridlines': {'visible': True, 'line': {'width': 1.25, 'dash_type': 'dash'}}})
            wb_general_info.insert_chart('G17', brain_chart)

        # GENOTYPE SECTION

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_phenotype = """SELECT
  Allele_symbol,
  Allele_type,
  CASE WHEN zygosity is null THEN 'NOT DECLARED' ELSE UPPER(zygosity) END AS zygosity,
  genotype,
  Phenotype
FROM phenotype WHERE Target_id='%s'
ORDER BY Allele_id,zygosity,genotype""" % target_id
        res_phenotype = dbase.get(query_phenotype)
        dbase.close()
        if res_phenotype:
            from itertools import groupby
            row = 0
            col_allele = 0
            col_zyg = 0
            col_gen = 0
            row_data = 4
            row_with_phen = []
            lethal_phen = re.compile('(lethal)', re.IGNORECASE)
            normal_phen = re.compile('(no abnormal phenotype detected)', re.IGNORECASE)
            for allele, data in groupby(res_phenotype, key=itemgetter('Allele_symbol')):
                tmp_row_with_phen = []
                for zygosity, d2 in groupby(data, key=itemgetter('zygosity')):
                    for genotype, d3 in groupby(d2, key=itemgetter('genotype')):
                        lethal = False
                        normal = False
                        for phen in d3:
                            if lethal_phen.search(phen['Phenotype']):
                                wb_genotypes.write(row_data, col_gen, phen['Phenotype'], red)
                                lethal = True
                            elif normal_phen.search(phen['Phenotype']):
                                wb_genotypes.write(row_data, col_gen, phen['Phenotype'], green)
                                normal = True
                            else:
                                wb_genotypes.write(row_data, col_gen, phen['Phenotype'])
                            row_data += 1
                            allele_type = phen['Allele_type']
                        tmp_row_with_phen.append(row_data)
                        row_data = row + 4
                        if lethal and normal:
                            wb_genotypes.write(row + 3, col_gen, genotype, col_header_orangebg)
                        elif lethal:
                            wb_genotypes.write(row + 3, col_gen, genotype, col_header_redbg)
                        elif normal:
                            wb_genotypes.write(row + 3, col_gen, genotype, col_header_greenbg)
                        else:
                            wb_genotypes.write(row + 3, col_gen, genotype, col_header)

                        col_gen += 1
                    if col_gen - col_zyg == 1:
                        wb_genotypes.write(row + 2, col_zyg, zygosity, col_header)
                    else:
                        wb_genotypes.merge_range(row + 2, col_zyg, row + 2, col_gen - 1, zygosity, col_header)
                    col_zyg = col_gen
                if col_zyg - col_allele == 1:
                    wb_genotypes.write(row + 1, col_allele, allele_type, col_header)
                    wb_genotypes.write(row, col_allele, allele, col_header)
                else:
                    wb_genotypes.merge_range(row + 1, col_allele, row + 1, col_zyg - 1, allele_type, col_header)
                    wb_genotypes.merge_range(row, col_allele, row, col_zyg - 1, allele, col_header)
                col_allele = 0
                col_zyg = 0
                col_gen = 0
                max_row_data = max(tmp_row_with_phen)
                row_with_phen.append((row + 4, max_row_data))
                row = max_row_data
                row += 1
                row_data = row + 4
            for row1, row2 in row_with_phen:
                for i in range(row1, row2):
                    wb_genotypes.set_row(i, None, None, {'level': 1, 'collapsed': True, 'hidden': True})

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_isoforms = """SELECT
  CONCAT(T.Gene_name,'-',I.Isoform_name) as isoform_name,
  I.Isoform_id,
  I.Sequence,
  I.n_residues,
  CASE WHEN I.Canonical = 1 THEN 'Yes' ELSE 'No' END AS is_canonical,
  I.Identity AS similarity
FROM Isoforms I
LEFT JOIN Targets T
  ON I.Target_id = T.Target_id
WHERE I.Target_id='%s' ORDER BY I.Canonical DESC""" % target_id
        query_isoforms_mod = """SELECT
  IM.isoform_id,
  M.start,
  M.stop,
  M.previous AS previous_seq,
  M.action AS modification_type,
  M.new AS new_seq,
  M.domains AS in_domains,
  M.comment AS comments
FROM isoform_modifications IM
LEFT JOIN modifications M
  on IM.mod_id = M.Unique_modID
WHERE IM.isoform_id in (SELECT I.Isoform_id FROM Isoforms I WHERE I.Target_id='%s')""" % target_id
        res_isoforms = dbase.get(query_isoforms)
        res_isoforms_mod = dbase.get(query_isoforms_mod)
        dbase.close()
        if res_isoforms:
            for iso in res_isoforms:
                iso['mod'] = []
                if res_isoforms_mod:
                    for mod in res_isoforms_mod:
                        if iso['Isoform_id'] == mod['isoform_id']:
                            iso['mod'].append(mod)
            row = 0
            mod_header = {'start': 0, 'stop': 1, 'previous_seq': 2, 'modification_type': 3, 'new_seq': 4,
                          'in_domains': 5, 'comments': 6}
            row_to_hide = []
            for iso in res_isoforms:
                wb_isoforms.merge_range(row, 0, row, 6, iso['isoform_name'], col_header)
                row += 1
                wb_isoforms.write(row, 0, 'Is Canonical', col_header)
                wb_isoforms.write(row, 1, iso['is_canonical'], bold_center)
                wb_isoforms.merge_range(row, 2, row, 3, 'Similarity', col_header)
                wb_isoforms.write(row, 4, iso['similarity'], bold_center)
                wb_isoforms.write(row, 5, 'number of residues', col_header)
                wb_isoforms.write(row, 6, iso['n_residues'], bold_center)
                row += 1
                to_hide_start = row
                wb_isoforms.write(row, 0, 'SEQUENCE', col_header)
                wb_isoforms.merge_range(row, 1, row, 6, iso['Sequence'], wrap)
                wb_isoforms.set_row(row, 50)
                row += 1
                if iso['mod']:
                    for header, col in mod_header.items():
                        wb_isoforms.write(row, col, header, col_header)
                    for mod in iso['mod']:
                        row += 1
                        for key, value in mod.items():
                            if key == 'isoform_id':
                                continue
                            elif key == 'previous_seq' or key == 'new_seq':
                                wb_isoforms.write(row, mod_header[key], value, wrap)
                            else:
                                wb_isoforms.write(row, mod_header[key], value)
                to_hide_stop = row + 1
                row_to_hide.append((to_hide_start, to_hide_stop))
                row += 2
            for row1, row2 in row_to_hide:
                for i in range(row1, row2):
                    wb_isoforms.set_row(i, None, None, {'level': 1, 'collapsed': True, 'hidden': True})

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_var = """SELECT
  M.start,
  M.stop,
  M.previous AS previous_seq,
  M.action AS modification_type,
  M.new AS new_seq,
  M.domains AS in_domains,
  M.comment AS comments
FROM modifications M
WHERE M.mod_type = 'VAR' AND M.Target_id='%s'""" % target_id
        query_mut = """SELECT
  M.start,
  M.stop,
  M.previous AS previous_seq,
  M.action AS modification_type,
  M.new AS new_seq,
  M.domains AS in_domains,
  M.comment AS comments
FROM modifications M
WHERE M.mod_type = 'MUTAGEN' AND M.Target_id='%s'""" % target_id
        res_var = dbase.get(query_var)
        res_mut = dbase.get(query_mut)

        mod_header = {'start': 0, 'stop': 1, 'previous_seq': 3, 'modification_type': 2, 'new_seq': 4, 'in_domains': 5,
                      'comments': 6}
        dbase.close()

        if res_var and res_mut:
            col_var = 0
            col_mut = 0
        elif res_var or res_mut:
            col_var = 0
            col_mut = 0
        if res_var:
            row = 0
            wb_var_mut.merge_range(row, col_var, row, col_var + 6, 'VARIANTS', col_header)
            row += 1
            for header, col in mod_header.items():
                wb_var_mut.write(row, col + col_var, header, col_header)
            for var in res_var:
                row += 1
                for key, value in var.items():
                    if key == 'previous_seq' or key == 'new_seq':
                        wb_var_mut.write(row, col_var + mod_header[key], value, wrap)
                    else:
                        wb_var_mut.write(row, col_var + mod_header[key], value)
        if res_mut:
            if res_var:
                row += 2
            else:
                row = 0
            wb_var_mut.merge_range(row, col_mut, row, col_mut + 6, 'MUTANTS', col_header)
            row += 1
            for header, col in mod_header.items():
                wb_var_mut.write(row, col + col_mut, header, col_header)
            for mut in res_mut:
                row += 1
                for key, value in mut.items():
                    if key == 'previous_seq' or key == 'new_seq':
                        wb_var_mut.write(row, col_mut + mod_header[key], value, wrap)
                    else:
                        wb_var_mut.write(row, col_mut + mod_header[key], value)

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_domains = """SELECT
  Domain_name,
  Domain_start as start,
  Domain_stop as stop,
  length,
  source_name as source
FROM Domain_targets
WHERE Target_id='%s'""" % target_id
        query_pdb_blast = """SELECT
  Hit_PDB_code as PDB_code,
  Chain_Letter as Chain,
  similarity,
  Hit_gene_name as gene,
  Hit_gene_species as species
FROM 3D_Blast
WHERE Query_target_id='%s'
ORDER BY similarity DESC""" % target_id
        query_pdb = """SELECT
  C.PDB_code,
  P.Technique,
  P.Resolution,
  GROUP_CONCAT(DISTINCT C.Chain SEPARATOR ',') AS Chain,
  C.n_residues,
  C.start_stop,
  GROUP_CONCAT(DISTINCT D.Domain_name SEPARATOR ',') AS Domain_name
FROM PDB_Chains C
  LEFT JOIN PDB P
    ON C.PDB_code = P.PDB_code
  LEFT JOIN PDBChain_Domain Domain
    ON C.Chain_id = Domain.Chain_id
  LEFT JOIN Domain_targets D
    ON Domain.Domain_id = D.domain_id
WHERE C.Target_id='%s'
GROUP BY C.PDB_code,P.Technique,P.Resolution,C.n_residues,C.start_stop""" % target_id
        res_domains = dbase.get(query_domains)
        res_pdb_blast = dbase.get(query_pdb_blast)
        res_pdb = dbase.get(query_pdb)
        dbase.close()

        row = 0
        col_orig = 0
        if sequence:
            wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'Total length', col_header)
            row += 1
            wb_struct.merge_range(row, col_orig, row, col_orig + 4, len(sequence), bold_center)
            row += 2
        if res_domains:
            wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'DOMAINS', col_header)
            row += 1
            dom = pd.DataFrame.from_records(res_domains)
            col_order = ['Domain_name', 'start', 'stop', 'length', 'source']
            dom = dom[col_order]
            dom.to_excel(writer, sheet_name='Structure', startrow=row, index=False)
            row += len(dom) + 1

        if res_pdb_blast:
            wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'PDB BLAST', col_header)
            row += 1
            blast = pd.DataFrame.from_records(res_pdb_blast)
            col_order = ['PDB_code', 'Chain', 'similarity', 'gene', 'species']
            blast = blast[col_order]
            blast.to_excel(writer, sheet_name='Structure', startrow=row, index=False)
        col_orig = 6
        if res_pdb:
            row = 0
            wb_struct.merge_range(row, col_orig, row, col_orig + 7, 'PDB', col_header)
            row += 1

            pdb = pd.DataFrame.from_records(res_pdb)
            pdb['% of full protein'] = round((pdb['n_residues'] / len(sequence)) * 100, 1)
            col_order = ['PDB_code', 'Technique', 'Resolution', 'Chain', 'Domain_name',
                         'n_residues', '% of full protein', 'start_stop']
            pdb = pdb[col_order]
            pdb.to_excel(writer, sheet_name='Structure', startrow=row, startcol=col_orig, index=False)

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_pockets = """SELECT
  F.PDB_code,
  F.DrugScore as druggability_score,
  round(F.total_sasa,1) as area,
  round(F.volume,1) as volume,
  round((F.apolar_sasa/F.total_sasa)*100,1) as fraction_apolar,
  F.Pocket_number as pocket_number,
  F.Score as pocket_score,
  GROUP_CONCAT(CONCAT(D.Domain_name,' (',Domain.Coverage,'%)') SEPARATOR ',') as domains
FROM fPockets F
  LEFT JOIN fPockets_Domain Domain
    ON F.Pocket_id = Domain.Pocket_id
  LEFT JOIN Domain_targets D
    ON Domain.Domain_id=D.domain_id
WHERE F.Target_id='{target}'
AND F.druggable='TRUE' AND F.blast='FALSE'
GROUP BY F.PDB_code,F.DrugScore,F.total_sasa,F.volume,fraction_apolar,pocket_number,pocket_score""".format(
            target=target_id)
        query_alt_pockets = """SELECT
  F.PDB_code,
  F.DrugScore as druggability_score,
  round(F.total_sasa,1) as area,
  round(F.volume,1) as volume,
  round((F.apolar_sasa/F.total_sasa)*100,1) as fraction_apolar,
  F.Pocket_number as pocket_number,
  F.Score as pocket_score,
  B.Hit_gene_name as gene,
  B.Hit_gene_species as species,
  B.similarity
FROM fPockets F
  LEFT JOIN `3D_Blast` B
    ON F.Target_id = B.Query_target_id AND F.PDB_code = B.Hit_PDB_code

WHERE F.Target_id='%s'
AND F.druggable='TRUE' AND F.blast='TRUE'
ORDER BY B.similarity DESC""" % target_id
        res_pockets = dbase.get(query_pockets)
        res_alt_pockets = dbase.get(query_alt_pockets)
        dbase.close()

        col_alt_pocket = 0

        if res_pockets:
            col_alt_pocket = 9
            pock = pd.DataFrame.from_records(res_pockets)
            col_order = ['PDB_code', 'druggability_score', 'pocket_score', 'pocket_number',
                         'volume', 'area', 'fraction_apolar', 'domains']
            pock = pock[col_order]
            pock.to_excel(writer, sheet_name='Pockets', startrow=1, index=False)
            wb_pockets = writer.sheets['Pockets']
            wb_pockets.merge_range(0, 0, 0, 7, 'DRUGGABLE POCKETS', col_header)

        if res_alt_pockets:
            alt_pock = pd.DataFrame.from_records(res_alt_pockets)
            col_order = ['PDB_code', 'druggability_score', 'pocket_score', 'pocket_number',
                         'volume', 'area', 'fraction_apolar', 'gene', 'species', 'similarity']
            alt_pock = alt_pock[col_order]
            alt_pock.to_excel(writer, sheet_name='Pockets', startrow=1, index=False, startcol=col_alt_pocket)
            wb_pockets = writer.sheets['Pockets']
            wb_pockets.merge_range(0, 0 + col_alt_pocket, 0, 9 + col_alt_pocket,
                                   'ALTERNATE DRUGGABLE POCKETS (PDB from blast)', col_header)

        dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
        query_bioactives = """SELECT
B.lig_id,
  B.assay_id,
  B.standard_type,
  B.operator,
  B.value_num,
  B.units,
  B.activity_comment,
  B.data_validity_comment,
  B.doi as ref_bio,
  B.pchembl_value as pX,
  L.mol_name,
  L.max_phase,
  L.oral,
  L.indication_class,
  L.class_def,
  L.alogp as aLogP,
  L.acd_logd as LogD,
  L.acd_logp as LogP,
  L.acd_most_apka as apKa,
  L.acd_most_bpka as bpKa,
  L.HBA,
  L.HBD,
  L.TPSA,
  L.molecularWeight as MW,
  L.rotatableBonds as rotB,
  L.n_Ar_rings as nAr,
  L.n_alerts as n_alerts,
  L.molecular_species,
  L.num_ro5_violations as ro5_violations,
  L.ro3_pass as pass_ro3,
  L.canonical_smiles as SMILES,
  A.assay_description,
  A.doi as assay_ref,
  A.species as assay_species,
  A.bioactivity_type,
  A.confidence_score
FROM Crossref C
  LEFT JOIN bioactivities B
  ON C.Chembl_id=B.Target_id
  LEFT JOIN ligands L
  ON B.lig_id=L.lig_id
  LEFT JOIN assays A
  ON B.assay_id=A.assay_id
WHERE C.target_id='%s'
AND B.operator!='>' AND B.operator!='<'
AND A.confidence_score>=8""" % target_id
        res_bio = dbase.get(query_bioactives)
        dbase.close()
        if res_bio:
            conc = re.compile(r'(?:of|at)\s(\d+\.*\d*)\s?((?:u|n)M)')
            bioactivity_types = ['Binding', 'Functionnal']
            percent = ['Activity', 'Residual activity', 'Residual_activity', 'Inhibition']
            percent_invert = ['Activity', 'Residual activity', 'Residual_activity']
            dose_response_type = ['IC50', 'Ki', 'EC50', 'Kd', 'Potency']
            top = ['Bio', 'Bio', 'Bio', 'Bio', 'Bio', 'Bio', 'Bio', 'Bio', 'Bio', 'Assay', 'Assay', 'Assay', 'Assay',
                   'Structure',
                   'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop', 'Prop',
                   'Prop', 'Prop', 'Prop', 'Lig_info', 'Lig_info', 'Lig_info', 'Lig_info', 'Lig_info', 'Lig_info',
                   'Ref', 'Ref', 'ID', 'ID', 'Structure']
            col = ['standard_type', 'operator', 'value_num', 'units', 'pX', 'Conc', 'Conc_units', 'activity_comment',
                   'data_validity_comment',
                   'bioactivity_type', 'assay_species', 'assay_description', 'confidence_score', 'SMILES', 'HBA', 'HBD',
                   'LogD', 'LogP', 'MW', 'TPSA', 'aLogP', 'apKa', 'bpKa', 'nAr', 'n_alerts', 'pass_ro3',
                   'ro5_violations', 'rotB', 'CNS_MPO', 'mol_name', 'molecular_species', 'indication_class',
                   'class_def', 'max_phase', 'oral', 'assay_ref', 'ref_bio', 'assay_id', 'lig_id']
            tup = tuple(zip(top, col))
            multi = pd.MultiIndex.from_tuples(tup)

            bioactives = pd.DataFrame.from_records(res_bio)

            bioactives[['Conc', 'Conc_units']] = bioactives[
                (bioactives['units'] == '%') & (bioactives['standard_type'].isin(percent)) & (
                    bioactives['bioactivity_type'].isin(bioactivity_types))].assay_description.str.extract(conc,
                                                                                                           expand=False)
            bioactives.bpKa = bioactives.bpKa.fillna(0)
            bioactives['CNS_MPO'] = mpo.calc_mpo_score(bpka=bioactives['bpKa'], logP=bioactives['LogP'],
                                                       logD=bioactives['LogD'], MW=bioactives['MW'],
                                                       HBD=bioactives['HBD'], TPSA=bioactives['TPSA'])
            bioactives = bioactives[col]
            bioactives.operator = ' ' + bioactives.operator
            bioactives.columns = multi

            percent_bio = bioactives[
                (bioactives.Bio['units'] == '%') & (bioactives.Bio['standard_type'].isin(percent)) & (
                    bioactives.Assay['bioactivity_type'].isin(bioactivity_types))].copy()

            for key in percent_invert:
                percent_bio.loc[percent_bio.Bio['standard_type'] == key, [('Bio', 'value_num')]] = 100 - \
                                                                                                   percent_bio.Bio[
                                                                                                       'value_num']
                percent_bio.loc[percent_bio.Bio['standard_type'] == key, [('Bio', 'standard_type')]] = '100 - ' + \
                                                                                                       percent_bio.Bio[
                                                                                                           'standard_type']
            percent_bio = percent_bio[(percent_bio.Bio['value_num'] > 50)]
            percent_bio.sort_values(by=[('Bio', 'value_num')], ascending=False, inplace=True)
            efficacy_bio = bioactives[
                (bioactives.Bio['units'] == '%') & (bioactives.Bio['standard_type'] == 'Efficacy')].copy()
            efficacy_bio = efficacy_bio[efficacy_bio.Bio.value_num >= 50]
            efficacy_bio.sort_values(by=[('Bio', 'value_num')], ascending=False, inplace=True)
            emax = bioactives[(bioactives.Bio['units'] == '%') & (bioactives.Bio['standard_type'] == 'Emax') & (
                bioactives.Assay['bioactivity_type'].isin(bioactivity_types))].copy()
            emax = emax[emax.Bio.value_num >= 50]
            emax.sort_values(by=[('Bio', 'value_num')], ascending=False, inplace=True)
            ADME = bioactives[(bioactives.Assay['bioactivity_type'] == 'ADME')].copy()
            ADME.sort_values(by=[('Assay', 'assay_description')], inplace=True)
            other = bioactives[~(bioactives.Bio['standard_type'].isin(
                ['Emax', 'Efficacy', 'Activity', 'Residual activity', 'Residual_activity', 'Inhibition', 'IC50', 'Ki',
                 'EC50', 'Kd', 'Potency'])) & ~(bioactives.Assay['bioactivity_type'] == 'ADME')].copy()
            other.sort_values(by=[('Bio', 'standard_type'), ('Assay', 'assay_description')], inplace=True)
            dose_response = bioactives[
                (bioactives.Bio['units'] == 'nM') & (bioactives.Bio['standard_type'].isin(dose_response_type)) & (
                    bioactives.Assay['bioactivity_type'].isin(bioactivity_types))].copy()
            dose_response = dose_response[dose_response.Bio.value_num <= 1000]
            dose_response.sort_values(by=[('Bio', 'standard_type'), ('Bio', 'value_num')], inplace=True)
            dose_response.loc[:, ('Bio', 'pX')].fillna(-np.log10(dose_response.Bio.value_num / 1000000000),
                                                       inplace=True)
            CNS_MPO_criteria = [{'criteria': '>=', 'type': 'number', 'value': 4.5},
                                {'criteria': '>=', 'type': 'number', 'value': 3.5},
                                {'criteria': '<', 'type': 'number', 'value': 3}]

            if not dose_response.empty:
                dose_response.to_excel(writer, sheet_name='Dose_response')
                w_dr = writer.sheets['Dose_response']
                w_dr.conditional_format('AD1:AD' + (str(len(dose_response) + 3)),
                                        {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                            , 'icons': CNS_MPO_criteria})

            if not percent_bio.empty:
                percent_bio.to_excel(writer, sheet_name='Percent_inhibition')
                w_per = writer.sheets['Percent_inhibition']
                w_per.conditional_format('AD1:AD' + (str(len(percent_bio) + 3)),
                                         {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                             , 'icons': [{'criteria': '>=', 'type': 'number', 'value': 4.5},
                                                         {'criteria': '>=', 'type': 'number', 'value': 3.5},
                                                         {'criteria': '<', 'type': 'number', 'value': 3}]})

            if not efficacy_bio.empty:
                efficacy_bio.to_excel(writer, sheet_name='Emax_Efficacy')
                w_eff = writer.sheets['Emax_Efficacy']
                w_eff.conditional_format('AD1:AD' + (str(len(efficacy_bio) + 3)),
                                         {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                             , 'icons': [{'criteria': '>=', 'type': 'number', 'value': 4.5},
                                                         {'criteria': '>=', 'type': 'number', 'value': 3.5},
                                                         {'criteria': '<', 'type': 'number', 'value': 3}]})
                row_efficacy = len(efficacy_bio) + len(efficacy_bio.columns.levels) + 1
            else:
                row_efficacy = 0
            if not emax.empty:
                emax.to_excel(writer, sheet_name='Emax_Efficacy', startrow=row_efficacy)
                w_eff = writer.sheets['Emax_Efficacy']
                w_eff.conditional_format('AD1:AD' + (str(len(emax) + row_efficacy + 3)),
                                         {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                             , 'icons': CNS_MPO_criteria})

            if not ADME.empty:
                ADME.to_excel(writer, sheet_name='ADME')
                w_adme = writer.sheets['ADME']
                w_adme.conditional_format('AD1:AD' + (str(len(ADME) + 3)),
                                          {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                              , 'icons': CNS_MPO_criteria})

            if not other.empty:
                other.to_excel(writer, sheet_name='Other_bioactivities')
                w_other = writer.sheets['Other_bioactivities']
                w_other.conditional_format('AD1:AD' + (str(len(other) + 3)),
                                           {'type': 'icon_set', 'icon_style': '3_traffic_lights'
                                               , 'icons': CNS_MPO_criteria})

        workbook.close()
    else:
        return print("Gene with ID [", target_id, '] not present in the database. Run the command without the -R flag '
                                                  'to run the analysis on it')


def get_list_excel(list_targets):
    list_to_do = {}
    not_in_db = []
    for g in list_targets:
        for entry in gene_in_db:
            if g.upper() == entry['name'].upper():
                list_to_do[entry['name']] = entry['ID']
                continue
            elif g.upper() in entry['syn']:
                list_to_do[entry['name']] = entry['ID']
                continue
        if g not in list_to_do.keys():
            not_in_db.append(g)
    if not list_to_do:
        return print("No genes that you entered are in the Database")

    query_list = '('
    for k, v in list_to_do.items():
        query_list += "'" + v + "',"
    query_list = query_list.rstrip(',') + ')'
    dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
    query = """SELECT T.Gene_name
      ,T.Target_id as Uniprot_id
      ,T.Species
      ,(CASE WHEN T.Number_isoforms=0 THEN 1 ELSE T.Number_isoforms END) Number_isoforms
      ,T.Protein_class_desc
      ,T.Protein_class_short
      ,T.Synonyms
      ,LENGTH(T.Sequence) as number_of_residues
      ,D.domain
      ,M.MUTANT
      ,V.VARIANT
      ,P.PDB
      ,B.protein_blast
      ,PDB_BLAST.pdb_blast
      ,POCK.pockets
      ,BIO.Number_of_ligands
      ,BIO.Max_phase
      ,A.BIO_TYPE AS Assay_types
      ,G.gwas
      ,ALT_POCK.alt_pockets
      ,D_EXP.disease_up as upregulated_in_disease
      ,D_EXP.disease_down as downregulated_in_disease
      ,T_EXP.tissue_up as overexpressed_in
      ,T_EXP.tissue_down as underexpressed_in
      ,PATH.pathways
      ,PHEN.genotypes
      ,PHEN.lethal_phenotype
      ,PHEN.normal_genotype AS normal_phenotype_for
      ,DIS.disease
      ,PROT_SEL.max_organ
      ,ROUND(PROT_SEL.Selectivity_entropy,3) AS expression_selectivity
      ,PROT_ATLAS.summary AS protein_atlas_expression
FROM Targets T
  LEFT JOIN (SELECT
    D.Target_id,
 GROUP_CONCAT(D.domain SEPARATOR '\n') AS domain
  FROM
(SELECT
   D.Target_id
  ,CONCAT(D.source_name,'\n',GROUP_CONCAT(CONCAT('\t',D.Domain_name,' (',D.Domain_start,'-',D.Domain_stop,')') ORDER BY D.Domain_start SEPARATOR '\n')) as domain
FROM Domain_targets D
    GROUP BY D.Target_id,D.source_name) D
GROUP BY D.Target_id) D
  ON T.Target_id=D.Target_id
  LEFT JOIN (SELECT
  Target_id,
  GROUP_CONCAT(CONCAT('(',start,') ',previous,'-->',new,' comment: ',SUBSTRING_INDEX(comment,'.',1),'; in domains: ',domains) ORDER BY start SEPARATOR '\n') as MUTANT
  FROM modifications
    WHERE mod_type ='MUTAGEN'
GROUP BY Target_id) M
    ON T.Target_id=M.Target_id
    LEFT JOIN (SELECT
  Target_id,
  concat(substring_index(GROUP_CONCAT(CONCAT('(',start,') ',previous,'-->',new,' comment: ',SUBSTRING_INDEX(comment,'.',1),'; in domains: ',domains) ORDER BY start SEPARATOR '\n'),'\n',15),case when count(comment) > 15 THEN  concat('\n+ ',count(comment)-15,' others') ELSE '' END)  as VARIANT
  FROM modifications
    WHERE mod_type = 'VAR'
GROUP BY Target_id) V
    ON T.Target_id=V.Target_id
  LEFT JOIN (SELECT
  T.Target_id,
  concat(substring_index(GROUP_CONCAT(CONCAT(T.PDB_code,': ',T.n_residues,' residues (',T.start_stop,', ',T.P_seq,'%) Chain(s): ',T.Chain,' Domain(s): ',CASE WHEN T.domain is NULL THEN '' ELSE T.domain END,' (',T.Technique,': ',T.Resolution,')') ORDER BY T.P_seq DESC SEPARATOR '\n'),'\n',15),case when count(T.PDB_code) > 15 THEN  concat('\n+ ',count(T.PDB_code)-15,' others') ELSE '' END) AS PDB
  FROM
  (SELECT
  C.Target_id,
  C.PDB_code,
  C.Chain_id,
  C.n_residues,
  C.start_stop,
  ROUND(C.n_residues/LENGTH(T.Sequence)*100) AS "P_seq",
  GROUP_CONCAT(DISTINCT C.Chain ORDER BY C.Chain) AS Chain,
  GROUP_CONCAT(DISTINCT DT.Domain_name ORDER BY DT.Domain_start) AS domain,
  DT.Domain_name,
  P.Technique,
  P.Resolution
  FROM PDB_Chains C
LEFT JOIN PDB P
    ON C.PDB_code=P.PDB_code
LEFT JOIN PDBChain_Domain D
    ON C.Chain_id=D.Chain_id
LEFT JOIN Domain_targets DT
    ON D.Domain_id=DT.domain_id
LEFT JOIN Targets T
    ON C.Target_id = T.Target_id
GROUP BY C.Target_id,C.PDB_code
  )T
    GROUP BY T.Target_id) P
  ON T.Target_id=P.Target_id
  LEFT JOIN (SELECT
  Query_target_id as target_id,
concat(substring_index(GROUP_CONCAT(CONCAT(Hit_gene_name,'_',Hit_gene_species,' (',similarity,'%)') ORDER BY similarity DESC SEPARATOR '\n'),'\n',10),case when count(Hit_gene_name) > 10 THEN  concat('\n+ ',count(Hit_gene_name)-10,' others') ELSE '' END) as protein_blast
  FROM protein_blast
GROUP BY Query_target_id) B
  ON T.Target_id=B.Target_id
  LEFT JOIN (SELECT
  Query_target_id as target_id,
concat(substring_index(GROUP_CONCAT(CONCAT(Hit_PDB_code,' Chain: ',Chain_Letter,' (',Hit_gene_name,'_',Hit_gene_species,' - ',similarity,'%)') ORDER BY similarity DESC SEPARATOR '\n'),'\n',10),case when count(Hit_gene_name) > 10 THEN  concat('\n+ ',count(Hit_gene_name)-10,' others') ELSE '' END) as pdb_blast
  FROM `3D_Blast`
GROUP BY Query_target_id) PDB_BLAST
  ON T.Target_id=PDB_BLAST.Target_id
  LEFT JOIN (SELECT
  P.Target_id,
  GROUP_CONCAT(P.pockets_domains SEPARATOR '\n') pockets
FROM
(SELECT
  POCK.Target_id,
  CONCAT(POCK.Domain_name,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',POCK.PDB_code,': Druggability_score=',POCK.DrugScore,' Volume=',POCK.volume,' Area=',POCK.total_sasa,' (',POCK.Fraction_apolar,'% apolar)(',POCK.Pocket_number,')')ORDER BY POCK.DrugScore DESC SEPARATOR '\n'),'\n',3),case when count(POCK.Pocket_number) > 3 THEN  concat('\n\t+ ',count(POCK.Pocket_number)-3,' others') ELSE '' END)) AS pockets_domains
  FROM
(SELECT T1.*,
  (CASE WHEN T1.Domain_id='other' THEN 'other' WHEN T1.Domain_id is NULL THEN 'other' ELSE D.Domain_name END) Domain_name
  FROM
(SELECT
  FP.Target_id,
  FP.PDB_code,
  FP.Pocket_number,
  FP.Score,
  FP.DrugScore,
  FP.total_sasa,
  ROUND((FP.apolar_sasa/FP.total_sasa)*100,1) AS Fraction_apolar,
  FP.volume,
  FPD.Domain_id
  FROM (SELECT * FROM fPockets WHERE druggable='TRUE'
  AND blast='FALSE')FP
  LEFT JOIN fPockets_Domain FPD
    ON FP.Pocket_id=FPD.Pocket_id
) T1
LEFT JOIN Domain_targets D
    ON T1.Domain_id=D.domain_id) POCK
GROUP BY POCK.Target_id,POCK.Domain_name) P
GROUP BY P.Target_id) POCK
  ON T.Target_id=POCK.Target_id
  LEFT JOIN (SELECT
  ALT_POCK.Target_id,
  GROUP_CONCAT(ALT_POCK.pocket_gene_name SEPARATOR '\n') alt_pockets
FROM
(SELECT POCK.Target_id,
    CONCAT(POCK.Hit_gene_name,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',POCK.PDB_code,'(Similarity=',ROUND(POCK.similarity),'%): Druggability_score=',ROUND(POCK.DrugScore,2),' Volume=',ROUND(POCK.volume,1),' Area=',ROUND(POCK.total_sasa,1),' (',ROUND(POCK.Fraction_apolar),'% apolar)(',POCK.Pocket_number,')')ORDER BY POCK.similarity DESC,POCK.DrugScore DESC SEPARATOR '\n'),'\n',3),case when count(POCK.Pocket_number) > 3 THEN  concat('\n\t+ ',count(POCK.Pocket_number)-3,' others') ELSE '' END)) AS pocket_gene_name
FROM
(SELECT
  FP.Target_id,
  FP.PDB_code,
  FP.Pocket_number,
  FP.Score,
  FP.DrugScore,
  FP.total_sasa,
  ROUND((FP.apolar_sasa/FP.total_sasa)*100,1) AS Fraction_apolar,
  FP.volume,
  3D.Hit_PDB_code,
  3D.Hit_gene_name,
  3D.similarity
  FROM fPockets FP
    LEFT JOIN 3D_Blast 3D
    ON 3D.Query_target_id=FP.Target_id AND 3D.Hit_PDB_code=FP.PDB_code
  WHERE
   FP.druggable='TRUE'
  AND blast='TRUE'
  AND 3D.similarity>=70) POCK
GROUP BY POCK.Target_id,POCK.Hit_gene_name) ALT_POCK
GROUP BY ALT_POCK.Target_id) ALT_POCK
  ON T.Target_id=ALT_POCK.Target_id
  LEFT JOIN (SELECT
  T1.Target_id,
  concat(substring_index(GROUP_CONCAT(CASE WHEN T1.up is null THEN null ELSE CONCAT(T1.up,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END))END ORDER BY T1.t_stat DESC SEPARATOR '\n'),'\n',20),case when count(T1.up) > 20 THEN  concat('\n+ ',count(T1.up)-20,' others') ELSE '' END) AS disease_up
,  concat(substring_index(GROUP_CONCAT(CASE WHEN T1.down is null THEN null ELSE CONCAT(T1.down,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END))END ORDER BY T1.t_stat SEPARATOR '\n'),'\n',20),case when count(T1.down) > 20 THEN  concat('\n+ ',count(T1.down)-20,' others') ELSE '' END) AS disease_down

FROM
    (SELECT *,
       (CASE WHEN t_stat<0 THEN disease END) AS down,
      (CASE WHEN t_stat>=0 THEN disease END) AS up
      FROM
( SELECT
  disease,
  avg(t_stat) as t_stat,
  stddev(t_stat) as std_dev_t,
  count(t_stat) as n_number,
  Target_id
  FROM diff_exp_disease
    WHERE t_stat > 5 or t_stat < -5
  GROUP BY Target_id,disease
  ) T1
    ) T1
GROUP BY T1.Target_id) D_EXP
  ON T.Target_id=D_EXP.Target_id
  LEFT JOIN (SELECT
  T1.Target_id,
  GROUP_CONCAT(CASE WHEN T1.up is null THEN null ELSE CONCAT(T1.up,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END)) END ORDER BY T1.t_stat DESC SEPARATOR '\n') AS tissue_up,
  GROUP_CONCAT(CASE WHEN T1.down is null THEN null ELSE CONCAT(T1.down,' (T-stat= ',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END)) END ORDER BY T1.t_stat SEPARATOR '\n') AS tissue_down
  FROM
( SELECT
  Tissue,
  avg(t_stat) as t_stat,
  stddev(t_stat) as std_dev_t,
  count(t_stat) as n_number,
  Target_id,
  (CASE WHEN t_stat<0 THEN Tissue END) AS down,
  (CASE WHEN t_stat>=0 THEN Tissue END) AS up
FROM diff_exp_tissue
WHERE t_stat > 5 or t_stat < -5
GROUP BY Target_id,Tissue)T1
GROUP BY T1.Target_id) T_EXP
  ON T.Target_id=T_EXP.Target_id
  LEFT JOIN (SELECT
  P.Target_id,
  GROUP_CONCAT(P.pathways SEPARATOR '\n') AS pathways
  FROM
(SELECT
  Target_id,
  CONCAT(pathway_dataset,'\n',GROUP_CONCAT(CONCAT('\t',pathway_name) ORDER BY pathway_name SEPARATOR '\n')) AS pathways
  FROM pathways
    WHERE pathway_dataset='KEGG pathways data set'
GROUP BY Target_id,pathway_dataset) P
GROUP BY P.Target_id) PATH
  ON T.Target_id=PATH.Target_id
  LEFT JOIN (SELECT
  T1.Target_id,
  GROUP_CONCAT(genotype_list ORDER BY T1.zygosity SEPARATOR '\n') AS genotypes,
  GROUP_CONCAT(T1.lethal_phenotype SEPARATOR '\n') lethal_phenotype
  ,GROUP_CONCAT(T1.normal_genotype SEPARATOR '\n') normal_genotype
  FROM
(SELECT
 T1.Target_id,
  T1.zygosity,
  CONCAT(' [',T1.zygosity,']\n',GROUP_CONCAT(CONCAT('\t',T1.normal_genotype) SEPARATOR '\n')) as normal_genotype,
  CONCAT(' [',T1.genotype,']\n',GROUP_CONCAT(DISTINCT CONCAT('\t',T1.lethal_phen) SEPARATOR '\n')) as lethal_phenotype,
  CONCAT(T1.zygosity,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',T1.genotype,' [',T1.organism,']',(CASE WHEN T1.phen_list like 'no abnormal phenotype detected' THEN ' [NORMAL PHENOTYPE]' WHEN T1.phen_list like '%lethal%' THEN ' [LETHAL PHENOTYPE OBSERVED]' ELSE '[P]' END)) ORDER BY T1.genotype SEPARATOR '\n'),'\n',5),case when count(T1.genotype) > 5 THEN  concat('\n\t+ ',count(T1.genotype)-5,' others') ELSE '' END)) as genotype_list
  FROM
(SELECT
  PHEN.Target_id,
  PHEN.genotype,
  (CASE WHEN PHEN.zygosity is NULL THEN 'not declared' ELSE PHEN.zygosity END) zygosity,
  PHEN.organism,
  GROUP_CONCAT(DISTINCT PHEN.Phenotype SEPARATOR ' ; ') AS phen_list,
  GROUP_CONCAT(DISTINCT (CASE WHEN PHEN.Phenotype like '%lethal%' THEN PHEN.Phenotype END) SEPARATOR '\n\t') AS lethal_phen,
  GROUP_CONCAT(DISTINCT (CASE WHEN PHEN.Phenotype like 'no abnormal phenotype detected' THEN PHEN.genotype END) SEPARATOR '\n') AS normal_genotype
 FROM phenotype PHEN
    GROUP BY PHEN.Target_id,PHEN.genotype,PHEN.zygosity)T1
GROUP BY T1.Target_id,T1.zygosity)T1
GROUP BY T1.Target_id) PHEN
  ON T.Target_id=PHEN.Target_id
  LEFT JOIN (SELECT
  Target_id,
  GROUP_CONCAT(CONCAT(disease_name,' [',disease_id,']') ORDER BY disease_name SEPARATOR '\n') AS disease
  FROM disease
GROUP BY Target_id) DIS
  ON T.Target_id=DIS.Target_id
  LEFT JOIN protein_expression_selectivity PROT_SEL
  ON T.Target_id=PROT_SEL.Target_id
  LEFT JOIN (SELECT
  T1.Target_id,
  GROUP_CONCAT(CONCAT(T1.level_graph,(CASE
                       WHEN 15-T1.n_cell = 0 THEN ''
                       WHEN 15-T1.n_cell = 1 THEN '      '
                       WHEN 15-T1.n_cell = 2 THEN '            '
                       WHEN 15-T1.n_cell = 3 THEN '                  '
                       WHEN 15-T1.n_cell = 4 THEN '                        '
                       WHEN 15-T1.n_cell = 5 THEN '                              '
                       WHEN 15-T1.n_cell = 6 THEN '                                    '
                       WHEN 15-T1.n_cell = 7 THEN '                                          '
                       WHEN 15-T1.n_cell = 8 THEN '                                                '
                       WHEN 15-T1.n_cell = 9 THEN '                                                      '
                       WHEN 15-T1.n_cell = 10 THEN '                                                            '
                       WHEN 15-T1.n_cell = 11 THEN '                                                                  '
                       WHEN 15-T1.n_cell = 12 THEN '                                                                        '
                       WHEN 15-T1.n_cell = 13 THEN '                                                                              '
                       WHEN 15-T1.n_cell = 14 THEN '                                                                                    '
                       WHEN 15-T1.n_cell = 15 THEN '                                                                                          '

                       END),'\t',T1.organ,' (',ROUND(T1.avg_level_num,1),')') ORDER BY T1.avg_level_num DESC SEPARATOR '\n') AS summary
  FROM
(SELECT
  T1.Target_id,
  T1.organ,
  GROUP_CONCAT(CONCAT(T1.level_graph) ORDER BY T1.level_numeric DESC SEPARATOR '') level_graph,
  SUM(T1.level_numeric) level_num,
  AVG(T1.level_numeric) avg_level_num,
  COUNT(T1.cell) n_cell
  FROM
(SELECT
    T1.Target_id,
    T1.organ,
    T1.tissue,
    T1.cell,
    CASE WHEN value=0 THEN '[ -]' WHEN value=1 THEN '[1]' WHEN value=2 THEN '[2]' WHEN value=3 THEN '[3]' END AS level_graph,
    T1.value AS level_numeric
  FROM protein_expression_levels T1) T1
GROUP BY T1.Target_id,T1.organ)T1
GROUP BY T1.Target_id) PROT_ATLAS
  ON T.Target_id=PROT_ATLAS.Target_id
  LEFT JOIN (SELECT
    C.target_id,
    COUNT(DISTINCT L.lig_id) AS Number_of_ligands,
    MAX(L.max_phase) AS Max_phase
    FROM bioactivities B
    INNER JOIN (SELECT * FROM Crossref C WHERE C.target_id in """ + query_list + """) C
    on B.Target_id=C.Chembl_id
    LEFT JOIN ligands L
    on B.lig_id=L.lig_id
    GROUP BY C.target_id) BIO
  ON T.Target_id=BIO.target_id
  LEFT JOIN (SELECT
  AT.target_id,
  GROUP_CONCAT(DISTINCT A.bioactivity_type ORDER BY A.bioactivity_type SEPARATOR '\n') AS BIO_TYPE
  FROM assay_target AT
  LEFT JOIN assays A
    ON A.assay_id=AT.assay_id
GROUP BY AT.target_id) A
  ON T.Target_id=A.target_id
  LEFT JOIN(SELECT
   G.Target_id
  ,GROUP_CONCAT(CONCAT('Phenotype: ',G.phenotype,' Organism: ',G.organism,' (',G.first_author,'-',G.publication_year,') doi:',G.doi,' PID:',G.pubmed_id) ORDER BY G.phenotype,G.publication_year DESC SEPARATOR '\n') as gwas
FROM gwas G
    GROUP BY G.Target_id) G
  ON T.Target_id=G.Target_id
WHERE T.Target_id in """ + query_list
    res = dbase.get(query)
    dbase.close()
    header = ['Gene_name', 'Uniprot_id', 'Synonyms', 'Species', 'pathways', 'disease', 'upregulated_in_disease',
              'downregulated_in_disease', 'gwas',
              'genotypes', 'lethal_phenotype', 'normal_phenotype_for', 'overexpressed_in', 'underexpressed_in',
              'protein_atlas_expression', 'max_organ',
              'expression_selectivity', 'Number_isoforms', 'Protein_class_desc', 'Protein_class_short',
              'number_of_residues', 'domain', 'MUTANT', 'VARIANT', 'PDB', 'pdb_blast', 'protein_blast', 'pockets',
              'alt_pockets', 'Number_of_ligands', 'Max_phase', 'Assay_types']
    header_index = {header[i]: i for i in range(len(header))}

    t = time.strftime("%d%b%Y_%H%M%S")
    workbook = xlsxwriter.Workbook(output_lists_path + 'Export_' + str(len(res)) + '_entries_' + t + '.xlsx')
    worksheet = workbook.add_worksheet('Druggability_list')
    bold = workbook.add_format({'bold': True})
    for k, v in header_index.items():
        worksheet.write(0, v, k, bold)
    row_counter = 1
    for entry in res:
        for k, v in entry.items():
            worksheet.write(row_counter, header_index[k], v)
        row_counter += 1
    workbook.close()
    print("[EXPORT]: Excel file: ", '[Export_' + t + '.xlsx]', ' successfully generated')
    return "Success"


def get_target_info(gene, gene_id):
    t = Target(gene, gene_id)
    if t.record is None:
        return 'Fail'
    else:
        return 'Target completed'


class Target:
    def __init__(self, gene, gene_id):
        script_start = time.time()
        # =====# INITIATING THE OBJECT #========#
        if '_' in gene:
            self.gene = gene.split('_')[0]
        else:
            self.gene = gene
        self.swissprotID = gene_id
        self.variants = ''
        self.modifications = ''
        self.assay = ''
        self.pockets = None
        self.record = None
        self.pathways = []
        self.disease = []
        self.gwas = []
        self.differential_exp_tissues = []
        self.differential_exp_disease = []
        self.phenotypes = []
        self.do_info = True
        self.do_lig = True
        self.target_in = False
        self.bioactivities = None
        self.alternate_pockets = None
        self.druggable_pockets = None

        while True:
            # ==============# IF NO UNIPROT ID --> SKIP #====================#
            if self.swissprotID == 'No Match':
                if args.verbose:
                    print('[GENE SKIPPED]: No UniprotID found for ' + self.gene)
                break

            if args.verbose:
                print('[BEGINNING OF GENE]: ' + self.gene + ' (' + str(self.swissprotID) + ')')

            # =============# Checking if entry already in DB #================#

            if self.swissprotID in list_of_entries:

                # =============== IF IN DB --> reconstruct the domain and variants variable for further use ===========#
                if args.verbose:
                    print("[GENE INFO]: Target already present in database, updating informations")
                dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
                query_variants = "SELECT iso.*,group_concat(modi.mod_id SEPARATOR ',') as seq_mod FROM Isoforms iso, " \
                                 "modifications modi,isoform_modifications isomod WHERE iso.Target_id ='%s' AND " \
                                 "iso.Isoform_id=isomod.isoform_id AND isomod.mod_id=modi.Unique_modID GROUP BY " \
                                 "iso.Isoform_id" % self.swissprotID
                query_domains = "SELECT * FROM Domain_targets WHERE Target_id ='%s'" % self.swissprotID
                query_sequence = "SELECT * FROM Targets WHERE Target_id ='%s'" % self.swissprotID
                res_variants = dbase.get(query_variants)
                res_domains = dbase.get(query_domains)
                res_protinfo = dbase.get(query_sequence)
                dbase.close()
                for i in res_variants:
                    i['isoid'] = i['Isoform_id']
                    i['name'] = i['Isoform_name']
                    i['sequence'] = i['Sequence']
                    i['seq_mod'] = i['seq_mod'].split(',')
                    i['seq_list'] = list(i['Sequence'])
                    i['alignment'] = {'gaps': i['Gaps'], 'score': i['Score'], 'identity': i['Identity']}
                for j in res_domains:
                    j['Start'] = j['Domain_start']
                    j['Stop'] = j['Domain_stop']
                    j['name'] = j['Domain_name']
                self.variants = res_variants
                self.domain = res_domains
                self.sequence = res_protinfo[0]['Sequence']
                self.target_in = True
                self.do_info = False

            # ============# COLLECTING THE UNIPROT RECORD #===============#
            self.record = get_uniprot(self.swissprotID)
            if self.record is None:
                if args.verbose:
                    print('[GENE SKIPPED]: No Uniprot data or wrong uniprotID for ' + self.gene)
                break

            # ===========# CREATING THE FILE PATH #==============#
            self.path = dbase_file_path + 'DB_files/'
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            # ===========# GET ALL THE CROSSREFERENCES (source: Uniprot)#==============#

            self.CrossRef, self.pdb, self.go, self.chembl_id = get_crossref_pdb_go_chembl(self.record)

            # ===========# GET PROTEIN EXPRESSION LEVELS (source: ProteinAtlas)#==============#

            if 'OpenTargets' in self.CrossRef.keys():
                self.ensembl_id = self.CrossRef['OpenTargets']
            elif 'Bgee' in self.CrossRef.keys():
                self.ensembl_id = self.CrossRef['Bgee']
            elif 'Ensembl' in self.CrossRef.keys():
                self.ensembl_id = self.CrossRef['Ensembl']
            else:
                self.ensembl_id = None

            self.protein_expression = patlas.ProteinExpression(self.gene, id=self.ensembl_id)
            if self.protein_expression.protein_lvl is None:
                self.protein_expression = None

            # ===========# IF THE ENTRY WASN'T IN THE DB IT WILL RUN #==============#

            if self.do_info:

                # ===========# GET INFO FROM HUMANMINE.ORG (disease, phenotypes, differential_exp_diseases, differential_exp_tissues, gwas,pathways) #=============#

                self.disease, self.phenotypes, self.differential_exp_disease, self.differential_exp_tissues, self.gwas, self.pathways = get_humanmine_data(
                    self.gene)

                # ==========# GET DOMAIN INFORMATION FROM BOTH CHEMBL AND UNIPROT #===========#

                self.domain = get_domains(record=self.record, chembl_id=self.chembl_id)

                # ==========# GET PROTEIN CLASS AND SYNONYMS FROM CHEMBL #===========#

                self.prot_info = get_chembl_info(self.chembl_id)
                if self.prot_info:
                    desc = str(self.prot_info[0]['protein_class_desc']).replace('  ', ' -> ')
                    self.ProteinClass = {'Class_name': self.prot_info[0]['pref_name'], 'Class_description': desc,
                                         'Class_short': self.prot_info[0]['pref_name']}
                    self.synonyms = self.prot_info[0]['Synonym']
                else:
                    self.ProteinClass = {'Class_name': '', 'Class_description': '', 'Class_short': ''}
                    self.synonyms = synonyms(self.record)

                # ==========# GET DISEASE ASSOCIATION (Source: Uniprot)#===========#

                self.comments = get_comments(self.record)

                # ==========# GET SEQUENCE INFORMATION (Source: Uniprot)#===========#

                self.sequence = self.record.sequence
                self.seq_list = list(self.sequence)

                # ==========# GET ISOFORMS INFORMATION (Source: Uniprot) #===========#

                self.variants, self.modifications = get_variants(self.record, self.domain)

            # ============================================================================#
            # ======================# END OF THE INFO SECTION #===========================#
            # ============================================================================#

            # ============================================================================#
            # ======================# START OF THE PDB SECTION #==========================#
            # ============================================================================#

            # ======# RETRIEVING LIST OF PDB ASSOCIATED TO TARGET ALREADY IN DB #=========#

            dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)
            query_pdb = "SELECT PDB_code FROM PDB_Chains WHERE Target_id ='%s' GROUP BY PDB_code" % self.swissprotID
            query_pockets = "SELECT Pocket_id FROM fPockets WHERE Target_id ='%s' AND druggable='TRUE' AND blast='FALSE'" % self.swissprotID
            res_pdb = dbase.get(query_pdb)
            res_pockets = dbase.get(query_pockets)
            dbase.close()
            list_pdb = [i['PDB_code'] for i in res_pdb]
            # =============# REMOVING PDB CODE FROM CrossRef if already done #============#
            if self.pdb:
                [self.pdb.pop(key) for key in list_pdb]

            if self.pdb:
                # ========================# DOWNLOADING THE PDBs #============================#

                get_pdb(self.pdb, self.path)

                # ========================# GET PDB SEQUENCE INFO #===========================#

                get_pdb_seq_info(self.pdb, self.domain, self.variants)

                # =====================# GET POCKETS (source: fpockets) ======================#

                self.pockets = pocket.get_pockets(self.path, sphere_size=args.sphere_size,
                                                  pdb_info=self.pdb, domain=self.domain,
                                                  uniprot_id=self.swissprotID)
                self.druggable_pockets = pocket.get_druggable_pockets(self.pockets)
            pdb_super_list = list_pdb + list(self.pdb.keys())

            # ============================================================================#
            # ======================# END OF THE PDB SECTION #============================#
            # ============================================================================#

            # ============================================================================#
            # ====================# START OF THE LIGAND SECTION #=========================#
            # ============================================================================#

            if self.chembl_id:

                # ===# GET LIST OF ASSAYS ASSOCIATED TO THE TARGET #=========#

                self.assay = get_assays(self.chembl_id)

                # ==============# GET LIST OF LIGAND TO ADD #================#

                ligands_to_do = get_ligands_to_do(self.chembl_id)

                if len(ligands_to_do) != 0:
                    # ========# GET ALL BIOACTIVITIES OF THESE LIGANDS #=========#
                    self.bioactivities = get_bioactivity(ligands_to_do)

            # ============================================================================#
            # =====================# END OF THE LIGAND SECTION #==========================#
            # ============================================================================#

            # ============================================================================#
            # ====================# START OF THE BLAST SECTION #==========================#
            # ============================================================================#

            self.alternate_pdb = pdb_blast(self.sequence, self.path, self.swissprotID, pdb_list=pdb_super_list,
                                           gene=self.gene)
            if self.alternate_pdb:
                very_close_pdb = {key: value for key, value in self.alternate_pdb.items() if
                                  value['similarity'] >= 90.0}
            else:
                very_close_pdb = None

            if not res_pockets and not self.druggable_pockets:
                self.alternate_pockets = pocket.get_pockets(self.path, sphere_size=args.sphere_size, alternate=True,
                                                            alternate_pdb=self.alternate_pdb,
                                                            uniprot_id=self.swissprotID)
            else:
                if very_close_pdb:
                    self.alternate_pockets = pocket.get_pockets(self.path, sphere_size=args.sphere_size, alternate=True,
                                                                alternate_pdb=very_close_pdb,
                                                                uniprot_id=self.swissprotID)

            self.neighbours = proteins_blast(self.sequence, self.swissprotID, self.gene, self.path)

            # ============================================================================#
            # ====================# END OF THE BLAST SECTION #============================#
            # ============================================================================#

            break
        # ============================================================================#
        # =================# END OF DATA GATTERING SECTION #==========================#
        # ============================================================================#

        # ======================# WRITING TO THE DATABASE #===========================#
        write_to_db(self, druggability_db)

        if args.verbose:
            script_stop = time.time()
            print('[END OF GENE]: ' + self.gene + ' (in ' + str(round(script_stop - script_start)) + ' sec.)')
            print('=======================================================================')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene', help='enter a single gene name', metavar='')
    parser.add_argument('-i', '--in_file', help='Name of the input file with a list of genes (.txt - 1 gene per line)',
                        metavar='')
    parser.add_argument('-U', '--Uniprot_ID', help='Enter the uniprot accession number followed by the name separated '
                                                   'by a comma(eg. Q9NZC2,TREM2)', metavar='')
    parser.add_argument('-l', '--list_genes', help='Enter a list of gene name separated by a ","', metavar='')
    parser.add_argument('-s', '--sphere_size', help='enter a value for the probe size of the pocket finder tool ('
                                                    'default = 3.0)', metavar='', type=float, default=3.0)
    parser.add_argument('-v', '--verbose', help="Print information", action='store_true', default=False)
    parser.add_argument('-R', '--report_only', help="Report mode ONLY, will not do any search (by default list "
                                                    "output, make sure to add -rs to get a single output)",
                        action='store_true',
                        default=False)
    parser.add_argument('-rl', '--report_list', help="produce an output file at the end of the analysis with a list "
                                                     "of genes entered", action='store_true',
                        default=False)
    parser.add_argument('-rs', '--report_single', help="produce an output file for each target listed (more detailed "
                                                       "information available) Caution : 1 File created per target",
                        action='store_true',
                        default=False)
    parser.add_argument('-update', '--update', help="Update record if already in database (default: No)",
                        action='store_true',
                        default=False)
    parser.add_argument('-blastcore', '--num_core',
                        help='Enter the value of processor core to use for the blast search (default=8)', metavar='',
                        type=int, default=8)
    parser.add_argument('-db_pwd', '--db_password',
                        help='enter the password of your database', metavar='',
                        type=str)
    parser.add_argument('-db_user', '--db_username',
                        help='enter the username of your database', metavar='',
                        type=str)
    args = parser.parse_args()

    while True:
        if args.in_file:
            if os.path.exists(args.in_file):
                with open(args.in_file, 'r') as gene_list:
                    list_of_genes = gene_list.readlines()
                gene_dict = gene_to_uniprotid(list_of_genes)
                break
            else:
                print('ERROR : file inputed as argument [-i] does not exist')
                sys.exit()
        if args.gene:
            list_of_genes = [args.gene]
            gene_dict = gene_to_uniprotid(list_of_genes)
            break
        elif args.list_genes:
            list_of_genes = args.list_genes.split(',')
            gene_dict = gene_to_uniprotid(list_of_genes)
            break
        elif args.Uniprot_ID:
            try:
                ID_name = args.Uniprot_ID.split(',')
                gene_dict = {ID_name[1]: ID_name[0]}
            except IndexError:
                print(
                    'ERROR : Please enter the uniprot ID followed by the gene name separated by a comma (eg. Q9NZC2,TREM2)')
                sys.exit()
            break
        else:
            print('Please use one of the optional input options :')
            parser.print_help()
            sys.exit()

    list_of_entries, gene_in_db = get_list_entries()

    if args.report_only:
        if args.report_single:
            for gene_id in gene_dict.values():
                get_single_excel(gene_id)
        else:
            get_list_excel(gene_dict.keys())
    else:
        targets_list = {}
        for key, value in gene_dict.items():
            if value in list_of_entries:
                if args.update:
                    status = get_target_info(key, value)
                    targets_list[key] = status
                else:
                    if args.verbose:
                        print('[GENE SKIPPED]: Already present in the database: ' + key)
                        print('=======================================================================')
                    targets_list[key] = 'Already present'

            else:
                status = get_target_info(key, value)
                targets_list[key] = status
        if args.report_list or args.report_single:
            list_of_entries, gene_in_db = get_list_entries()
            if args.report_list:
                export_list = []
                for key in targets_list:
                    if targets_list[key] == 'Target completed' or targets_list[key] == 'Already present':
                        export_list.append(key)
                get_list_excel(export_list)
            if args.report_single:
                for key in targets_list:
                    if targets_list[key] == 'Target completed' or targets_list[key] == 'Already present':
                        get_single_excel(key)
