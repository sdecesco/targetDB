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

import sys, os, mygene, requests, subprocess, re
from Bio import ExPASy, Entrez, Medline
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
from sqlalchemy import create_engine
from druggability3 import cns_mpo as mpo
from druggability3 import drugg_errors
from druggability3 import target_descriptors as td
from druggability3 import target_features as tf
import pkg_resources as pkg

# ===================# SETTING UP PATHS #============================#

# TODO: Path to be stored in a config file in a more generic location (e.g. ~/druggability_data/...)

dbase_file_path = '/data/sdecesco/databases/druggability/'
output_lists_path = '/data/sdecesco/databases/druggability/outputs/lists/'
output_single_path = '/data/sdecesco/databases/druggability/outputs/single_targets/'


if not os.path.exists(dbase_file_path):
	os.makedirs(dbase_file_path)
if not os.path.exists(output_lists_path):
	os.makedirs(output_lists_path)
if not os.path.exists(output_single_path):
	os.makedirs(output_single_path)

# ===================# ChembL DATABASE VERSION IN USE #============================#

chembL_db_version = 'chembl_23'
druggability_db = 'druggability'

# =============================# PATH TO SQLITE DB #============================#

targetDB = pkg.resource_filename(__name__, 'data/TargetDB_v1.db')
chembl_24 = pkg.resource_filename(__name__, 'data/chembl_24.db')
tcrd = pkg.resource_filename(__name__, 'data/tcrd_v5.2.0.db')


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


def retryer_pubmed(max_retries=10, timeout=5):
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
				return pd.DataFrame(data=None)

		return inner

	return wraps


@retryer(max_retries=10, timeout=10)
def gene_to_uniprotid(list_of_gene_name):
	mg = mygene.MyGeneInfo()
	gene_id = {}
	request = mg.querymany(list_of_gene_name, scopes="symbol", fields=['uniprot'], species=9606, as_dataframe=True)
	request.dropna(axis=0, subset=['uniprot'], inplace=True)
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

def gene_to_uniprotid_local(list_of_genes):
	engine = create_engine("mysql://" + args.db_username + ":" + args.db_password + "@localhost")
	connector = engine.connect()
	gene_list = []
	for gene in list_of_genes:
		gene = gene.rstrip('\n')
		gene = gene.rstrip('\r')
		gene_list.append(gene)
	gene_ids = "','".join(gene_list)
	gene_id_query = """SELECT * FROM druggability.hgnc as hgn WHERE hgn.hgnc_id in (SELECT hg.hgnc_id FROM druggability.hgnc as hg WHERE hg.xref_value in ('%s'))""" % gene_ids
	gene_xref = pd.read_sql(gene_id_query, con=connector)
	connector.close()
	output = {}
	for gene in gene_list:
		gene_id = gene_xref[(gene_xref.xref_value == gene) & (gene_xref.xref_name == 'symbol')]['hgnc_id'].values
		if gene_id.size == 0:
			gene_id = gene_xref[(gene_xref.xref_value == gene) & ((gene_xref.xref_name == 'prev_symbol') | (gene_xref.xref_name == 'alias_symbol'))]['hgnc_id'].values
		if gene_id.size == 0:
			output[gene] = "No Match"
			continue
		elif gene_id.size > 1:
			for g_id in gene_id:
				gene_name = gene_xref[(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'symbol')].xref_value.values[0]
				gene_uniprot = gene_xref[(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'uniprot_ids')].xref_value.values
				if gene_uniprot.size > 1:
					count = 0
					for uniprot_id in gene_uniprot:
						name = gene_name + '_' + str(count)
						output[name] = uniprot_id
						count += 1
				elif gene_uniprot.size == 0:
					output[gene] = "No Match"
				else:
					output[gene_name] = gene_uniprot[0]
		else:
			gene_name = gene_xref[(gene_xref.hgnc_id == gene_id[0]) & (gene_xref.xref_name == 'symbol')].xref_value.values[0]
			gene_uniprot = gene_xref[(gene_xref.hgnc_id == gene_id[0]) & (gene_xref.xref_name == 'uniprot_ids')].xref_value.values
			if gene_uniprot.size > 1:
				count = 0
				for uniprot_id in gene_uniprot:
					name = gene_name + '_' + str(count)
					output[name] = uniprot_id
					count += 1
			elif gene_uniprot.size == 0:
				output[gene] = "No Match"
			else:
				output[gene_name] = gene_uniprot[0]
	return output


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
	service = Service("http://www.humanmine.org/humanmine/service")
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

	if len(seq1) > 5000 or len(seq2) > 5000:
		return {'score': 0, 'identity': 'too long', 'gaps': 0}

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
		# chain_done = []
		# for chain in seq.keys():
		#     if seq[chain]['chain_name'] not in pdb_list[keys]['Chain']:
		#         continue
		#     for i in seq[chain]['equal']:
		#         if i in chain_done:
		#             seq[chain]['aligned'] = seq[i]['aligned']
		#             chain_done.append(chain)
		#             break
		#     if seq[chain]['aligned']:
		#         continue
		#     if isoforms:
		#         max_score = 0
		#         isoform_match = []
		#         list_of_align = []
		#         for variants in isoforms:
		#             partial_seq = ''.join(variants['seq_list'][seq[chain]['start'] - 1:seq[chain]['stop'] - 1])
		#             seq_align = align(partial_seq, seq[chain]['sequence'], end_gaps=False)
		#             if seq_align is None:
		#                 continue
		#             if seq_align['score'] > max_score:
		#                 max_score = seq_align['score']
		#             list_of_align.append((seq_align, variants['isoid']))
		#         for item in list_of_align:
		#             if item[0]['score'] == max_score:
		#                 isoform_match.append({'isoform': item[1], 'identity': item[0]['identity']})
		#         seq[chain]['aligned'] = isoform_match
		#         chain_done.append(chain)
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
	if counter == 0:
		pass
	elif counter < 1000:
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


@retryer_pubmed(max_retries=10, timeout=5)
def pubmed_search(gene_name, email, return_number=False, mesh_term=None):
	dict_medline = {"AB": "Abstract", "CI": "Copyright Information", "AD": "Affiliation", "AUID": "Author ID",
					"IRAD": "Investigator Affiliation", "AID": "Article Identifier", "AU": "Author",
					"FAU": "Full Author", "CN": "Corporate Author", "DCOM": "Date Completed", "DA": "Date Created",
					"LR": "Date Last Revised", "DEP": "Date of Electronic Publication", "DP": "Date of Publication",
					"EDAT": "Entrez Date", "GS": "Gene Symbol", "GN": "General Note", "GR": "Grant Number",
					"IR": "Investigator Name", "FIR": "Full Investigator Name", "IS": "ISSN", "IP": "Issue",
					"TA": "Journal Title Abbreviation", "JT": "Journal Title", "LA": "Language",
					"LID": "Location Identifier", "MID": "Manuscript Identifier", "MHDA": "MeSH Date",
					"MH": "MeSH Terms", "JID": "NLM Unique ID", "RF": "Number of References", "OAB": "Other Abstract",
					"OCI": "Other Copyright Information", "OID": "Other ID", "OT": "Other Term",
					"OTO": "Other Term Owner", "OWN": "Owner", "PG": "Pagination", "PS": "Personal Name as Subject",
					"FPS": "Full Personal Name as Subject", "PL": "Place of Publication",
					"PHST": "Publication History Status", "PST": "Publication Status", "PT": "Publication Type",
					"PUBM": "Publishing Model", "PMC": "PubMed Central Identifier", "PMID": "PMID",
					"RN": "Registry Number/EC Number", "NM": "Substance Name", "SI": "Secondary Source ID",
					"SO": "Source", "SFM": "Space Flight Mission", "STAT": "Status", "SB": "Subset", "TI": "Title",
					"TT": "Transliterated Title", "VI": "Volume", "CON": "Comment on", "CIN": "Comment in",
					"EIN": "Erratum in", "EFR": "Erratum for", "CRI": "Corrected and Republished in",
					"CRF": "Corrected and Republished from", "PRIN": "Partial retraction in",
					"PROF": "Partial retraction of", "RPI": "Republished in", "RPF": "Republished from",
					"RIN": "Retraction in", "ROF": "Retraction of", "UIN": "Update in", "UOF": "Update of",
					"SPIN": "Summary for patients in", "ORI": "Original report in"}
	pubmed_url = 'https://www.ncbi.nlm.nih.gov/pubmed/'

	Entrez.email = email
	if mesh_term:
		mesh = '"' + mesh_term + '"[Mesh] AND '
		search_term = mesh + gene_name
	else:
		search_term = gene_name
	protein_id = Entrez.esearch(db='pubmed', term=search_term, retmax=500)
	pid = Entrez.read(protein_id)
	if return_number:
		return pid['Count']
	if pid['Count'] == '0':
		return pd.DataFrame()
	# handle = Entrez.elink(db='pubmed', dbfrom="gene", id=pid['IdList'][0], linkname="gene_pubmed")
	# rec = Entrez.read(handle)
	# pub_id = [i['Id'] for i in rec[0]['LinkSetDb'][0]['Link']]
	info_pub = Entrez.efetch(db='pubmed', id=pid['IdList'], rettype='medline', retmode='text')
	data = [i for i in Medline.parse(info_pub)]

	df = pd.DataFrame.from_records(data)
	df.rename(index=str, columns=dict_medline, inplace=True)
	pub_type_list = ['Journal Article', 'Case Reports', 'Clinical Trial', 'Comparative Study', 'Letter',
					 'Meta-Analysis', 'Review']
	for pub_type in pub_type_list:
		df[pub_type] = [pub_type in i for i in df['Publication Type'].values]
	columns_to_keep = ['Abstract', 'Affiliation', 'Author', 'Date of Publication',
					   'Journal Title', 'MeSH Terms', 'Other Term',
					   'Other Term Owner', 'Place of Publication', 'PMID',
					   'Subset', 'Source', 'Journal Title Abbreviation', 'Title', 'Volume',
					   'Journal Article', 'Case Reports', 'Clinical Trial',
					   'Comparative Study', 'Letter', 'Meta-Analysis', 'Review']
	for i in columns_to_keep:
		if i not in df.columns:
			df[i] = ''
	df = df[columns_to_keep]
	df['Year of Publication'] = df['Date of Publication'].str.split(' ', expand=True)[0]

	df['PMID'] = pubmed_url + df.PMID + '/'
	neurodeg = []
	chem = []
	major_keywords = []
	for i in df['MeSH Terms'].values:
		if type(i) == float:
			neurodeg.append(False)
			chem.append(False)
			major_keywords.append([])
		else:
			major = []
			neuro = False
			chemistry = False
			for k in i:
				if 'neurodege' in k.lower() or 'alzheimer' in k.lower() or 'dementia' in k.lower() or 'parkinson' in k.lower():
					neuro = True
				if '*' in k:
					major.append(k)
				if '*chemistry' in k or '*Chemistry' in k:
					chemistry = True
			major_keywords.append(' / '.join(major))
			if neuro:
				neurodeg.append(True)
			else:
				neurodeg.append(False)
			if chemistry:
				chem.append(True)
			else:
				chem.append(False)

	df['Neurodegeneration'] = neurodeg
	df['Major Keywords'] = major_keywords
	df['Chemistry'] = chem
	return df


@retryer(max_retries=10, timeout=10)
def open_target_association(gene_id):
	from opentargets import OpenTargetsClient
	ot = OpenTargetsClient()
	associations = ot.get_associations_for_target(gene_id)

	df = associations.to_dataframe()
	df = df[(df['is_direct'] == True)]
	cols = ['target.gene_info.symbol', 'disease.efo_info.therapeutic_area.labels', 'disease.efo_info.label',
			'association_score.overall',
			'association_score.datatypes.genetic_association',
			'association_score.datatypes.known_drug',
			'association_score.datatypes.literature',
			'association_score.datatypes.animal_model',
			'association_score.datatypes.affected_pathway',
			'association_score.datatypes.rna_expression',
			'association_score.datatypes.somatic_mutation']
	rename = {'association_score.datatypes.affected_pathway': 'affected_pathway',
			  'association_score.datatypes.animal_model': 'animal_model',
			  'association_score.datatypes.genetic_association': 'genetic_association',
			  'association_score.datatypes.known_drug': 'known_drug',
			  'association_score.datatypes.literature': 'litterature_mining',
			  'association_score.datatypes.rna_expression': 'rna_expression',
			  'association_score.datatypes.somatic_mutation': 'somatic_mutation',
			  'association_score.overall': 'overall_score', 'disease.efo_info.label': 'disease_name',
			  'disease.efo_info.therapeutic_area.labels': 'disease_area', 'target.gene_info.symbol': 'gene_symbol'}
	df = df[cols]
	df = df.round(2)
	df = df[df['association_score.overall'] > 0.05]
	df.rename(columns=rename, inplace=True)
	return df


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


def write_excel_header(header_dict, worksheet, format):
	for head in header_dict.keys():
		if len(header_dict[head]) == 2:
			row, col = header_dict[head]
			worksheet.write(row, col, head, format)
		elif len(header_dict[head]) == 4:
			row, col, last_row, last_col = header_dict[head]
			worksheet.merge_range(row, col, last_row, last_col, head, format)


def get_single_excel(target_id):
	if target_id in list_of_entries:
		writer = pd.ExcelWriter(output_single_path + list_of_entries[target_id] + '_' + target_id + '.xlsx',
								engine='xlsxwriter',options={'nan_inf_to_errors': True})

		workbook = writer.book

		# ============================ STYLES ===============================#

		bold_center = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center'})
		red = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'red'})
		green = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'green'})
		col_header = workbook.add_format({'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'vcenter'})
		vert_col_header = workbook.add_format({'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'bottom','rotation':90})

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

		# ================= DIFFERENT TAB CREATION ==========================#

		wb_general_info = workbook.add_worksheet('General info')
		writer.sheets['General info'] = wb_general_info
		wb_references = workbook.add_worksheet('References')
		writer.sheets['References'] = wb_references
		wb_disease = workbook.add_worksheet('diseases')
		writer.sheets['diseases'] = wb_disease
		wb_opentarget = workbook.add_worksheet('open_target_association')
		writer.sheets['open_target_association'] = wb_opentarget
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

		# ================= GET THE DIFFERENT DATAFRAMES ====================#

		res = tf.get_single_features(target_id, user=db_user, pwd=db_pwd)

		# =================== FILLING THE WORKSHEETS ========================#

		# =============== GETTING OPENTARGETS + PUBMED DATAFRAMES ======================#
		sequence = None
		pubmed = pd.DataFrame(data=None)
		opentarget = pd.DataFrame(data=None)

		if not res['general_info'].empty:
			if args.email:
				pubmed = pubmed_search(res['general_info'].iloc[0]['Gene_name'], args.email)
			search_term = res['general_info'].iloc[0]['Gene_name']
			opentarget = open_target_association(search_term)

		# ============================ PUBMED TAB ===============================#
		if not pubmed.empty:
			col_order = ['Title', 'Journal Title', 'Year of Publication', 'Journal Article', 'Case Reports',
						 'Clinical Trial', 'Comparative Study', 'Letter', 'Meta-Analysis', 'Review',
						 'Neurodegeneration', 'Chemistry', 'Major Keywords', 'Abstract', 'Author', 'Affiliation',
						 'PMID', 'MeSH Terms', 'Other Term']
			pubmed = pubmed[col_order]
			pubmed.sort_values(by='Year of Publication', ascending=False, inplace=True)
			pubmed.to_excel(writer, sheet_name='Pubmed_search', index=False)
			for col_num, value in enumerate(pubmed.columns.values):
				writer.sheets['Pubmed_search'].write(0, col_num, value, vert_col_header)

		# ========================== OPENTARGET TAB ===============================#
		if not opentarget.empty:
			opentarget.to_excel(writer, sheet_name='open_target_association', index=False)
			for col_num, value in enumerate(opentarget.columns.values):
				writer.sheets['open_target_association'].write(0, col_num, value, vert_col_header)

		# ============================ GENERAL TAB ===============================#

		# GENERAL INFO HEADER WRITING
		header_index = {'Gene_name': (0, 0), 'Synonyms': (1, 0), 'Target_id': (2, 0), 'Protein_class': (3, 0),
						'Protein_class_desc': (4, 0), 'Species': (5, 0), 'Number_isoforms': (6, 0),
						'Sequence': (0, 3), 'Cell_location': (0, 4), 'DISEASE': (8, 0, 8, 1), 'disease_id': (9, 0),
						'disease_name': (9, 1), 'PATHWAYS': (8, 3, 8, 4), 'Reactome': (9, 3), 'KEGG': (9, 4)}
		write_excel_header(header_index, wb_general_info, col_header)

		if not res['general_info'].empty:
			sequence = ''
			for k, v in res['general_info'].iloc[0].items():
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
			target_desc = td.get_descriptors(target_id, user=args.db_username, pwd=args.db_password)
			target_score = td.make_score(target_desc)
			spider_plot = td.make_spider_plot(target_score.loc[0].values, target_score.columns,
											  target_name=res['general_info'].iloc[0]['Gene_name'])
			wb_general_info.insert_image('G1', 'spider_plot', {'image_data': spider_plot})

		if not res['disease'].empty:
			for i in range(len(res['disease'])):
				for k, v in res['disease'].iloc[i].items():
					row, col = header_index[k]
					row = row + i + 1
					if k == 'disease_id':
						wb_general_info.write_url(row, col, 'https://omim.org/entry/' + v.split(':')[1], link)
						wb_general_info.write(row, col, v, link)
					else:
						wb_general_info.write(row, col, v, v_center)

		if not res['reactome'].empty:
			row, col = header_index['Reactome']
			row = row + 1
			wb_general_info.write_column(row, col, list(res['reactome']['pathway_name']), v_center)

		if not res['kegg'].empty:
			row, col = header_index['KEGG']
			row = row + 1
			wb_general_info.write_column(row, col, list(res['kegg']['pathway_name']), v_center)

		# ============================ DISEASE TAB ===============================#

		# DISEASE HEADER WRITING
		dis_header_index = {'DISEASE REGULATION': (0, 0, 0, 4), 'GWAS': (0, 6, 0, 11), 'disease': (1, 0),
							't_stat': (1, 1),
							'std_dev_t': (1, 2), 'n': (1, 3), 'direction': (1, 4), 'phenotype': (1, 6),
							'organism': (1, 7), 'author': (1, 8), 'year': (1, 9), 'p_value': (1, 10),
							'pubmed_id': (1, 11)}
		write_excel_header(dis_header_index, wb_disease, col_header)

		if not res['disease_exp'].empty:
			for i in range(len(res['disease_exp'])):
				for k, v in res['disease_exp'].iloc[i].items():
					row, col = dis_header_index[k]
					row = row + i + 1
					wb_disease.write(row, col, v)
			wb_disease.conditional_format(1, 1, row, 1, {'type': 'data_bar'})
			wb_disease.conditional_format(1, 2, row, 2,
										  {'type': 'icon_set', 'reverse_icons': True, 'icon_style': '3_traffic_lights'})

		if not res['gwas'].empty:
			for i in range(len(res['gwas'])):
				for k, v in res['gwas'].iloc[i].items():
					row, col = dis_header_index[k]
					row = row + i + 1
					if k == 'pubmed_id':
						wb_disease.write_url(row, col, 'https://www.ncbi.nlm.nih.gov/pubmed/' + v, link)
						wb_disease.write(row, col, v, link)
					else:
						wb_disease.write(row, col, v)

		# ============================ EXPRESSION TAB ===============================#

		# EXPRESSION HEADER
		expression_header_index = {'Tissue Expression': (0, 0, 0, 3), 'Tissue': (1, 0), 't_stat': (1, 1),
								   'std_dev_t': (1, 2), 'n': (1, 3), 'Selectivity': (0, 5, 0, 6),
								   'ORGANS': (1, 5, 1, 8), 'organ_name': (2, 5), 'Total_value': (2, 6),
								   'n_tissues': (2, 7), 'avg_value': (2, 8)}
		write_excel_header(expression_header_index, wb_expression, col_header)

		if not res['tissue'].empty:
			for i in range(len(res['tissue'])):
				for k, v in res['tissue'].iloc[i].items():
					row, col = expression_header_index[k]
					row = row + i + 1
					wb_expression.write(row, col, v)
			wb_expression.conditional_format(1, 1, row, 1, {'type': 'data_bar'})
			wb_expression.conditional_format(1, 2, row, 2, {'type': 'icon_set', 'reverse_icons': True,
															'icon_style': '3_traffic_lights'})

		if not res['selectivity'].empty:
			wb_expression.merge_range(0, 7, 0, 8, res['selectivity'].iloc[0]['Selectivity_entropy'], col_header)

		if not res['organ_expression'].empty:
			for i in range(len(res['organ_expression'])):
				for k, v in res['organ_expression'].iloc[i].items():
					row, col = expression_header_index[k]
					row = row + i + 1
					wb_expression.write(row, col, v)
			wb_expression.conditional_format(3, 8, row, 8, {'type': 'data_bar'})

			# organ_chart = workbook.add_chart({'type': 'bar'})
			# organ_chart.add_series({'values': '=expression!$I$4:$I$16',
			#                         'categories': '=expression!$F$4:$F$16',
			#                         'name': 'Organ Expression'})
			# organ_chart.set_legend({'none': True})
			# organ_chart.set_x_axis({'min': 0, 'max': 3, 'major_unit': 1, 'minor_unit_type': 'level',
			#                         'major_gridlines': {'visible': True, 'line': {'width': 1.25, 'dash_type': 'dash'}}})
			# wb_general_info.insert_chart('G1', organ_chart)

		if not res['tissue_expression'].empty:
			previous_organ = ''
			row = 0
			col = 10
			organ_count = 0
			for i in range(len(res['tissue_expression'])):
				if res['tissue_expression'].iloc[i]['organ'] != previous_organ:
					if row >= 65:
						col += 5
						row = 0
					if row == 0:
						pass
					else:
						row += 2
					wb_expression.merge_range(row, col, row, col + 3, res['tissue_expression'].iloc[i]['organ'].upper(),
											  col_header)
					row += 1
					wb_expression.write(row, col, 'tissue name', col_header)
					wb_expression.merge_range(row, col + 1, row, col + 2, 'Cell type', col_header)
					wb_expression.write(row, col + 3, 'Value', col_header)
					previous_organ = res['tissue_expression'].iloc[i]['organ']
					organ_count += 1
				row += 1
				wb_expression.write(row, col, res['tissue_expression'].iloc[i]['tissue'])
				wb_expression.merge_range(row, col + 1, row, col + 2, res['tissue_expression'].iloc[i]['cell'])
				wb_expression.write(row, col + 3, res['tissue_expression'].iloc[i]['value'])
			# brain_chart = workbook.add_chart({'type': 'bar'})
			# brain_chart.add_series({'values': '=expression!$N$3:$N$13',
			#                         'categories': '=expression!$K$3:$L$13',
			#                         'name': 'Brain Expression'})
			# brain_chart.set_legend({'none': True})
			# brain_chart.set_x_axis({'min': 0, 'max': 3, 'major_unit': 1, 'minor_unit_type': 'level',
			#                         'major_gridlines': {'visible': True, 'line': {'width': 1.25, 'dash_type': 'dash'}}})
			# wb_general_info.insert_chart('G17', brain_chart)

		# ============================ PHENOTYPE TAB ===============================#

		if not res['phenotype'].empty:
			row = 0
			col_allele = 0
			col_zyg = 0
			col_gen = 0
			row_data = 4
			row_with_phen = []
			lethal_phen = re.compile('(lethal)|(death)', re.IGNORECASE)
			normal_phen = re.compile('(no abnormal phenotype detected)', re.IGNORECASE)

			for allele, data in res['phenotype'].groupby(['Allele_symbol']):
				tmp_row_with_phen = []
				for zygosity, d2 in data.groupby(['zygosity']):
					for genotype, d3 in d2.groupby(['genotype']):
						lethal = False
						normal = False
						for phen in list(d3['Phenotype'].values):
							if lethal_phen.search(phen):
								wb_genotypes.write(row_data, col_gen, phen, red)
								lethal = True
							elif normal_phen.search(phen):
								wb_genotypes.write(row_data, col_gen, phen, green)
								normal = True
							else:
								wb_genotypes.write(row_data, col_gen, phen)
							row_data += 1
							allele_type = d3.iloc[0]['Allele_type']
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

		# ============================ ISOFORM TAB ===============================#

		if not res['isoforms'].empty:
			row = 0
			mod_header = {'start': 0, 'stop': 1, 'previous_seq': 2, 'modification_type': 3, 'new_seq': 4,
						  'in_domains': 5, 'comments': 6}
			row_to_hide = []
			res['isoforms'] = res['isoforms'].sort_values(by='similarity', ascending=False)
			for iso in res['isoforms'].to_dict(orient='records'):
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
				if not res['isoforms_mod'].empty:
					for header, col in mod_header.items():
						wb_isoforms.write(row, col, header, col_header)
					row+=1
					to_print = res['isoforms_mod'][res['isoforms_mod']['isoform_id']==iso['Isoform_id']]
					to_print = to_print[['start', 'stop', 'previous_seq', 'modification_type', 'new_seq','in_domains', 'comments']]
					to_print.to_excel(writer, sheet_name='isoforms', startrow=row, index=False,header=False)
					row = row + len(to_print)
				to_hide_stop = row + 1
				row_to_hide.append((to_hide_start, to_hide_stop))
				row += 2
			for row1, row2 in row_to_hide:
				for i in range(row1, row2):
					wb_isoforms.set_row(i, None, None, {'level': 1, 'collapsed': True, 'hidden': True})

		# =================== VARIANTS AND MUTANTS TAB ============================#

		if not res['var'].empty and not res['mut'].empty:
			col_var = 0
			col_mut = 0
		elif not res['var'].empty or not res['mut'].empty:
			col_var = 0
			col_mut = 0
		if not res['var'].empty:
			mod_header = {'start': 0, 'stop': 1, 'previous_seq': 3, 'modification_type': 2, 'new_seq': 4,
						  'in_domains': 5,
						  'comments': 6}
			row = 0
			wb_var_mut.merge_range(row, col_var, row, col_var + 6, 'VARIANTS', col_header)
			row += 1
			for header, col in mod_header.items():
				wb_var_mut.write(row, col + col_var, header, col_header)
			for i in range(len(res['var'])):
				row += 1
				for key, value in res['var'].iloc[i].items():
					if key == 'previous_seq' or key == 'new_seq':
						wb_var_mut.write(row, col_var + mod_header[key], value, wrap)
					else:
						wb_var_mut.write(row, col_var + mod_header[key], value)
		if not res['mut'].empty:
			mod_header = {'start': 0, 'stop': 1, 'previous_seq': 3, 'modification_type': 2, 'new_seq': 4,
						  'in_domains': 5,
						  'comments': 6}
			if not res['var'].empty:
				row += 2
			else:
				row = 0
			wb_var_mut.merge_range(row, col_mut, row, col_mut + 6, 'MUTANTS', col_header)
			row += 1
			for header, col in mod_header.items():
				wb_var_mut.write(row, col + col_mut, header, col_header)
			for i in range(len(res['mut'])):
				row += 1
				for key, value in res['mut'].iloc[i].items():
					if key == 'previous_seq' or key == 'new_seq':
						wb_var_mut.write(row, col_mut + mod_header[key], value, wrap)
					else:
						wb_var_mut.write(row, col_mut + mod_header[key], value)

		# ======================== STRUCTURE TAB ==================================#
		row = 2
		col_orig = 0

		if sequence:
			wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'Total length', col_header)
			row += 1
			wb_struct.merge_range(row, col_orig, row, col_orig + 4, len(sequence), bold_center)
			row += 2

		if not res['domains'].empty:
			wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'DOMAINS', col_header)
			row += 1
			col_order = ['Domain_name', 'start', 'stop', 'length', 'source']
			res['domains'] = res['domains'][col_order]
			res['domains'].to_excel(writer, sheet_name='Structure', startrow=row, index=False)
			row += len(res['domains']) + 2

		if not res['domain_drugE'].empty:
			wb_struct.merge_range(row, col_orig, row, col_orig + 4, 'DOMAINS - DrugEbillity', col_header)
			row += 1
			col_order = ['pdb_list', 'domain_fold', 'domain_superfamily', 'tractable', 'druggable']
			res['domain_drugE'] = res['domain_drugE'][col_order]
			res['domain_drugE'].to_excel(writer, sheet_name='Structure', startrow=row, index=False)
			row += len(res['domain_drugE']) + 2

		if not res['pdb_blast'].empty:
			wb_struct.merge_range(row, col_orig, row, col_orig + 6, 'PDB BLAST', col_header)
			row += 1
			col_order = ['PDB_code', 'Chain', 'similarity', 'gene', 'species', 'SITES_tractable', 'SITES_druggable']
			res['pdb_blast'] = res['pdb_blast'][col_order]
			res['pdb_blast'].to_excel(writer, sheet_name='Structure', startrow=row, index=False)

		if not res['pdb'].empty:
			col_orig = 6
			if not res['pdb_blast'].empty:
				col_orig = col_orig + 2
			row = 0
			wb_struct.merge_range(row, col_orig, row, col_orig + 7, 'PDB', col_header)
			wb_struct.merge_range(row, col_orig + 8, row, col_orig + 14, 'PDB: Ligand', col_header)
			wb_struct.merge_range(row, col_orig + 15, row, col_orig + 16, 'ChEMBL - DruggEbillity', col_header)
			row += 1

			pdb = res['pdb'].copy()
			pdb['% of full protein'] = round((pdb['n_residues'] / len(sequence)) * 100, 1)
			pdb.operator = ' ' + pdb.operator
			col_order = ['PDB_code', 'Technique', 'Resolution', 'Chain', 'Domain_name',
						 'n_residues', '% of full protein', 'start_stop', 'type_of_binder', 'binding_type', 'operator',
						 'value', 'units',
						 'Ligand_name', 'publication_year', 'SITES_tractable', 'SITES_druggable']
			pdb = pdb[col_order]
			pdb.to_excel(writer, sheet_name='Structure', startrow=row, startcol=col_orig, index=False)
			for col_num, value in enumerate(pdb.columns.values):
				wb_struct.write(row, col_num + col_orig, value, vert_col_header)

		# ======================== POCKETS TAB ==================================#
		col_alt_pocket = 0

		if not res['pockets'].empty:
			col_alt_pocket = 9
			col_order = ['PDB_code', 'druggability_score', 'pocket_score', 'pocket_number',
						 'volume', 'area', 'fraction_apolar', 'domains']
			res['pockets'] = res['pockets'][col_order]
			res['pockets'].to_excel(writer, sheet_name='Pockets', startrow=1, index=False)
			wb_pockets = writer.sheets['Pockets']
			wb_pockets.merge_range(0, 0, 0, 7, 'DRUGGABLE POCKETS', col_header)
		
		if not res['alt_pockets'].empty:
			col_order = ['PDB_code', 'druggability_score', 'pocket_score', 'pocket_number',
						 'volume', 'area', 'fraction_apolar', 'gene', 'species', 'similarity']
			res['alt_pockets'] = res['alt_pockets'][col_order]
			res['alt_pockets'].to_excel(writer, sheet_name='Pockets', startrow=1, index=False, startcol=col_alt_pocket)
			wb_pockets = writer.sheets['Pockets']
			wb_pockets.merge_range(0, 0 + col_alt_pocket, 0, 9 + col_alt_pocket,
								   'ALTERNATE DRUGGABLE POCKETS (PDB from blast)', col_header)

		# =================== CHEMBL BIOACTIVITIES TABS ===========================#
		header_groups = {'Bioactivity info': (0, 9), 'Assay info': (10, 14), 'Structure': (15, 15),
		                 'Ligand properties': (16, 30), 'Ligand info': (31, 36), 'References': (37, 38)}
		CNS_MPO_criteria = [{'criteria': '>=', 'type': 'number', 'value': 4.5},
							{'criteria': '>=', 'type': 'number', 'value': 3.5},
							{'criteria': '<', 'type': 'number', 'value': 3}]
		CNS_MPO_col = 'AE1:AE'

		if not res['binding'].empty:
			res['binding'].to_excel(writer, sheet_name='Binding', startrow=1, index=False)
			w_bd = writer.sheets['Binding']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_bd.write(0, span[0], head, col_header)
				else:
					w_bd.merge_range(0, span[0], 0, span[1], head, col_header)
			w_bd.conditional_format(CNS_MPO_col + (str(len(res['binding']) + 3)),
									{'type': 'icon_set', 'icon_style': '3_traffic_lights'
										, 'icons': CNS_MPO_criteria})
			for col_num, value in enumerate(res['binding'].columns.values):
				writer.sheets['Binding'].write(1, col_num, value, vert_col_header)

		if not res['dose_response'].empty:
			res['dose_response'].to_excel(writer, sheet_name='Dose_response', startrow=1, index=False)
			w_dr = writer.sheets['Dose_response']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_dr.write(0, span[0], head, col_header)
				else:
					w_dr.merge_range(0, span[0], 0, span[1], head, col_header)
			w_dr.conditional_format(CNS_MPO_col + (str(len(res['dose_response']) + 3)),
									{'type': 'icon_set', 'icon_style': '3_traffic_lights'
										, 'icons': CNS_MPO_criteria})
			for col_num, value in enumerate(res['dose_response'].columns.values):
				writer.sheets['Dose_response'].write(1, col_num, value, vert_col_header)

		if not res['percent_inhibition'].empty:
			res['percent_inhibition'].to_excel(writer, sheet_name='Percent_inhibition', startrow=1, index=False)
			w_per = writer.sheets['Percent_inhibition']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_per.write(0, span[0], head, col_header)
				else:
					w_per.merge_range(0, span[0], 0, span[1], head, col_header)
			for col_num, value in enumerate(res['percent_inhibition'].columns.values):
				writer.sheets['Percent_inhibition'].write(1, col_num, value, vert_col_header)
			w_per.conditional_format(CNS_MPO_col + (str(len(res['percent_inhibition']) + 3)),
									 {'type': 'icon_set', 'icon_style': '3_traffic_lights'
										 , 'icons': CNS_MPO_criteria})

		if not res['efficacy_bio'].empty:
			res['efficacy_bio'].to_excel(writer, sheet_name='Emax_Efficacy', startrow=1, index=False)
			w_eff = writer.sheets['Emax_Efficacy']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_eff.write(0, span[0], head, col_header)
				else:
					w_eff.merge_range(0, span[0], 0, span[1], head, col_header)
			for col_num, value in enumerate(res['efficacy_bio'].columns.values):
				writer.sheets['Emax_Efficacy'].write(1, col_num, value, vert_col_header)
			w_eff.conditional_format(CNS_MPO_col + (str(len(res['efficacy_bio']) + 3)),
									 {'type': 'icon_set', 'icon_style': '3_traffic_lights'
										 , 'icons': CNS_MPO_criteria})
			row_efficacy = len(res['efficacy_bio']) + len(res['efficacy_bio'].columns.levels) + 1
		else:
			row_efficacy = 0
		if not res['emax'].empty:
			res['emax'].to_excel(writer, sheet_name='Emax_Efficacy', startrow=row_efficacy + 1, index=False)
			w_eff = writer.sheets['Emax_Efficacy']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_eff.write(row_efficacy, span[0], head, col_header)
				else:
					w_eff.merge_range(row_efficacy, span[0], row_efficacy, span[1], head, col_header)
			for col_num, value in enumerate(res['emax'].columns.values):
				writer.sheets['Emax_Efficacy'].write(row_efficacy+1, col_num, value, vert_col_header)
			w_eff.conditional_format(CNS_MPO_col + (str(len(res['emax']) + row_efficacy + 3)),
									 {'type': 'icon_set', 'icon_style': '3_traffic_lights'
										 , 'icons': CNS_MPO_criteria})

		if not res['ADME'].empty:
			res['ADME'].to_excel(writer, sheet_name='ADME', startrow=1, index=False)
			w_adme = writer.sheets['ADME']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_adme.write(0, span[0], head, col_header)
				else:
					w_adme.merge_range(0, span[0], 0, span[1], head, col_header)
			for col_num, value in enumerate(res['ADME'].columns.values):
				writer.sheets['ADME'].write(1, col_num, value, vert_col_header)
			w_adme.conditional_format(CNS_MPO_col + (str(len(res['ADME']) + 3)),
									  {'type': 'icon_set', 'icon_style': '3_traffic_lights'
										  , 'icons': CNS_MPO_criteria})

		if not res['other'].empty:
			res['other'].to_excel(writer, sheet_name='Other_bioactivities', startrow=1, index=False)
			w_other = writer.sheets['Other_bioactivities']
			for head, span in header_groups.items():
				if span[0] == span[1]:
					w_other.write(0, span[0], head, col_header)
				else:
					w_other.merge_range(0, span[0], 0, span[1], head, col_header)
			for col_num, value in enumerate(res['other'].columns.values):
				writer.sheets['Other_bioactivities'].write(1, col_num, value, vert_col_header)
			w_other.conditional_format(CNS_MPO_col + (str(len(res['other']) + 3)),
									   {'type': 'icon_set', 'icon_style': '3_traffic_lights'
										   , 'icons': CNS_MPO_criteria})

		# ======================== BINDING DB TAB ==================================#

		if not res['bindingDB'].empty:
			bdb = res['bindingDB'].copy()
			columns = ['ZincID', 'IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)', 'kon(M-1s-1)', 'koff(s-1)', 'pH', 'Temp',
					   'Source', 'DOI', 'patent_number', 'institution', 'ligand_name', 'SMILES', 'HBA', 'HBD', 'LogD',
					   'LogP', 'MW', 'TPSA', 'aLogP', 'apKa', 'bpKa', 'nAr', 'n_alerts', 'pass_ro3', 'ro5_violations',
					   'rotB', 'CNS_MPO', 'mol_name', 'molecular_species', 'indication_class', 'class_def', 'max_phase',
					   'oral']
			bdb = bdb[((bdb[['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)', 'kon(M-1s-1)', 'koff(s-1)']] <= 10000) & (
				bdb[['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)', 'kon(M-1s-1)', 'koff(s-1)']].notna())).any(axis=1)]
			bdb.loc[(bdb[['LogP', 'LogD', 'MW', 'HBD', 'TPSA']].notnull().all(axis=1)) & (
				bdb.bpKa.isnull()), 'bpKa'] = 0
			bdb['CNS_MPO'] = mpo.calc_mpo_score(bpka=bdb['bpKa'], logP=bdb['LogP'], logD=bdb['LogD'], MW=bdb['MW'],
												HBD=bdb['HBD'], TPSA=bdb['TPSA'])
			if not bdb.empty:
				bdb = bdb[columns]
				bdb.to_excel(writer, sheet_name='BindingDB', index=False)
				w_bindingDB = writer.sheets['BindingDB']
				for col_num, value in enumerate(bdb.columns.values):
					w_bindingDB.write(0, col_num, value, vert_col_header)
				CNS_MPO_criteria = [{'criteria': '>=', 'type': 'number', 'value': 4.5},
									{'criteria': '>=', 'type': 'number', 'value': 3.5},
									{'criteria': '<', 'type': 'number', 'value': 3}]
				CNS_MPO_col = 'AD1:AD'
				w_bindingDB.conditional_format(CNS_MPO_col + (str(len(bdb) + 3)),
											   {'type': 'icon_set', 'icon_style': '3_traffic_lights'
												   , 'icons': CNS_MPO_criteria})

		# ======================== COMMERCIAL CPDS TAB ==================================#

		if not res['commercials'].empty:
			col_order = ['smiles', 'affinity_type', 'op', 'affinity_value', 'affinity_unit', 'price', 'website']
			comm = res['commercials'].copy()
			comm.sort_values(by='affinity_value', inplace=True)
			comm = comm[col_order]
			comm.to_excel(writer, sheet_name='Commercial compounds', index=False)
			for col_num, value in enumerate(comm.columns.values):
				writer.sheets['Commercial compounds'].write(0, col_num, value, vert_col_header)

		# ======================== CLOSING FILE ==================================#
		workbook.close()
	else:
		return print("Gene with ID [", target_id, '] not present in the database. Run the command without the -R flag '
												  'to run the analysis on it')


def get_list_excel(list_targets):
	list_to_do = {}
	not_in_db = {'Not present in DB': []}
	for g in list_targets:
		in_db = False
		if '_' in g:
			g = g.split('_')[0]
		for entry in gene_in_db:
			if g.upper() == entry['name'].upper():
				list_to_do[entry['ID']] = g
				in_db = True
				continue
			elif g.upper() in entry['syn']:
				list_to_do[entry['ID']] = g
				in_db = True
				continue
		if not in_db:
			not_in_db['Not present in DB'].append(g)
	if not list_to_do:
		return print("No genes that you entered are in the Database")

	not_in_db = pd.DataFrame.from_dict(not_in_db)
	pubmed = {'ID': [], 'total # publications': [], 'number of Dementia publications': []}
	for id, name in list_to_do.items():
		pubmed['ID'].append(id)
		pubmed['number of Dementia publications'].append(
			pubmed_search(name, args.email, return_number=True, mesh_term='Dementia'))
		pubmed['total # publications'].append(pubmed_search(name, args.email, return_number=True))
		# TODO: Mesh term is hard-coded here
	pubmed = pd.DataFrame.from_dict(pubmed)

	dbase = db.open_db(druggability_db, pwd=args.db_password, user=args.db_username)

	gene_ids = "','".join(list_to_do.keys())

	list_queries = {'gen': """SELECT
T.Gene_name
,T.Target_id as Uniprot_id
,T.Target_id as ID
,T.Species
,(CASE WHEN T.Number_isoforms=0 THEN 1 ELSE T.Number_isoforms END) Number_isoforms
,T.Protein_class_desc
,T.Protein_class_short
,T.Synonyms
,LENGTH(T.Sequence) as number_of_residues
FROM Targets T
WHERE Target_id in ('%s')""" % gene_ids,
					'domains': """SELECT
 D.Target_id ID,
 GROUP_CONCAT(D.domain SEPARATOR '\n') AS domain
  FROM
(SELECT
   D.Target_id
  ,CONCAT(D.source_name,'\n',GROUP_CONCAT(CONCAT('\t',D.Domain_name,' (',D.Domain_start,'-',D.Domain_stop,')') ORDER BY D.Domain_start SEPARATOR '\n')) as domain
FROM Domain_targets D
  WHERE D.Target_id in ('%s')
	GROUP BY D.Target_id,D.source_name) D
GROUP BY D.Target_id""" % gene_ids,
					'mutant': """SELECT
  Target_id ID,
  GROUP_CONCAT(CONCAT('(',start,') ',previous,'-->',new,' comment: ',SUBSTRING_INDEX(comment,'.',1),'; in domains: ',domains) ORDER BY start SEPARATOR '\n') as MUTANT
  FROM modifications
	WHERE mod_type ='MUTAGEN'
	AND Target_id in ('%s')
GROUP BY Target_id""" % gene_ids,
					'variant': """SELECT
  Target_id ID,
  concat(substring_index(GROUP_CONCAT(CONCAT('(',start,') ',previous,'-->',new,' comment: ',SUBSTRING_INDEX(comment,'.',1),'; in domains: ',domains) ORDER BY start SEPARATOR '\n'),'\n',15),case when count(comment) > 15 THEN  concat('\n+ ',count(comment)-15,' others') ELSE '' END)  as VARIANT
  FROM modifications
	WHERE mod_type = 'VAR'
	AND Target_id in ('%s')
GROUP BY Target_id""" % gene_ids,
					'pdb': """SELECT
  T.Target_id ID,
  concat(substring_index(GROUP_CONCAT(CONCAT(T.PDB_code,': ',T.n_residues,' residues (',T.start_stop,', ',T.P_seq,'%%) Chain(s): ',T.Chain,' Domain(s): ',CASE WHEN T.domain is NULL THEN '' ELSE T.domain END,' (',T.Technique,': ',T.Resolution,')') ORDER BY T.P_seq DESC SEPARATOR '\n'),'\n',15),case when count(T.PDB_code) > 15 THEN  concat('\n+ ',count(T.PDB_code)-15,' others') ELSE '' END) AS PDB
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
  FROM (SELECT * FROM PDB_Chains C WHERE C.Target_id in ('%s')) C
LEFT JOIN PDB P
	ON C.PDB_code=P.PDB_code
LEFT JOIN PDBChain_Domain D
	ON C.Chain_id=D.Chain_id
LEFT JOIN Domain_targets DT
	ON D.Domain_id=DT.domain_id
LEFT JOIN Targets T
	ON C.Target_id = T.Target_id
GROUP BY C.Target_id,C.PDB_code)T
	GROUP BY T.Target_id""" % gene_ids,
					'blast': """SELECT
  Query_target_id as ID,
concat(substring_index(GROUP_CONCAT(CONCAT(Hit_gene_name,'_',Hit_gene_species,' (',similarity,'%%)') ORDER BY similarity DESC SEPARATOR '\n'),'\n',10),case when count(Hit_gene_name) > 10 THEN  concat('\n+ ',count(Hit_gene_name)-10,' others') ELSE '' END) as protein_blast
  FROM protein_blast
	WHERE Query_target_id in ('%s')
GROUP BY Query_target_id""" % gene_ids,
					'pdbblast': """SELECT
  Query_target_id as ID,
concat(substring_index(GROUP_CONCAT(CONCAT(Hit_PDB_code,' Chain: ',Chain_Letter,' (',Hit_gene_name,'_',Hit_gene_species,' - ',similarity,'%%)') ORDER BY similarity DESC SEPARATOR '\n'),'\n',10),case when count(Hit_gene_name) > 10 THEN  concat('\n+ ',count(Hit_gene_name)-10,' others') ELSE '' END) as pdb_blast
  FROM `3D_Blast`
	WHERE Query_target_id in ('%s')
GROUP BY Query_target_id""" % gene_ids,
					'pockets': """SELECT
  P.Target_id ID,
  GROUP_CONCAT(P.pockets_domains SEPARATOR '\n') pockets
FROM
(SELECT
  POCK.Target_id,
  CONCAT(POCK.Domain_name,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',POCK.PDB_code,': Druggability_score=',POCK.DrugScore,' Volume=',POCK.volume,' Area=',POCK.total_sasa,' (',POCK.Fraction_apolar,'%% apolar)(',POCK.Pocket_number,')')ORDER BY POCK.DrugScore DESC SEPARATOR '\n'),'\n',3),case when count(POCK.Pocket_number) > 3 THEN  concat('\n\t+ ',count(POCK.Pocket_number)-3,' others') ELSE '' END)) AS pockets_domains
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
  AND blast='FALSE' AND Target_id in ('%s'))FP
  LEFT JOIN fPockets_Domain FPD
	ON FP.Pocket_id=FPD.Pocket_id
) T1
LEFT JOIN Domain_targets D
	ON T1.Domain_id=D.domain_id) POCK
GROUP BY POCK.Target_id,POCK.Domain_name) P
GROUP BY P.Target_id""" % gene_ids,
					'altpockets': """SELECT
  ALT_POCK.Target_id ID,
  GROUP_CONCAT(ALT_POCK.pocket_gene_name SEPARATOR '\n') alt_pockets
FROM
(SELECT POCK.Target_id,
	CONCAT(POCK.Hit_gene_name,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',POCK.PDB_code,'(Similarity=',ROUND(POCK.similarity),'%%): Druggability_score=',ROUND(POCK.DrugScore,2),' Volume=',ROUND(POCK.volume,1),' Area=',ROUND(POCK.total_sasa,1),' (',ROUND(POCK.Fraction_apolar),'%% apolar)(',POCK.Pocket_number,')')ORDER BY POCK.similarity DESC,POCK.DrugScore DESC SEPARATOR '\n'),'\n',3),case when count(POCK.Pocket_number) > 3 THEN  concat('\n\t+ ',count(POCK.Pocket_number)-3,' others') ELSE '' END)) AS pocket_gene_name
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
  FROM (SELECT * FROM fPockets WHERE Target_id in ('%s') AND druggable='TRUE' AND blast='TRUE') FP
	LEFT JOIN 3D_Blast 3D
	ON 3D.Query_target_id=FP.Target_id AND 3D.Hit_PDB_code=FP.PDB_code
  WHERE 3D.similarity>=70) POCK
GROUP BY POCK.Target_id,POCK.Hit_gene_name) ALT_POCK
GROUP BY ALT_POCK.Target_id""" % gene_ids,
					'disease_expression': """SELECT
  T1.Target_id ID,
  concat(substring_index(GROUP_CONCAT(CASE WHEN T1.up is null THEN null ELSE CONCAT(T1.up,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END))END ORDER BY T1.t_stat DESC SEPARATOR '\n'),'\n',20),case when count(T1.up) > 20 THEN  concat('\n+ ',count(T1.up)-20,' others') ELSE '' END) AS upregulated_in_disease
,  concat(substring_index(GROUP_CONCAT(CASE WHEN T1.down is null THEN null ELSE CONCAT(T1.down,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END))END ORDER BY T1.t_stat SEPARATOR '\n'),'\n',20),case when count(T1.down) > 20 THEN  concat('\n+ ',count(T1.down)-20,' others') ELSE '' END) AS downregulated_in_disease

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
	AND Target_id in ('%s')
  GROUP BY Target_id,disease
  ) T1
	) T1
GROUP BY T1.Target_id""" % gene_ids,
					'tissue_expression': """SELECT
  T1.Target_id ID,
  GROUP_CONCAT(CASE WHEN T1.up is null THEN null ELSE CONCAT(T1.up,' (T-stat=',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END)) END ORDER BY T1.t_stat DESC SEPARATOR '\n') AS overexpressed_in,
  GROUP_CONCAT(CASE WHEN T1.down is null THEN null ELSE CONCAT(T1.down,' (T-stat= ',round(T1.t_stat,1),CASE WHEN T1.n_number>1 THEN CONCAT(' +/- ',ROUND(T1.std_dev_t,2)) ELSE '' END,')',(CASE WHEN T1.n_number>1 THEN CONCAT('(n=',T1.n_number,')') ELSE '' END)) END ORDER BY T1.t_stat SEPARATOR '\n') AS underexpressed_in
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
  AND Target_id in ('%s')
GROUP BY Target_id,Tissue)T1
GROUP BY T1.Target_id""" % gene_ids,
					'pathways': """SELECT
  P.Target_id ID,
  GROUP_CONCAT(P.pathways SEPARATOR '\n') AS pathways
  FROM
(SELECT
  Target_id,
  CONCAT(pathway_dataset,'\n',GROUP_CONCAT(CONCAT('\t',pathway_name) ORDER BY pathway_name SEPARATOR '\n')) AS pathways
  FROM pathways
	WHERE pathway_dataset='KEGG pathways data set' AND Target_id in ('%s')
GROUP BY Target_id,pathway_dataset) P
GROUP BY P.Target_id""" % gene_ids,
					'phenotypes': """SELECT
  T1.Target_id ID,
  GROUP_CONCAT(genotype_list ORDER BY T1.zygosity SEPARATOR '\n') AS genotypes,
  GROUP_CONCAT(T1.lethal_phenotype SEPARATOR '\n') lethal_phenotype
  ,GROUP_CONCAT(T1.normal_genotype SEPARATOR '\n') normal_phenotype_for
  FROM
(SELECT
 T1.Target_id,
  T1.zygosity,
  CONCAT(' [',T1.zygosity,']\n',GROUP_CONCAT(CONCAT('\t',T1.normal_genotype) SEPARATOR '\n')) as normal_genotype,
  CONCAT(' [',T1.genotype,']\n',GROUP_CONCAT(DISTINCT CONCAT('\t',T1.lethal_phen) SEPARATOR '\n')) as lethal_phenotype,
  CONCAT(T1.zygosity,'\n',concat(substring_index(GROUP_CONCAT(CONCAT('\t',T1.genotype,' [',T1.organism,']',(CASE WHEN T1.phen_list like 'no abnormal phenotype detected' THEN ' [NORMAL PHENOTYPE]' WHEN T1.phen_list like '%%lethal%%' THEN ' [LETHAL PHENOTYPE OBSERVED]' ELSE '[P]' END)) ORDER BY T1.genotype SEPARATOR '\n'),'\n',5),case when count(T1.genotype) > 5 THEN  concat('\n\t+ ',count(T1.genotype)-5,' others') ELSE '' END)) as genotype_list
  FROM
(SELECT
  PHEN.Target_id,
  PHEN.genotype,
  (CASE WHEN PHEN.zygosity is NULL THEN 'not declared' ELSE PHEN.zygosity END) zygosity,
  PHEN.organism,
  GROUP_CONCAT(DISTINCT PHEN.Phenotype SEPARATOR ' ; ') AS phen_list,
  GROUP_CONCAT(DISTINCT (CASE WHEN PHEN.Phenotype like '%%lethal%%' THEN PHEN.Phenotype END) SEPARATOR '\n\t') AS lethal_phen,
  GROUP_CONCAT(DISTINCT (CASE WHEN PHEN.Phenotype like 'no abnormal phenotype detected' THEN PHEN.genotype END) SEPARATOR '\n') AS normal_genotype
 FROM phenotype PHEN
   WHERE Target_id in ('%s')
	GROUP BY PHEN.Target_id,PHEN.genotype,PHEN.zygosity)T1
GROUP BY T1.Target_id,T1.zygosity)T1
GROUP BY T1.Target_id""" % gene_ids,
					'diseases': """SELECT
  Target_id ID,
  GROUP_CONCAT(CONCAT(disease_name,' [',disease_id,']') ORDER BY disease_name SEPARATOR '\n') AS disease
  FROM disease
	WHERE Target_id in ('%s')
GROUP BY Target_id""" % gene_ids,
					'protexpression_sel': """SELECT PROT_SEL.Target_id ID
	  ,PROT_SEL.max_organ
	  ,ROUND(PROT_SEL.Selectivity_entropy,3) AS expression_selectivity
  FROM protein_expression_selectivity PROT_SEL
WHERE PROT_SEL.Target_id in ('%s')""" % gene_ids,
					'protAtlas': """SELECT
  T1.Target_id ID,
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

					   END),'\t',T1.organ,' (',ROUND(T1.avg_level_num,1),')') ORDER BY T1.avg_level_num DESC SEPARATOR '\n') AS protein_atlas_expression
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
  FROM protein_expression_levels T1
  WHERE T1.Target_id in ('%s')) T1
GROUP BY T1.Target_id,T1.organ)T1
GROUP BY T1.Target_id""" % gene_ids,
					'bioactivities': """SELECT
	B.ID,
	COUNT(DISTINCT L.lig_id) AS Number_of_ligands,
	MAX(L.max_phase) AS Max_phase
	FROM
	(SELECT B.*,C.target_id as ID FROM bioactivities B
	INNER JOIN (SELECT * FROM Crossref C WHERE C.target_id in ('%s')) C
	on B.Target_id=C.Chembl_id) B
	LEFT JOIN ligands L
	on B.lig_id=L.lig_id
	GROUP BY B.target_id""" % gene_ids,
					'assays': """SELECT
  AT.target_id ID,
  GROUP_CONCAT(DISTINCT A.bioactivity_type ORDER BY A.bioactivity_type SEPARATOR '\n') AS Assay_types
  FROM (SELECT * FROM assay_target AT WHERE AT.target_id in ('%s')) AT
  LEFT JOIN assays A
	ON A.assay_id=AT.assay_id
GROUP BY AT.target_id""" % gene_ids,
					'gwas': """SELECT
   G.Target_id ID
  ,GROUP_CONCAT(CONCAT('Phenotype: ',G.phenotype,' Organism: ',G.organism,' (',G.first_author,'-',G.publication_year,') doi:',G.doi,' PID:',G.pubmed_id) ORDER BY G.phenotype,G.publication_year DESC SEPARATOR '\n') as gwas
FROM gwas G
  WHERE G.Target_id in ('%s')
	GROUP BY G.Target_id""" % gene_ids,
					'commercial': """SELECT
   target_id ID,
   GROUP_CONCAT(CONCAT('Affinity: ',affinity_type,': ',affinity_value,affinity_unit,' (price: ',price,') (website: ',website,')') SEPARATOR '\n') as commercially_available
FROM purchasable_compounds
WHERE target_id in ('%s') AND affinity_value <= 500
GROUP BY target_id""" % gene_ids
					}

	results = {qname: pd.read_sql(query, con=dbase.db) for qname, query in list_queries.items()}

	if not results['gen'].empty:
		all = pd.DataFrame.from_records(results['gen'])
	else:
		print(
			"[EXPORT ERROR]: Something went wrong during the export process \nYour requested genes might not be present in the database, or the database is not available at the moment\nPlease try again later")
		return "Failure"

	for name, res in results.items():
		if name == 'gen':
			continue
		else:
			all = all.merge(res, on='ID', how='left')
	all = all.merge(pubmed, on='ID', how='left')

	dbase.close()
	header = ['Gene_name', 'Uniprot_id', 'Synonyms', 'Species', 'pathways', 'disease', 'upregulated_in_disease',
			  'downregulated_in_disease', 'gwas', 'genotypes', 'lethal_phenotype', 'normal_phenotype_for',
			  'overexpressed_in', 'underexpressed_in', 'protein_atlas_expression', 'max_organ',
			  'expression_selectivity',
			  'Number_isoforms', 'Protein_class_desc', 'Protein_class_short', 'number_of_residues', 'domain', 'MUTANT',
			  'VARIANT', 'PDB', 'pdb_blast', 'protein_blast', 'pockets', 'alt_pockets', 'Number_of_ligands',
			  'commercially_available', 'Max_phase', 'Assay_types', 'total # publications'
		, 'number of Dementia publications']
	all = all[header]
	t = time.strftime("%d%b%Y_%H%M%S")
	writer = pd.ExcelWriter(output_lists_path + 'Export_' + str(len(all)) + '_entries_' + t + '.xlsx',
							engine='xlsxwriter')
	all.to_excel(writer, 'Druggability_list', index=False)
	not_in_db.to_excel(writer, 'Not in DB', index=False)
	writer.save()
	print("[EXPORT]: Excel file: ", '[Export_' + str(len(all)) + '_entries_' + t + '.xlsx]', ' successfully generated')
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
		self.pubmed = None
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
	parser.add_argument('-pubmed_email', '--email',
						help='enter your email address (required for the pubmed search', metavar='',
						type=str)
	args = parser.parse_args()

	if args.db_password and args.db_username:
		db_pwd = args.db_password
		db_user = args.db_username
	else:
		print("[ERROR]: Please provide username and password for the MySQL database (-db_user / -db_pwd)")
		parser.print_help()
		sys.exit()

	while True:
		if args.in_file:
			if os.path.exists(args.in_file):
				with open(args.in_file, 'r') as gene_list:
					if args.Uniprot_ID:
						gene_dict = {}
						for i in gene_list.readlines():
							i = i.rstrip('\n').split(',')
							gene_dict[i[1]] = i[0]
					else:
						list_of_genes = gene_list.readlines()
						gene_dict = gene_to_uniprotid_local(list_of_genes)
				break
			else:
				print('ERROR : file inputed as argument [-i] does not exist')
				sys.exit()
		if args.gene:
			list_of_genes = [args.gene]
			gene_dict = gene_to_uniprotid_local(list_of_genes)
			break
		elif args.list_genes:
			list_of_genes = args.list_genes.split(',')
			gene_dict = gene_to_uniprotid_local(list_of_genes)
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

		# if args.report_list or args.report_single:
		# 	list_of_entries, gene_in_db = get_list_entries()
		# 	if args.report_list:
		# 		export_list = []
		# 		for key in targets_list:
		# 			if targets_list[key] == 'Target completed' or targets_list[key] == 'Already present':
		# 				export_list.append(key)
		# 		get_list_excel(export_list)
		# 	if args.report_single:
		# 		for key in targets_list:
		# 			if targets_list[key] == 'Target completed' or targets_list[key] == 'Already present':
		# 				get_single_excel(gene_dict[key])
