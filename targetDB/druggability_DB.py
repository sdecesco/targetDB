#!/usr/bin/env python
import argparse
import configparser
import subprocess
import sys
import os
import time
import sqlite3
import urllib.request as urllib
from urllib.error import *

import pandas as pd
from Bio import ExPASy, Entrez, Medline
from Bio import SwissProt
from Bio import pairwise2
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from intermine.webservice import Service
from opentargets import OpenTargetsClient
from pathlib import Path

from targetDB import drugg_errors
from targetDB.protein_atlas_api import proteinatlas as patlas
from targetDB.utils import config as cf, pocket_finder as pocket
from targetDB.utils import pdb_parser
from targetDB.utils import retryers as ret
from targetDB.utils import gene2id as g2id
from targetDB.utils import targetDB_init as tinit


def get_list_entries(target_db_path=None):
	connector = sqlite3.connect(target_db_path)
	query = "SELECT Target_id,Gene_name FROM Targets"
	entries_list = pd.read_sql(query, con=connector, index_col='Target_id')
	connector.close()
	return entries_list


@ret.retryer(max_retries=10, timeout=10)
def get_uniprot(gene_id):
	try:
		handle = ExPASy.get_sprot_raw(gene_id)
	except:
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
	PDB_list = pd.DataFrame(columns=['PDB_code', 'Technique', 'Resolution','Chain', 'Domain', 'path'])
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
			PDB_list.loc[ref[1]] = {'PDB_code': ref[1], 'Technique': ref[2], 'Resolution': ref[3],
			                    'Chain': chain, 'Domain': [], 'path': ''}
		elif ref[0] == 'GO':
			GO = ref[2].split(':')
			if GO[0] == 'P':
				GO_process.append(GO[1])
			if GO[0] == 'F':
				GO_Function.append(GO[1])
			if GO[0] == 'C':
				GO_Component.append(GO[1])

		elif ref[0] == 'ChEMBL':
			chembl_id = ref[1]
	go = pd.DataFrame.from_dict(
		{'P': ' / '.join(GO_process), 'F': ' / '.join(GO_Function), 'C': ' / '.join(GO_Component)},
		orient='index').rename(columns={0: 'value'})

	return PDB_list, go, chembl_id


@ret.retryer(max_retries=10, timeout=10)
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

	if query.count() == 1:
		primary_id = list(query.rows())[0]["id"]
	else:
		pathways_df = pd.DataFrame(columns=['pathway_dataset', 'pathway_name'])
		disease_df = pd.DataFrame(columns=['disease', 'disease_id'])
		phenotypes_df = pd.DataFrame(
			columns=['gene', 'organism', 'Allele_id', 'Phenotype', 'Phenotype_desc', 'genotype', 'zygosity',
			         'Allele_symbol', 'Allele_type'])
		differential_exp_diseases_df = pd.DataFrame(
			columns=['Condition', 'T_statistic', 'expression_status', 'p_value'])
		differential_exp_tissues_df = pd.DataFrame(columns=['T_statistic', 'Tissue', 'expression_status', 'p_value'])
		gwas_df = pd.DataFrame(columns=['doi', 'first_author', 'organism', 'p_value', 'phenotype',
		                                'publication_year', 'pubmed_id'])
		return disease_df, phenotypes_df, differential_exp_diseases_df, differential_exp_tissues_df, gwas_df, pathways_df

	query = service.new_query("Gene")
	query.add_view("diseases.primaryIdentifier", "diseases.name")
	query.add_constraint("id", "=", primary_id, code="A")

	for row in query.rows():
		disease.append({'disease': row["diseases.name"], 'disease_id': row["diseases.primaryIdentifier"]})

	disease_df = pd.DataFrame.from_records(disease)
	if disease_df.empty:
		disease_df = pd.DataFrame(columns=['disease', 'disease_id'])

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
	phenotypes_df = pd.DataFrame.from_records(phenotypes)
	if phenotypes_df.empty:
		phenotypes_df = pd.DataFrame(
			columns=['gene', 'organism', 'Allele_id', 'Phenotype', 'Phenotype_desc', 'genotype', 'zygosity',
			         'Allele_symbol', 'Allele_type'])

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
	differential_exp_diseases_df = pd.DataFrame.from_records(differential_exp_diseases)
	if differential_exp_diseases_df.empty:
		differential_exp_diseases_df = pd.DataFrame(
			columns=['Condition', 'T_statistic', 'expression_status', 'p_value'])
	differential_exp_tissues_df = pd.DataFrame.from_records(differential_exp_tissues)
	if differential_exp_tissues_df.empty:
		differential_exp_tissues_df = pd.DataFrame(columns=['T_statistic', 'Tissue', 'expression_status', 'p_value'])

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

	gwas_df = pd.DataFrame.from_records(gwas)
	if gwas_df.empty:
		gwas_df = pd.DataFrame(columns=['doi', 'first_author', 'organism', 'p_value', 'phenotype',
		                                'publication_year', 'pubmed_id'])

	query = service.new_query("Gene")
	query.add_view("pathways.name", "pathways.dataSets.name")
	query.add_constraint("id", "=", primary_id, code="A")
	query.add_constraint("organism.name", "=", "Homo sapiens", code="B")

	for row in query.rows():
		pathways.append({'pathway_name': row["pathways.name"], 'pathway_dataset': row["pathways.dataSets.name"]})

	pathways_df = pd.DataFrame.from_records(pathways)
	if pathways_df.empty:
		pathways_df = pd.DataFrame(columns=['pathway_dataset', 'pathway_name'])

	return disease_df, phenotypes_df, differential_exp_diseases_df, differential_exp_tissues_df, gwas_df, pathways_df


def get_domains(record=None, gene_id=None, chembl_id=None):
	# ------------ domain and binding sites ----------#
	domain = []
	if not record:
		if gene_id:
			record = get_uniprot(gene_id)
			if record is None:
				return pd.DataFrame()
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
		connector = sqlite3.connect(chembl_24)
		query = "SELECT CD.start_position,CD.end_position,DOM.domain_name,DOM.source_domain_id,DOM.domain_type FROM target_dictionary " \
		        "TD,target_components TC,domains DOM,component_domains CD WHERE TD.chembl_id='%s' AND TD.tid=TC.tid AND " \
		        "TC.component_id=CD.component_id AND DOM.domain_id=CD.domain_id GROUP BY TD.chembl_id,CD.domain_id" % chembl_id
		domains_df = pd.read_sql(query, con=connector)
		connector.close()
		for i in domains_df.index:
			domain_id = str(record.accessions[0]) + str(domains_df.loc[i].start_position) + str(
				domains_df.loc[i].end_position) + '_' + domains_df.loc[i].source_domain_id
			domain.append({'Start': int(domains_df.loc[i].start_position), 'Stop': int(domains_df.loc[i].end_position),
			               'name': domains_df.loc[i].domain_name,
			               'length': int(domains_df.end_position.loc[i]) - int(domains_df.start_position.loc[i]),
			               'domain_id': domain_id, 'source_name': domains_df.domain_type.loc[i],
			               'Source_id': domains_df.source_domain_id.loc[i]})

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

	domain_df = pd.DataFrame.from_records(domain)
	return domain_df


def get_chembl_info(chembl_id):
	connector = sqlite3.connect(chembl_24)
	query = "SELECT PC.pref_name ,PC.short_name ,PC.protein_class_desc ,group_concat(DISTINCT(" \
	        "CSYN.component_synonym)) AS Synonym FROM target_dictionary TD, target_components TC, " \
	        "component_sequences CS, component_class CC, protein_classification PC, " \
	        "component_synonyms CSYN WHERE TD.chembl_id='%s' AND TD.tid=TC.tid AND " \
	        "TC.component_id=CS.component_id AND TC.component_id=CC.component_id AND " \
	        "TC.component_id=CSYN.component_id AND CC.protein_class_id=PC.protein_class_id GROUP BY " \
	        "TD.chembl_id" % chembl_id
	entry_info = pd.read_sql(query, con=connector)
	connector.close()
	entry_info.protein_class_desc = entry_info.protein_class_desc.str.replace('  ', ' -> ')
	entry_info.rename(
		columns={'short_name': 'Class_short', 'pref_name': 'Class_name', 'protein_class_desc': 'Class_description'},
		inplace=True)
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
	modifications = pd.DataFrame(
		columns=['mod_id', 'mod_type', 'start', 'stop', 'previous', 'new', 'action', 'comment', 'domains'])
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
			modifications.loc[feat[4]] = {'mod_id': feat[4], 'mod_type': str(feat[4]).split('_')[0], 'start': feat[1],
			                              'stop': feat[2],
			                              'previous': previous,
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
			modifications.loc['MUTAGEN_' + str(feat[1])] = {'mod_id': 'MUTAGEN_' + str(feat[1]), 'mod_type': 'MUTAGEN',
			                                                'start': feat[1], 'stop': feat[2], 'previous': previous,
			                                                'new': new, 'action': action, 'comment': comment,
			                                                'domains': []}
	for i in domains.index:
		for m in modifications.index:
			if modifications.loc[m]['start'] >= domains.Start.loc[i] and modifications.loc[m]['stop'] <= \
					domains.Stop.loc[i]:
				modifications.loc[m]['domains'].append(str(domains.name.loc[i]))
	for m in modifications.index:
		modifications.loc[m]['domains'] = ','.join(modifications.loc[m]['domains'])

	isoforms = pd.DataFrame(
		columns=['name', 'isoid', 'seq_mod', 'Canonical', 'sequence', 'seq_list', 'n_residues', 'Score', 'Identity',
		         'Gaps'])
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
				isoforms.loc[name + '_' + isoid] = {'name': name, 'isoid': isoid, 'seq_mod': seq_mod,'Canonical':'', 'sequence':'', 'seq_list':'','n_residues':'', 'Score':'', 'Identity':'','Gaps':''}
	if verbose:
		print("[ISOFORMS SEQUENCE ALIGNMENT IN PROGRESS]")
	for isoform in isoforms.index:
		temp_seq = list(record.sequence)
		if 'Displayed' in isoforms.loc[isoform]['seq_mod']:
			isoforms.loc[isoform]['Canonical'] = 1
			isoforms.loc[isoform]['sequence'] = record.sequence
			isoforms.loc[isoform]['seq_list'] = list(record.sequence)
			isoforms.loc[isoform]['n_residues'] = len(record.sequence)
		elif 'External' in isoforms.loc[isoform]['seq_mod']:
			isoforms.loc[isoform]['Canonical'] = 0
			isoforms.loc[isoform]['sequence'] = get_uniprot(isoforms.loc[isoform]['isoid']).sequence
			isoforms.loc[isoform]['seq_list'] = list(isoforms.loc[isoform]['sequence'])
			isoforms.loc[isoform]['n_residues'] = len(isoforms.loc[isoform]['sequence'])
		else:
			for changes in isoforms.loc[isoform]['seq_mod']:
				if changes == 'Notdescribed':
					temp_seq = ['no_sequence']
					break
				if modifications.loc[changes]['action'] == 'remove':
					start = modifications.loc[changes]['start'] - 1
					stop = modifications.loc[changes]['stop']
					for i in range(start, stop):
						temp_seq[i] = ''
				if modifications.loc[changes]['action'] == 'replace':
					start = modifications.loc[changes]['start'] - 1
					stop = modifications.loc[changes]['stop']
					temp_seq[start] = modifications.loc[changes]['new']
					for i in range(start + 1, stop):
						temp_seq[i] = ''

			isoforms.loc[isoform]['sequence'] = ''.join(temp_seq)
			isoforms.loc[isoform]['seq_list'] = temp_seq
			isoforms.loc[isoform]['Canonical'] = 0
			isoforms.loc[isoform]['n_residues'] = len(isoforms.loc[isoform]['sequence'])
		seq_object = Seq(isoforms.loc[isoform]['sequence'])
		alignment = align(record.sequence, seq_object)
		if alignment:
			isoforms.loc[isoform]['Score'] = alignment['score']
			isoforms.loc[isoform]['Identity'] = alignment['identity']
			isoforms.loc[isoform]['Gaps'] = alignment['gaps']
		else:
			isoforms.loc[isoform]['Score'] = ''
			isoforms.loc[isoform]['Identity'] = ''
			isoforms.loc[isoform]['Gaps'] = ''
	if verbose:
		print("[ISOFORMS SEQUENCE ALIGNMENT DONE]")
	return isoforms, modifications


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


@ret.retryer(max_retries=10, timeout=10)
def get_pdb(list_of_pdb, path):
	file_path = path.joinpath('PDB')
	if not file_path.is_dir():
		file_path.mkdir(parents=True,exist_ok=True)
	for k in list_of_pdb.index:
		saved_pdb = file_path.joinpath(str(list_of_pdb.loc[k]['PDB_code']) + '.pdb')
		if saved_pdb.is_file():
			if verbose:
				print('[FILE ALREADY THERE]: ',saved_pdb.name)
			list_of_pdb.at[k, 'path'] = str(saved_pdb)
		else:
			try:
				urllib.urlretrieve('http://files.rcsb.org/download/' + saved_pdb.name, str(saved_pdb))
				list_of_pdb.at[k, 'path'] = str(saved_pdb)
				if verbose:
					print('[PDB DOWNLOAD]: ' ,saved_pdb.name)
			except HTTPError:
				pass
				if verbose:
					print('[ERROR]: ' + str(k) + '.pdb file not found (do not exist or is '
					                             'too big)')
	if verbose:
		print('[PDB DOWNLOAD DONE]')


def get_pdb_seq_info(pdb_list, domain):
	if verbose:
		print("[PDB INFO]: Extracting PDB information")
	seq = []
	for keys in pdb_list.index:
		if pdb_list.path.loc[keys] == '':
			continue
		seq.append(pdb_parser.get_sequence(keys, pdb_list.loc[keys]['path'], pdb_list.loc[keys]['Chain'], domain))
	if verbose:
		print("[PDB INFO]: Done")
	if seq:
		seq_df = pd.concat(seq)
	else:
		seq_df = pd.DataFrame(columns=['PDB_code','chain_name', 'sequence', 'equal','start', 'stop','length','start_stop_pairs','domain','domain_id','seq_list'])
	return seq_df


def get_ligands_to_do(chembl_code):
	# ====================# GETTING THE LIST OF LIGAND WITH BIOACTIVITIES IN THE DB #=========================#
	lig_to_do = []
	if verbose:
		print("[LIGANDS]: Extracting ligand informations")
	# =====================# GETTING THE LIST OF LIGAND ASSOCIATED TO  THE TARGET #==========================#
	connector = sqlite3.connect(chembl_24)
	query_lig_target = "SELECT TD.chembl_id AS target_id,MOL.chembl_id AS lig_id,version.name as " \
	                   "chembl_version FROM target_dictionary TD,activities BIO,assays AC," \
	                   "molecule_dictionary MOL,version WHERE TD.chembl_id='%s' AND TD.tid=AC.tid AND " \
	                   "AC.assay_id=BIO.assay_id AND BIO.published_value is not null AND " \
	                   "BIO.molregno=MOL.molregno" % chembl_code

	res_lig_target = pd.read_sql(query_lig_target, con=connector)
	res_lig_target.drop_duplicates(['lig_id'], inplace=True)
	connector.close()
	if res_lig_target.empty:
		return lig_to_do

	count = 0
	limit = 500
	if len(res_lig_target) <= limit:
		limit = len(res_lig_target)
	n_loop = 0
	excess = len(res_lig_target) % limit
	connector2 = sqlite3.connect(targetDB)
	lig_ids = '('
	query_lig_bioact_db_base = """SELECT DISTINCT lig_id,chembl_version FROM bioactivities WHERE lig_id IN """
	res_lig_in_db = pd.DataFrame()
	for i in res_lig_target.index:
		count += 1
		lig_ids += "'" + res_lig_target.loc[i]['lig_id'] + "',"
		if count % limit == 0:
			lig_ids = lig_ids.rstrip(',') + ')'
			full_query = query_lig_bioact_db_base + lig_ids
			tmp_res = pd.read_sql(full_query, con=connector2)
			lig_ids = '('
			count = 0
			n_loop += 1
			res_lig_in_db = res_lig_in_db.append(tmp_res)
	if excess != 0:
		lig_ids = lig_ids.rstrip(',') + ')'
		full_query = query_lig_bioact_db_base + lig_ids
		tmp_res = pd.read_sql(full_query, con=connector2)
		res_lig_in_db = res_lig_in_db.append(tmp_res)
	res_lig_in_db.index = res_lig_in_db.lig_id
	connector2.close()

	# =========# ADDING IN THE TO-DO LIST ONLY LIGAND WITH NO BIOACTIVITY IN THE DB #========================#

	for i in res_lig_target.index:
		if res_lig_target.loc[i]['lig_id'] in res_lig_in_db.index:
			if int(str(res_lig_target.loc[i]['chembl_version']).split('_')[1]) <= int(
					str(res_lig_in_db.loc[res_lig_target.loc[i]['lig_id']].chembl_version).split('_')[1]):
				pass
			else:
				lig_to_do.append(res_lig_target.loc[i]['lig_id'])
		else:
			lig_to_do.append(res_lig_target.loc[i]['lig_id'])
	return lig_to_do


def get_assays():
	if verbose:
		print("[ASSAYS]: Extracting assays information")
	connector = sqlite3.connect(chembl_24)

	query = """SELECT
		  TD.chembl_id AS target_chembl_id,
		  AC.chembl_id AS assay_id,
		  AC.description AS assay_description,
		  t.relationship_desc AS relationship_desc,
		  DOC.doc_type AS ref_type,
		  DOC.doi AS doi,
		  DOC.patent_id AS patent_id,
		  AC.assay_organism AS species,
		  AT.assay_desc AS bioactivity_type,
		  AC.curated_by AS curated_by,
		  AC.confidence_score AS confidence_score,
		  lookup.target_mapping AS confidence_txt,
		  dictionary.cell_name AS cell_type,
		  dictionary.cell_source_organism AS cell_organism,
		  v.accession AS variant_uniprot,
		  v.mutation AS mutant,
		  VER.name AS chembl_version
		
		FROM
	    assays AC 
	      LEFT JOIN target_dictionary TD ON TD.tid = AC.tid
		  LEFT JOIN assay_type AT ON AC.assay_type = AT.assay_type
		  LEFT JOIN docs DOC ON AC.doc_id = DOC.doc_id
		  LEFT JOIN relationship_type t ON AC.relationship_type = t.relationship_type
		  LEFT JOIN confidence_score_lookup lookup ON AC.confidence_score = lookup.confidence_score
		  LEFT JOIN cell_dictionary dictionary ON AC.cell_id = dictionary.cell_id
		  LEFT JOIN variant_sequences v ON AC.variant_id = v.variant_id,
		    version VER"""
	res = pd.read_sql(query, con=connector)
	res.drop_duplicates(['assay_id'], inplace=True)
	connector.close()

	tdb_connector = sqlite3.connect(targetDB)
	res.to_sql('assays', con=tdb_connector, if_exists='append', index=False)
	tdb_connector.close()


def get_ligands_info():
	if verbose:
		print("[LIGANDS]: Extracting ligands information")
	connector = sqlite3.connect(chembl_24)
	sql_lig = """SELECT
	  MOL.pref_name AS mol_name,
	  MOL.chembl_id AS lig_id,
	  MOL.max_phase,
	  MOL.oral,
	  MOL.black_box_warning,
	  MOL.indication_class,
	  MOL.usan_stem_definition AS class_def,
	  PROP.alogp,
	  PROP.acd_logd,
	  PROP.acd_logp,
	  PROP.acd_most_bpka,
	  PROP.acd_most_apka,
	  PROP.hbd AS HBD,
	  PROP.hba AS HBA,
	  PROP.psa AS TPSA,
	  PROP.heavy_atoms AS n_heavy_atoms,
	  PROP.full_mwt AS molecularWeight,
	  PROP.rtb AS rotatableBonds,
	  PROP.aromatic_rings AS n_Ar_rings,
	  PROP.molecular_species,
	  PROP.num_ro5_violations,
	  PROP.ro3_pass,
	  PROP.full_molformula AS mol_formula,
	  STRUCT.canonical_smiles,
	  STRUCT.standard_inchi AS std_inchi,
	  STRUCT.standard_inchi_key AS std_inchi_key,
	  version.name AS chembl_version
	FROM
	  molecule_dictionary MOL
	     LEFT JOIN compound_structures STRUCT ON MOL.molregno = STRUCT.molregno
	     LEFT JOIN compound_properties PROP ON MOL.molregno = PROP.molregno
	  ,version"""
	df = pd.read_sql(sql_lig, con=connector)
	connector.close()

	tdb_con = sqlite3.connect(targetDB)
	df.to_sql('ligands', con=tdb_con, if_exists='append', index=False)
	tdb_con.close()


def get_bioactivity(lig_to_do):
	# =========# CREATING THE STRING OF LIGANDS ('CHEMBLXXXX','CHEMBLXXXX',....) #========================#
	if verbose:
		print("[BIOACTIVITIES]: Extracting bioactivities from chembl - CAN TAKE A FEW MINUTES ")

	connector = sqlite3.connect(chembl_24)

	# ========# FETCHING ALL BIOACTIVITIES (FOR THE LIST OF LIGANDS: ALL TARGETS) #=======================#
	res_bioactivity = pd.DataFrame()
	query = """SELECT EXP.chembl_id AS assay_id,
	  TD.chembl_id AS Target_id,
	  TD.target_type AS target_type,
	  TD.pref_name AS target_name,
	  TD.organism AS target_organism,
	  MOL.chembl_id AS lig_id,
	  ACT.standard_units AS units,
	  ACT.standard_relation AS operator,
	  ACT.standard_value AS value_num,
	  DOC.doi AS doi,
	  ((CASE WHEN ACT.standard_relation!='=' THEN ACT.standard_relation ELSE '' END)||ROUND(ACT.standard_value,2)||' '||
	  CASE WHEN ACT.standard_units IS NULL THEN '' ELSE ACT.standard_units END) AS value_text,
	  VER.name AS chembl_version,
	  ACT.standard_type AS standard_type,
	  ACT.activity_comment,
	  ACT.data_validity_comment,
	  ACT.pchembl_value
	FROM molecule_dictionary MOL
	  LEFT JOIN activities ACT ON MOL.molregno = ACT.molregno
	  LEFT JOIN assays EXP ON ACT.assay_id = EXP.assay_id
	  LEFT JOIN target_dictionary TD ON EXP.tid = TD.tid
	  LEFT JOIN docs DOC ON ACT.doc_id = DOC.doc_id,
	  version VER
	WHERE
	  MOL.chembl_id IN (%s)
	  AND ACT.published_value IS NOT NULL"""
	lig_str = ''
	counter = 0
	for i in lig_to_do:
		counter += 1
		lig_str += "'" + str(i) + "',"
		if counter == 1000:
			lig_str = lig_str.rstrip(',')
			query_bioactivity = query % lig_str
			res_bioactivity = res_bioactivity.append(pd.read_sql(query_bioactivity, con=connector))
			counter = 0
			lig_str = ''
	if counter == 0:
		pass
	elif counter < 1000:
		lig_str = lig_str.rstrip(',')
		query_bioactivity = query % lig_str
		res_bioactivity = res_bioactivity.append(pd.read_sql(query_bioactivity, con=connector))

	connector.close()
	return res_bioactivity


def blast_launcher(sequence, seq_file, db, output_name, num_core=8):
	if seq_file.is_file():
		subprocess.check_output(
			[blast_exe, '-db', str(db), '-query', str(seq_file), '-out', str(output_name),
			 '-num_threads', str(num_core), '-outfmt', str(5), '-max_target_seqs', str(100)],
			env={'BLASTDB': blast_db})
	else:
		seq_file.write_text(str(sequence))
		subprocess.check_output(
			[blast_exe, '-db', str(db), '-query', str(seq_file), '-out', str(output_name),
			 '-num_threads', str(num_core), '-outfmt', str(5), '-max_target_seqs', str(100)],
			env={'BLASTDB': blast_db})


def pdb_blast(sequence, path, gene_id, gene='', pdb_list=None):
	if pdb_list is None:
		pdb_list = []
	file_path = path.joinpath('PDB_BLAST')
	blast_file = file_path.joinpath(gene + '_' + gene_id + '.xml')
	seq_file = file_path.joinpath(gene + '_' + gene_id + '.seq')
	if not file_path.is_dir():
		file_path.mkdir(parents=True,exist_ok=True)
	if verbose:
		print('[3D BLAST]:' + gene + '(' + gene_id + ')')
	columns = ['PDB_code', 'seq_percent', 'similarity', 'Chain_id', 'chain_letter', 'gene', 'organism', 'uniprot_id',
	           'path']
	alternate_pdb = pd.DataFrame(columns=columns)
	if blast_file.is_file():
		if ((((time.time() - blast_file.stat().st_mtime) / 60) / 60) / 24) <= 15:
			if verbose:
				print('[3D BLAST FILE FOUND]:' + gene)
			result_handle = blast_file.open()
		else:
			if verbose:
				print('[3D BLAST FILE FOUND]: File older than 2 weeks, a new blast will be performed (' + gene + ')')
			blast_file.unlink()
			blast_launcher(sequence, seq_file, 'pdbaa', blast_file, num_core=bcore)
			if blast_file.is_file():
				result_handle = blast_file.open()
			else:
				print("[3D BLAST][ERROR]: Something went wrong, no blast result generated")
				return alternate_pdb
	else:
		blast_launcher(sequence, seq_file, 'pdbaa', blast_file, num_core=bcore)
		if blast_file.is_file():
			result_handle = blast_file.open()
		else:
			print("[3D BLAST][ERROR]: Something went wrong, no blast result generated")
			return alternate_pdb

	if verbose:
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
				data = {'PDB_code': pdb_code, 'seq_percent': percent_seq,'similarity': similarity,'Chain_id': chain_id, 'chain_letter': chain_id.split('_')[1], 'gene': '','organism': '', 'uniprot_id': '', 'path': ''}
				alternate_pdb = alternate_pdb.append(data,ignore_index=True)
	alternate_pdb.drop_duplicates(subset=['PDB_code', 'Chain_id'], inplace=True)
	get_pdb(alternate_pdb, path)
	for id in alternate_pdb.index:
		if Path(alternate_pdb.path.loc[id]).is_file():
			result = pdb_to_uniprot(alternate_pdb.loc[id]['chain_letter'], alternate_pdb.loc[id]['path'])
			if result:
				alternate_pdb.at[id, 'gene'] = result['gene']
				alternate_pdb.at[id, 'organism'] = result['organism']
				alternate_pdb.at[id, 'uniprot_id'] = result['uniprot_id']
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

	file_path = path.joinpath('PROTEIN_BLAST')
	blast_file = file_path.joinpath(gene + '_' + gene_id + '.xml')
	seq_file = file_path.joinpath(gene + '_' + gene_id + '.seq')

	if not file_path.is_dir():
		file_path.mkdir(parents=True,exist_ok=True)
	if verbose:
		print('[PROTEIN BLAST] ' + gene + '(' + gene_id + ')')
	if blast_file.is_file():
		if ((((time.time() - blast_file.stat().st_mtime) / 60) / 60) / 24) <= 15:
			if verbose:
				print('[PROTEIN BLAST FILE FOUND]:' + gene)
			result_handle = blast_file.open()
		else:
			if verbose:
				print(
					'[PROTEIN BLAST FILE FOUND]: File older than 2 weeks, a new blast will be performed (' + gene + ')')
			blast_file.unlink()
			blast_launcher(sequence, seq_file, 'swissprot', blast_file, num_core=bcore)
			if blast_file.is_file():
				result_handle = blast_file.open()
			else:
				print("[PROTEIN BLAST][ERROR]: Something went wrong, no blast result generated")
				return []
	else:
		blast_launcher(sequence, seq_file, 'swissprot', blast_file, num_core=bcore)
		if blast_file.is_file():
			result_handle = blast_file.open()
		else:
			print("[PROTEIN BLAST][ERROR]: Something went wrong, no blast result generated")
			return []
	blast_record = NCBIXML.read(result_handle)
	result_handle.close()
	if verbose:
		print('[PROTEIN BLAST DONE]: Now parsing the data - ' + gene + '(' + gene_id + ')')
	e_value_treshold = 0.0001
	query_length = float(blast_record.query_length)
	list_of_neighbours = pd.DataFrame(
		columns=['Query_target_id', 'Hit_gene_id', 'sequence_in_common', 'similarity', 'Hit_gene_name',
		         'Hit_gene_species', 'Unique_ID'])
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
			if percent_seq > 100:
				percent_seq = 100
			similarity = round((float(hsp.positives) / length) * 100, 1)
			if hsp.score > 200 and hsp.expect < e_value_treshold and (
					length / query_length) > 0.5 and similarity >= 40 and neighbour_accession_code != gene_id:
				list_of_neighbours.loc[gene_id + '_' + neighbour_accession_code] = {'Query_target_id': gene_id,
				                                                                    'Hit_gene_id': neighbour_accession_code,
				                                                                    'sequence_in_common': percent_seq,
				                                                                    'similarity': similarity,
				                                                                    'Hit_gene_name': neighbour_gene_name,
				                                                                    'Hit_gene_species': neighbour_gene_species,
				                                                                    'Unique_ID': gene_id + '_' + neighbour_accession_code}

	return list_of_neighbours


@ret.retryer_pubmed(max_retries=10, timeout=5)
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


@ret.retryer(max_retries=10, timeout=10)
def open_target_association(ensembl_id):
	ot = OpenTargetsClient()
	if ensembl_id == '':
		return pd.DataFrame(columns=['affected_pathway','animal_model','genetic_association', 'known_drug', 'litterature_mining', 'rna_expression', 'somatic_mutation', 'overall_score','disease_name','disease_area','gene_symbol'])
	associations = ot.get_associations_for_target(ensembl_id)

	df = associations.to_dataframe()

	if df.empty:
		return pd.DataFrame(columns=['affected_pathway', 'animal_model', 'genetic_association', 'known_drug', 'litterature_mining',
				         'rna_expression', 'somatic_mutation', 'overall_score', 'disease_name', 'disease_area',
				         'gene_symbol'])

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


def write_to_db(target, db_path):
	if target.record is None:
		return None

	# ========# OPEN DATABASE #========#
	if verbose:
		print("[DATABASE]: Start to write info into the database")
	connector = sqlite3.connect(db_path)

	# ========# FILLING THE TARGETS TABLE #=========#
	target.prot_info.to_sql('Targets', con=connector, index=False, if_exists='append')

	if verbose:
		print("[DATABASE]: Target table populated/updated")

	# ========# FILLING THE DOMAIN TABLE #=========#

	target.domain['Target_id'] = target.swissprotID
	target.domain.rename(
		columns={'domain_id': 'Domain_id', 'Start': 'Domain_start', 'Stop': 'Domain_stop', 'name': 'Domain_name'},
		inplace=True)
	target.domain.to_sql('Domain_targets', con=connector, index=False, if_exists='append')

	if verbose:
		print("[DATABASE]: Domain table populated/updated")

	# ========# FILLING THE 3D BLAST TABLE #=========#

	target.alternate_pdb['Query_target_id'] = target.swissprotID
	target.alternate_pdb.rename(
		columns={'PDB_code': 'Hit_PDB_code', 'Chain_id': 'Hit_Chain_id', 'chain_letter': 'Chain_Letter',
		         'gene': 'Hit_gene_name', 'organism': 'Hit_gene_species', 'uniprot_id': 'Hit_gene_id'}, inplace=True)
	target.alternate_pdb.drop(columns=['seq_percent', 'path'], inplace=True)
	target.alternate_pdb[
		'Unique_ID'] = target.alternate_pdb.Query_target_id + '_' + target.alternate_pdb.Hit_gene_id + '_' + target.alternate_pdb.Hit_Chain_id
	target.alternate_pdb.to_sql('3D_Blast', con=connector, index=False, if_exists='append')

	if verbose:
		print("[DATABASE]: PDB Blast table populated/updated")

	# ========# FILLING THE PROTEIN BLAST TABLE #=========#

	target.neighbours.to_sql('protein_blast', con=connector, index=False, if_exists='append')

	if verbose:
		print("[DATABASE]: Blast table populated/updated")

	# ========# FILLING THE MODIFICATIONS TABLE #=========#

	target.modifications['Target_id'] = target.swissprotID
	target.modifications['Unique_modID'] = target.swissprotID + '_' + target.modifications.mod_id
	target.modifications.to_sql('modifications', con=connector, index=False, if_exists='append')

	# ========# FILLING THE ISOFORMS + ISOFORMS MODIFICATIONS TABLES #=========#

	isoforms_mods = []
	target.isoforms.isoid = target.swissprotID+'_'+target.isoforms.isoid
	for i in target.isoforms.index:
		for j in target.isoforms.seq_mod.loc[i]:
			if j in ('Displayed', 'External'):
				continue
			else:
				isoforms_mods.append(
					{'mod_id': target.swissprotID + '_' + j, 'isoform_id': target.isoforms.isoid.loc[i]})
	isoform_mods = pd.DataFrame.from_records(isoforms_mods)
	target.isoforms['Target_id'] = target.swissprotID
	target.isoforms.rename(columns={'isoid': 'Isoform_id', 'name': 'Isoform_name', 'sequence': 'Sequence'},
	                       inplace=True)
	target.isoforms.drop(columns=['seq_list', 'seq_mod'], inplace=True)
	target.isoforms.to_sql('Isoforms', con=connector, if_exists='append', index=False)
	isoform_mods.to_sql('isoform_modifications', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Isoforms tables populated/updated")

	# ========# FILLING THE PDB RELATED TABLES #=========#

	pdb_in = pd.read_sql('SELECT PDB_code FROM PDB', con=connector)
	target.pdb.drop(pdb_in.PDB_code, errors='ignore', inplace=True)

	pdb_domains = pd.DataFrame(columns=['Chain_id', 'Domain_id'])

	target.pdb.drop(columns=['Chain', 'path', 'Domain'], inplace=True)
	target.pdb.to_sql('PDB', con=connector, if_exists='append', index=False)

	if not target.pdb_info.empty:

		target.pdb_info.start_stop_pairs = target.pdb_info.start_stop_pairs.apply(
			lambda x: ' | '.join(str(i[0]) + '-' + str(i[1]) for i in x))
		target.pdb_info['Target_id'] = target.swissprotID
		target.pdb_info.equal = target.pdb_info.equal.apply(lambda x: ','.join(list(set(x))))
		target.pdb_info.index = target.swissprotID+'_'+target.pdb_info.index
		for chain in target.pdb_info.index:
			for domain in target.pdb_info.domain_id.loc[chain]:
				pdb_domains.loc[chain] = {'Chain_id': chain, 'Domain_id': domain}
		target.pdb_info.reset_index(inplace=True)
		target.pdb_info.rename(columns={'index': 'Chain_id', 'chain_name': 'Chain', 'equal': 'equal_chains',
		                                'length': 'n_residues', 'sequence': 'Sequence', 'start': 'Start',
		                                'start_stop_pairs': 'start_stop', 'stop': 'Stop'}, inplace=True)
		target.pdb_info.drop(columns=['domain', 'domain_id', 'seq_list'], inplace=True)

		target.pdb_info.to_sql('PDB_Chains', con=connector, if_exists='append', index=False)

	pdb_domains.to_sql('PDBChain_Domain', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: PDB tables populated/updated")

	# ========# FILLING THE fPOCKET RELATED TABLES #=========#

	if not target.pockets['pockets_chain'].empty:
		target.pockets['pockets_chain'].List_of_contacts = target.pockets['pockets_chain'].List_of_contacts.apply(lambda x: ','.join(x))
		target.pockets['pockets_chain'].to_sql('fPockets_Chain', con=connector, if_exists='append', index=False)

	target.pockets['pockets'].to_sql('fPockets', con=connector, if_exists='append', index=False)
	target.pockets['pockets_domain'].to_sql('fPockets_Domain', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: fPockets tables populated/updated")

	target.alternate_pockets['pockets'].to_sql('fPockets', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: alternate fPockets tables populated/updated")

	# ========# FILLING THE BIOACTIVITIES TABLE #=========#

	target.bioactivities.to_sql('bioactivities', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Bioactivities table populated")

	# ========# FILLING THE PROTEIN EXPRESSION TABLE #=========#

	if target.protein_expression is not None:
		protein_level = pd.DataFrame.from_records(target.protein_expression.protein_lvl['Cells'])
		protein_selectivity = pd.DataFrame({'Selectivity_entropy': [target.protein_expression.selective_entropy],
		                                    'max_organ': [target.protein_expression.max_organ],
		                                    'Target_id': [target.swissprotID]})

		protein_level.cell_id = target.swissprotID + '_' + protein_level.cell_id
		protein_level['Target_id'] = target.swissprotID
		protein_level.rename(columns={'cell_id': 'Entry_id', 'level': 'value', 'cell_type': 'cell'}, inplace=True)

		protein_level.to_sql('protein_expression_levels', con=connector, if_exists='append', index=False)
		protein_selectivity.to_sql('protein_expression_selectivity', con=connector, if_exists='append', index=False)

		if verbose:
			print("[DATABASE]: Protein expression levels tables populated/updated")

	# ========# FILLING THE HUMAN MINE DATA TABLES #=========#

	target.disease['Target_id'] = target.swissprotID
	target.disease['Unique_id'] = target.swissprotID + '_' + target.disease.disease_id
	target.disease.rename(columns={'disease': 'disease_name'}, inplace=True)
	target.disease.to_sql('disease', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Disease table populated/updated")

	target.differential_exp_tissues['Target_id'] = target.swissprotID
	target.differential_exp_tissues.rename(columns={'T_statistic': 't_stat'}, inplace=True)
	target.differential_exp_tissues.to_sql('diff_exp_tissue', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Differential expression (tissues) table populated/updated")

	target.differential_exp_disease['Target_id'] = target.swissprotID
	target.differential_exp_disease.rename(columns={'T_statistic': 't_stat', 'Condition': 'disease'}, inplace=True)
	target.differential_exp_disease.to_sql('diff_exp_disease', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Differential expression (disease) table populated/updated")

	target.gwas['Target_id'] = target.swissprotID
	target.gwas.to_sql('gwas', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: GWAS table populated/updated")

	target.phenotypes['Target_id'] = target.swissprotID
	target.phenotypes.drop(columns=['gene'], inplace=True)
	target.phenotypes.to_sql('phenotype', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: phenotype table populated/updated")

	target.pathways['Target_id'] = target.swissprotID
	target.pathways.to_sql('pathways', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Pathway table populated/updated")

	# ========# FILLING THE CROSSREF TABLE #=========#

	crossref = pd.DataFrame(
		{'target_id': [target.swissprotID], 'Chembl_id': [target.chembl_id], 'hgnc_id': [target.hgnc_id],
		 'ensembl_id': [target.ensembl_id]})
	crossref.to_sql('Crossref', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Cross-references table populated/updated")

	# ========# FILLING THE CROSSREF TABLE #=========#

	target.open_targets['target_id'] = target.swissprotID
	target.open_targets.disease_area = target.open_targets.disease_area.apply(lambda x: ','.join(x))
	target.open_targets.drop(columns=['gene_symbol'],inplace=True)
	target.open_targets.to_sql('opentarget_association', con=connector, if_exists='append', index=False)

	if verbose:
		print("[DATABASE]: Open-targets table populated/updated")

	connector.close()


class Target:
	def __init__(self, gname, uniprot_id=None, ensembl_id=None, hgnc_id=None, v=False,
	             db_files_path=None, chembl=None, target_db=None, blast_cores=8):
		global verbose
		global chembl_24
		global targetDB
		global bcore
		bcore = blast_cores
		targetDB = target_db
		chembl_24 = chembl
		verbose = v
		script_start = time.time()
		# =====# INITIATING THE OBJECT #========#
		if '_' in gname:
			self.gene = gname.split('_')[0]
		else:
			self.gene = gname
		self.swissprotID = uniprot_id
		self.ensembl_id = ensembl_id
		self.hgnc_id = hgnc_id
		self.isoforms = ''
		self.modifications = ''
		self.pockets = {'pockets':pd.DataFrame(),'pockets_chain':pd.DataFrame(),'pockets_domain':pd.DataFrame()}
		self.record = None
		self.pubmed = None
		self.pathways = []
		self.disease = []
		self.gwas = []
		self.differential_exp_tissues = []
		self.differential_exp_disease = []
		self.phenotypes = []
		self.bioactivities = pd.DataFrame()
		self.alternate_pdb = pd.DataFrame()
		self.alternate_pockets = {'pockets':pd.DataFrame(),'pockets_chain':pd.DataFrame(),'pockets_domain':pd.DataFrame()}
		self.druggable_pockets = pd.DataFrame()
		self.pdb_info = pd.DataFrame()

		self.path = Path(db_files_path).joinpath('DB_files')
		if not self.path.is_dir():
			self.path.mkdir(parents=True,exist_ok=True)

		while True:
			# ==============# IF NO UNIPROT ID --> SKIP #====================#
			if self.swissprotID == '' or self.swissprotID is None:
				if verbose:
					print('[GENE SKIPPED]: No UniprotID found for ' + self.gene)
				break

			if verbose:
				print('[BEGINNING OF GENE]: ' + self.gene + ' (' + str(self.swissprotID) + ')')

			# ============# COLLECTING THE UNIPROT RECORD #===============#
			self.record = get_uniprot(self.swissprotID)
			if self.record is None:
				if verbose:
					print('[GENE SKIPPED]: No Uniprot data or wrong uniprotID for ' + self.gene)
				break

			# ===========# GET ALL THE CROSSREFERENCES (source: Uniprot)#==============#

			self.pdb, self.go, self.chembl_id = get_crossref_pdb_go_chembl(self.record)

			# ===========# GET PROTEIN EXPRESSION LEVELS (source: ProteinAtlas)#==============#

			self.protein_expression = patlas.ProteinExpression(self.gene, id=self.ensembl_id)
			if self.protein_expression.protein_lvl is None:
				self.protein_expression = None

			# ===========# GET INFO FROM HUMANMINE.ORG (disease, phenotypes, differential_exp_diseases,
			# differential_exp_tissues, gwas,pathways) #=============#

			self.disease, self.phenotypes, self.differential_exp_disease, self.differential_exp_tissues, self.gwas, self.pathways = get_humanmine_data(
				self.gene)

			# ==========# GET DOMAIN INFORMATION FROM BOTH CHEMBL AND UNIPROT #===========#

			self.domain = get_domains(record=self.record, chembl_id=self.chembl_id)

			# ==========# GET DISEASE ASSOCIATION (Source: OpenTargets)#===========#

			self.open_targets = open_target_association(self.ensembl_id)

			# ==========# GET SEQUENCE INFORMATION (Source: Uniprot)#===========#

			self.sequence = self.record.sequence
			self.seq_list = list(self.sequence)

			# ==========# GET ISOFORMS INFORMATION (Source: Uniprot) #===========#

			self.isoforms, self.modifications = get_variants(self.record, self.domain)

			# ==========# GET PROTEIN CLASS AND SYNONYMS FROM CHEMBL #===========#

			self.prot_info = get_chembl_info(self.chembl_id)
			if self.prot_info.empty:
				self.prot_info.loc[0] = [None, None, None, synonyms(self.record)]
			self.prot_info['Target_id'] = self.swissprotID
			self.prot_info.index = self.prot_info.Target_id
			self.prot_info.rename(
				columns={'Synonym': 'Synonyms', 'Class_name': 'Protein_class', 'Class_short': 'Protein_class_short',
				         'Class_description': 'Protein_class_desc'}, inplace=True)
			self.prot_info['Gene_name'] = self.gene
			self.prot_info['Species'] = self.record.organism
			self.prot_info['species_id'] = self.record.taxonomy_id[0]
			self.prot_info['Sequence'] = self.sequence
			self.prot_info['Cell_location'] = self.go.loc['C'].value
			self.prot_info['Process'] = self.go.loc['P'].value
			self.prot_info['Function'] = self.go.loc['F'].value
			self.prot_info['Number_isoforms'] = len(self.isoforms)
			self.prot_info['chembl_id'] = self.chembl_id

			# ============================================================================#
			# ======================# END OF THE INFO SECTION #===========================#
			# ============================================================================#

			# ============================================================================#
			# ======================# START OF THE PDB SECTION #==========================#
			# ============================================================================#

			# ======# RETRIEVING LIST OF PDB ASSOCIATED TO TARGET ALREADY IN DB #=========#

			connector = sqlite3.connect(targetDB)
			query_pdb = "SELECT PDB_code FROM PDB_Chains WHERE Target_id ='%s' GROUP BY PDB_code" % self.swissprotID
			query_pockets = "SELECT Pocket_id FROM fPockets WHERE Target_id ='%s' AND druggable='TRUE' AND blast='FALSE'" % self.swissprotID
			res_pdb = pd.read_sql(query_pdb, con=connector)
			res_pockets = pd.read_sql(query_pockets, con=connector)
			connector.close()

			# =============# REMOVING PDB CODE FROM CrossRef if already done #============#

			self.pdb.drop(res_pdb.PDB_code, inplace=True)

			if not self.pdb.empty:
				# ========================# DOWNLOADING THE PDBs #============================#

				get_pdb(self.pdb, self.path)

				# ========================# GET PDB SEQUENCE INFO #===========================#

				self.pdb_info = get_pdb_seq_info(self.pdb, self.domain)

				# =====================# GET POCKETS (source: fpockets) ======================#

				self.pockets = pocket.get_pockets(self.path, sphere_size=3, pdb_info=self.pdb, domain=self.domain,
				                                  uniprot_id=self.swissprotID, fpocket_exe=fpocket_exe,verbose=verbose)
				self.druggable_pockets = self.pockets['pockets'][self.pockets['pockets'].druggable == 'TRUE'].copy()

			pdb_super_list = list(res_pdb.PDB_code) + list(self.pdb.index)

			# ============================================================================#
			# ======================# END OF THE PDB SECTION #============================#
			# ============================================================================#

			# ============================================================================#
			# ====================# START OF THE LIGAND SECTION #=========================#
			# ============================================================================#

			# ==============# FILL THE ASSAY TABLE IF EMPTY #================#

			connector = sqlite3.connect(targetDB)
			if pd.read_sql("SELECT assay_id FROM assays", con=connector).empty:
				if verbose:
					print("[INFO]: Assay table is being filled, it may take a few minutes")
				get_assays()
				if verbose:
					print("[DATABASE]: Assay table populated")
			if pd.read_sql("SELECT lig_id FROM ligands", con=connector).empty:
				if verbose:
					print("[INFO]: Ligands table is being filled, it may take a few minutes")
				get_ligands_info()
				if verbose:
					print("[DATABASE]: Ligand table populated")
			connector.close()

			if self.chembl_id:
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
			very_close_pdb = self.alternate_pdb[self.alternate_pdb.similarity >= 90.0].copy()

			if res_pockets.empty and self.druggable_pockets.empty:
				self.alternate_pockets = pocket.get_pockets(self.path, sphere_size=3, alternate=True,
				                                            alternate_pdb=self.alternate_pdb,
				                                            uniprot_id=self.swissprotID, fpocket_exe=fpocket_exe,verbose=verbose)
			else:
				if not very_close_pdb.empty:
					self.alternate_pockets = pocket.get_pockets(self.path, sphere_size=3, alternate=True,
					                                            alternate_pdb=very_close_pdb,
					                                            uniprot_id=self.swissprotID, fpocket_exe=fpocket_exe,verbose=verbose)

			self.neighbours = proteins_blast(self.sequence, self.swissprotID, self.gene, self.path)

			# ============================================================================#
			# ====================# END OF THE BLAST SECTION #============================#
			# ============================================================================#

			# ======================# WRITING TO THE DATABASE #===========================#
			write_to_db(self, targetDB)
			break
		# ============================================================================#
		# =================# END OF DATA GATHERING SECTION #==========================#
		# ============================================================================#



		script_stop = time.time()
		if verbose:
			print('[END OF GENE]: ' + self.gene + ' (in ' + str(round(script_stop - script_start)) + ' sec)')
			print('=======================================================================')


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gene', help='enter a single gene name', metavar='')
	parser.add_argument('-i', '--in_file', help='Name of the input file with a list of genes (.txt - 1 gene per line)',
	                    metavar='')
	parser.add_argument('-l', '--list_genes', help='Enter a list of gene name separated by a ","', metavar='')
	parser.add_argument('-s', '--sphere_size', help='enter a value for the probe size of the pocket finder tool ('
	                                                'default = 3.0)', metavar='', type=float, default=3.0)
	parser.add_argument('-v', '--verbose', help="Print information", action='store_true', default=False)

	parser.add_argument('-update', '--update', help="Update record if already in database (default: No)",
	                    action='store_true', default=False)
	parser.add_argument('-blastcore', '--num_core',
	                    help='Enter the value of processor core to use for the blast search (default=8)', metavar='',
	                    type=int, default=8)
	parser.add_argument('-update_config', '--update_config', help="use this if you want to update the config file",
	                    action='store_true', default=False)
	parser.add_argument('-create_db', '--create_db', help='Use this option to create an empty targetDB database',
	                    action='store_true', default=False)

	arguments = parser.parse_args()
	if not arguments.gene and not arguments.in_file and not arguments.list_genes:
		print('[ERROR]: Please use one of the optional input options : -g / -i / -l ')
		parser.print_help()
		sys.exit()
	return arguments


def main():
	global args, targetDB, chembl_24, pubmed_email, dbase_file_path, entries_in_db, blast_exe, blast_db, fpocket_exe
	args = parse_args()
	update_config = args.update_config
	if args.create_db:
		createdDB_path = tinit.create_db()
	while True:
		config = configparser.ConfigParser()
		config_file_path = Path('~/.druggability/config.ini').expanduser()
		config_file_path.parent.mkdir(exist_ok=True, parents=True)

		if config_file_path.is_file() and not update_config:
			config.read(str(config_file_path))
			if args.create_db:
				config['database_path']['targetdb'] = str(createdDB_path)
			todo = []
			for var_name in config['database_path']:
				if not Path(config['database_path'][var_name]).is_file():
					todo.append(var_name)
			if not Path(config['executables']['fpocket']).is_file() and os.name != 'nt':
				todo.append('fpocket')
			if not Path(config['executables']['blast']).is_file() or not Path(config['executables']['blastdb_path']).is_dir():
				todo.append('blast')
			if not Path(config['output_path']['db_files']).is_dir() or config['output_path']['db_files'] == '':
				todo.append('db_files')
			if not cf.is_email(config['pubmed_email']['email']):
				todo.append('email')
			if todo:
				config = cf.get_config_from_user(config, todo=todo)
				with config_file_path.open(mode='w') as cfile:
					config.write(cfile)
			else:
				if args.create_db:
					with config_file_path.open(mode='w') as cfile:
						config.write(cfile)
				# =============================# PATH TO SQLITE DB #============================#

				targetDB = config['database_path']['targetdb']
				chembl_24 = config['database_path']['chembl']
				pubmed_email = config['pubmed_email']['email']
				dbase_file_path = config['output_path']['db_files']
				blast_exe = config['executables']['blast']
				blast_db = config['executables']['blastdb_path']
				fpocket_exe = config['executables']['fpocket']
				break
		else:
			todo = ['targetdb', 'chembl', 'email', 'db_files', 'blast','fpocket']
			if args.create_db:
				todo.remove('targetdb')
			if os.name == 'nt':
				todo.remove('fpocket')

			config = cf.get_config_from_user(config,todo=todo, new=True)
			if args.create_db:
				config['database_path']['targetdb'] = str(createdDB_path)
			with config_file_path.open(mode='w') as cfile:
				config.write(cfile)
			update_config = False

	while True:
		if args.in_file:
			if Path(args.in_file).is_file():
				with open(args.in_file, 'r') as gene_list:
					list_of_genes = gene_list.readlines()
					gene_df = g2id.gene_to_id(list_of_genes, targetDB_path=targetDB)
				break
			else:
				print('ERROR : file inputed as argument [-i] does not exist')
				sys.exit()
		if args.gene:
			list_of_genes = [args.gene]
			gene_df = g2id.gene_to_id(list_of_genes, targetDB_path=targetDB)
			break
		elif args.list_genes:
			list_of_genes = args.list_genes.split(',')
			gene_df = g2id.gene_to_id(list_of_genes, targetDB_path=targetDB)
			break

	for g_id in gene_df.index:
		if len(gene_df.uniprot_ids.loc[g_id]) == 0:
			if args.verbose:
				print('[GENE SKIPPED]: No uniprot id was found for the entered gene name: ',gene_df.symbol.loc[g_id])
		for uniprot_id in gene_df.uniprot_ids.loc[g_id]:
			entries_in_db = get_list_entries(target_db_path=targetDB)
			if uniprot_id in entries_in_db.index:
				if args.update:
					if args.verbose:
						print("[GENE INFO]: Target already present in database, updating informations")
					activate_fk_sql = """PRAGMA foreign_keys = 1"""
					delete_sql = """DELETE FROM Targets WHERE Target_id='%s'""" % uniprot_id
					tdb = sqlite3.connect(targetDB)
					tdb.execute(activate_fk_sql)
					tdb.execute(delete_sql)
					tdb.commit()
					tdb.close()

					Target(gene_df.symbol.loc[g_id], uniprot_id=uniprot_id,
					       ensembl_id=gene_df.ensembl_gene_id.loc[g_id], hgnc_id=g_id,db_files_path=dbase_file_path, chembl=chembl_24, target_db=targetDB,
					       blast_cores=args.num_core, v=args.verbose)
				else:
					if args.verbose:
						print('[GENE SKIPPED]: Already present in the database: ' + gene_df.symbol.loc[g_id])
						print('=======================================================================')

			else:
				Target(gene_df.symbol.loc[g_id], uniprot_id=uniprot_id,
				       ensembl_id=gene_df.ensembl_gene_id.loc[g_id], hgnc_id=g_id, db_files_path=dbase_file_path, chembl=chembl_24, target_db=targetDB,
				       blast_cores=args.num_core, v=args.verbose)


def entry_point():
	main()


if __name__ == "__main__":
	entry_point()
