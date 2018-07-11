#!/usr/bin/env python

try:
	import urllib.request as urllib
	from urllib.error import *
except ImportError:
	import urllib
	from urllib2 import *

import argparse, sys, os, sqlite3, re,requests,time
from Bio import Entrez, Medline
from druggability3 import target_descriptors as td
from druggability3 import target_features as tf
from druggability3 import cns_mpo as mpo
import pandas as pd
import pkg_resources as pkg

# =============================# PATH TO SAVE REPORT FILES #============================#

output_lists_path = '/data/sdecesco/databases/druggability/outputs/lists/'
output_single_path = '/data/sdecesco/databases/druggability/outputs/single_targets/'

# =============================# PATH TO SQLITE DB #============================#

targetDB = pkg.resource_filename(__name__, 'data/TargetDB_v1.db')
chembl_24 = pkg.resource_filename(__name__, 'data/chembl_24.db')
tcrd = pkg.resource_filename(__name__, 'data/tcrd_v5.2.0.db')


def get_list_entries():
	connector = sqlite3.connect(targetDB)
	query = "SELECT Target_id,Gene_name,Synonyms FROM Targets"
	entries_list = pd.read_sql(query, con=connector)
	connector.close()
	return entries_list


def gene_to_uniprotid_local(list_of_genes):
	connector = sqlite3.connect(targetDB)
	gene_list = []
	for gene in list_of_genes:
		gene = gene.rstrip('\n')
		gene = gene.rstrip('\r')
		gene_list.append(gene)
	gene_ids = "','".join(gene_list)
	gene_id_query = """SELECT * FROM hgnc as hgn WHERE hgn.hgnc_id in (SELECT hg.hgnc_id FROM hgnc as hg WHERE hg.xref_value in ('%s'))""" % gene_ids
	gene_xref = pd.read_sql(gene_id_query, con=connector)
	connector.close()
	output = {}
	for gene in gene_list:
		gene_id = gene_xref[(gene_xref.xref_value == gene) & (gene_xref.xref_name == 'symbol')]['hgnc_id'].values
		if gene_id.size == 0:
			gene_id = gene_xref[(gene_xref.xref_value == gene) & (
					(gene_xref.xref_name == 'prev_symbol') | (gene_xref.xref_name == 'alias_symbol'))][
				'hgnc_id'].values
		if gene_id.size == 0:
			output[gene] = "No Match"
			continue
		elif gene_id.size > 1:
			for g_id in gene_id:
				gene_name = \
					gene_xref[(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'symbol')].xref_value.values[0]
				gene_uniprot = gene_xref[
					(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'uniprot_ids')].xref_value.values
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
			gene_name = \
				gene_xref[(gene_xref.hgnc_id == gene_id[0]) & (gene_xref.xref_name == 'symbol')].xref_value.values[0]
			gene_uniprot = gene_xref[
				(gene_xref.hgnc_id == gene_id[0]) & (gene_xref.xref_name == 'uniprot_ids')].xref_value.values
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


class NetworkError(RuntimeError):
	pass


def retryer(max_retries=10, timeout=5):
	def wraps(func):
		request_exceptions = (
			requests.exceptions.Timeout,
			requests.exceptions.ConnectionError,
			requests.exceptions.HTTPError, urllib.URLError, urllib.HTTPError
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
			requests.exceptions.HTTPError, urllib.URLError, urllib.HTTPError
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


def write_excel_header(header_dict, worksheet, format):
	for head in header_dict.keys():
		if len(header_dict[head]) == 2:
			row, col = header_dict[head]
			worksheet.write(row, col, head, format)
		elif len(header_dict[head]) == 4:
			row, col, last_row, last_col = header_dict[head]
			worksheet.merge_range(row, col, last_row, last_col, head, format)


def get_single_excel(target_id):
	if (list_of_entries.Target_id == target_id).any():
		writer = pd.ExcelWriter(output_single_path + list_of_entries[list_of_entries.Target_id == target_id].to_string(
			columns=['Gene_name'], header=False, index=False) + '_' + target_id + '.xlsx',
		                        engine='xlsxwriter', options={'nan_inf_to_errors': True})

		workbook = writer.book

		# ============================ STYLES ===============================#

		bold_center = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center'})
		red = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'red'})
		green = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'green'})
		col_header = workbook.add_format({'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'vcenter'})
		vert_col_header = workbook.add_format(
			{'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'bottom', 'rotation': 90})

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

		res = tf.get_single_features(target_id,dbase=targetDB)

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
			target_desc = td.get_descriptors(target_id, targetdb=targetDB,tcrdDB=tcrd)
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
					row += 1
					to_print = res['isoforms_mod'][res['isoforms_mod']['isoform_id'] == iso['Isoform_id']]
					to_print = to_print[
						['start', 'stop', 'previous_seq', 'modification_type', 'new_seq', 'in_domains', 'comments']]
					to_print.to_excel(writer, sheet_name='isoforms', startrow=row, index=False, header=False)
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
				writer.sheets['Emax_Efficacy'].write(row_efficacy + 1, col_num, value, vert_col_header)
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
			           'Source', 'DOI', 'Patent_number', 'Institution', 'ligand_name', 'SMILES', 'HBA', 'HBD', 'LogD',
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
		return print("Gene with ID [", target_id, '] not present in the database. Run the druggability_DB command '
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


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gene', help='enter a single gene name', metavar='')
	parser.add_argument('-i', '--in_file',
	                    help='Name of the input file with a list of genes (.txt - 1 gene per line)',
	                    metavar='')
	parser.add_argument('-U', '--Uniprot_ID',
	                    help='Enter the uniprot accession number followed by the name separated '
	                         'by a comma(eg. Q9NZC2,TREM2)', metavar='')
	parser.add_argument('-l', '--list_genes', help='Enter a list of gene name separated by a ","', metavar='')

	parser.add_argument('-v', '--verbose', help="Print information", action='store_true', default=False)

	parser.add_argument('-rl', '--report_list',
	                    help="produce an output file at the end of the analysis with a list "
	                         "of genes entered", action='store_true',
	                    default=False)
	parser.add_argument('-rs', '--report_single',
	                    help="produce an output file for each target listed (more detailed "
	                         "information available) Caution : 1 File created per target",
	                    action='store_true',
	                    default=False)

	parser.add_argument('-pubmed_email', '--email',
	                    help='enter your email address (required for the pubmed search', metavar='',
	                    type=str)
	args = parser.parse_args()

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

	list_of_entries = get_list_entries()

	if args.report_single:
		for gene_id in gene_dict.values():
			get_single_excel(gene_id)
	elif args.report_list:
		get_list_excel(gene_dict.keys())
	else:
		print('ERROR: Please provide a report type "-rl" for a list report, "-rs" for a single target report')
		parser.print_help()
		sys.exit()
