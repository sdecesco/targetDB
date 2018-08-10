#!/usr/bin/env python

import sqlite3
import pandas as pd


def gene_to_id(list_of_genes,targetDB_path=None):
	connector = sqlite3.connect(targetDB_path)
	gene_list = []
	for gene in list_of_genes:
		gene = gene.rstrip('\n')
		gene = gene.rstrip('\r')
		gene_list.append(gene)
	gene_ids = "','".join(gene_list)
	gene_id_query = """SELECT * FROM hgnc as hgn WHERE hgn.hgnc_id in (SELECT hg.hgnc_id FROM hgnc as hg WHERE hg.xref_value in ('%s'))""" % gene_ids
	gene_xref = pd.read_sql(gene_id_query, con=connector)
	connector.close()
	output = pd.DataFrame(columns=['uniprot_id','hgnc_id','ensembl_id'])
	for gene in gene_list:
		gene_id = gene_xref[(gene_xref.xref_value == gene) & (gene_xref.xref_name == 'symbol')]['hgnc_id']
		if gene_id.size == 0:
			gene_id = gene_xref[(gene_xref.xref_value == gene) & ((gene_xref.xref_name == 'prev_symbol') | (gene_xref.xref_name == 'alias_symbol'))]['hgnc_id']
		if gene_id.size == 0:
			output.loc[gene] = [None, None, None]
			continue
		elif gene_id.size > 1:
			for g_id in gene_id:
				gene_name = gene_xref.xref_value.loc[gene_xref[(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'symbol')].first_valid_index()]
				gene_uniprot = gene_xref[(gene_xref.hgnc_id == g_id) & (gene_xref.xref_name == 'uniprot_ids')].xref_value.values
				gene_ensembl = gene_xref.at[gene_xref[(gene_xref.hgnc_id == g_id) &(gene_xref.xref_name=='ensembl_gene_id')].first_valid_index(),'xref_value']
				if gene_uniprot.size > 1:
					count = 0
					for uniprot_id in gene_uniprot:
						name = gene_name + '_' + str(count)
						output.loc[name] = [uniprot_id,g_id,gene_ensembl]
						count += 1
				elif gene_uniprot.size == 0:
					output.loc[gene_name] = [None, g_id, gene_ensembl]
				else:
					output.loc[gene_name] = [gene_uniprot.loc[gene_uniprot.first_valid_index()], g_id, gene_ensembl]
		else:
			gene_name = gene_xref.xref_value.loc[gene_xref[(gene_xref.hgnc_id == gene_id.loc[gene_id.first_valid_index()]) & (gene_xref.xref_name == 'symbol')].first_valid_index()]
			gene_uniprot = gene_xref[(gene_xref.hgnc_id == gene_id.loc[gene_id.first_valid_index()]) & (gene_xref.xref_name == 'uniprot_ids')].xref_value
			gene_ensembl = gene_xref.at[gene_xref[(gene_xref.hgnc_id == gene_id.loc[gene_id.first_valid_index()]) & (gene_xref.xref_name == 'ensembl_gene_id')].first_valid_index(), 'xref_value']
			if gene_uniprot.size > 1:
				count = 0
				for uniprot_id in gene_uniprot:
					name = gene_name + '_' + str(count)
					output.loc[name] = [uniprot_id, gene_id.loc[gene_id.first_valid_index()], gene_ensembl]
					count += 1
			elif gene_uniprot.size == 0:
				output.loc[gene_name] = [None, gene_id.loc[gene_id.first_valid_index()], gene_ensembl]
			else:
				output.loc[gene_name] = [gene_uniprot.loc[gene_uniprot.first_valid_index()], gene_id.loc[gene_id.first_valid_index()], gene_ensembl]
	return output