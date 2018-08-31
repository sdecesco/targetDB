#!/usr/bin/env python

import sqlite3
import pandas as pd


def gene_to_id(list_of_genes,targetDB_path=None):
	print('[NAME CONVERSION]: Converting gene names into IDs (Uniprot,Ensembl,HGNC)')
	connector = sqlite3.connect(targetDB_path)
	gene_list = []
	required_columns = ['alias_symbol', 'ensembl_gene_id', 'prev_symbol', 'symbol','uniprot_ids']
	for gene in list_of_genes:
		gene = gene.rstrip('\n')
		gene = gene.rstrip('\r')
		gene_list.append(gene)
	gene_ids = "','".join(gene_list)
	gene_id_query = """SELECT * FROM hgnc as hgn WHERE hgn.hgnc_id in (SELECT hg.hgnc_id FROM hgnc as hg WHERE hg.xref_value in ('%s'))""" % gene_ids
	gene_xref = pd.read_sql(gene_id_query, con=connector)
	xref_df = gene_xref[gene_xref.xref_name.isin(['symbol', 'uniprot_ids', 'ensembl_gene_id', 'prev_symbol', 'alias_symbol'])]
	xref_piv = xref_df.pivot_table(index='hgnc_id', values='xref_value', columns='xref_name',
	                               aggfunc=lambda x: ';'.join(x), fill_value='')
	for col in required_columns:
		if col not in xref_piv.columns:
			xref_piv[col] = ''
	xref_piv.uniprot_ids = xref_piv.uniprot_ids.apply(lambda x: list(filter(None, x.split(';'))))
	connector.close()
	list_of_ids = []
	for gene in gene_list:
		gene_id = xref_piv[xref_piv.symbol == gene].index.tolist()
		if len(gene_id) == 0:
			gene_id = xref_df[(xref_df.xref_value == gene) & (
						(xref_df.xref_name == 'prev_symbol') | (xref_df.xref_name == 'alias_symbol'))][
				'hgnc_id'].tolist()
		if len(gene_id) != 0:
			list_of_ids.extend(gene_id)
	output = xref_piv.loc[list_of_ids]
	print('[NAME CONVERSION]: Conversion Done')
	output = output.reset_index()
	output.drop_duplicates(subset=['hgnc_id','ensembl_gene_id','symbol'],inplace=True)
	output.index = output['hgnc_id']
	return output