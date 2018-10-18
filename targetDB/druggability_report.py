#!/usr/bin/env python

import argparse
import configparser
import re
import sqlite3
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, Medline

from targetDB import cns_mpo as mpo
from targetDB import target_descriptors as td
from targetDB import target_features as tf
from targetDB.utils import config as cf
from targetDB.utils import retryers as ret
from targetDB.utils import gene2id as g2id
from targetDB.utils import druggability_ml as dml

ml_model = dml.generate_model()


def get_list_entries():
    connector = sqlite3.connect(targetDB)
    query = "SELECT Target_id,Gene_name FROM Targets"
    entries_list = pd.read_sql(query, con=connector, index_col='Target_id')
    connector.close()
    return entries_list


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


def write_excel_header(header_dict, worksheet, format):
    for head in header_dict.keys():
        if len(header_dict[head]) == 2:
            row, col = header_dict[head]
            worksheet.write(row, col, head, format)
        elif len(header_dict[head]) == 4:
            row, col, last_row, last_col = header_dict[head]
            worksheet.merge_range(row, col, last_row, last_col, head, format)


def get_single_excel(target):
    for uniprot_id in target.uniprot_ids:
        if uniprot_id in list_of_entries.index:
            output_name = Path(output_single_path).joinpath(target.symbol + '_' + uniprot_id + '.xlsx')
            writer = pd.ExcelWriter(str(output_name), engine='xlsxwriter', options={'nan_inf_to_errors': True})

            workbook = writer.book

            # ============================ STYLES ===============================#

            bold_center = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center'})
            red = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'red'})
            green = workbook.add_format({'bold': True, 'valign': 'vcenter', 'color': 'green'})
            col_header = workbook.add_format(
                {'bold': True, 'bg_color': '#D9D9D9', 'align': 'center', 'valign': 'vcenter'})
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
            hv_center = workbook.add_format({'valign': 'vcenter','align': 'center'})
            left_v_center = workbook.add_format({'valign': 'vcenter', 'align': 'left'})
            link = workbook.add_format(
                {'bold': True, 'valign': 'vcenter', 'align': 'center', 'color': 'blue', 'underline': True})

            # ================= DIFFERENT TAB CREATION ==========================#

            wb_general_info = workbook.add_worksheet('General info')
            writer.sheets['General info'] = wb_general_info
            wb_references = workbook.add_worksheet('Pubmed_search')
            writer.sheets['Pubmed_search'] = wb_references
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

            res = tf.get_single_features(uniprot_id, dbase=targetDB)

            # =================== FILLING THE WORKSHEETS ========================#

            # =============== GETTING PUBMED DATAFRAMES ======================#
            sequence = None
            pubmed = pd.DataFrame(data=None)

            if not res['general_info'].empty:
                if pubmed_email:
                    pubmed = pubmed_search(target.symbol, pubmed_email)

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
            if not res['open_target'].empty:
                res['open_target'].to_excel(writer, sheet_name='open_target_association', index=False)
                for col_num, value in enumerate(res['open_target'].columns.values):
                    writer.sheets['open_target_association'].write(0, col_num, value, vert_col_header)

            # ============================ GENERAL TAB ===============================#

            # GENERAL INFO HEADER WRITING
            header_index = {'Gene_name': (0, 0), 'Synonyms': (1, 0), 'Target_id': (2, 0), 'Protein_class': (3, 0),
                            'Protein_class_desc': (4, 0), 'Species': (5, 0), 'Number_isoforms': (6, 0),
                            'Tractable': (0, 3), 'Tractability_probability': (0, 4),
                            'In Training data set ?': (2, 3, 2, 4), 'DISEASE': (8, 0, 8, 1), 'disease_id': (9, 0),
                            'disease_name': (9, 1), 'PATHWAYS': (8, 3, 8, 4), 'Reactome': (9, 3), 'KEGG': (9, 4)}
            write_excel_header(header_index, wb_general_info, col_header)

            if not res['general_info'].empty:
                sequence = ''
                for k, v in res['general_info'].iloc[0].items():
                    if k in header_index:
                        row, col = header_index[k]
                        col = col + 1
                        wb_general_info.write(row, col, v, left_v_center)

                target_desc = td.get_descriptors_list(uniprot_id, targetdb=targetDB)
                tscore = td.target_scores(target_desc, mode='single')
                druggability_pred = dml.predict(ml_model, tscore.score_components)
                drug_proba = pd.DataFrame(dml.predict_prob(ml_model, tscore.score_components),
                                          columns=ml_model.classes_)
                tscore.scores['Tractable'] = druggability_pred
                tscore.scores['Tractability_probability'] = round(drug_proba[1] * 100, 2)
                tscore.scores['Tractable'] = tscore.scores['Tractable'].replace({0: 'False', 1: 'True'})
                tscore.scores['In_training_set'] = dml.in_training_set(tscore.score_components)
                target_desc = target_desc.merge(tscore.scores, on='Target_id', how='left')
                score_col = ['structure_info_score', 'structural_drug_score', 'chemistry_score', 'biology_score',
                             'disease_score', 'genetic_score', 'information_score', 'safety_score']
                target_score = target_desc[score_col] * 10
                target_score.index = target_desc.Target_id
                target_score.fillna(0, inplace=True)
                target_score = target_score.rename(columns={'structure_info_score': 'Structural information',
                                                            'structural_drug_score': 'Structural Druggability',
                                                            'chemistry_score': 'Chemistry', 'biology_score': 'Biology',
                                                            'disease_score': 'Diseases',
                                                            'genetic_score': 'Genetic Association',
                                                            'information_score': 'Literature',
                                                            'safety_score': 'Safety'})
                spider_plot = td.make_spider_plot(target_score.loc[uniprot_id].values, target_score.columns,
                                                  target_name=res['general_info'].iloc[0]['Gene_name'])
                wb_general_info.insert_image('G1', 'spider_plot', {'image_data': spider_plot})

                wb_general_info.write(1, 3, target_desc['Tractable'].iloc[0], hv_center)
                wb_general_info.write(1, 4, target_desc['Tractability_probability'].iloc[0], hv_center)
                wb_general_info.merge_range(3, 3, 3, 4, target_desc['In_training_set'].iloc[0], hv_center)

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
                                              {'type': 'icon_set', 'reverse_icons': True,
                                               'icon_style': '3_traffic_lights'})

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
                        wb_expression.merge_range(row, col, row, col + 3,
                                                  res['tissue_expression'].iloc[i]['organ'].upper(),
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
                             'n_residues', '% of full protein', 'start_stop', 'type_of_binder', 'binding_type',
                             'operator',
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
                res['alt_pockets'].to_excel(writer, sheet_name='Pockets', startrow=1, index=False,
                                            startcol=col_alt_pocket)
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
                row_efficacy = len(res['efficacy_bio']) + len(res['efficacy_bio'].columns) + 1
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
                columns = ['ZincID', 'IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)', 'kon(M-1s-1)', 'koff(s-1)', 'pH',
                           'Temp',
                           'Source', 'DOI', 'Patent_number', 'Institution', 'ligand_name', 'SMILES', 'HBA', 'HBD',
                           'LogD',
                           'LogP', 'MW', 'TPSA', 'aLogP', 'apKa', 'bpKa', 'nAr', 'pass_ro3', 'ro5_violations',
                           'rotB', 'CNS_MPO', 'mol_name', 'molecular_species', 'indication_class', 'class_def',
                           'max_phase',
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
            print('[FILE]: File for ', target.symbol, ' generated [', output_name, ']')
        else:
            return print("Gene named ", target.symbol, " with ID [", uniprot_id,
                         '] not present in the database. Run the druggability_DB command to run the analysis on it')


def get_list_excel(list_targets):
    not_in_db = {'Not present in DB': []}

    if list_targets.empty:
        return print("No genes that you entered are in the Database")

    if pubmed_email:
        pubmed = {'Target_id': [], 'total # publications': [], 'number of Dementia publications': []}
        for gene_symbol in list_targets.index:
            pubmed['Target_id'].extend(list_targets.uniprot_ids.loc[gene_symbol])
            pubmed['number of Dementia publications'].extend([pubmed_search(list_targets.symbol.loc[gene_symbol],
                                                                            pubmed_email, return_number=True,
                                                                            mesh_term='Dementia')] * len(
                list_targets.uniprot_ids.loc[gene_symbol]))
            pubmed['total # publications'].extend(
                [pubmed_search(list_targets.symbol.loc[gene_symbol], pubmed_email, return_number=True)] * len(
                    list_targets.uniprot_ids.loc[gene_symbol]))
        # TODO: Mesh term is hard-coded here
        pubmed = pd.DataFrame.from_dict(pubmed)
    else:
        pubmed = pd.DataFrame(columns=['Target_id', 'total # publications', 'number of Dementia publications'])

    gene_ids = "','".join(list_targets.uniprot_ids.astype(str))
    gene_ids = gene_ids.replace('[\'', '').replace('\']', '')

    data = td.get_descriptors_list(gene_ids, targetdb=targetDB)
    tscore = td.target_scores(data)
    druggability_pred = dml.predict(ml_model, tscore.score_components)
    drug_proba = pd.DataFrame(dml.predict_prob(ml_model, tscore.score_components), columns=ml_model.classes_)
    tscore.scores['Tractable'] = druggability_pred
    tscore.scores['Tractability_probability'] = round(drug_proba[1] * 100, 2)
    tscore.scores['Tractable'] = tscore.scores['Tractable'].replace({0: False, 1: True})
    tscore.scores['In_training_set'] = dml.in_training_set(tscore.score_components)
    data = data.merge(tscore.scores, on='Target_id', how='left')
    data = data.merge(pubmed, on='Target_id', how='left')
    list_done = data.Target_id.values.tolist()

    for gene_symbol in list_targets.index:
        for tid in list_targets.uniprot_ids.loc[gene_symbol]:
            if tid not in list_done:
                not_in_db['Not present in DB'].append(list_targets.symbol.loc[gene_symbol])
    not_in_db = pd.DataFrame.from_dict(not_in_db)

    t = time.strftime("%d%b%Y_%H%M%S")
    output_file_name = Path(output_lists_path).joinpath('Export_' + str(len(data)) + '_entries_' + t + '.xlsx')

    writer = pd.ExcelWriter(str(output_file_name), engine='xlsxwriter')
    workbook = writer.book

    col_order = ["Target_id", "Gene_name", "Pharos_class", "protein_family", "protein_family_detail", "Number_isoforms",
                 'mpo_score', 'Tractable', 'Tractability_probability', 'In_training_set', 'structure_info_score',
                 'structural_drug_score',
                 'chemistry_score', 'biology_score', 'disease_score', 'genetic_score',
                 'information_score', 'safety_score',
                 "EBI Total Patent Count", "JensenLab PubMed Score", "NCBI Gene PubMed Count", "PubTator Score",
                 "total_patent_count", "year_max_patents", "count_patents_max_year", "novelty_score",
                 "total # publications", "number of Dementia publications", "Brain", "Endocrine_tissue",
                 "Female_tissue", "Immune", "Kidney", "Liver_gallbladder", "Lung", "Male_tissue", "Muscle_tissue",
                 "Pancreas", "Skin", "Soft_tissue", "gitract", "Expression_Selectivity", "tissue_max_expression",
                 "expression_max_tissue", 'EXP_LVL_AVG', 'EXP_LVL_STDDEV', 'Heart_alert', 'Heart_value', 'Liver_alert',
                 'Liver_value', 'Kidney_alert', 'Kidney_value',
                 "variants_count", "mutants_count", "gwas_count", 'number_of_genotypes',
                 "phenotypes_heterozygotes_lethal_count", "phenotypes_homozygotes_lethal_count",
                 "phenotypes_heterozygotes_normal_count", "phenotypes_homozygotes_normal_count", "Ab Count",
                 "MAb Count", "kegg_list", "kegg_count", "reactome_list", "reactome_count", "disease_count_uniprot",
                 "disease_list_tcrd", "disease_count_tcrd", "max_disease_score", "name_max_disease",
                 "OT_number_of_associations", "OT_number_of_disease_areas", "OT_list_max_disease_area",
                 "OT_max_association_diseaseArea_score", "OT_list_max_diseases", "OT_TOP10_diseases",
                 "OT_max_association_score", "OT_%_genetic_association", "OT_%_known_drug", "OT_%_litterature_mining",
                 "OT_%_animal_model", "OT_%_affected_pathway", "OT_%_rna_expression", "OT_%_somatic_mutation",
                 "OT_MAX_VAL_genetic_association", "OT_NUM_MAX_genetic_association", "OT_MAX_VAL_known_drug",
                 "OT_NUM_MAX_known_drug", "OT_MAX_VAL_litterature_mining", "OT_NUM_MAX_litterature_mining",
                 "OT_MAX_VAL_animal_model", "OT_NUM_MAX_animal_model", "OT_MAX_VAL_affected_pathway",
                 "OT_NUM_MAX_affected_pathway", "OT_MAX_VAL_rna_expression", "OT_NUM_MAX_rna_expression",
                 "OT_MAX_VAL_somatic_mutation", "OT_NUM_MAX_somatic_mutation", "PDB_total_count",
                 "PDB_with_Ligand_count", '%_sequence_covered', '%_domain_covered', "PDB_sites_tractable_count",
                 "PDB_sites_druggable_count",
                 "PDB_blast_close_count", "PDB_blast_max_similarity", "domains_count", "domain_tractable",
                 "domain_druggable", "mean_druggability_score", "stddev_druggability_score", "mean_area", "mean_volume",
                 "mean_fraction_apolar", "mean_pocket_score", "pdb_with_druggable_pocket", "druggable_pockets_total",
                 "mean_alt_druggability_score", "alt_stddev_druggability_score", "mean_alt_area", "mean_alt_volume",
                 "mean_alt_fraction_apolar", "mean_alt_pocket_score", "mean_alt_similarity", "max_alt_similarity",
                 "alt_pdb_with_druggable_pocket", "alt_druggable_pockets_total", "BindingDB_count",
                 "BindingDB_potent_count", "BindingDB_potent_phase2_count", "ChEMBL_bioactives_count",
                 "ChEMBL_bioactives_potent_count", "ChEMBL_bioactives_moderate_selectivity_count",
                 "ChEMBL_bioactives_good_selectivity_count", "ChEMBL_bioactives_great_selectivity_count",
                 "commercial_total", "commercial_potent_total"]

    data = data[col_order]
    data.to_excel(writer, sheet_name='Druggability_list', index=False, startrow=1)

    gen_info_len = 6
    score_len = 12
    litt_len = 10
    bio_len = 34
    pathways_len = 37
    structure_len = 29
    chemistry_len = 10

    header_groups = {'GENERAL INFO': (0, gen_info_len - 1),
                     'SCORES': (gen_info_len, gen_info_len + score_len - 1),
                     'LITTERATURE/PATENT INFORMATION': (
                     gen_info_len + score_len, gen_info_len + score_len + litt_len - 1),
                     'BIOLOGY': (
                     gen_info_len + score_len + litt_len, gen_info_len + score_len + litt_len + bio_len - 1),
                     'PATHWAYS AND DISEASES': (gen_info_len + score_len + litt_len + bio_len,
                                               gen_info_len + score_len + litt_len + bio_len + pathways_len - 1),
                     'STRUCTURAL INFORMATION': (gen_info_len + score_len + litt_len + bio_len + pathways_len,
                                                gen_info_len + score_len + litt_len + bio_len + pathways_len + structure_len - 1),
                     'CHEMISTRY': (gen_info_len + score_len + litt_len + bio_len + pathways_len + structure_len,
                                   gen_info_len + score_len + litt_len + bio_len + pathways_len + structure_len + chemistry_len - 1)}
    color_dict = {'GENERAL INFO': '#fde9d9', 'SCORES': '#ffff99', 'LITTERATURE/PATENT INFORMATION': '#d9d9d9',
                  'BIOLOGY': '#ebf1de',
                  'PATHWAYS AND DISEASES': '#f2dcdb', 'STRUCTURAL INFORMATION': '#dce6f1', 'CHEMISTRY': '#e4dfec'}
    for head, span in header_groups.items():
        col_header = workbook.add_format(
            {'bold': True, 'bg_color': color_dict[head], 'align': 'center', 'valign': 'vcenter'})
        if span[0] == span[1]:
            writer.sheets['Druggability_list'].write(0, span[0], head, col_header)
        else:
            writer.sheets['Druggability_list'].merge_range(0, span[0], 0, span[1], head, col_header)

    for col_num, value in enumerate(data.columns.values):
        vert_col_header = workbook.add_format(
            {'bold': True, 'bg_color': '#d9d9d9', 'align': 'center', 'valign': 'bottom', 'rotation': 90})
        for head, span in header_groups.items():
            if span[0] <= col_num <= span[1]:
                vert_col_header = workbook.add_format(
                    {'bold': True, 'bg_color': color_dict[head], 'align': 'center', 'valign': 'bottom', 'rotation': 90})
                break
        writer.sheets['Druggability_list'].write(1, col_num, value, vert_col_header)

    not_in_db.to_excel(writer, 'Not in DB', index=False)
    writer.save()
    print("[EXPORT]: Excel file: ", '[Export_' + str(len(data)) + '_entries_' + t + '.xlsx]', ' successfully generated')
    return "Success"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene', help='enter a single gene name', metavar='')
    parser.add_argument('-i', '--in_file',
                        help='Name of the input file with a list of genes (.txt - 1 gene per line)',
                        metavar='')
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
    parser.add_argument('-update_config', '--update_config',
                        help="use this if you want to update the config file",
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    if not args.gene and not args.in_file and not args.list_genes:
        print('Please use one of the optional input options :')
        parser.print_help()
        sys.exit()
    if not args.report_list and not args.report_single:
        print('ERROR: Please provide a report type "-rl" for a list report, "-rs" for a single target report')
        parser.print_help()
        sys.exit()
    return args


def main():
    global output_lists_path, output_single_path, targetDB, pubmed_email, list_of_entries
    args = parse_args()
    update_config = args.update_config
    while True:
        config = configparser.ConfigParser()
        config_file_path = Path('~/.druggability/config.ini').expanduser()
        config_file_path.parent.mkdir(exist_ok=True, parents=True)

        if config_file_path.is_file() and not update_config:
            config.read(str(config_file_path))
            todo = []
            for var_name in config['database_path']:
                if var_name in ['chembl']:
                    continue
                if not Path(config['database_path'][var_name]).is_file():
                    todo.append(var_name)
            for var_name in config['output_path']:
                if var_name in ['db_files']:
                    continue
                if not Path(config['output_path'][var_name]).is_dir() or config['output_path'][var_name] == '':
                    todo.append(var_name)
            if not cf.is_email(config['pubmed_email']['email']):
                todo.append('email')
            if todo:
                config = cf.get_config_from_user(config, todo=todo)
                with config_file_path.open(mode='w') as cfile:
                    config.write(cfile)
            else:
                # =============================# PATH TO SAVE REPORT FILES #============================#
                output_lists_path = config['output_path']['list']
                output_single_path = config['output_path']['single']
                # =============================# PATH TO SQLITE DB #============================#

                targetDB = config['database_path']['targetdb']
                pubmed_email = config['pubmed_email']['email']
                break
        else:
            config = cf.get_config_from_user(config, todo=['list', 'single', 'targetdb', 'email'], new=True)
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

    list_of_entries = get_list_entries()

    if args.report_single:
        for gene_name in gene_df.index:
            get_single_excel(gene_df.loc[gene_name])
    elif args.report_list:
        get_list_excel(gene_df)


def entry_point():
    main()


if __name__ == "__main__":
    entry_point()
