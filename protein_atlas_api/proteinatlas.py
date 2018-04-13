#!/usr/bin/env python
from __future__ import print_function

try:
    input = raw_input
except NameError:
    pass
import requests, sys, os

try:
    import urllib.request as urllib
    from urllib.error import *
except ImportError:
    import urllib
    from urllib2 import *
import xmltodict
import numpy as np
import progressbar
from operator import attrgetter
from operator import itemgetter
import xlsxwriter
import argparse
import pandas as pd

# FOR PROGRESSBAR
widgets = [' [', progressbar.Timer(), '] ', progressbar.Bar(), ' [', progressbar.Percentage(), '] ', ' (',
           progressbar.ETA(), ') ']

# GENERAL DICTIONARIES
sub_class = {'Brain': ['caudate', 'cerebellum', 'cerebral cortex', 'hippocampus'], 'Skin': ['skin','skin 1','skin 2'],
             'Soft_tissue': ['soft tissue'],
             'Female_tissue': ['breast', 'cervix, uterine', 'endometrium','endometrium 1','endometrium 2','fallopian tube',
                               'ovary', 'placenta', 'vagina'],
             'Endocrine_tissue': ['adrenal gland', 'parathyroid gland', 'thyroid gland'],
             'Immune': ['appendix', 'bone marrow', 'lymph node', 'spleen', 'tonsil'],
             'Muscle_tissue': ['heart muscle', 'skeletal muscle', 'smooth muscle'],
             'Lung': ['bronchus', 'lung', 'nasopharynx'],
             'Liver_gallbladder': ['liver', 'gallbladder'],
             'Pancreas': ['pancreas'],
             'gitract': ['colon', 'duodenum', 'esophagus', 'oral mucosa', 'rectum', 'salivary gland',
                         'small intestine', 'stomach','stomach 1','stomach 2'],
             'Kidney': ['kidney', 'urinary bladder'],
             'Male_tissue': ['epididymis', 'prostate', 'seminal vesicle', 'testis']}
conversion_dict = {'low': 1, 'medium': 2, 'high': 3, 'not detected': 0}
tissue_to_organ_dict = {}
all_cells=[]
for organ, tissues in sub_class.items():
    for k in tissues:
        tissue_to_organ_dict[k] = organ


def gene_name_to_ensemblid(gene_name, alt=False):
    base_url = "http://rest.ensembl.org/lookup/symbol/homo_sapiens/"
    if alt:
        base_url = "http://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/"
    constructed_request = base_url + gene_name + '?content-type=application/json;format=condensed'
    try:
        request = requests.get(constructed_request)
        # request = urllib.urlopen(constructed_request).read().decode('utf-8')
    except HTTPError:
        return None
    if not request.ok:
        return None
    # print(request)
    data = request.json()
    # print(data)
    return data['id']


def get_xml(ensemblid):
    base_url = 'http://www.proteinatlas.org/'
    constructed_url = base_url + ensemblid + '.xml'
    return _parse_xml(constructed_url)


def _parse_xml(xml_URL):
    try:
        xml_online = urllib.urlopen(xml_URL)
        xml_dict = xmltodict.parse(xml_online.read())
        return xml_dict
    except (HTTPError, URLError):
        return None


def get_protein_level_tissue(xml_dict):
    protein_level = {'Cells':[],'Tissues':{},'Organs':{}}
    tissue_count={}
    for sub_key in sub_class.keys():
        protein_level['Organs'][sub_key] = 0
    try:
        for i in xml_dict['proteinAtlas']['entry']['tissueExpression']['data']:
            organ = tissue_to_organ_dict[i['tissue']]
            if i['tissue'] in protein_level['Tissues'].keys():
                if i['tissue'] not in tissue_count.keys():
                    tissue_count[i['tissue']]=1
                else:
                    tissue_count[i['tissue']]+=1
                i['tissue']=i['tissue']+'-'+str(tissue_count[i['tissue']])
            if isinstance(i['tissueCell'],list):
                for j in i['tissueCell']:
                    cell_id=organ+'_'+i['tissue']+'_'+j['cellType']
                    protein_level['Cells'].append({'level':int(conversion_dict[j['level']['#text']]),'tissue':i['tissue'],'cell_type':j['cellType'],'organ':organ,'cell_id':cell_id})
                    all_cells.append({'tissue':i['tissue'],'cell_type':j['cellType'],'organ':organ,'cell_id':cell_id})
            else:
                cell_id = organ + '_' + i['tissue'] + '_' + i['tissueCell']['cellType']
                protein_level['Cells'].append({'level': int(conversion_dict[i['tissueCell']['level']['#text']]),'tissue': i['tissue'], 'cell_type': i['tissueCell']['cellType'],'organ':organ,'cell_id':cell_id})
                all_cells.append({'tissue': i['tissue'], 'cell_type': i['tissueCell']['cellType'], 'organ': organ, 'cell_id': cell_id})
            protein_level['Tissues'][i['tissue']] = int(conversion_dict[i['level']['#text']])
            protein_level['Organs'][organ] += protein_level['Tissues'][i['tissue']]
    except KeyError:
        return None
    return protein_level


def get_level_tissue_class(tissue_level):

    tissue_class_level = {}
    for sub_key in sub_class.keys():
        tissue_class_level[sub_key] = 0
    for key in tissue_level.keys():
        for sub_key in sub_class.keys():
            if key in sub_class[sub_key]:
                tissue_class_level[sub_key] += tissue_level[key]
    return tissue_class_level


def selectivity(general_tissue_level):
    ratio = []
    sum_all_levels = np.sum(list(general_tissue_level.values()))
    if sum_all_levels == 0:
        return 10
    for keys in general_tissue_level.keys():
        part = general_tissue_level[keys] / sum_all_levels
        ratio.append(part)
    sum_ratioln = 0
    for i in ratio:
        if i == 0:
            continue
        sum_ratioln += i * np.log(i)
    entrop_sel = -sum_ratioln
    return entrop_sel


def writeXL(list_of_proteinexpression, no_data, output=None):
    workbook = xlsxwriter.Workbook(file_path.rstrip('.txt') + '_out.xlsx')
    if output:
        workbook = xlsxwriter.Workbook(output + '_out.xlsx')
    worksheet_Organs = workbook.add_worksheet('Expression_organs')
    worksheet_Tissues = workbook.add_worksheet('Expression_tissues')
    worksheet_Cells = workbook.add_worksheet('Expression_cells')
    bold = workbook.add_format({'bold': True})
    cells_set = sorted([dict(t) for t in set([tuple(d.items()) for d in all_cells])],key=itemgetter('organ', 'tissue', 'cell_type'))
    tissue_set = sorted([dict(t) for t in set([(('tissue',i['tissue']),('organ',i['organ'])) for i in cells_set])],key=itemgetter('organ', 'tissue'))
    organ_set = sorted(list(set([i['organ'] for i in cells_set])))
    cells_rowval={cells_set[i]['cell_id']:i+1 for i in range(len(cells_set))}
    tissue_rowval={tissue_set[i]['tissue']:i+1 for i in range(len(tissue_set))}
    organ_rowval={organ_set[i]:i+1 for i in range(len(organ_set))}

    cells_header=['Organ','Tissue','Cell type']
    tissue_header=['Organ','Tissue']
    organ_header=['Organ']
    for i in list_of_proteinexpression:
        cells_header.append(i.gene)
        tissue_header.append(i.gene)
        organ_header.append(i.gene)

    for x in range(len(cells_header)):
        worksheet_Cells.write(0, x, cells_header[x], bold)
    for x in range(len(tissue_header)):
        worksheet_Tissues.write(0, x, tissue_header[x], bold)
    for x in range(len(organ_header)):
        worksheet_Organs.write(0, x, organ_header[x], bold)
    for i in cells_set:
        worksheet_Cells.write(cells_rowval[i['cell_id']],0,i['organ'])
        worksheet_Cells.write(cells_rowval[i['cell_id']], 1, i['tissue'])
        worksheet_Cells.write(cells_rowval[i['cell_id']], 2, i['cell_type'],bold)
    for i in tissue_set:
        worksheet_Tissues.write(tissue_rowval[i['tissue']],0,i['organ'])
        worksheet_Tissues.write(tissue_rowval[i['tissue']], 1, i['tissue'],bold)
    for i in organ_set:
        worksheet_Organs.write(organ_rowval[i],0,i,bold)
    worksheet_Organs.write(len(organ_rowval)+1,0,'Selectivity Score',bold)

    for i in list_of_proteinexpression:
        for c in i.protein_lvl['Cells']:
            worksheet_Cells.write(cells_rowval[c['cell_id']],cells_header.index(i.gene),c['level'])
        for t,val in i.protein_lvl['Tissues'].items():
            worksheet_Tissues.write(tissue_rowval[t],tissue_header.index(i.gene),val)
        for o,val in i.protein_lvl['Organs'].items():
            worksheet_Organs.write(organ_rowval[o],organ_header.index(i.gene), val)
            worksheet_Organs.write(len(organ_rowval)+1,organ_header.index(i.gene),i.selective_entropy)
    worksheet_nodata = workbook.add_worksheet('No Expression Data')
    x = 1
    worksheet_nodata.write(0, 0, 'List of Gene with no expression data')
    no_data=sorted(no_data)
    for j in no_data:
        worksheet_nodata.write(x, 0, str(j))
        x += 1

    workbook.close()


class ProteinExpression:
    def __init__(self, gene_name, id=None):
        self.gene = gene_name
        if not id:
            id = gene_name_to_ensemblid(gene_name)
        if not id:
            id = gene_name_to_ensemblid(gene_name, alt=True)
        if id:
            doc = get_xml(id)
            if not doc:
                id = gene_name_to_ensemblid(gene_name, alt=True)
                if id:
                    doc = get_xml(id)
            if doc:
                self.protein_lvl = get_protein_level_tissue(doc)
                if self.protein_lvl:
                    self.protein_lvl['Cells'] = sorted(self.protein_lvl['Cells'],
                                                       key=itemgetter('organ', 'tissue', 'cell_type'))
                    self.max_organ = max(self.protein_lvl['Organs'], key=self.protein_lvl['Organs'].get)
                    self.selective_entropy = selectivity(self.protein_lvl['Organs'])
                    self.sort_selective_entropy = -self.selective_entropy
                else:
                    pass  # No protein expression data
            else:
                self.protein_lvl = None
                print(self.gene, ': No Human Protein Atlas data found  (', id, ')')
        else:
            self.protein_lvl = None
            print(self.gene, ': No ensembl id found')


    def show(self):
        return print(self.gene, ': ', self.max_organ, ' (Tissue Max = ',
                     self.protein_lvl['Organs'][self.max_organ],
                     '|| Ssel = ', self.selective_entropy, ')')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', help='Name of the input file', metavar='')
    parser.add_argument('-g', '--gene', help='enter a single gene name', metavar='')
    parser.add_argument('-e', '--ensembl_ID', help='Enter the ensembl accession number', metavar='')
    parser.add_argument('-l', '--list_genes', help='Enter a list of gene name separated by a ","', metavar='')
    parser.add_argument('-X', '--excel', help='Write an excel spreadsheet output (default terminal output, use excel '
                                              'for more detailed results)', action='store_true', default=False)
    args = parser.parse_args()

    while True:
        gene_name = []
        id = None

        if args.in_file:
            if os.path.exists(args.in_file):
                file_path = args.in_file
                with open(file_path, 'r') as gene_list:
                    list_of_genes = gene_list.readlines()
                    gene_name = list(set(list_of_genes))
                break
            else:
                print("The inputed file does not exist ... program exit")
                sys.exit()
        if args.gene:
            file_path = os.getcwd() + '/' + args.gene
            gene_name = [args.gene]
            break
        if args.list_genes:
            list_of_genes = args.list_genes.split(',')
            gene_name = list(set(list_of_genes))
            file_path = os.getcwd() + '/' + gene_name[0] + '_and_' + str(len(gene_name)) + '_others'
            break
        if args.ensembl_ID:
            id = args.ensembl_ID
            file_path = os.getcwd() + '/' + id
            gene_name = [input('Please input a gene name for the corresponding ensembl ID: ')]
            break
        print('No input data found :')
        parser.print_help()
        sys.exit()

    list_expression = []
    bar = progressbar.ProgressBar(maxval=len(gene_name), widgets=widgets)
    bar.start()
    i = 0
    if gene_name:
        for gene in gene_name:
            bar.update(i)
            gene = gene.rstrip('\n')
            if id:
                list_expression.append(ProteinExpression(gene, id=id))
            else:
                list_expression.append(ProteinExpression(gene))
            i += 1
    bar.finish()

    no_expression_data = []
    for i in reversed(list_expression):
        if not i.protein_lvl:
            no_expression_data.append(str(i.gene))
            list_expression.remove(i)
    if args.excel:
        writeXL(list_expression, no_expression_data)

    for i in list_expression:
        i.show()
