#!/usr/bin/env python

import requests
import urllib.request as urllib
from urllib.error import *
import xmltodict
import numpy as np
from operator import itemgetter
import targetDB.utils.retryers as ret


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

@ret.retryer(max_retries=10, timeout=10)
def gene_name_to_ensemblid(gene_name, alt=False):
    base_url = "http://rest.ensembl.org/lookup/symbol/homo_sapiens/"
    if alt:
        base_url = "http://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/"
    constructed_request = base_url + gene_name + '?content-type=application/json;format=condensed'
    try:
        request = requests.get(constructed_request)
    except HTTPError:
        return None
    if not request.ok:
        return None
    data = request.json()
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