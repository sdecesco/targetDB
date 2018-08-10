#!/usr/bin/env python
import pandas as pd
from targetDB import cns_mpo as mpo
import numpy as np
import io, sqlite3, math
import scipy.stats as sc
import matplotlib.pyplot as plt


class StdevFunc:
	def __init__(self):
		self.M = 0.0
		self.S = 0.0
		self.k = 1

	def step(self, value):
		if value is None:
			return
		tM = self.M
		self.M += (value - tM) / self.k
		self.S += (value - tM) * (value - self.M)
		self.k += 1

	def finalize(self):
		if self.k == 2:
			return 0
		elif self.k < 3:
			return None
		return math.sqrt(self.S / (self.k - 2))


def get_descriptors(target_id, targetdb=None):
	connector_targetDB = sqlite3.connect(targetdb)
	connector_targetDB.create_aggregate('stddev', 1, StdevFunc)

	list_queries = {'gen_info': """SELECT * FROM Targets WHERE Target_id='%s'""" % target_id,
					'disease': """SELECT Target_id,disease_name,disease_id FROM disease WHERE Target_id='%s'""" % target_id,
					'reactome': """SELECT pathway_name FROM pathways WHERE pathway_dataset='Reactome pathways data set' AND Target_id='%s'""" % target_id,
					'kegg': """SELECT pathway_name FROM pathways WHERE pathway_dataset='KEGG pathways data set' AND Target_id='%s'""" % target_id,
					'disease_exp': """SELECT
	  disease,
	  round(avg(t_stat),1) as t_stat,
	  round(stddev(t_stat),1) as std_dev_t,
	  count(t_stat) as n,
	  max(expression_status) as direction
	  FROM diff_exp_disease
	  WHERE Target_id='%s'
	  GROUP BY Target_id,disease
	  ORDER BY t_stat DESC""" % target_id,
					'gwas': """SELECT
	  phenotype,
	  organism,
	  p_value,
	  first_author as author,
	  publication_year as 'year',
	  pubmed_id
	FROM gwas
	WHERE Target_id='%s'
	ORDER BY phenotype""" % target_id,
					'tissue': """SELECT
	  Tissue,
	  round(avg(t_stat),1) as t_stat,
	  round(stddev(t_stat),1) as std_dev_t,
	  count(t_stat) as n
	  FROM diff_exp_tissue
	  WHERE Target_id='%s'
	  GROUP BY Tissue
	  ORDER BY t_stat DESC""" % target_id,
					'selectivity': """SELECT
	Selectivity_entropy
	FROM protein_expression_selectivity
	WHERE Target_id='%s'""" % target_id,
					'tissue_expression': """SELECT
	  organ,
	  tissue,
	  cell,
	  value
	  FROM protein_expression_levels
	  WHERE Target_id='%s'""" % target_id,
					'phenotype': """SELECT
	  Allele_symbol,
	  Allele_type,
	  CASE WHEN zygosity is null THEN 'NOT DECLARED' ELSE UPPER(zygosity) END AS zygosity,
	  genotype,
	  Phenotype
	FROM phenotype WHERE Target_id='%s'
	ORDER BY Allele_id,zygosity,genotype""" % target_id,
					'isoforms': """SELECT
	  (T.Gene_name||'-'||I.Isoform_name) as isoform_name,
	  I.Isoform_id,
	  I.Sequence,
	  I.n_residues,
	  CASE WHEN I.Canonical = 1 THEN 'Yes' ELSE 'No' END AS is_canonical,
	  I.Identity AS similarity
	FROM Isoforms I
	LEFT JOIN Targets T
	  ON I.Target_id = T.Target_id
	WHERE I.Target_id='%s' ORDER BY I.Canonical DESC""" % target_id,
					'isoforms_mod': """SELECT
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
	WHERE IM.isoform_id in (SELECT I.Isoform_id FROM Isoforms I WHERE I.Target_id='%s')""" % target_id,
					'var': """SELECT
	  M.start,
	  M.stop,
	  M.previous AS previous_seq,
	  M.action AS modification_type,
	  M.new AS new_seq,
	  M.domains AS in_domains,
	  M.comment AS comments
	FROM modifications M
	WHERE M.mod_type = 'VAR' AND M.Target_id='%s'""" % target_id,
					'mut': """SELECT
	  M.start,
	  M.stop,
	  M.previous AS previous_seq,
	  M.action AS modification_type,
	  M.new AS new_seq,
	  M.domains AS in_domains,
	  M.comment AS comments
	FROM modifications M
	WHERE M.mod_type = 'MUTAGEN' AND M.Target_id='%s'""" % target_id,
					'domains': """SELECT
	  Domain_name,
	  Domain_start as start,
	  Domain_stop as stop,
	  length,
	  source_name as source
	FROM Domain_targets
	WHERE Target_id='%s'""" % target_id,
					'pdb_blast': """SELECT
	  Hit_PDB_code as PDB_code,
	  Chain_Letter as Chain,
	  similarity,
	  Hit_gene_name as gene,
	  Hit_gene_species as species,
	  max(tractable) SITES_tractable,
	  max(druggable) SITES_druggable
	FROM `3D_Blast`
	  LEFT JOIN drugEbility_sites DS
	  ON DS.pdb_code=Hit_PDB_code
	WHERE Query_target_id='%s'
	GROUP BY Hit_PDB_code
	ORDER BY similarity DESC""" % target_id,
					'pdb': """SELECT
	  C.PDB_code,
	  P.Technique,
	  P.Resolution,
	  GROUP_CONCAT(DISTINCT C.Chain) AS Chain,
	  C.n_residues,
	  C.start_stop,
	  GROUP_CONCAT(DISTINCT D.Domain_name) AS Domain_name,
	  B.type type_of_binder,
	  B.binding_type,
	  B.binding_operator operator,
	  B.binding_value 'value',
	  B.binding_units units,
	  B.lig_name Ligand_name,
	  B.pub_year publication_year,
	  max(DS.tractable) SITES_tractable,
	  max(DS.druggable) SITES_druggable
	FROM PDB_Chains C
	  LEFT JOIN PDB P
		ON C.PDB_code = P.PDB_code
	  LEFT JOIN PDBChain_Domain Domain
		ON C.Chain_id = Domain.Chain_id
	  LEFT JOIN Domain_targets D
		ON Domain.Domain_id = D.domain_id
	  LEFT JOIN pdb_bind B
		ON B.pdb_code = C.PDB_code
	  LEFT JOIN drugEbility_sites DS
		ON DS.pdb_code = LOWER(C.PDB_code)
	WHERE C.Target_id='%s'
	GROUP BY C.PDB_code,P.Technique,P.Resolution,C.n_residues,C.start_stop""" % target_id,
					'pockets': """SELECT
	  F.PDB_code,
	  F.DrugScore as druggability_score,
	  round(F.total_sasa,1) as area,
	  round(F.volume,1) as volume,
	  round((F.apolar_sasa/F.total_sasa)*100,1) as fraction_apolar,
	  F.Pocket_number as pocket_number,
	  F.Score as pocket_score,
	  GROUP_CONCAT((D.Domain_name||' ('||Domain.Coverage||'%)') ) as domains
	FROM fPockets F
	  LEFT JOIN fPockets_Domain Domain
		ON F.Pocket_id = Domain.Pocket_id
	  LEFT JOIN Domain_targets D
		ON Domain.Domain_id=D.domain_id
	WHERE F.Target_id='{target}'
	AND F.druggable='TRUE' AND F.blast='FALSE'
	GROUP BY F.PDB_code,F.DrugScore,F.total_sasa,F.volume,fraction_apolar,pocket_number,pocket_score""".format(
						target=target_id),
					'alt_pockets': """SELECT
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
	ORDER BY B.similarity DESC""" % target_id,
					'bioactives': """SELECT
	B.lig_id,
	  B.assay_id,
	  B.target_id,
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
	AND A.confidence_score>=8""" % target_id,
					'commercials': """SELECT
	smiles,
	affinity_type,
	op,
	affinity_value,
	affinity_unit,
	price,
	website
	FROM purchasable_compounds
	WHERE target_id='%s'""" % target_id,
					'bindingDB': """SELECT
	  B.ligand_name,
	  B.ZincID,
	  B.`IC50(nM)`,
	  B.`EC50(nM)`,
	  B.`Ki(nM)`,
	  B.`Kd(nM)`,
	  B.`kon(M-1s-1)`,
	  B.`koff(s-1)`,
	  B.pH,
	  B.`Temp`,
	  B.Source,
	  B.DOI,
	  B.institution,
	  B.patent_number,
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
	  L.molecular_species,
	  L.num_ro5_violations as ro5_violations,
	  L.ro3_pass as pass_ro3,
	  B.ligand_smiles as SMILES
	FROM BindingDB B
	  LEFT JOIN ligands L
	  ON B.inchi_key = L.std_inchi_key
	WHERE target_id = '%s'""" % target_id,
					'domain_drugE': """SELECT
	GROUP_CONCAT(DISTINCT UPPER(pdb_code) ) pdb_list,
	domain_fold,
	domain_superfamily,
	max(tractable) tractable,
	max(druggable) druggable
	FROM drugEbility_domains
	WHERE pdb_code in (SELECT DISTINCT LOWER(PDB_code)
	FROM PDB_Chains
	WHERE target_id = '%s')
	GROUP BY domain_fold""" % target_id}
	results = {qname: pd.read_sql(query, con=connector_targetDB) for qname, query in list_queries.items()}

	data = results['gen_info'].drop(['Species', 'species_id', 'Sequence',
									 'Cell_location', 'Process', 'Function', 'Synonyms',
									 'Date_modified', 'Protein_class', 'Protein_class_desc',
									 'Protein_class_short', 'chembl_id'], axis=1)

	if not results['kegg'].empty:
		data['pathways_kegg'] = len(results['kegg'])
	else:
		data['pathways_kegg'] = [0]

	if not results['pockets'].empty:
		if results['pdb'].empty:
			data['number_druggable_pockets_NORM'] = 1
		else:
			data['number_druggable_pockets_NORM'] = results['pockets']['PDB_code'].nunique() / len(results['pdb'])
		data['druggable_pockets_total'] = len(results['pockets'])
		data['pockets_mean_area'] = results['pockets']['area'].mean()
		data['pockets_mean_volume'] = results['pockets']['volume'].mean()
		data['pockets_mean_apolarfrac'] = results['pockets']['fraction_apolar'].mean()
		data['pockets_mean_druggability_score'] = results['pockets']['druggability_score'].mean()
		data['pockets_stddev_druggability_score'] = results['pockets']['druggability_score'].std()
		data['alternate_pockets'] = [False]
	elif not results['alt_pockets'].empty:
		data['number_druggable_pockets_NORM'] = results['alt_pockets']['PDB_code'].nunique() / len(results['pdb_blast'])
		data['druggable_pockets_total'] = len(results['alt_pockets'])
		data['pockets_mean_area'] = results['alt_pockets']['area'].mean()
		data['pockets_mean_volume'] = results['alt_pockets']['volume'].mean()
		data['pockets_mean_apolarfrac'] = results['alt_pockets']['fraction_apolar'].mean()
		data['pockets_mean_druggability_score'] = results['alt_pockets']['druggability_score'].mean()
		data['pockets_stddev_druggability_score'] = results['alt_pockets']['druggability_score'].std()
		data['alternate_pockets'] = [True]
	else:
		data['number_druggable_pockets_NORM'] = [0]
		data['druggable_pockets_total'] = [0]
		data['pockets_mean_area'] = [None]
		data['pockets_mean_volume'] = [None]
		data['pockets_mean_apolarfrac'] = [None]
		data['pockets_mean_druggability_score'] = [None]
		data['pockets_stddev_druggability_score'] = [None]
		data['alternate_pockets'] = [None]

	data['bindingDB_count'] = len(results['bindingDB'])
	potent = results['bindingDB'][((results['bindingDB'][['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)']] <= 10) &
								   (results['bindingDB'][['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)']].notna())).any(
		axis=1)]
	data['bindingDB_potent_count'] = len(potent)
	data['bindingDB_potent_phase2_count'] = len(potent[potent['max_phase'] >= 2])

	data['isoforms_count'] = len(results['isoforms'])

	if not results['tissue_expression'].empty:
		tissue = results['tissue_expression'].groupby('tissue').max().groupby('organ').mean().round(1).transpose()
		tissue['Target_id'] = target_id
		data = pd.merge(data, tissue)
		data['Expression_selectivity'] = results['selectivity']
		data['tissue_max_expression'] = results['tissue_expression'].groupby('tissue').max().groupby(
			'organ').mean().round(1).idxmax().value
	else:
		data['Expression_selectivity'] = np.nan

	data['variants_count'] = len(results['var'])
	data['disease_count'] = len(results['disease'])
	data['gwas_count'] = len(results['gwas'])
	data['mutant_count'] = len(results['mut'])

	results['bioactives']['CNS_MPO'] = mpo.calc_mpo_score(bpka=results['bioactives']['bpKa'],
														  logP=results['bioactives']['LogP'],
														  logD=results['bioactives']['LogD'],
														  MW=results['bioactives']['MW'],
														  HBD=results['bioactives']['HBD'],
														  TPSA=results['bioactives']['TPSA'])
	best = results['bioactives'][results['bioactives']['pX'].notnull()]
	total_bioact = len(best)
	best = best[best['pX'] >= 8]
	total_potent = len(best)
	best = best[best['standard_type'].isin(['Ki', 'Kd'])]

	if not best.empty:
		best['pX'].fillna(-np.log10(best.value_num / 1000000000), inplace=True)
		query_lig = "','".join(best.lig_id.unique())
		query = """SELECT
				B.lig_id,
				B.Target_id,
				B.target_name,
				ROUND(AVG(B.value_num),2) avg_value,
				ROUND(STDDEV(B.value_num),2) sttdev,
				COUNT(*) n_values
				FROM bioactivities B
				  LEFT JOIN assays A
				  ON B.assay_id=A.assay_id
				WHERE B.operator='=' 
				  AND B.lig_id in ('%s')
				  AND A.bioactivity_type='Binding'
				  AND UPPER(B.standard_type) in ('KD','KI')
				  AND B.data_validity_comment is NULL
				  AND A.confidence_score>=8
		GROUP BY B.lig_id,B.Target_id""" % query_lig
		entropies = []
		binding_data = pd.read_sql(query, con=connector_targetDB)
		best_target_id = best.iloc[0]['Target_id']
		if not binding_data.empty:
			for name, group in binding_data.groupby('lig_id'):
				best_target = True
				group = group[(group['sttdev'] < group['avg_value'])].copy()
				group['association'] = (1 / group.avg_value)
				group['association_prob'] = group.association / group.association.sum()
				if len(group) > 1:
					if group.loc[group['association_prob'].idxmax()]['Target_id'] == best_target_id:
						best_target = True
					else:
						best_target = False
				entropies.append({'Selectivity': round(sc.entropy(group.association_prob), 2), 'lig_id': name,
								  'number of other targets': len(group),
								  'targets name': ' / '.join(np.unique(group['target_name'].values)),
								  'best_target': best_target})

			entropy = pd.DataFrame(data=entropies)

			best = pd.merge(best, entropy, on='lig_id')

			total_moderate_selectivity = len(best[(best['Selectivity'] <= 2) & (best['best_target'] == True) & (
					best['number of other targets'] > 1)])
			total_good_selectivity = len(best[(best['Selectivity'] <= 1.5) & (best['best_target'] == True) & (
					best['number of other targets'] > 1)])
			total_great_selectivity = len(best[(best['Selectivity'] <= 1) & (best['best_target'] == True) & (
					best['number of other targets'] > 1)])
		else:
			total_moderate_selectivity = 0
			total_good_selectivity = 0
			total_great_selectivity = 0
	else:
		total_moderate_selectivity = 0
		total_good_selectivity = 0
		total_great_selectivity = 0

	data['ChEMBL_bioactives_count'] = total_bioact
	data['ChEMBL_bioactives_potent_count'] = total_potent
	data['ChEMBL_bioactives_moderate_selectivity_count'] = total_moderate_selectivity
	data['ChEMBL_bioactives_good_selectivity_count'] = total_good_selectivity
	data['ChEMBL_bioactives_great_selectivity_count'] = total_great_selectivity

	results['phenotype']['lethal'] = results['phenotype']['Phenotype'].str.contains('lethal|death', case=False)
	results['phenotype']['normal'] = results['phenotype']['Phenotype'].str.contains('no abnormal phenotype detected',
																					case=False)
	data['phenotypes_count'] = len(results['phenotype'])
	data['phenotypes_heterozygotes_lethal_count'] = len(results['phenotype'][
															(results['phenotype']['lethal'] == True) & (
																	results['phenotype'][
																		'zygosity'] == 'HETEROZYGOTE')])
	data['phenotypes_homozygotes_lethal_count'] = len(results['phenotype'][(results['phenotype']['lethal'] == True) & (
			results['phenotype']['zygosity'] == 'HOMOZYGOTE')])
	data['phenotypes_heterozygotes_normal_count'] = len(results['phenotype'][
															(results['phenotype']['normal'] == True) & (
																	results['phenotype'][
																		'zygosity'] == 'HETEROZYGOTE')])
	data['phenotypes_homozygotes_normal_count'] = len(results['phenotype'][(results['phenotype']['normal'] == True) & (
			results['phenotype']['zygosity'] == 'HOMOZYGOTE')])

	data['PDB_total_count'] = len(results['pdb'])
	data['PDB_with_Ligand_count'] = len(results['pdb'][results['pdb']['type_of_binder'].notnull()])
	data['PDB_sites_tractable'] = len(results['pdb'][results['pdb']['SITES_tractable'] == 1])
	data['PDB_sites_druggable'] = len(results['pdb'][results['pdb']['SITES_druggable'] == 1])
	data['PDB_blast_close_count'] = len(results['pdb_blast'])

	data['domains_count'] = len(results['domains'])
	data['domain_druggable'] = len(results['domain_drugE'][results['domain_drugE']['druggable'] == 1])
	data['domain_tractable'] = len(results['domain_drugE'][results['domain_drugE']['tractable'] == 1])

	data['commercials_total'] = len(results['commercials'])
	data['commercials_potent_total'] = len(results['commercials'][results['commercials']['affinity_value'] <= 100])


	query_id = """SELECT tcrd_id FROM tcrd_id WHERE Target_id= '%s'""" % target_id
	tcrd_id = pd.read_sql(query_id, con=connector_targetDB)
	if tcrd_id.empty:
		tcrd_id = 'None'
	else:
		tcrd_id = tcrd_id.iloc[0]['tcrd_id']

	tcrd_queries = {'target': """SELECT * FROM tcrd_target WHERE tcrd_id = '%s' """ % tcrd_id,
					'tdl_info': """SELECT * FROM tcrd_info WHERE protein_id = '%s'""" % tcrd_id,
					'disease': """SELECT protein_id,disease_id,doid,score,name,parent
	FROM
	tcrd_disease
	WHERE protein_id='%s'
	ORDER BY score DESC""" % tcrd_id,
					'novelty': """SELECT score FROM tcrd_novelty WHERE protein_id = '%s'""" % tcrd_id}
	tcrd_res = {qname: pd.read_sql(query, con=connector_targetDB) for qname, query in tcrd_queries.items()}

	tcrd_data = tcrd_res['target']
	tcrd_data['Target_id'] = target_id
	tcrd_data.set_index('tcrd_id', inplace=True)

	tcrd_res['tdl_info'].index = tcrd_res['tdl_info'].protein_id
	tcrd_data = tcrd_data.join(tcrd_res['tdl_info'])

	tcrd_data['tcrd_disease_count'] = len(tcrd_res['disease'][(~tcrd_res['disease']['doid'].isin(
		tcrd_res['disease']['parent'].unique())) & (tcrd_res['disease']['score'] > 1)])

	if len(tcrd_res['novelty']) != 0:
		tcrd_data['tcrd_novelty_score'] = tcrd_res['novelty']['score'].iloc[0]
	else:
		tcrd_data['tcrd_novelty_score'] = np.nan

	data = pd.merge(data, tcrd_data, on='Target_id')
	connector_targetDB.close()
	return data


def get_descriptors_list(target_id, targetdb=None):
	connector_targetDB = sqlite3.connect(targetdb)
	connector_targetDB.create_aggregate('stddev', 1, StdevFunc)
	list_queries = {'gen_info': """SELECT * FROM Targets WHERE Target_id in ('%s')""" % target_id,
					'disease': """SELECT Target_id,disease_name,disease_id FROM disease WHERE Target_id in ('%s')""" % target_id,
					'reactome': """SELECT pathway_name,Target_id FROM pathways WHERE pathway_dataset='Reactome pathways data set' AND Target_id in ('%s')""" % target_id,
					'kegg': """SELECT Target_id,pathway_name FROM pathways WHERE pathway_dataset='KEGG pathways data set' AND Target_id in ('%s')""" % target_id,
					'gwas': """SELECT
	  Target_id,
	  phenotype,
	  organism,
	  p_value,
	  first_author as author,
	  publication_year as 'year',
	  pubmed_id
	FROM gwas
	WHERE Target_id in ('%s')
	ORDER BY phenotype""" % target_id,
					'selectivity': """SELECT
	Target_id,
	Selectivity_entropy
	FROM protein_expression_selectivity
	WHERE Target_id in ('%s')""" % target_id,
					'tissue_expression': """SELECT Target_id,
	  organ,
	  tissue,
	  cell,
	  value
	  FROM protein_expression_levels
	  WHERE Target_id in ('%s')""" % target_id,
					'phenotype': """SELECT Target_id,
	  Allele_symbol,
	  Allele_type,
	  CASE WHEN zygosity is null THEN 'NOT DECLARED' ELSE UPPER(zygosity) END AS zygosity,
	  genotype,
	  Phenotype
	FROM phenotype WHERE Target_id in ('%s')
	ORDER BY Allele_id,zygosity,genotype""" % target_id,
					'isoforms': """SELECT I.Target_id,
	  (T.Gene_name||'-'||I.Isoform_name) as isoform_name,
	  I.Isoform_id,
	  I.Sequence,
	  I.n_residues,
	  CASE WHEN I.Canonical = 1 THEN 'Yes' ELSE 'No' END AS is_canonical,
	  I.Identity AS similarity
	FROM Isoforms I
	LEFT JOIN Targets T
	  ON I.Target_id = T.Target_id
	WHERE I.Target_id in ('%s') ORDER BY I.Canonical DESC""" % target_id,
					'var': """SELECT M.Target_id,
	  M.start,
	  M.stop,
	  M.previous AS previous_seq,
	  M.action AS modification_type,
	  M.new AS new_seq,
	  M.domains AS in_domains,
	  M.comment AS comments
	FROM modifications M
	WHERE M.mod_type = 'VAR' AND M.Target_id in ('%s')""" % target_id,
					'mut': """SELECT M.Target_id,
	  M.start,
	  M.stop,
	  M.previous AS previous_seq,
	  M.action AS modification_type,
	  M.new AS new_seq,
	  M.domains AS in_domains,
	  M.comment AS comments
	FROM modifications M
	WHERE M.mod_type = 'MUTAGEN' AND M.Target_id in ('%s')""" % target_id,
					'domains': """SELECT
	  Target_id,
	  Domain_name,
	  Domain_start as start,
	  Domain_stop as stop,
	  length,
	  source_name as source
	FROM Domain_targets
	WHERE Target_id in ('%s')""" % target_id,
					'pdb_blast': """SELECT Query_target_id,
	  Hit_PDB_code as PDB_code,
	  Chain_Letter as Chain,
	  similarity,
	  Hit_gene_name as gene,
	  Hit_gene_species as species,
	  max(tractable) SITES_tractable,
	  max(druggable) SITES_druggable
	FROM `3D_Blast`
	  LEFT JOIN drugEbility_sites DS
	  ON DS.pdb_code=Hit_PDB_code
	WHERE Query_target_id in ('%s')
	GROUP BY Query_target_id,Hit_PDB_code
	ORDER BY similarity DESC""" % target_id,
					'pdb': """SELECT
	  C.Target_id,
	  C.PDB_code,
	  P.Technique,
	  P.Resolution,
	  GROUP_CONCAT(DISTINCT C.Chain) AS Chain,
	  C.n_residues,
	  C.start_stop,
	  GROUP_CONCAT(DISTINCT D.Domain_name) AS Domain_name,
	  B.type type_of_binder,
	  B.binding_type,
	  B.binding_operator operator,
	  B.binding_value 'value',
	  B.binding_units units,
	  B.lig_name Ligand_name,
	  B.pub_year publication_year,
	  max(DS.tractable) SITES_tractable,
	  max(DS.druggable) SITES_druggable
	FROM PDB_Chains C
	  LEFT JOIN PDB P
		ON C.PDB_code = P.PDB_code
	  LEFT JOIN PDBChain_Domain Domain
		ON C.Chain_id = Domain.Chain_id
	  LEFT JOIN Domain_targets D
		ON Domain.Domain_id = D.domain_id
	  LEFT JOIN pdb_bind B
		ON B.pdb_code = C.PDB_code
	  LEFT JOIN drugEbility_sites DS
		ON DS.pdb_code = LOWER(C.PDB_code)
	WHERE C.Target_id in ('%s')
	GROUP BY C.Target_id,C.PDB_code,P.Technique,P.Resolution,C.n_residues,C.start_stop""" % target_id,
					'pockets': """SELECT
	  F.Target_id,
	  F.PDB_code,
	  F.DrugScore as druggability_score,
	  round(F.total_sasa,1) as area,
	  round(F.volume,1) as volume,
	  round((F.apolar_sasa/F.total_sasa)*100,1) as fraction_apolar,
	  F.Pocket_number as pocket_number,
	  F.Score as pocket_score,
	  GROUP_CONCAT((D.Domain_name||' ('||Domain.Coverage||'%)') ) as domains
	FROM fPockets F
	  LEFT JOIN fPockets_Domain Domain
		ON F.Pocket_id = Domain.Pocket_id
	  LEFT JOIN Domain_targets D
		ON Domain.Domain_id=D.domain_id
	WHERE F.Target_id in ('{target}')
	AND F.druggable='TRUE' AND F.blast='FALSE'
	GROUP BY F.Target_id,F.PDB_code,F.DrugScore,F.total_sasa,F.volume,fraction_apolar,pocket_number,pocket_score""".format(
						target=target_id),
					'alt_pockets': """SELECT F.Target_id,
	  F.PDB_code,
	  F.DrugScore as alt_druggability_score,
	  round(F.total_sasa,1) as alt_area,
	  round(F.volume,1) as alt_volume,
	  round((F.apolar_sasa/F.total_sasa)*100,1) as alt_fraction_apolar,
	  F.Pocket_number as alt_pocket_number,
	  F.Score as alt_pocket_score,
	  B.Hit_gene_name as alt_gene,
	  B.Hit_gene_species as alt_species,
	  B.similarity alt_similarity
	FROM fPockets F
	  LEFT JOIN `3D_Blast` B
		ON F.Target_id = B.Query_target_id AND F.PDB_code = B.Hit_PDB_code

	WHERE F.Target_id in ('%s')
	AND F.druggable='TRUE' AND F.blast='TRUE'
	ORDER BY B.similarity DESC""" % target_id,
					'bioactives': """SELECT C.target_id as Target_id,
	B.lig_id,
	  B.assay_id,
	  B.target_id as chembl_target_id,
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
	WHERE C.target_id in ('%s')
	AND B.operator!='>' AND B.operator!='<'
	AND A.confidence_score>=8""" % target_id,
					'commercials': """SELECT target_id,
	smiles,
	affinity_type,
	' =' as op,
	affinity_value,
	affinity_unit,
	price,
	website
	FROM purchasable_compounds
	WHERE target_id in ('%s')""" % target_id,
					'bindingDB': """SELECT B.target_id,
	  B.ligand_name,
	  B.ZincID,
	  B.`IC50(nM)`,
	  B.`EC50(nM)`,
	  B.`Ki(nM)`,
	  B.`Kd(nM)`,
	  B.`kon(M-1s-1)`,
	  B.`koff(s-1)`,
	  B.pH,
	  B.`Temp`,
	  B.Source,
	  B.DOI,
	  B.institution,
	  B.patent_number,
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
	  L.molecular_species,
	  L.num_ro5_violations as ro5_violations,
	  L.ro3_pass as pass_ro3,
	  B.ligand_smiles as SMILES
	FROM BindingDB B
	  LEFT JOIN ligands L
	  ON B.inchi_key = L.std_inchi_key
	WHERE target_id in ('%s')""" % target_id,
					'domain_drugE': """SELECT Target_id,
						max(tractable) tractable,
						max(druggable) druggable
					FROM
					(SELECT DISTINCT LOWER(PDB_code) as pdb_codes,Target_id
						FROM PDB_Chains
						WHERE target_id in ('%s')) T
					LEFT JOIN drugEbility_domains
					ON pdb_code = T.pdb_codes
					GROUP BY Target_id""" % target_id}
	results = {qname: pd.read_sql(query, con=connector_targetDB) for qname, query in list_queries.items()}


	data = results['gen_info'].drop(['Species', 'species_id', 'Sequence',
									 'Cell_location', 'Process', 'Function', 'Synonyms',
									 'Date_modified', 'Protein_class', 'Protein_class_desc',
									 'Protein_class_short', 'chembl_id'], axis=1)

	data = data.merge(results['kegg'].groupby('Target_id')['pathway_name'].apply('\n'.join).reset_index().rename(
		columns={'pathway_name': 'kegg_list'}), on='Target_id', how='left')

	data = data.merge(results['kegg'].groupby('Target_id')['pathway_name'].count().reset_index().rename(
		columns={'pathway_name': 'kegg_count'}), on='Target_id', how='left')

	data = data.merge(results['reactome'].groupby('Target_id')['pathway_name'].apply('\n'.join).reset_index().rename(
		columns={'pathway_name': 'reactome_list'}), on='Target_id', how='left')

	data = data.merge(results['reactome'].groupby('Target_id')['pathway_name'].count().reset_index().rename(
		columns={'pathway_name': 'reactome_count'}), on='Target_id', how='left')

	data = data.merge(results['disease'].groupby('Target_id')['disease_id'].count().reset_index().rename(
		columns={'disease_id': 'disease_count_uniprot'}), on='Target_id', how='left')
	if not results['pockets'].empty:
		data = data.merge(results['pockets'].groupby('Target_id').mean().add_prefix('mean_').reset_index().round(2), on='Target_id',
						  how='left')
		data = data.merge(results['pockets'].groupby('Target_id')['druggability_score'].std().reset_index().rename(
			columns={'druggability_score': 'stddev_druggability_score'}).round(2), on='Target_id', how='left')
	data = data.merge(results['pockets'].groupby('Target_id')['PDB_code'].nunique().reset_index().rename(
		columns={'PDB_code': 'pdb_with_druggable_pocket'}), on='Target_id', how='left')
	data = data.merge(results['pockets'].groupby('Target_id')['PDB_code'].count().reset_index().rename(
		columns={'PDB_code': 'druggable_pockets_total'}), on='Target_id', how='left')
	if not results['alt_pockets'].empty:
		data = data.merge(results['alt_pockets'].groupby('Target_id').mean().add_prefix('mean_').reset_index().round(2),
						  on='Target_id',
						  how='left')
		data = data.merge(results['alt_pockets'].groupby('Target_id')['alt_druggability_score'].std().reset_index().rename(
			columns={'alt_druggability_score': 'alt_stddev_druggability_score'}).round(2), on='Target_id', how='left')
	data = data.merge(results['alt_pockets'].groupby('Target_id')['alt_similarity'].max().reset_index().rename(
		columns={'alt_similarity': 'max_alt_similarity'}), on='Target_id', how='left')

	data = data.merge(results['alt_pockets'].groupby('Target_id')['PDB_code'].nunique().reset_index().rename(
		columns={'PDB_code': 'alt_pdb_with_druggable_pocket'}), on='Target_id', how='left')
	data = data.merge(results['alt_pockets'].groupby('Target_id')['PDB_code'].count().reset_index().rename(
		columns={'PDB_code': 'alt_druggable_pockets_total'}), on='Target_id', how='left')

	potent = results['bindingDB'][((results['bindingDB'][['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)']] <= 10) & (
		results['bindingDB'][['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)']].notna())).any(axis=1)]
	potent_max_phase = potent[potent['max_phase'] >= 2]
	data = data.merge(results['bindingDB'].groupby('target_id')['ligand_name'].count().reset_index().rename(
		columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_count'}), on='Target_id', how='left')
	data = data.merge(potent.groupby('target_id')['ligand_name'].count().reset_index().rename(
		columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_potent_count'}), on='Target_id', how='left')
	data = data.merge(potent_max_phase.groupby('target_id')['ligand_name'].count().reset_index().rename(
		columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_potent_phase2_count'}), on='Target_id', how='left')

	if not results['tissue_expression'].empty:
		tissue_grouped = results['tissue_expression'].groupby(['Target_id', 'tissue']).max().groupby(
			['Target_id', 'organ']).mean().round(1).reset_index()
		tissue = tissue_grouped.pivot(index='Target_id', columns='organ', values='value').reset_index()
		tissue_max = tissue_grouped.loc[tissue_grouped.groupby('Target_id').idxmax()['value']].rename(
				columns={'organ': 'tissue_max_expression', 'value': 'expression_max_tissue'})
	else:
		tissue = pd.DataFrame(columns=['Target_id','Brain','Endocrine_tissue','Female_tissue','Immune','Kidney','Liver_gallbladder','Lung','Male_tissue','Muscle_tissue','Pancreas','Skin','Soft_tissue','gitract'])
		tissue_max = pd.DataFrame(columns=['tissue_max_expression','expression_max_tissue','Target_id'])

	data = data.merge(tissue, on='Target_id', how='left')
	data = data.merge(results['selectivity'].round(2).rename(columns={'Selectivity_entropy': 'Expression_Selectivity'}),
					  on='Target_id', how='left')
	data = data.merge(tissue_max, on='Target_id', how='left')

	data = data.merge(
		results['var'].groupby('Target_id')['start'].count().reset_index().rename(columns={'start': 'variants_count'}),
		on='Target_id', how='left')
	data = data.merge(
		results['mut'].groupby('Target_id')['start'].count().reset_index().rename(columns={'start': 'mutants_count'}),
		on='Target_id', how='left')
	data = data.merge(
		results['gwas'].groupby('Target_id')['phenotype'].count().reset_index().rename(
			columns={'phenotype': 'gwas_count'}),
		on='Target_id', how='left')

	best = results['bioactives'][results['bioactives']['pX'].notnull()]
	total_bioact = best.groupby('Target_id')['lig_id'].nunique().reset_index().rename(
		columns={'lig_id': 'ChEMBL_bioactives_count'})
	best = best[best['pX'] >= 8]
	total_potent = best.groupby('Target_id')['lig_id'].nunique().reset_index().rename(
		columns={'lig_id': 'ChEMBL_bioactives_potent_count'})
	best = best[best['standard_type'].isin(['Ki', 'Kd'])]
	if not best.empty:
		best['pX'].fillna(-np.log10(best.value_num / 1000000000), inplace=True)
		query_lig = "','".join(best.lig_id.unique())
		query = """SELECT
					B.lig_id,
					B.Target_id,
					B.target_name,
					ROUND(AVG(B.value_num),2) avg_value,
					ROUND(STDDEV(B.value_num),2) sttdev,
					COUNT(*) n_values
					FROM bioactivities B
					  LEFT JOIN assays A
					  ON B.assay_id=A.assay_id
					WHERE B.operator='=' 
					  AND B.lig_id in ('%s')
					  AND A.bioactivity_type='Binding'
					  AND UPPER(B.standard_type) in ('KD','KI')
					  AND B.data_validity_comment is NULL
					  AND A.confidence_score>=8
			GROUP BY B.lig_id,B.Target_id""" % query_lig
		entropies = []
		binding_data = pd.read_sql(query, con=connector_targetDB)
		if not binding_data.empty:
			for name, group in binding_data.groupby('lig_id'):
				group = group[(group['sttdev'] < group['avg_value'])].copy()
				if group.empty:
					continue
				group['association'] = (1 / group.avg_value)
				group['association_prob'] = group.association / group.association.sum()
				best_target_id = group.loc[group['association_prob'].idxmax()]['Target_id']
				entropies.append({'Selectivity': round(sc.entropy(group.association_prob), 2), 'lig_id': name,
								  'number of other targets': len(group),
								  'targets name': ' / '.join(np.unique(group['target_name'].values)),
								  'best_target_id': best_target_id})

			entropy = pd.DataFrame(data=entropies)

			best = pd.merge(best, entropy, on='lig_id')

			total_moderate_selectivity = best[
				(best['Selectivity'] <= 2) & (best['best_target_id'] == best['chembl_target_id']) & (
							best['number of other targets'] > 1)].groupby(
				'Target_id')['lig_id'].nunique().reset_index().rename(
				columns={'lig_id': 'ChEMBL_bioactives_moderate_selectivity_count'})
			total_good_selectivity = best[(best['Selectivity'] <= 1.5) & (best['best_target_id'] == best['chembl_target_id']) & (
					best['number of other targets'] > 1)].groupby('Target_id')['lig_id'].nunique().reset_index().rename(
				columns={'lig_id': 'ChEMBL_bioactives_good_selectivity_count'})
			total_great_selectivity = best[
				(best['Selectivity'] <= 1) & (best['best_target_id'] == best['chembl_target_id']) & (
							best['number of other targets'] > 1)].groupby(
				'Target_id')['lig_id'].nunique().reset_index().rename(
				columns={'lig_id': 'ChEMBL_bioactives_great_selectivity_count'})
		else:
			total_moderate_selectivity = pd.DataFrame(
				{'Target_id': [], 'ChEMBL_bioactives_moderate_selectivity_count': []})
			total_good_selectivity = pd.DataFrame({'Target_id': [], 'ChEMBL_bioactives_good_selectivity_count': []})
			total_great_selectivity = pd.DataFrame({'Target_id': [], 'ChEMBL_bioactives_great_selectivity_count': []})
	else:
		total_moderate_selectivity = pd.DataFrame({'Target_id': [], 'ChEMBL_bioactives_moderate_selectivity_count': []})
		total_good_selectivity = pd.DataFrame({'Target_id': [], 'ChEMBL_bioactives_good_selectivity_count': []})
		total_great_selectivity = pd.DataFrame({'Target_id': [], 'ChEMBL_bioactives_great_selectivity_count': []})

	data = data.merge(total_bioact, on='Target_id', how='left')
	data = data.merge(total_potent, on='Target_id', how='left')
	data = data.merge(total_moderate_selectivity, on='Target_id', how='left')
	data = data.merge(total_good_selectivity, on='Target_id', how='left')
	data = data.merge(total_great_selectivity, on='Target_id', how='left')

	results['phenotype']['lethal'] = results['phenotype']['Phenotype'].str.contains('lethal|death', case=False)
	results['phenotype']['normal'] = results['phenotype']['Phenotype'].str.contains('no abnormal phenotype detected',
																					case=False)
	pheno_hetero_lethal = results['phenotype'][
		(results['phenotype']['lethal'] == True) & (results['phenotype']['zygosity'] == 'HETEROZYGOTE')].groupby(
		'Target_id')['Allele_symbol'].count().reset_index().rename(
		columns={'Allele_symbol': 'phenotypes_heterozygotes_lethal_count'})
	pheno_homo_lethal = results['phenotype'][
		(results['phenotype']['lethal'] == True) & (results['phenotype']['zygosity'] == 'HOMOZYGOTE')].groupby(
		'Target_id')[
		'Allele_symbol'].count().reset_index().rename(columns={'Allele_symbol': 'phenotypes_homozygotes_lethal_count'})
	pheno_hetero_normal = results['phenotype'][
		(results['phenotype']['normal'] == True) & (results['phenotype']['zygosity'] == 'HETEROZYGOTE')].groupby(
		'Target_id')['Allele_symbol'].count().reset_index().rename(
		columns={'Allele_symbol': 'phenotypes_heterozygotes_normal_count'})
	pheno_homo_normal = results['phenotype'][
		(results['phenotype']['normal'] == True) & (results['phenotype']['zygosity'] == 'HOMOZYGOTE')].groupby(
		'Target_id')[
		'Allele_symbol'].count().reset_index().rename(columns={'Allele_symbol': 'phenotypes_homozygotes_normal_count'})

	data = data.merge(pheno_hetero_lethal, on='Target_id', how='left')
	data = data.merge(pheno_homo_lethal, on='Target_id', how='left')
	data = data.merge(pheno_hetero_normal, on='Target_id', how='left')
	data = data.merge(pheno_homo_normal, on='Target_id', how='left')

	total_pdb = results['pdb'].groupby('Target_id')['PDB_code'].nunique().reset_index().rename(
		columns={'PDB_code': 'PDB_total_count'})
	pdb_w_lig = results['pdb'][results['pdb']['type_of_binder'].notnull()].groupby('Target_id')[
		'PDB_code'].nunique().reset_index().rename(columns={'PDB_code': 'PDB_with_Ligand_count'})
	pdb_sites_tractable = results['pdb'][results['pdb']['SITES_tractable'] == 1].groupby('Target_id')[
		'PDB_code'].nunique().reset_index().rename(columns={'PDB_code': 'PDB_sites_tractable_count'})
	pdb_sites_druggable = results['pdb'][results['pdb']['SITES_druggable'] == 1].groupby('Target_id')[
		'PDB_code'].nunique().reset_index().rename(columns={'PDB_code': 'PDB_sites_druggable_count'})
	pdb_blast_close = results['pdb_blast'].groupby('Query_target_id')['PDB_code'].nunique().reset_index().rename(
		columns={'PDB_code': 'PDB_blast_close_count', 'Query_target_id': 'Target_id'})
	pdb_blast_max_simil = results['pdb_blast'].groupby('Query_target_id')['similarity'].max().reset_index().rename(
		columns={'similarity': 'PDB_blast_max_similarity', 'Query_target_id': 'Target_id'})

	data = data.merge(total_pdb, on='Target_id', how='left')
	data = data.merge(pdb_w_lig, on='Target_id', how='left')
	data = data.merge(pdb_sites_tractable, on='Target_id', how='left')
	data = data.merge(pdb_sites_druggable, on='Target_id', how='left')
	data = data.merge(pdb_blast_close, on='Target_id', how='left')
	data = data.merge(pdb_blast_max_simil, on='Target_id', how='left')

	data = data.merge(
		results['domains'].groupby('Target_id')['start'].count().reset_index().rename(
			columns={'start': 'domains_count'}),
		on='Target_id', how='left')
	data = data.merge(
		results['domain_drugE'].rename(columns={'tractable': 'domain_tractable', 'druggable': 'domain_druggable'}),
		on='Target_id', how='left')

	data = data.merge(results['commercials'].groupby('target_id')['smiles'].count().reset_index().rename(
		columns={'smiles': 'commercial_total', 'target_id': 'Target_id'}), on='Target_id', how='left')
	data = data.merge(results['commercials'][results['commercials']['affinity_value'] <= 100].groupby('target_id')[
		'smiles'].count().reset_index().rename(
		columns={'smiles': 'commercial_potent_total', 'target_id': 'Target_id'}), on='Target_id', how='left')


	query_id = """SELECT * FROM tcrd_id WHERE Target_id in ('%s')""" % target_id
	tcrd_id = pd.read_sql(query_id, con=connector_targetDB)
	tcrd_id_list = "','".join([str(i) for i in tcrd_id['tcrd_id'].values.tolist()])


	tcrd_queries = {'target': """SELECT * FROM tcrd_target WHERE tcrd_id in ('%s') """ % tcrd_id_list,
				'tdl_info': """SELECT * FROM tcrd_info WHERE protein_id in ('%s')""" % tcrd_id_list,
				'patent': """SELECT * FROM tcrd_patent where protein_id in ('%s')""" % tcrd_id_list,
				'disease': """SELECT protein_id,disease_id,doid,score,name,parent
					FROM
					tcrd_disease
	WHERE protein_id in ('%s')
	ORDER BY score DESC""" % tcrd_id_list,
				'novelty': """SELECT score,protein_id as tcrd_id FROM tcrd_novelty WHERE protein_id in ('%s')""" % tcrd_id_list}
	tcrd_res = {qname: pd.read_sql(query, con=connector_targetDB) for qname, query in tcrd_queries.items()}

	tcrd_data = tcrd_id.copy()

	tcrd_data = tcrd_data.merge(tcrd_res['target'], on='tcrd_id', how='left')

	tcrd_res['tdl_info'].rename(columns={'protein_id': 'tcrd_id'},inplace=True)
	tcrd_res['tdl_info'] = tcrd_res['tdl_info'].round(2)

	tcrd_data = tcrd_data.merge(tcrd_res['tdl_info'], on='tcrd_id', how='left')

	tcrd_data = tcrd_data.merge(tcrd_res['patent'].groupby('protein_id')['count'].sum().reset_index().rename(
		columns={'count': 'total_patent_count', 'protein_id': 'tcrd_id'}), on='tcrd_id', how='left')
	tcrd_data = tcrd_data.merge(
		tcrd_res['patent'].iloc[tcrd_res['patent'].groupby('protein_id')['count'].idxmax()].rename(
			columns={'protein_id': 'tcrd_id', 'year': 'year_max_patents', 'count': 'count_patents_max_year'}), on='tcrd_id',
		how='left')

	disease_clean = tcrd_res['disease'][
		(tcrd_res['disease']['score'] > 1) & (~tcrd_res['disease']['doid'].isin(tcrd_res['disease']['parent'].unique()))]
	disease_list = disease_clean.groupby('protein_id')['name'].unique().apply('\n'.join).reset_index().rename(
		columns={'protein_id': 'tcrd_id', 'name': 'disease_list_tcrd'})
	disease_count = disease_clean.groupby('protein_id')['name'].nunique().reset_index().rename(
		columns={'protein_id': 'tcrd_id', 'name': 'disease_count_tcrd'})
	disease_max = disease_clean.loc[disease_clean.groupby('protein_id')['score'].idxmax().values].drop(
		['parent', 'doid', 'disease_id'], axis=1).rename(
		columns={'protein_id': 'tcrd_id', 'score': 'max_disease_score', 'name': 'name_max_disease'}).round(2)

	tcrd_data = tcrd_data.merge(disease_list, on='tcrd_id', how='left')
	tcrd_data = tcrd_data.merge(disease_count, on='tcrd_id', how='left')
	tcrd_data = tcrd_data.merge(disease_max, on='tcrd_id', how='left')

	tcrd_data = tcrd_data.merge(tcrd_res['novelty'].rename(columns={'score': 'novelty_score'}), on='tcrd_id', how='left')
	tcrd_data = tcrd_data.drop(['tcrd_id'], axis=1)

	data = data.merge(tcrd_data, on='Target_id',how='left')
	connector_targetDB.close()
	return data


def cap_score(df, cap_value):
	df2 = df.copy()
	df2[df2 > cap_value] = cap_value
	df2[df2 <= cap_value] = (df2 / cap_value) * 10
	return df2


def make_spider_plot(data, labels, target_name=''):
	fig = plt.figure(figsize=(8, 8))
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
	N = len(data)
	theta = np.arange(0, 2 * np.pi, 2 * np.pi / N)
	ax.bar(0, 1, bottom=9, width=2 * np.pi, color='r', linewidth=0, alpha=0.3)
	ax.bar(0, 5, bottom=4, width=2 * np.pi, color='lime', linewidth=0, alpha=0.2)
	ax.bar(0, 3, bottom=1, width=2 * np.pi, color='gold', linewidth=0, alpha=0.2)
	ax.bar(0, 1, width=2 * np.pi, color='r', linewidth=0)
	bars = ax.bar(theta, data, width=2 * np.pi / N, align='center')
	plt.title(target_name, y=1.08)
	ax.set_xticks(theta)
	ax.set_xticklabels(labels)
	ax.yaxis.grid(False)
	ax.set_yticks([])
	counter = 0
	for bar in bars:
		counter += 1
		if counter <= 3:
			bar.set_facecolor('mediumorchid')
		elif counter <= 7:
			bar.set_facecolor('cornflowerblue')
		elif counter <= 11:
			bar.set_facecolor('forestgreen')
		else:
			bar.set_facecolor('powderblue')

	buf = io.BytesIO()
	fig.savefig(buf, format='png')
	buf.seek(0)
	return buf


def make_score(df):
	col_to_keep = ['Ab Count', 'MAb Count', 'tcrd_disease_count', 'disease_count',
				   'gwas_count', 'variants_count', 'mutant_count', 'isoforms_count', 'domains_count',
				   'ChEMBL_bioactives_count', 'ChEMBL_bioactives_good_selectivity_count',
				   'ChEMBL_bioactives_great_selectivity_count',
				   'ChEMBL_bioactives_moderate_selectivity_count',
				   'ChEMBL_bioactives_potent_count', 'bindingDB_count',
				   'bindingDB_potent_count', 'commercials_potent_total',
				   'commercials_total',
				   'JensenLab PubMed Score', 'NCBI Gene PubMed Count', 'PubTator Score',
				   'PDB_blast_close_count',
				   'PDB_sites_druggable', 'PDB_sites_tractable', 'PDB_total_count',
				   'PDB_with_Ligand_count', 'druggable_pockets_total',
				   'number_druggable_pockets_NORM', 'phenotypes_count',
				   'phenotypes_heterozygotes_lethal_count',
				   'phenotypes_heterozygotes_normal_count',
				   'phenotypes_homozygotes_lethal_count',
				   'phenotypes_homozygotes_normal_count']
	df2 = df[col_to_keep]
	dflog2 = df2.copy()
	dfsqrt = df2.copy()
	dfsqrt = np.sqrt(dfsqrt)
	with np.errstate(divide='ignore'):
		dflog2 = np.log2(dflog2).replace(-np.inf, 0)
	spider_score = df2.copy()
	spider_score.drop(spider_score.columns, axis=1, inplace=True)
	spider_score['ChEMBL'] = (dflog2['ChEMBL_bioactives_count'] + dflog2['ChEMBL_bioactives_great_selectivity_count'] +
							  dflog2['ChEMBL_bioactives_potent_count']) / 2
	spider_score['ChEMBL'] = cap_score(spider_score['ChEMBL'], 10)
	spider_score['BindingDB'] = cap_score(dfsqrt['bindingDB_potent_count'], 7)
	spider_score['Commercial'] = cap_score(df2['commercials_potent_total'], 10)
	# spider_score['Screening'] = cap_score(dflog2['screens_compounds_total'], 7)
	spider_score['PDB_total'] = cap_score(dflog2['PDB_total_count'], 5)
	spider_score['PDB_with_Lig'] = (df2['PDB_with_Ligand_count'] / df2['PDB_total_count']) * 10
	spider_score['PDB_druggable'] = cap_score(
		(df2['druggable_pockets_total'] / (df2['PDB_total_count'] + df2['PDB_blast_close_count'])) * 10, 10)
	spider_score['PDB_blast'] = cap_score(dfsqrt['PDB_blast_close_count'], 6)
	spider_score['MAb'] = cap_score(dflog2['MAb Count'], 10)
	spider_score['Diseases'] = cap_score((dflog2['tcrd_disease_count'] + 2 * dflog2['disease_count']), 10)
	spider_score['GWAS'] = cap_score(df2['gwas_count'], 10)
	spider_score['phenotypes'] = cap_score(dflog2['phenotypes_count'], 5)
	spider_score['Pubmed Score'] = cap_score(dfsqrt['JensenLab PubMed Score'], 32)
	spider_score['Pubmed count'] = cap_score(dfsqrt['NCBI Gene PubMed Count'] / 2, 10)
	# spider_score['Grants'] = cap_score(dflog2['grants_total_count'], 10)
	# spider_score['Grants money'] = cap_score(
	#     ((np.sqrt(df2['grants_amount_total'] / df2['grants_total_count']) - 200) / 500) * 10, 10)
	# spider_score.loc[spider_score['Grants money'] < 0, 'Grants money'] = 0
	spider_score = spider_score.fillna(0)
	return spider_score


if __name__ == '__main__':
	pass
