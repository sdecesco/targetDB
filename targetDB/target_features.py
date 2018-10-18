#!/usr/bin/env python
import re, sqlite3, math
import pandas as pd
from targetDB import cns_mpo as mpo
import numpy as np
import scipy.stats as sc


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


def get_single_features(target_id, dbase=None):
	single_queries = {'general_info': "SELECT * FROM Targets WHERE Target_id='" + target_id + "'",
	                  'disease': "SELECT disease_name,disease_id FROM disease WHERE Target_id='" + target_id + "'",
	                  'reactome': "SELECT pathway_name FROM pathways WHERE pathway_dataset='Reactome pathways data set' AND Target_id='" + target_id + "'",
	                  'kegg': "SELECT pathway_name FROM pathways WHERE pathway_dataset='KEGG pathways data set' AND Target_id='" + target_id + "'",
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
	                  'organ_expression': """SELECT
              organ as organ_name,
              sum(value) as Total_value,
              count(value)as n_tissues,
              avg(value) as avg_value
              FROM protein_expression_levels
              WHERE Target_id='%s'
              GROUP BY organ
              ORDER BY avg_value DESC""" % target_id,
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
              (T.Gene_name || '-' || I.Isoform_name) as isoform_name,
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
              GROUP_CONCAT((D.Domain_name || ' (' || Domain.Coverage || '%)')) as domains
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
              B.target_id as target_id,
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
       ' =' as op,
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
      GROUP_CONCAT(DISTINCT UPPER(pdb_code)) pdb_list,
      domain_fold,
      domain_superfamily,
    max(tractable) tractable,
      max(druggable) druggable
    FROM drugEbility_domains
    WHERE pdb_code in (SELECT DISTINCT LOWER(PDB_code)
      FROM PDB_Chains
    WHERE target_id = '%s')
    GROUP BY domain_fold""" % target_id,
                      'open_target':"""SELECT * FROM opentarget_association WHERE target_id='%s'""" % target_id}
	connector = sqlite3.connect(dbase)
	connector.create_aggregate('stddev', 1, StdevFunc)
	results = {qname: pd.read_sql(query, con=connector) for qname, query in single_queries.items()}
	results.update(transform_bioactivities(results['bioactives'], connector))
	connector.close()
	return results


def transform_bioactivities(results, dbase):
	if results.empty:
		return {'binding': pd.DataFrame(), 'dose_response': pd.DataFrame(), 'other': pd.DataFrame(),
		        'ADME': pd.DataFrame(), 'emax': pd.DataFrame(),
		        'efficacy_bio': pd.DataFrame(), 'percent_inhibition': pd.DataFrame()}

	conc = re.compile(r'(?:of|at)\s(\d+\.*\d*)\s?((?:u|n)M)')
	bioactivity_types = ['Binding', 'Functionnal']
	percent = ['Activity', 'Residual activity', 'Residual_activity', 'Residual Activity', 'Inhibition']
	percent_invert = ['Activity', 'Residual activity', 'Residual Activity', 'Residual_activity']
	binding_affinity = ['Ki', 'Kd']
	dose_response_type = ['IC50', 'EC50', 'Potency']

	col = ['lig_id', 'standard_type', 'operator', 'value_num', 'units', 'pX', 'Conc', 'Conc_units',
	       'activity_comment', 'data_validity_comment', 'bioactivity_type', 'assay_species',
	       'assay_description',
	       'confidence_score', 'assay_id', 'SMILES', 'HBA', 'HBD', 'LogD', 'LogP', 'MW', 'TPSA', 'aLogP',
	       'apKa', 'bpKa', 'nAr', 'pass_ro3', 'ro5_violations', 'rotB', 'CNS_MPO', 'mol_name',
	       'molecular_species', 'indication_class', 'class_def', 'max_phase', 'oral', 'assay_ref', 'ref_bio',
	       'target_id']

	bioactives = results.copy()

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

	percent_bio = bioactives[
		(bioactives['units'] == '%') & (bioactives['standard_type'].isin(percent)) & (
			bioactives['bioactivity_type'].isin(bioactivity_types))].copy()

	for key in percent_invert:
		percent_bio.loc[percent_bio['standard_type'] == key, 'value_num'] = 100 - percent_bio['value_num']
		percent_bio.loc[percent_bio['standard_type'] == key, 'standard_type'] = '100 - ' + percent_bio[
			'standard_type']
	percent_bio = percent_bio[(percent_bio['value_num'] > 50)]
	percent_bio.sort_values(by='value_num', ascending=False, inplace=True)
	efficacy_bio = bioactives[
		(bioactives['units'] == '%') & (bioactives['standard_type'] == 'Efficacy')].copy()
	efficacy_bio = efficacy_bio[efficacy_bio.value_num >= 50]
	efficacy_bio.sort_values(by='value_num', ascending=False, inplace=True)
	emax = bioactives[(bioactives['units'] == '%') & (bioactives['standard_type'] == 'Emax') & (
		bioactives['bioactivity_type'].isin(bioactivity_types))].copy()
	emax = emax[emax.value_num >= 50]
	emax.sort_values(by='value_num', ascending=False, inplace=True)
	ADME = bioactives[(bioactives['bioactivity_type'] == 'ADME')].copy()
	ADME.sort_values(by='assay_description', inplace=True)
	other = bioactives[~(bioactives['standard_type'].isin(
		['Emax', 'Efficacy', 'Activity', 'Residual activity', 'Residual_activity', 'Residual Activity',
		 'Inhibition', 'IC50', 'Ki',
		 'EC50', 'Kd', 'Potency'])) & ~(bioactives['bioactivity_type'] == 'ADME')].copy()
	other.sort_values(by=['standard_type', 'assay_description'], inplace=True)
	dose_response = bioactives[
		(bioactives['units'] == 'nM') & (bioactives['standard_type'].isin(dose_response_type)) & (
			bioactives['bioactivity_type'].isin(bioactivity_types))].copy()
	dose_response = dose_response[dose_response.value_num <= 1000]
	dose_response.sort_values(by=['standard_type', 'value_num'], inplace=True)
	dose_response['pX'].fillna(-np.log10(dose_response.value_num / 1000000000),
	                           inplace=True)

	binding = bioactives[
		(bioactives['units'] == 'nM') & (bioactives['standard_type'].isin(binding_affinity)) & (
			bioactives['bioactivity_type'].isin(bioactivity_types))].copy()
	binding = binding[binding.value_num <= 1000]

	if not binding.empty:
		binding.sort_values(by=['standard_type', 'value_num'], inplace=True)
		binding['pX'].fillna(-np.log10(binding.value_num / 1000000000), inplace=True)

		query_lig = "','".join(binding.lig_id.unique())
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
		binding_data = pd.read_sql(query, con=dbase)
		best_target_id = binding.iloc[0]['target_id']
		if not binding_data.empty:
			for name, group in binding_data.groupby('lig_id'):
				best_target = True
				group = group[(group['sttdev'] < group['avg_value'])].copy()
				if group.empty:
					continue
				group['association'] = (1 / group.avg_value)
				group['association_prob'] = group.association / group.association.sum()
				if len(group) > 1:
					if group.loc[group['association_prob'].idxmax()]['Target_id'] == best_target_id:
						best_target = True
						best_target_name = group.loc[group['association_prob'].idxmax()]['target_name']
					else:
						best_target = False
						best_target_name = group.loc[group['association_prob'].idxmax()]['target_name']
				else:
					best_target_name = group.iloc[0]['target_name']
				entropies.append({'Selectivity': round(sc.entropy(group.association_prob), 2), 'lig_id': name,
				                  'number of other targets': len(group),
				                  'targets name': ' / '.join(np.unique(group['target_name'].values)),
				                  'best_target': best_target, 'best_target_name': best_target_name})

			entropy = pd.DataFrame(data=entropies)

			binding = pd.merge(binding, entropy, on='lig_id')

			col_order = ['lig_id', 'standard_type', 'operator', 'value_num', 'units', 'pX', 'Selectivity',
			             'number of other targets',
			             'best_target_name', 'activity_comment', 'bioactivity_type', 'assay_species',
			             'assay_description',
			             'confidence_score', 'assay_id', 'SMILES', 'HBA', 'HBD', 'LogD', 'LogP', 'MW',
			             'TPSA', 'aLogP', 'apKa', 'bpKa', 'nAr', 'pass_ro3',
			             'ro5_violations', 'rotB', 'CNS_MPO', 'mol_name', 'molecular_species',
			             'indication_class', 'class_def', 'max_phase', 'oral', 'assay_ref',
			             'ref_bio']
			binding = binding[col_order]
	return {'binding': binding, 'dose_response': dose_response, 'other': other, 'ADME': ADME, 'emax': emax,
	        'efficacy_bio': efficacy_bio, 'percent_inhibition': percent_bio}
