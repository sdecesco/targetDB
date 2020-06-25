#!/usr/bin/env python
import pandas as pd
import numpy as np
import io, sqlite3, math
import scipy.stats as sc
# start fix for macOS - as per psterk comment on github
from sys import platform

if platform == 'darwin':
    import matplotlib

    matplotlib.use('TkAgg')
# end fix for macOS
import matplotlib.pyplot as plt
from operator import itemgetter
from targetDB.utils import targetDB_gui as tgui


def norm_min_max(df, df_min=0.0, df_max=1.0):
    return ((df - df_min) * 1.0 / (df_max - df_min)).clip(0, 1)


def norm_min_median(df, median):
    return ((df - 0) * 1.0 / (median - 0)).clip(0, 1)


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


class target_scores:
    def __init__(self, data, mode='list'):
        self.mode = mode
        self.data = data
        self.coeff = {'sbio': 1, 'sdrug': 1, 'chem': 1, 'bio': 1, 'dise': 1, 'gen': 1, 'info': 1, 'safe': 1}
        self.scores = self.data.Target_id.reset_index().drop(['index'], axis=1)
        self.scores_quality = self.data.Target_id.reset_index().drop(['index'], axis=1)
        self.score_components = self.data.Target_id.reset_index().drop(['index'], axis=1)

        # ==============================================================================================#
        # ==================================== SCORE DATAFRAME =========================================#
        # ==============================================================================================#

        self.get_structure_info_score()
        self.get_structure_druggability_score()
        self.get_chem_score()
        self.get_bio_score()
        self.get_disease_score()
        self.get_genetic_score()
        self.get_info_score()
        self.get_safety_score()
        self.scores = self.scores.fillna(0)
        self.get_mpo_score()
        self.scores = self.scores.round(2)

    def get_mpo_score(self):
        if self.mode == 'list':
            values = tgui.get_mpo_coeff_gui()
            coeff_df = values.coeff
            self.coeff = coeff_df
            self.scores_mpo = self.scores.copy()
            self.scores_mpo.set_index('Target_id', inplace=True)
            self.scores_mpo = self.scores_mpo - 0.5
            self.scores_mpo = self.scores_mpo * 2
            self.scores_mpo['mpo_score'] = (self.scores_mpo.safety_score * coeff_df[
                'safe'] + self.scores_mpo.information_score * coeff_df['info'] + self.scores_mpo.genetic_score *
                                            coeff_df['gen'] + self.scores_mpo.disease_score * coeff_df[
                                                'dise'] + self.scores_mpo.biology_score * coeff_df[
                                                'bio'] + self.scores_mpo.chemistry_score * coeff_df[
                                                'chem'] + self.scores_mpo.structural_drug_score * coeff_df[
                                                'sdrug'] + self.scores_mpo.structure_info_score * coeff_df[
                                                'sbio']) / sum([abs(v) for v in coeff_df.values()])
            self.scores_mpo.mpo_score = self.scores_mpo.mpo_score / 2
            self.scores_mpo.mpo_score = self.scores_mpo.mpo_score + 0.5
            self.scores_mpo.reset_index(inplace=True)
            self.scores['mpo_score'] = self.scores_mpo.mpo_score
        else:
            self.scores['mpo_score'] = np.nan

    def get_structure_info_score(self):
        # ==============================================================================================#
        # ================================= STRUCTURE INFO SCORE =======================================#
        # ==============================================================================================#
        structural_info_cols = ["Target_id", "PDB_total_count", '%_sequence_covered', '%_domain_covered',
                                "PDB_blast_close_count", "PDB_blast_max_similarity"]
        structural_info_df = self.data[structural_info_cols].copy()
        structural_info_df['pdb_count_score'] = np.where(structural_info_df.PDB_total_count >= 3, 1, 0)
        structural_info_df['pdb_alt_count_score'] = np.where(((structural_info_df.PDB_blast_close_count * (
                structural_info_df.PDB_blast_max_similarity / 100)) / 2) >= 3, 1, 0)
        structural_info_df['structure_info_score'] = structural_info_df[
            ['%_sequence_covered', '%_domain_covered', 'pdb_count_score', 'pdb_alt_count_score']].mean(axis=1)
        struct_info_components = structural_info_df[
            ['Target_id', '%_sequence_covered', '%_domain_covered', 'pdb_count_score', 'pdb_alt_count_score']]
        struct_info_score = structural_info_df[['Target_id', "structure_info_score"]]
        self.scores = self.scores.merge(struct_info_score, on='Target_id', how='left')
        self.score_components = self.score_components.merge(struct_info_components, on='Target_id', how='left')

    def get_structure_druggability_score(self):
        structural_drug_cols = ["Target_id", "domain_tractable", "domain_druggable", "mean_druggability_score",
                                "mean_alt_druggability_score", "mean_alt_similarity"]
        structural_drug_df = self.data[structural_drug_cols].copy()
        structural_drug_df['domain_drug_score'] = structural_drug_df[["domain_tractable", "domain_druggable"]].mean(
            axis=1)
        structural_drug_df['alt_drug_score'] = structural_drug_df.mean_alt_druggability_score * (
                structural_drug_df.mean_alt_similarity / 100)
        structural_drug_df['structural_drug_score'] = structural_drug_df[
            ['mean_druggability_score', 'domain_drug_score', 'alt_drug_score']].mean(axis=1)
        self.scores_quality = self.scores_quality.merge(structural_drug_df[['Target_id', 'structural_drug_score']],
                                                        on='Target_id',
                                                        how='left')
        self.scores = self.scores.merge(structural_drug_df[['Target_id', 'structural_drug_score']], on='Target_id',
                                        how='left')
        self.score_components = self.score_components.merge(structural_drug_df[
                                                                ['Target_id', 'mean_druggability_score',
                                                                 'domain_drug_score', 'alt_drug_score']],
                                                            on='Target_id', how='left')

    def get_chem_score(self):
        chemistry_cols = ["Target_id", "BindingDB_potent_count", "BindingDB_potent_phase2_count",
                          "ChEMBL_bioactives_potent_count", "ChEMBL_bioactives_moderate_selectivity_count",
                          "ChEMBL_bioactives_good_selectivity_count", "ChEMBL_bioactives_great_selectivity_count",
                          "commercial_potent_total"]
        chemistry_df = self.data[chemistry_cols].copy()
        chemistry_df['bindingDB_potent'] = norm_min_max(np.log(chemistry_df.BindingDB_potent_count), df_max=8)
        chemistry_df['bindingDB_potent_log'] = norm_min_median(np.log(chemistry_df.BindingDB_potent_count), median=3.13)
        chemistry_df['bindingDB_phase2'] = norm_min_max(np.log(chemistry_df.BindingDB_potent_phase2_count), df_max=4.7)
        chemistry_df['chembl_potent'] = norm_min_max(np.log(chemistry_df.ChEMBL_bioactives_potent_count), df_max=8)
        chemistry_df['chembl_potent_log'] = norm_min_median(np.log(chemistry_df.ChEMBL_bioactives_potent_count),
                                                            median=2.64)
        chemistry_df['chembl_selective_M'] = norm_min_max(
            np.log(chemistry_df.ChEMBL_bioactives_moderate_selectivity_count), df_max=7)
        chemistry_df['chembl_selective_G'] = norm_min_max(np.log(chemistry_df.ChEMBL_bioactives_good_selectivity_count),
                                                          df_max=7)
        chemistry_df['chembl_selective_E'] = norm_min_max(
            np.log(chemistry_df.ChEMBL_bioactives_great_selectivity_count), df_max=7)
        chemistry_df['commercial_potent'] = norm_min_max(np.log(chemistry_df.commercial_potent_total), df_max=5)

        chemistry_df['chemistry_score'] = chemistry_df[['chembl_potent_log', 'bindingDB_potent_log']].mean(axis=1)
        chemistry_df['chemistry_qual_score'] = np.where(chemistry_df.BindingDB_potent_phase2_count > 0, 1,
                                                        np.where(
                                                            chemistry_df.ChEMBL_bioactives_great_selectivity_count > 0,
                                                            0.8,
                                                            np.where(
                                                                chemistry_df.ChEMBL_bioactives_good_selectivity_count > 0,
                                                                0.7,
                                                                np.where(
                                                                    chemistry_df.ChEMBL_bioactives_moderate_selectivity_count > 0,
                                                                    0.6,
                                                                    np.where(
                                                                        chemistry_df.ChEMBL_bioactives_potent_count > 0,
                                                                        0.3,
                                                                        np.where(
                                                                            chemistry_df.BindingDB_potent_count > 0,
                                                                            0.3,
                                                                            0))))))
        self.score_components = self.score_components.merge(chemistry_df[
                                                                ['Target_id', 'bindingDB_potent', 'bindingDB_phase2',
                                                                 'chembl_potent', 'chembl_selective_M',
                                                                 'chembl_selective_G', 'chembl_selective_E',
                                                                 'commercial_potent']], on='Target_id', how='left')
        self.scores = self.scores.merge(chemistry_df[['Target_id', 'chemistry_score']], on='Target_id', how='left')
        self.scores_quality = self.scores_quality.merge(chemistry_df[['Target_id', 'chemistry_qual_score']],
                                                        on='Target_id', how='left')

    def get_bio_score(self):
        biology_cols = ['Target_id', 'EXP_LVL_AVG', 'Ab Count', 'variants_count', 'mutants_count',
                        'number_of_genotypes',
                        'kegg_count', 'reactome_count']
        biology_df = self.data[biology_cols].copy()
        biology_df['bio_EScore'] = np.where(biology_df.EXP_LVL_AVG > 0, 1, 0)
        biology_df['bio_AScore'] = np.where(biology_df['Ab Count'] > 50, 1, 0)
        biology_df['bio_VScore'] = np.where(biology_df.variants_count > 0, 1, 0)
        biology_df['bio_MScore'] = np.where(biology_df.mutants_count > 0, 1, 0)
        biology_df['bio_GScore'] = np.where(biology_df.number_of_genotypes > 0, 1, 0)
        biology_df['bio_PScore'] = (np.where(biology_df.kegg_count > 0, 1, 0) + np.where(biology_df.reactome_count > 0,
                                                                                         1,
                                                                                         0)) / 2
        biology_df['biology_score'] = biology_df[
            ['bio_EScore', 'bio_AScore', 'bio_VScore', 'bio_MScore', 'bio_GScore', 'bio_PScore']].mean(
            axis=1)

        self.scores = self.scores.merge(biology_df[['Target_id', 'biology_score']], on='Target_id', how='left')
        self.score_components = self.score_components.merge(biology_df[
                                                                ['Target_id', 'bio_EScore', 'bio_AScore', 'bio_VScore',
                                                                 'bio_MScore', 'bio_GScore', 'bio_PScore']])

    def get_disease_score(self):
        disease_cols = ['Target_id', 'disease_count_uniprot', 'disease_count_tcrd', 'OT_number_of_disease_areas',
                        'OT_max_association_score']
        disease_df = self.data[disease_cols].copy()
        disease_df['dis_AScore'] = (np.where(disease_df.OT_number_of_disease_areas > 0, 1, 0) + np.where(
            disease_df.OT_number_of_disease_areas > 1, 1, 0)) / 2
        disease_df['dis_UScore'] = np.where(disease_df.disease_count_uniprot > 0, 1, 0)
        disease_df['dis_TScore'] = np.where(disease_df.disease_count_tcrd > 0, 1, 0)
        disease_df['disease_score'] = (disease_df[
                                           'OT_max_association_score'] * 2 + disease_df.dis_AScore + disease_df.dis_UScore + disease_df.dis_TScore) / 5
        disease_score_component = disease_df[
            ['Target_id', 'OT_max_association_score', 'dis_AScore', 'dis_UScore', 'dis_TScore']]
        self.scores = self.scores.merge(disease_df[['Target_id', 'disease_score']], on='Target_id', how='left')
        self.score_components = self.score_components.merge(disease_score_component, on='Target_id', how='left')

    def get_genetic_score(self):
        col_genetic = ['Target_id', 'gwas_count', 'gwas_sig', 'OT_MAX_VAL_genetic_association',
                       'OT_NUM_MAX_genetic_association', 'OT_TOP10_avg_score']
        genetic_df = self.data[col_genetic].copy()
        genetic_df['genetic_NORM'] = norm_min_max(np.log2(genetic_df.OT_NUM_MAX_genetic_association.replace(0, 0.0001)),
                                                  df_max=5)
        genetic_df['gen_GScore'] = norm_min_max(np.log10(genetic_df.gwas_count * 10), df_min=0, df_max=2)
        genetic_df['gen_GQualScore'] = genetic_df.gwas_sig / genetic_df.gwas_count
        genetic_df['gen_AScore'] = genetic_df.OT_MAX_VAL_genetic_association * genetic_df.genetic_NORM
        genetic_df['gen_AQualScore'] = genetic_df.OT_TOP10_avg_score
        genetic_df['genetic_score'] = genetic_df[['gen_GScore', 'gen_AScore']].mean(axis=1)
        genetic_df['genetic_score_qual'] = genetic_df[['gen_GQualScore', 'gen_AQualScore']].mean(axis=1)

        self.score_components = self.score_components.merge(
            genetic_df[['Target_id', 'gen_AScore', 'gen_GScore', 'gen_GQualScore', 'gen_AQualScore','genetic_NORM']],
            on='Target_id', how='left')
        self.scores = self.scores.merge(genetic_df[['Target_id', 'genetic_score']], on='Target_id', how='left')
        self.scores_quality = self.scores_quality.merge(genetic_df[['Target_id', 'genetic_score_qual']],
                                                        on='Target_id', how='left')

    def get_info_score(self):
        col_info = ['Target_id', 'JensenLab PubMed Score']
        info_df = self.data[col_info].copy()
        info_df['information_score'] = norm_min_max(np.log(info_df['JensenLab PubMed Score'].replace(0, 0.00001)),
                                                    df_max=12)

        self.scores = self.scores.merge(info_df[['Target_id', 'information_score']], on='Target_id', how='left')
        self.score_components = self.score_components.merge(info_df[['Target_id', 'information_score']], on='Target_id',
                                                            how='left')

    def get_safety_score(self):
        col_safety = ['Target_id', 'Heart_alert', 'Liver_alert', 'Kidney_alert', 'number_of_genotypes',
                      'phenotypes_heterozygotes_lethal_count', 'phenotypes_heterozygotes_normal_count',
                      'phenotypes_homozygotes_lethal_count', 'phenotypes_homozygotes_normal_count']
        safety_df = self.data[col_safety].copy()

        safety_df['safe_GScore'] = ((safety_df.phenotypes_homozygotes_lethal_count.fillna(
            0) * 2 + safety_df.phenotypes_heterozygotes_lethal_count.fillna(
            0) - safety_df.phenotypes_heterozygotes_normal_count.fillna(
            0) - 2 * safety_df.phenotypes_homozygotes_normal_count.fillna(0)) / safety_df.number_of_genotypes).clip(0,
                                                                                                                    1)

        safety_df['safe_EScore'] = np.where(safety_df[['Heart_alert', 'Liver_alert', 'Kidney_alert']].any(axis=1), 1,
                                            np.where(
                                                safety_df[['Heart_alert', 'Liver_alert', 'Kidney_alert']].isna().all(
                                                    axis=1),
                                                np.nan, 0))
        safety_df['safety_qual'] = safety_df[['safe_GScore', 'safe_EScore']].mean(axis=1)
        safety_df['n_genotypes'] = norm_min_max(np.log(safety_df.number_of_genotypes), df_max=6)
        safety_df['safety_score'] = (safety_df.n_genotypes + np.where(self.data.Brain.isna(), 0, 0.3)).clip(0, 1)
        safety_score_component = safety_df[
            ['Target_id', 'safe_GScore', 'safe_EScore','n_genotypes', 'Heart_alert', 'Liver_alert', 'Kidney_alert']].copy()
        rep_dict = {'True': 1, 'False': 0}
        safety_score_component.replace(rep_dict, inplace=True)

        self.score_components = self.score_components.merge(safety_score_component, on='Target_id', how='left')
        self.scores = self.scores.merge(safety_df[['Target_id', 'safety_score']], on='Target_id', how='left')
        self.scores_quality = self.scores_quality.merge(safety_df[['Target_id', 'safety_qual']], on='Target_id',
                                                        how='left')


def range_list(x):
    out = []
    for y in x:
        for i in y.split('|'):
            j = i.strip().split('-')
            if '' not in j:
                k = (int(j[0]), int(j[1]))
                out.append(k)
    longest = remove_overlap(out)
    return longest


def max_range(x, dx=0.1):
    x.index = x.target_id
    x = x.drop(['target_id'], axis=1)
    x = x[(x >= (x.max() - dx)) & (x != 0)]
    return x.count()


def domain_ranges(x):
    out = []
    for y in x:
        out.append(y)
    out_sorted = sorted(out, key=itemgetter(0))
    return out_sorted


def remove_overlap(ranges):
    ranges_sorted = sorted(ranges, key=itemgetter(0))
    merged = []

    for higher in ranges_sorted:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1] or higher[0] == lower[1] + 1:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


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
                    GROUP BY Target_id""" % target_id,
                    'opentargets': """SELECT * FROM opentarget_association WHERE target_id in ('%s')""" % target_id}

    results = {qname: pd.read_sql(query, con=connector_targetDB) for qname, query in list_queries.items()}

    # ================================================================================================================#
    # ========================================= GENERAL INFO =========================================================#
    # ================================================================================================================#

    data = results['gen_info'].drop(['Species', 'species_id', 'Sequence',
                                     'Cell_location', 'Process', 'Function',
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

    # ================================================================================================================#
    # ========================================= OPENTARGETS ==========================================================#
    # ================================================================================================================#

    OT_diseases = results['opentargets'][(results['opentargets']['disease_area'] != '')]
    if OT_diseases.empty:
        OT = pd.DataFrame(
            columns=["Target_id", "OT_number_of_associations", "OT_list_max_diseases", "OT_TOP10_diseases",
                     "OT_max_association_score", "OT_%_genetic_association", "OT_%_known_drug",
                     "OT_%_litterature_mining",
                     "OT_%_animal_model", "OT_%_affected_pathway", "OT_%_rna_expression", "OT_%_somatic_mutation",
                     "OT_MAX_VAL_genetic_association", "OT_NUM_MAX_genetic_association", "OT_MAX_VAL_known_drug",
                     "OT_NUM_MAX_known_drug", "OT_MAX_VAL_litterature_mining", "OT_NUM_MAX_litterature_mining",
                     "OT_MAX_VAL_animal_model", "OT_NUM_MAX_animal_model", "OT_MAX_VAL_affected_pathway",
                     "OT_NUM_MAX_affected_pathway", "OT_MAX_VAL_rna_expression", "OT_NUM_MAX_rna_expression",
                     "OT_MAX_VAL_somatic_mutation", "OT_NUM_MAX_somatic_mutation"])
    else:
        OT_n_associations = OT_diseases.groupby('target_id')['disease_area'].count().reset_index().rename(
            columns={'target_id': 'Target_id', 'disease_area': 'OT_number_of_associations'})
        OT_disease_max = OT_diseases.loc[OT_diseases.groupby('target_id').idxmax().overall_score][
            ['target_id', 'overall_score']].rename(
            columns={'target_id': 'Target_id', 'overall_score': 'OT_max_association_score'})
        OT_max_disease_list = OT_diseases.groupby('target_id').apply(
            lambda x: x[x.overall_score == x.overall_score.max()]).reset_index(drop=True).groupby('target_id')[
            'disease_name'].apply(','.join).reset_index().rename(
            columns={'target_id': 'Target_id', 'disease_name': 'OT_list_max_diseases'})
        OT_disease_association_type_proportion = OT_diseases.groupby('target_id').apply(
            lambda x: x[x != 0].count() / x[x != 0]['target_id'].count()).drop(
            ['target_id', 'disease_area', 'disease_name', 'overall_score'], axis=1).reset_index().rename(
            columns={'target_id': 'Target_id',
                     'genetic_association': 'OT_%_genetic_association', 'known_drug': 'OT_%_known_drug',
                     'litterature_mining': 'OT_%_litterature_mining',
                     'animal_model': 'OT_%_animal_model', 'affected_pathway': 'OT_%_affected_pathway',
                     'rna_expression': 'OT_%_rna_expression',
                     'somatic_mutation': 'OT_%_somatic_mutation'}).round(2)
        OT_diseases_associations_max_count = OT_diseases.drop(['disease_area', 'disease_name', 'overall_score'],
                                                              axis=1).groupby('target_id').apply(
            lambda x: max_range(x, dx=0.1)).reset_index().rename(
            columns={'target_id': 'Target_id',
                     'genetic_association': 'OT_NUM_MAX_genetic_association', 'known_drug': 'OT_NUM_MAX_known_drug',
                     'litterature_mining': 'OT_NUM_MAX_litterature_mining',
                     'animal_model': 'OT_NUM_MAX_animal_model', 'affected_pathway': 'OT_NUM_MAX_affected_pathway',
                     'rna_expression': 'OT_NUM_MAX_rna_expression',
                     'somatic_mutation': 'OT_NUM_MAX_somatic_mutation'})
        OT_diseases_associations_max_val = OT_diseases.groupby('target_id').max().drop(
            ['disease_area', 'disease_name', 'overall_score'], axis=1).reset_index().rename(
            columns={'target_id': 'Target_id',
                     'genetic_association': 'OT_MAX_VAL_genetic_association', 'known_drug': 'OT_MAX_VAL_known_drug',
                     'litterature_mining': 'OT_MAX_VAL_litterature_mining',
                     'animal_model': 'OT_MAX_VAL_animal_model', 'affected_pathway': 'OT_MAX_VAL_affected_pathway',
                     'rna_expression': 'OT_MAX_VAL_rna_expression',
                     'somatic_mutation': 'OT_MAX_VAL_somatic_mutation'})
        OT_top10 = OT_diseases.groupby('target_id').apply(lambda x: x.loc[x.overall_score.nlargest(10).index])
        OT_top10_qual = OT_diseases.groupby('target_id').apply(lambda x: x.loc[x.overall_score.nlargest(10).index])
        OT_top10_qual = OT_top10_qual['overall_score'].reset_index().groupby('target_id')[
            'overall_score'].mean().reset_index().rename(
            columns={'target_id': 'Target_id', 'overall_score': 'OT_TOP10_avg_score'})
        OT_top10 = OT_top10[OT_top10.overall_score >= 0.1]
        OT_top10 = OT_top10['disease_name'].reset_index().groupby('target_id')['disease_name'].apply(
            ','.join).reset_index().rename(columns={'target_id': 'Target_id', 'disease_name': 'OT_TOP10_diseases'})
        OT = OT_n_associations
        OT = OT.merge(OT_disease_max, on='Target_id', how='left')
        OT = OT.merge(OT_max_disease_list, on='Target_id', how='left')
        OT = OT.merge(OT_disease_association_type_proportion, on='Target_id', how='left')
        OT = OT.merge(OT_diseases_associations_max_count, on='Target_id', how='left')
        OT = OT.merge(OT_diseases_associations_max_val, on='Target_id', how='left')
        OT = OT.merge(OT_top10, on='Target_id', how='left')
        OT = OT.merge(OT_top10_qual, on='Target_id', how='left')

    OT_disease_areas = results['opentargets'][(results['opentargets']['disease_area'] == '')]

    if OT_disease_areas.empty:
        OT = OT.merge(pd.DataFrame(columns=["Target_id", "OT_number_of_disease_areas", "OT_list_max_disease_area",
                                            "OT_max_association_diseaseArea_score"]), on='Target_id', how='left')
    else:
        OT_n_disease_areas = OT_disease_areas.groupby('target_id')['disease_area'].count().reset_index().rename(
            columns={'target_id': 'Target_id', 'disease_area': 'OT_number_of_disease_areas'})
        OT_disease_area_max = OT_disease_areas.loc[OT_disease_areas.groupby('target_id').idxmax().overall_score][
            ['target_id', 'overall_score']].rename(
            columns={'target_id': 'Target_id', 'overall_score': 'OT_max_association_diseaseArea_score'})
        OT_max_disease_area_list = OT_disease_areas.groupby('target_id').apply(
            lambda x: x[x.overall_score == x.overall_score.max()]).reset_index(drop=True).groupby('target_id')[
            'disease_name'].apply(','.join).reset_index().rename(
            columns={'target_id': 'Target_id', 'disease_name': 'OT_list_max_disease_area'})
        OT = OT.merge(OT_n_disease_areas, on='Target_id', how='left')
        OT = OT.merge(OT_disease_area_max, on='Target_id', how='left')
        OT = OT.merge(OT_max_disease_area_list, on='Target_id', how='left')

    col_order = ["Target_id", "OT_number_of_associations", "OT_number_of_disease_areas", "OT_list_max_disease_area",
                 "OT_max_association_diseaseArea_score", "OT_list_max_diseases", "OT_TOP10_diseases",
                 'OT_TOP10_avg_score',
                 "OT_max_association_score", "OT_%_genetic_association", "OT_%_known_drug", "OT_%_litterature_mining",
                 "OT_%_animal_model", "OT_%_affected_pathway", "OT_%_rna_expression", "OT_%_somatic_mutation",
                 "OT_MAX_VAL_genetic_association", "OT_NUM_MAX_genetic_association", "OT_MAX_VAL_known_drug",
                 "OT_NUM_MAX_known_drug", "OT_MAX_VAL_litterature_mining", "OT_NUM_MAX_litterature_mining",
                 "OT_MAX_VAL_animal_model", "OT_NUM_MAX_animal_model", "OT_MAX_VAL_affected_pathway",
                 "OT_NUM_MAX_affected_pathway", "OT_MAX_VAL_rna_expression", "OT_NUM_MAX_rna_expression",
                 "OT_MAX_VAL_somatic_mutation", "OT_NUM_MAX_somatic_mutation"]
    OT = OT[col_order]

    data = data.merge(OT, on='Target_id', how='left')

    # ================================================================================================================#
    # ============================================= POCKETS ==========================================================#
    # ================================================================================================================#

    if not results['pockets'].empty:
        data = data.merge(results['pockets'].groupby('Target_id').mean().add_prefix('mean_').reset_index().round(2),
                          on='Target_id',
                          how='left')
        data = data.merge(results['pockets'].groupby('Target_id')['druggability_score'].std().reset_index().rename(
            columns={'druggability_score': 'stddev_druggability_score'}).round(2), on='Target_id', how='left')
        data = data.merge(results['pockets'].groupby('Target_id')['PDB_code'].nunique().reset_index().rename(
            columns={'PDB_code': 'pdb_with_druggable_pocket'}), on='Target_id', how='left')
        data = data.merge(results['pockets'].groupby('Target_id')['PDB_code'].count().reset_index().rename(
            columns={'PDB_code': 'druggable_pockets_total'}), on='Target_id', how='left')
    else:
        pocket = pd.DataFrame(
            columns=["Target_id", "mean_druggability_score", "stddev_druggability_score", "mean_area", "mean_volume",
                     "mean_fraction_apolar", "mean_pocket_score", "pdb_with_druggable_pocket", "druggable_pockets_total"
                     ])
        data = data.merge(pocket, on="Target_id", how='left')

    if not results['alt_pockets'].empty:
        data = data.merge(results['alt_pockets'].groupby('Target_id').mean().add_prefix('mean_').reset_index().round(2),
                          on='Target_id',
                          how='left')
        data = data.merge(
            results['alt_pockets'].groupby('Target_id')['alt_druggability_score'].std().reset_index().rename(
                columns={'alt_druggability_score': 'alt_stddev_druggability_score'}).round(2), on='Target_id',
            how='left')
        data = data.merge(results['alt_pockets'].groupby('Target_id')['alt_similarity'].max().reset_index().rename(
            columns={'alt_similarity': 'max_alt_similarity'}), on='Target_id', how='left')

        data = data.merge(results['alt_pockets'].groupby('Target_id')['PDB_code'].nunique().reset_index().rename(
            columns={'PDB_code': 'alt_pdb_with_druggable_pocket'}), on='Target_id', how='left')
        data = data.merge(results['alt_pockets'].groupby('Target_id')['PDB_code'].count().reset_index().rename(
            columns={'PDB_code': 'alt_druggable_pockets_total'}), on='Target_id', how='left')
    else:
        alt_pocket = pd.DataFrame(
            columns=["Target_id", "mean_alt_druggability_score", "alt_stddev_druggability_score", "mean_alt_area",
                     "mean_alt_volume", "mean_alt_fraction_apolar", "mean_alt_pocket_score", "mean_alt_similarity",
                     "max_alt_similarity", "alt_pdb_with_druggable_pocket", "alt_druggable_pockets_total"])
        data = data.merge(alt_pocket, on="Target_id", how='left')

    # ================================================================================================================#
    # ========================================= BINDING DB ===========================================================#
    # ================================================================================================================#

    results['bindingDB']['IC50(nM)'] = pd.to_numeric(results['bindingDB']['IC50(nM)'], errors='coerce')
    results['bindingDB']['EC50(nM)'] = pd.to_numeric(results['bindingDB']['EC50(nM)'], errors='coerce')
    results['bindingDB']['Ki(nM)'] = pd.to_numeric(results['bindingDB']['Ki(nM)'], errors='coerce')
    results['bindingDB']['Kd(nM)'] = pd.to_numeric(results['bindingDB']['Kd(nM)'], errors='coerce')

    potent = results['bindingDB'][((results['bindingDB'][['IC50(nM)', 'EC50(nM)', 'Kd(nM)', 'Ki(nM)']] <= 100) & (~
                                                                                                                  results[
                                                                                                                      'bindingDB'][
                                                                                                                      [
                                                                                                                          'IC50(nM)',
                                                                                                                          'EC50(nM)',
                                                                                                                          'Kd(nM)',
                                                                                                                          'Ki(nM)']].isnull())).any(
        axis=1)]
    potent_max_phase = potent[potent['max_phase'] >= 2]
    data = data.merge(results['bindingDB'].groupby('target_id')['ligand_name'].count().reset_index().rename(
        columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_count'}), on='Target_id', how='left')
    data = data.merge(potent.groupby('target_id')['ligand_name'].count().reset_index().rename(
        columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_potent_count'}), on='Target_id', how='left')
    data = data.merge(potent_max_phase.groupby('target_id')['ligand_name'].count().reset_index().rename(
        columns={'target_id': 'Target_id', 'ligand_name': 'BindingDB_potent_phase2_count'}), on='Target_id', how='left')

    # ================================================================================================================#
    # ========================================= TISSUE EXPRESSION ====================================================#
    # ================================================================================================================#

    if not results['tissue_expression'].empty:
        cell_values = results['tissue_expression'].groupby(['Target_id', 'organ', 'tissue']).max()
        cell_values.value = cell_values.value.replace(0, np.nan)
        tissue_grouped = cell_values.groupby(['Target_id', 'organ']).mean().round(1).reset_index().fillna(0)
        tissue = tissue_grouped.pivot(index='Target_id', columns='organ', values='value').reset_index()
        tissue_avg = tissue_grouped.groupby('Target_id')['value'].agg(['mean', 'std']).rename(
            columns={'mean': 'EXP_LVL_AVG', 'std': 'EXP_LVL_STDDEV'}).reset_index()
        tissue_avg.index = tissue_avg.Target_id
        tissue_avg['Heart_alert'] = False
        tissue_avg['Heart_value'] = np.nan
        tissue_avg['Liver_alert'] = False
        tissue_avg['Liver_value'] = np.nan
        tissue_avg['Kidney_alert'] = False
        tissue_avg['Kidney_value'] = np.nan
        for tid in tissue_avg.index:
            avg_exp_1stddev = tissue_avg.loc[tid].EXP_LVL_AVG + tissue_avg.loc[tid].EXP_LVL_STDDEV
            try:
                if cell_values.loc[(tid, 'Kidney', 'kidney')].value > avg_exp_1stddev or cell_values.loc[
                    (tid, 'Kidney', 'kidney')].value == 3:
                    tissue_avg.loc[tid, ['Kidney_alert']] = True
                    tissue_avg.loc[tid, ['Kidney_value']] = cell_values.loc[(tid, 'Kidney', 'kidney')].value
            except KeyError:
                pass
            try:
                if cell_values.loc[(tid, 'Liver_gallbladder', 'liver')].value > avg_exp_1stddev or cell_values.loc[
                    (tid, 'Liver_gallbladder', 'liver')].value == 3:
                    tissue_avg.loc[tid, ['Liver_alert']] = True
                    tissue_avg.loc[tid, ['Liver_value']] = cell_values.loc[(tid, 'Liver_gallbladder', 'liver')].value
            except KeyError:
                pass
            try:
                if cell_values.loc[(tid, 'Muscle_tissue', 'heart muscle')].value > avg_exp_1stddev or cell_values.loc[
                    (tid, 'Muscle_tissue', 'heart muscle')].value == 3:
                    tissue_avg.loc[tid, ['Heart_alert']] = True
                    tissue_avg.loc[tid, ['Heart_value']] = cell_values.loc[(tid, 'Muscle_tissue', 'heart muscle')].value
            except KeyError:
                pass
        tissue_avg.reset_index(drop=True, inplace=True)
        tissue_max = tissue_grouped.loc[tissue_grouped.groupby('Target_id').idxmax()['value']].rename(
            columns={'organ': 'tissue_max_expression', 'value': 'expression_max_tissue'})
    else:
        tissue = pd.DataFrame(
            columns=['Target_id','Brain', 'Adipose & soft tissue',
                 'Bone marrow & lymphoid tissues', 'Endocrine tissues',
                 'Female tissues', 'Gastrointestinal tract', 'Kidney & urinary bladder',
                 'Liver & gallbladder', 'Lung', 'Male tissues', 'Muscle tissues',
                 'Pancreas', 'Proximal digestive tract', 'Skin'])
        tissue_avg = pd.DataFrame(columns=['Target_id', 'EXP_LVL_AVG', 'EXP_LVL_STDDEV', 'Heart_alert',
                                           'Heart_value', 'Liver_alert', 'Liver_value', 'Kidney_alert',
                                           'Kidney_value'])
        tissue_max = pd.DataFrame(columns=['tissue_max_expression', 'expression_max_tissue', 'Target_id'])

    data = data.merge(tissue, on='Target_id', how='left')
    data = data.merge(results['selectivity'].round(2).rename(columns={'Selectivity_entropy': 'Expression_Selectivity'}),
                      on='Target_id', how='left')
    data = data.merge(tissue_avg, on='Target_id', how='left')
    data = data.merge(tissue_max, on='Target_id', how='left')

    # ================================================================================================================#
    # =================================== VARIANTS/MUTANTS/GWAS ======================================================#
    # ================================================================================================================#

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
    data = data.merge(
        results['gwas'][results['gwas']['p_value'] <= 0.000000005].groupby('Target_id')[
            'phenotype'].count().reset_index().rename(
            columns={'phenotype': 'gwas_sig'}),
        on='Target_id', how='left')

    # ================================================================================================================#
    # ========================================= CHEMBL BIOACTIVITIES =================================================#
    # ================================================================================================================#

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

            if not entropy.empty:
                best = pd.merge(best, entropy, on='lig_id')
            else:
                best['Selectivity'] = np.nan
                best['lig_id'] = np.nan
                best['number of other targets'] = np.nan
                best['targets name'] = np.nan
                best['best_target_id'] = np.nan

            total_moderate_selectivity = best[
                (best['Selectivity'] <= 2) & (best['best_target_id'] == best['chembl_target_id']) & (
                        best['number of other targets'] > 1)].groupby(
                'Target_id')['lig_id'].nunique().reset_index().rename(
                columns={'lig_id': 'ChEMBL_bioactives_moderate_selectivity_count'})
            total_good_selectivity = \
                best[(best['Selectivity'] <= 1.5) & (best['best_target_id'] == best['chembl_target_id']) & (
                        best['number of other targets'] > 1)].groupby('Target_id')[
                    'lig_id'].nunique().reset_index().rename(
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

    # ================================================================================================================#
    # ===================================== GENOTYPE/PHENOTYPES ======================================================#
    # ================================================================================================================#

    results['phenotype']['lethal'] = results['phenotype']['Phenotype'].str.contains('lethal|death', case=False)
    results['phenotype']['normal'] = results['phenotype']['Phenotype'].str.contains('no abnormal phenotype detected',
                                                                                    case=False)
    number_of_genotypes = results['phenotype'].groupby(['Target_id'])['genotype'].nunique().reset_index().rename(
        columns={'genotype': 'number_of_genotypes'})
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

    data = data.merge(number_of_genotypes, on='Target_id', how='left')
    data = data.merge(pheno_hetero_lethal, on='Target_id', how='left')
    data = data.merge(pheno_homo_lethal, on='Target_id', how='left')
    data = data.merge(pheno_hetero_normal, on='Target_id', how='left')
    data = data.merge(pheno_homo_normal, on='Target_id', how='left')

    # ================================================================================================================#
    # ================================================= PDB ==========================================================#
    # ================================================================================================================#

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

    general_data = pd.read_sql("""SELECT Target_id,Sequence FROM Targets WHERE Target_id in ('%s')""" % target_id,
                               con=connector_targetDB)
    general_data['length'] = general_data.Sequence.str.len()
    general_data.index = general_data.Target_id

    results['domains']['ranges'] = list(zip(results['domains'].start, results['domains'].stop))
    domain_range = results['domains'].groupby('Target_id')['ranges'].apply(lambda x: domain_ranges(x)).reset_index()
    domain_range.index = domain_range.Target_id
    pdb_range = results['pdb'].groupby('Target_id')['start_stop'].apply(lambda x: range_list(x)).reset_index()
    pdb_range.index = pdb_range.Target_id
    pdb_range['n_length'] = np.nan
    pdb_range['%_sequence_covered'] = np.nan
    pdb_range['%_domain_covered'] = np.nan
    for tid in pdb_range.index:
        n_residues = 0
        domain_coverages = []
        for i in pdb_range.loc[tid].start_stop:
            if i[0] > general_data.loc[tid].length:
                continue
            n_residues += (i[1] - i[0] + 1)
            if tid in domain_range.index:
                for j in domain_range.loc[tid]['ranges']:
                    d_range = range(j[0], j[1] + 1)
                    p_range = range(i[0], i[1] + 1)
                    overlap = len(set(d_range).intersection(p_range)) / len(d_range)
                    domain_coverages.append(overlap)
        if tid in domain_range.index:
            pdb_range.loc[tid, ['%_domain_covered']] = sum(domain_coverages) / len(
                domain_range.loc[tid]['ranges'])
        pdb_range.loc[tid, ['n_length']] = n_residues
        pdb_range.loc[tid, ['%_sequence_covered']] = np.where((n_residues / general_data.loc[tid].length) > 1, 1,
                                                              n_residues / general_data.loc[tid].length)
    pdb_range = pdb_range.drop(['start_stop', 'n_length'], axis=1)
    pdb_range.reset_index(drop=True, inplace=True)

    data = data.merge(pdb_range, on='Target_id', how='left')
    data = data.merge(total_pdb, on='Target_id', how='left')
    data = data.merge(pdb_w_lig, on='Target_id', how='left')
    data = data.merge(pdb_sites_tractable, on='Target_id', how='left')
    data = data.merge(pdb_sites_druggable, on='Target_id', how='left')
    data = data.merge(pdb_blast_close, on='Target_id', how='left')
    data = data.merge(pdb_blast_max_simil, on='Target_id', how='left')

    # ================================================================================================================#
    # ========================================= DOMAINS ==============================================================#
    # ================================================================================================================#

    data = data.merge(
        results['domains'].groupby('Target_id')['start'].count().reset_index().rename(
            columns={'start': 'domains_count'}),
        on='Target_id', how='left')
    data = data.merge(
        results['domain_drugE'].rename(columns={'tractable': 'domain_tractable', 'druggable': 'domain_druggable'}),
        on='Target_id', how='left')

    # ================================================================================================================#
    # =========================================== COMMERCIAL =========================================================#
    # ================================================================================================================#

    data = data.merge(results['commercials'].groupby('target_id')['smiles'].count().reset_index().rename(
        columns={'smiles': 'commercial_total', 'target_id': 'Target_id'}), on='Target_id', how='left')
    data = data.merge(results['commercials'][results['commercials']['affinity_value'] <= 100].groupby('target_id')[
        'smiles'].count().reset_index().rename(
        columns={'smiles': 'commercial_potent_total', 'target_id': 'Target_id'}), on='Target_id', how='left')

    # ================================================================================================================#
    # ========================================= TCRD INFO =========================================================#
    # ================================================================================================================#

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

    tcrd_res['tdl_info'].rename(columns={'protein_id': 'tcrd_id'}, inplace=True)
    tcrd_res['tdl_info'] = tcrd_res['tdl_info'].round(2)

    tcrd_data = tcrd_data.merge(tcrd_res['tdl_info'], on='tcrd_id', how='left')

    tcrd_data = tcrd_data.merge(tcrd_res['patent'].groupby('protein_id')['count'].sum().reset_index().rename(
        columns={'count': 'total_patent_count', 'protein_id': 'tcrd_id'}), on='tcrd_id', how='left')
    tcrd_data = tcrd_data.merge(
        tcrd_res['patent'].iloc[tcrd_res['patent'].groupby('protein_id')['count'].idxmax()].rename(
            columns={'protein_id': 'tcrd_id', 'year': 'year_max_patents', 'count': 'count_patents_max_year'}),
        on='tcrd_id',
        how='left')

    disease_clean = tcrd_res['disease'][
        (tcrd_res['disease']['score'] > 1) & (
            ~tcrd_res['disease']['doid'].isin(tcrd_res['disease']['parent'].unique()))]
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

    tcrd_data = tcrd_data.merge(tcrd_res['novelty'].rename(columns={'score': 'novelty_score'}), on='tcrd_id',
                                how='left')
    tcrd_data = tcrd_data.drop(['tcrd_id'], axis=1)

    data = data.merge(tcrd_data, on='Target_id', how='left')

    connector_targetDB.close()
    return data


def get_green_red_grad(number, v_min, v_max):
    middle = (v_min + v_max) / 2
    scale = 255 / (middle - v_min)
    if number <= v_min:
        color = "#FF0000"
    elif number >= v_max:
        color = "#00FF00"
    elif number < middle:
        color = "#FF%02X00" % int((number - v_min) * scale)
    else:
        color = "#%02XFF00" % (255 - int((number - middle) * scale))
    return color


def make_spider_plot_v3(data, data_qual, labels, druggability_val=0, target_name=''):
    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_axes([0, 0, 0.6, 0.8], polar=True)
    ax.spines['polar'].set_visible(False)
    N = len(data)
    bar_width = 2 * np.pi / N - 0.15
    theta = np.arange(0 * np.pi, 2 * np.pi, 2 * np.pi / N)
    theta = np.roll(theta, 4)
    colors = ['#5c88ed', '#b56dd6', '#6bcf67', '#f0ed37', '#f08522', '#7a7a7a', '#d41c1c']
    # ax.bar(0, 1, bottom=9, width=2 * np.pi, color='r', linewidth=0, alpha=0.3)
    # ax.bar(0, 5, bottom=4, width=2 * np.pi, color='lime', linewidth=0, alpha=0.2)
    # ax.bar(0, 3, bottom=1, width=2 * np.pi, color='gold', linewidth=0, alpha=0.2)
    # ax.bar(0, 1, bottom=0, width=2 * np.pi,alpha=0.2,color='red', linewidth=0)
    ax.bar(0, 0.3, bottom=13, width=2 * np.pi, alpha=1, color='black', linewidth=0)
    ax.bar(0, 3, width=2 * np.pi, alpha=1, color=get_green_red_grad(druggability_val, 0, 100), linewidth=0)
    value = '%2d%%' % druggability_val
    ax.text(1.1 * np.pi, 0.7 * np.pi, value, weight='bold', size='large')
    for i in reversed(range(len(data))):
        if list(labels)[i] in ['Safety', 'Structural Biology', 'Chemistry', 'Genetic Association']:
            color_to_use = get_green_red_grad(data_qual[list(labels)[i]], 0, 10)
        else:
            color_to_use = '#d1d1d1'
        # ax.bar(theta[i], data[i],bottom=3, width=bar_width, align='center', color=color_to_use)
        ax.bar(theta[i], data[i], bottom=3, width=bar_width, align='center', color=color_to_use, label=list(labels)[i],
               edgecolor=colors[i], linewidth=4)

    plt.title(target_name, weight='bold', fontsize=14)
    ax.set_xticks(theta)
    x_labels = [''] * len(theta)
    ax.set_xticklabels(x_labels)
    ax.yaxis.grid(True)
    ax.xaxis.grid(False)
    fig.legend(loc=7, fontsize=10, fancybox=True, markerscale=1.2)
    ax.set_yticks([])

    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)
    return buf


def make_spider_plot(data, labels, target_name=''):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0, 0, 0.6, 1], polar=True)
    ax.spines['polar'].set_visible(False)
    N = len(data)
    bar_width = 2 * np.pi / N
    theta = np.arange(0, 2 * np.pi, 2 * np.pi / N)
    colors = ['#95d0fc', '#0485d1', '#b790d4', '#87ae73', '#fec615', '#fb7d07', '#95a3a6', '#ccff33']
    ax.bar(0, 1, bottom=9, width=2 * np.pi, color='r', linewidth=0, alpha=0.3)
    ax.bar(0, 5, bottom=4, width=2 * np.pi, color='lime', linewidth=0, alpha=0.2)
    ax.bar(0, 3, bottom=1, width=2 * np.pi, color='gold', linewidth=0, alpha=0.2)
    ax.bar(0, 1, width=2 * np.pi, color='r', linewidth=0)
    for i in range(len(data)):
        ax.bar(theta[i], data[i], width=bar_width, align='center', label=list(labels)[i], color=colors[i],
               edgecolor='black', linewidth=1.5)
    plt.title(target_name, weight='bold', fontsize=20)
    ax.set_xticks(theta)
    x_labels = [''] * len(theta)
    ax.set_xticklabels(x_labels)
    ax.yaxis.grid(True)
    ax.xaxis.grid(False)
    fig.legend(loc=7, fontsize=14, fancybox=True, markerscale=1.2)
    ax.set_yticks([])

    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)
    return buf


if __name__ == '__main__':
    pass
