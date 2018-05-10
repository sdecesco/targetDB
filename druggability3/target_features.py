#!/usr/bin/env python
import re
import pandas as pd
from druggability3 import cns_mpo as mpo
from druggability3 import db_connection as db
import numpy as np
import scipy.stats as sc


def get_single_features(target_id, user=None, pwd=None):
	single_queries = {'general_info': "SELECT * FROM Targets WHERE Target_id='" + target_id + "'",
	                  'disease': "SELECT Target_id,disease_name,disease_id FROM disease WHERE Target_id='" + target_id + "'",
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
              CONCAT(T.Gene_name,'-',I.Isoform_name) as isoform_name,
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
            FROM 3D_Blast
              LEFT JOIN drugEbility_sites DS
              ON DS.pdb_code=Hit_PDB_code
            WHERE Query_target_id='%s'
            GROUP BY Hit_PDB_code
            ORDER BY similarity DESC""" % target_id,
	                  'pdb': """SELECT
              C.PDB_code,
              P.Technique,
              P.Resolution,
              GROUP_CONCAT(DISTINCT C.Chain SEPARATOR ',') AS Chain,
              C.n_residues,
              C.start_stop,
              GROUP_CONCAT(DISTINCT D.Domain_name SEPARATOR ',') AS Domain_name,
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
                ON DS.pdb_code = C.PDB_code
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
              GROUP_CONCAT(CONCAT(D.Domain_name,' (',Domain.Coverage,'%)') SEPARATOR ',') as domains
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
              L.n_alerts as n_alerts,
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
              L.n_alerts as n_alerts,
              L.molecular_species,
              L.num_ro5_violations as ro5_violations,
              L.ro3_pass as pass_ro3,
              B.ligand_smiles as SMILES
            FROM BindingDB B
              LEFT JOIN ligands L
              ON B.inchi_key = L.std_inchi_key
            WHERE target_id = '%s'""" % target_id,
	                  'domain_drugE': """SELECT
      GROUP_CONCAT(DISTINCT UPPER(pdb_code) SEPARATOR ',') pdb_list,
      domain_fold,
      domain_superfamily,
    max(tractable) tractable,
      max(druggable) druggable
    FROM drugEbility_domains
    WHERE pdb_code in (SELECT DISTINCT PDB_code
      FROM PDB_Chains
    WHERE target_id = '%s')
    GROUP BY domain_fold""" % target_id}

	dbase = db.open_db('druggability', user=user, pwd=pwd)
	results = {qname: pd.read_sql(query, con=dbase.db) for qname, query in single_queries.items()}

	if not results['bioactives'].empty:
		results.update(transform_bioactivities(results['bioactives'],dbase))

	return results
	dbase.close()


def get_list_features(gene_ids, user=None, pwd=None):
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
	dbase = db.open_db('druggability', user=user, pwd=pwd)

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


def transform_bioactivities(results, dbase):
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
	       'apKa', 'bpKa', 'nAr', 'n_alerts', 'pass_ro3', 'ro5_violations', 'rotB', 'CNS_MPO', 'mol_name',
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
		res_lig = dbase.get(query)

		entropies = []
		binding_data = pd.DataFrame.from_records(res_lig)
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
			             'TPSA', 'aLogP', 'apKa', 'bpKa', 'nAr', 'n_alerts', 'pass_ro3',
			             'ro5_violations', 'rotB', 'CNS_MPO', 'mol_name', 'molecular_species',
			             'indication_class', 'class_def', 'max_phase', 'oral', 'assay_ref',
			             'ref_bio']
			binding = binding[col_order]
	return {'binding': binding, 'dose_response': dose_response, 'other': other, 'ADME': ADME, 'emax': emax,
            'efficacy_bio': efficacy_bio, 'percent_inhibition': percent_bio}
