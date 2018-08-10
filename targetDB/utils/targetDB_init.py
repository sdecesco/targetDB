#!/usr/bin/env python

import sqlite3
import time
from tkinter.filedialog import askdirectory
from pathlib import Path

import pandas as pd
import pkg_resources

creation_sql = """DROP TABLE IF EXISTS BindingDB;
CREATE TABLE IF NOT EXISTS BindingDB
(
ligand_name   text,
ligand_smiles text,
Inchi         text,
Inchi_key     varchar(100) default NULL,
target_name   text,
target_id     varchar(25)  default NULL,
"IC50(nM)"    float        default NULL,
"Ki(nM)"      float        default NULL,
"Kd(nM)"      float        default NULL,
"EC50(nM)"    float        default NULL,
"kon(M-1s-1)" float        default NULL,
"koff(s-1)"   float        default NULL,
pH            float        default NULL,
Temp          varchar(20)  default NULL,
DOI           varchar(50)  default NULL,
Patent_number varchar(50)  default NULL,
Institution   text,
ZincID        varchar(50)  default NULL,
Source        text
);

create index idx_BindingDB_BindingDB_Inchi_key_index
on BindingDB (Inchi_key);

create index idx_BindingDB_BindingDB_ZincID_index
on BindingDB (ZincID);

create index idx_BindingDB_BindingDB_target_id_index
on BindingDB (target_id);

DROP TABLE IF EXISTS PDB;
CREATE TABLE IF NOT EXISTS PDB
(
PDB_code   varchar(5) not null,
Technique  varchar(20) default NULL,
Resolution varchar(20) default NULL,
Date       timestamp   default CURRENT_TIMESTAMP not null,
primary key (PDB_code)
);

create index idx_PDB_PDB_PDB_code_index
on PDB (PDB_code);

DROP TABLE IF EXISTS Targets;
CREATE TABLE IF NOT EXISTS Targets
(
Target_id           varchar(10) not null,
Gene_name           varchar(20) not null,
Species             varchar(40)  default NULL,
species_id          varchar(25)  default NULL,
Sequence            text,
Cell_location       text,
Process             text,
Function            text,
Number_isoforms     integer      default NULL,
Synonyms            text,
Date_modified       timestamp    default CURRENT_TIMESTAMP not null,
Protein_class       text,
Protein_class_desc  text,
Protein_class_short varchar(100) default NULL,
chembl_id           varchar(25)  default NULL,
primary key (Target_id)
);

DROP TABLE IF EXISTS "3D_Blast";
CREATE TABLE IF NOT EXISTS "3D_Blast"
(
Query_target_id  varchar(10) not null,
Hit_gene_id      varchar(10) not null,
Hit_PDB_code     varchar(5)  not null,
Hit_Chain_id     varchar(20) not null,
similarity       float       not null,
last_updated     timestamp   default CURRENT_TIMESTAMP not null,
Unique_ID        varchar(50) default '' not null,
Hit_gene_name    varchar(15) default NULL,
Hit_gene_species varchar(50) default NULL,
Chain_Letter     varchar(3)  default NULL,
primary key (Unique_ID),
constraint "3D_Blast_Targets_Target_id_fk"
foreign key (Query_target_id) references Targets
on delete cascade
);

create index idx_3D_Blast_3D_Blast_Hit_PDB_code_index
on "3D_Blast" (Hit_PDB_code);

create index idx_3D_Blast_3D_Blast_Targets_Target_id_fk
on "3D_Blast" (Query_target_id);

DROP TABLE IF EXISTS Crossref;
CREATE TABLE IF NOT EXISTS Crossref
(
target_id  varchar(15) default NULL,
Chembl_id  varchar(50) default NULL,
date       timestamp   default CURRENT_TIMESTAMP not null,
hgnc_id    varchar(50) default NULL,
ensembl_id varchar(50) default NULL,
constraint Crossrefunip_Targets_Target_id_fk
foreign key (target_id) references Targets
on update cascade
on delete cascade
);

create index idx_Crossref_Chembl_id_index
on Crossref (Chembl_id);

create index idx_Crossref_Targets_Ensembl_id
on Crossref (ensembl_id);

create index idx_Crossref_Targets_Target_id
on Crossref (target_id);

create index idx_Crossref_Targets_hgnc_id
on Crossref (hgnc_id);

DROP TABLE IF EXISTS Domain_targets;
CREATE TABLE IF NOT EXISTS Domain_targets
(
Target_id    varchar(10) not null,
Domain_start integer     default NULL,
Domain_stop  integer     default NULL,
Domain_name  text        not null,
domain_id    varchar(50) default NULL,
source_name  varchar(50) default NULL,
Source_id    varchar(20) default NULL,
length       integer     default NULL,
unique (domain_id),
constraint Domain_targets_Targets_Target_id_fk
foreign key (Target_id) references Targets
on delete cascade
);

create index idx_Domain_targets_Domain_targets_Targets_Target_id_fk
on Domain_targets (Target_id);

DROP TABLE IF EXISTS Isoforms;
CREATE TABLE IF NOT EXISTS Isoforms
(
Target_id    varchar(10) not null,
Isoform_name varchar(10) not null,
Isoform_id   varchar(20) not null,
Sequence     text,
n_residues   integer      default NULL,
Canonical    integer      default NULL,
Identity     float(10, 2) default NULL,
Gaps         float(10, 2) default NULL,
Date         timestamp    default CURRENT_TIMESTAMP not null,
Score        float        default NULL,
primary key (Isoform_id),
constraint Isoforms_Targets_Target_id_fk
foreign key (Target_id) references Targets
on delete cascade
);

create index idx_Isoforms_Isoforms_Targets_Target_id_fk
on Isoforms (Target_id);

DROP TABLE IF EXISTS PDB_Chains;
CREATE TABLE IF NOT EXISTS PDB_Chains
(
Chain_id     varchar(10) not null,
PDB_code     varchar(5)  not null,
Chain        varchar(2)  not null,
Sequence     text,
Target_id    varchar(10) not null,
Start        integer      default NULL,
Stop         integer      default NULL,
n_residues   integer      default NULL,
DATE         timestamp    default CURRENT_TIMESTAMP,
start_stop   varchar(100) default NULL,
equal_chains VARCHAR(50),
primary key (Chain_id),
constraint PDB_Chains_PDB_PDB_code_fk
foreign key (PDB_code) references PDB,
constraint PDB_Chains_Targets_Target_id_fk
foreign key (Target_id) references Targets
on delete cascade
);

DROP TABLE IF EXISTS PDBChain_Domain;
CREATE TABLE IF NOT EXISTS PDBChain_Domain
(
Chain_id  varchar(20) default NULL,
Domain_id varchar(50) default NULL,
constraint PDBChain_Domain_PDB_Chains_Chain_id_fk
foreign key (Chain_id) references PDB_Chains
on delete cascade
);

create index idx_PDBChain_Domain_PDBChain_Domain_Domain_targets_Domain_id_fk
on PDBChain_Domain (Domain_id);

create index idx_PDBChain_Domain_PDBChain_Domain_PDB_Chains_Chain_id_fk
on PDBChain_Domain (Chain_id);

create index idx_PDB_Chains_PDB_Chains_Chain_id_index
on PDB_Chains (Chain_id);

create index idx_PDB_Chains_PDB_Chains_PDB_PDB_code_fk
on PDB_Chains (PDB_code);

create index idx_PDB_Chains_PDB_Chains_Targets_Target_id_fk
on PDB_Chains (Target_id);

create index idx_Targets_Targets_Target_id_index
on Targets (Target_id);

DROP TABLE IF EXISTS assays;
CREATE TABLE IF NOT EXISTS assays
(
assay_id          varchar(30) not null,
assay_description text,
relationship_desc varchar(40)  default NULL,
date              timestamp    default CURRENT_TIMESTAMP not null,
doi               varchar(100) default NULL,
species           varchar(100) default NULL,
bioactivity_type  varchar(50)  default NULL,
chembl_version    varchar(50)  default NULL,
curated_by        varchar(50)  default NULL,
confidence_score  integer      default NULL,
ref_type          VARCHAR(20),
patent_id         VARCHAR(30),
confidence_txt    TEXT,
cell_type         VARCHAR(50),
cell_organism     VARCHAR(50),
variant_uniprot   VARCHAR(50),
mutant            VARCHAR(50),
target_chembl_id  VARCHAR(25),
primary key (assay_id)
);

create index assays_target_id_index
on assays (target_chembl_id);

create index idx_assays_assays_assay_id_index
on assays (assay_id);

create index idx_assays_assays_bioactivity_type_index
on assays (bioactivity_type);

DROP TABLE IF EXISTS bioactivities;
CREATE TABLE IF NOT EXISTS bioactivities
(
lig_id                VARCHAR(50) not null,
Target_id             VARCHAR(40) not null,
assay_id              VARCHAR(30) not null,
units                 VARCHAR(30)  default NULL,
operator              VARCHAR(10)  default NULL,
value_num             DOUBLE       default NULL,
doi                   TEXT,
date                  TIMESTAMP    default CURRENT_TIMESTAMP not null,
value_text            VARCHAR(50)  default NULL,
chembl_version        VARCHAR(20)  default NULL,
standard_type         VARCHAR(100) default NULL,
activity_comment      TEXT,
data_validity_comment TEXT,
pchembl_value         FLOAT        default NULL,
target_type           VARCHAR(100) default NULL,
target_name           TEXT,
target_organism       TEXT
);

create index idx_bioactivities_bioactivities_Target_id_index
on bioactivities (Target_id);

create index idx_bioactivities_bioactivities_assays_assay_id_fk
on bioactivities (assay_id);

create index idx_bioactivities_bioactivities_lig_id_index
on bioactivities (lig_id);

create index idx_bioactivities_bioactivities_standard_type_index
on bioactivities (standard_type);

DROP TABLE IF EXISTS diff_exp_disease;
CREATE TABLE IF NOT EXISTS diff_exp_disease
(
Target_id         varchar(50) default NULL,
disease           text,
t_stat            float       default NULL,
expression_status varchar(15) default NULL,
p_value           float       default NULL,
date              timestamp   default CURRENT_TIMESTAMP not null,
constraint diff_exp_disease_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_diff_exp_disease_diff_exp_disease_Targets_Target_id_fk
on diff_exp_disease (Target_id);

DROP TABLE IF EXISTS diff_exp_tissue;
CREATE TABLE IF NOT EXISTS diff_exp_tissue
(
Target_id         varchar(50) default NULL,
t_stat            float       default NULL,
Tissue            text,
expression_status varchar(20) default NULL,
p_value           float       default NULL,
date              timestamp   default CURRENT_TIMESTAMP,
constraint diff_exp_tissue_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_diff_exp_tissue_diff_exp_tissue_Targets_Target_id_fk
on diff_exp_tissue (Target_id);

DROP TABLE IF EXISTS disease;
CREATE TABLE IF NOT EXISTS disease
(
Target_id    varchar(50) default NULL,
disease_name text,
disease_id   varchar(50) default NULL,
Unique_id    varchar(100) not null,
date         timestamp   default CURRENT_TIMESTAMP,
primary key (Unique_id),
constraint disease_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_disease_disease_Targets_Target_id_fk
on disease (Target_id);

DROP TABLE IF EXISTS drugEbility_domains;
CREATE TABLE IF NOT EXISTS drugEbility_domains
(
PX_number          integer     default NULL,
pdb_code           varchar(10) default NULL,
domain_fold        text,
domain_superfamily text,
domain_family      text,
ensemble           float       default NULL,
tractable          integer     default NULL,
druggable          integer     default NULL
);

create index idx_drugEbility_domains_drugEbility_domains_PX_number_index
on drugEbility_domains (PX_number);

create index idx_drugEbility_domains_drugEbility_domains_pdb_code_index
on drugEbility_domains (pdb_code);

DROP TABLE IF EXISTS drugEbility_sites;
CREATE TABLE IF NOT EXISTS drugEbility_sites
(
site_id     integer     default NULL,
px_number   integer     default NULL,
pdb_code    varchar(10) default NULL,
site_number integer     default NULL,
ensemble    float       default NULL,
tractable   integer     default NULL,
druggable   integer     default NULL
);

create index idx_drugEbility_sites_drugEbility_sites_pdb_code_index
on drugEbility_sites (pdb_code);

create index idx_drugEbility_sites_drugEbility_sites_px_number_index
on drugEbility_sites (px_number);

DROP TABLE IF EXISTS ensembl;
CREATE TABLE IF NOT EXISTS ensembl
(
ensembl_gene varchar(50) default NULL,
uniprot_id   varchar(20) default NULL,
gene_name    varchar(50) default NULL,
gene_symbol  varchar(50) default NULL
);

create index idx_ensembl_ensembl_ensembl_gene_index
on ensembl (ensembl_gene);

create index idx_ensembl_ensembl_gene_symbol_index
on ensembl (gene_symbol);

create index idx_ensembl_ensembl_uniprot_id_index
on ensembl (uniprot_id);

DROP TABLE IF EXISTS fPockets;
CREATE TABLE IF NOT EXISTS fPockets
(
PDB_code      varchar(5)  not null,
Target_id     varchar(10) not null,
Pocket_number varchar(5)  not null,
Pocket_id     varchar(50) not null,
Score         float       default NULL,
DrugScore     float       default NULL,
apolar_sasa   float       default NULL,
polar_sasa    float       default NULL,
total_sasa    float       default NULL,
volume        float       default NULL,
last_updated  timestamp   default CURRENT_TIMESTAMP,
blast         varchar(10) default 'FALSE',
druggable     varchar(10) default NULL,
primary key (Pocket_id),
constraint fPockets_Targets_Target_id_fk
foreign key (Target_id) references Targets
on delete cascade
);

create index idx_fPockets_fPockets_PDB_code_index
on fPockets (PDB_code);

create index idx_fPockets_fPockets_Pocket_id_index
on fPockets (Pocket_id);

create index idx_fPockets_fPockets_Targets_Target_id_fk
on fPockets (Target_id);

DROP TABLE IF EXISTS fPockets_Chain;
CREATE TABLE IF NOT EXISTS fPockets_Chain
(
Pocket_id        varchar(50) default NULL,
Chain_id         varchar(10) default NULL,
List_of_contacts text,
constraint fPockets_Chain_fPockets_Pocket_id_fk
foreign key (Pocket_id) references fPockets
on delete cascade
);

create index idx_fPockets_Chain_fPockets_Chain_Chain_id_index
on fPockets_Chain (Chain_id);

create index idx_fPockets_Chain_fPockets_Chain_fPockets_Pocket_id_fk
on fPockets_Chain (Pocket_id);

DROP TABLE IF EXISTS fPockets_Domain;
CREATE TABLE IF NOT EXISTS fPockets_Domain
(
Pocket_id varchar(50) default NULL,
Domain_id varchar(75) default NULL,
Coverage  double      default NULL,
constraint fPockets_Domain_fPockets_Pocket_id_fk
foreign key (Pocket_id) references fPockets
on delete cascade
);

create index idx_fPockets_Domain_fPockets_Domain_Domain_targets_Domain_id_fk
on fPockets_Domain (Domain_id);

create index idx_fPockets_Domain_fPockets_Domain_fPockets_Pocket_id_fk
on fPockets_Domain (Pocket_id);

DROP TABLE IF EXISTS gwas;
CREATE TABLE IF NOT EXISTS gwas
(
Target_id        varchar(50) default NULL,
doi              text,
first_author     text,
organism         varchar(50) default NULL,
p_value          double      default NULL,
phenotype        text,
publication_year varchar(15) default NULL,
pubmed_id        varchar(50) default NULL,
date             timestamp   default CURRENT_TIMESTAMP not null,
id               integer not null,
primary key (id autoincrement),
constraint gwas_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_gwas_gwas_Targets_Target_id_fk
on gwas (Target_id);

DROP TABLE IF EXISTS hgnc;
CREATE TABLE IF NOT EXISTS hgnc
(
hgnc_id    varchar(25)  default NULL,
xref_name  varchar(50)  default NULL,
xref_value varchar(200) default NULL
);

create index idx_hgnc_hgnc_hgnc_id_index
on hgnc (hgnc_id);

create index idx_hgnc_hgnc_hgnc_id_xref_name_xref_value_index
on hgnc (hgnc_id, xref_name, xref_value);

create index idx_hgnc_hgnc_xref_name_index
on hgnc (xref_name);

create index idx_hgnc_hgnc_xref_value_index
on hgnc (xref_value);

DROP TABLE IF EXISTS ligands;
CREATE TABLE IF NOT EXISTS ligands
(
mol_name           varchar(100) default NULL,
lig_id             varchar(50) not null,
max_phase          integer      default NULL,
oral               integer      default NULL,
black_box_warning  integer      default NULL,
indication_class   varchar(100) default NULL,
class_def          text,
alogp              float        default NULL,
acd_logd           float        default NULL,
acd_logp           float        default NULL,
acd_most_bpka      float        default NULL,
acd_most_apka      float        default NULL,
HBD                integer      default NULL,
HBA                integer      default NULL,
TPSA               float        default NULL,
n_heavy_atoms      integer      default NULL,
molecularWeight    float        default NULL,
rotatableBonds     integer      default NULL,
n_Ar_rings         integer      default NULL,
CNS_MPO            float        default NULL,
molecular_species  varchar(20)  default NULL,
num_ro5_violations integer      default NULL,
ro3_pass           varchar(5)   default NULL,
mol_formula        varchar(50)  default NULL,
canonical_smiles   text,
std_inchi          text,
std_inchi_key      varchar(100) default NULL,
chembl_version     integer      default NULL,
date_updated       timestamp    default CURRENT_TIMESTAMP not null,
primary key (lig_id)
);

create index idx_ligands_ligands_lig_id_index
on ligands (lig_id);

create index idx_ligands_ligands_std_inchi_key_index
on ligands (std_inchi_key);

DROP TABLE IF EXISTS modifications;
CREATE TABLE IF NOT EXISTS modifications
(
mod_id       varchar(20) not null,
mod_type     varchar(10) default NULL,
action       varchar(20) default NULL,
comment      text,
previous     text,
new          text,
start        integer     default NULL,
stop         integer     default NULL,
date_mod     timestamp   default CURRENT_TIMESTAMP not null,
Unique_modID varchar(50) default '' not null,
Target_id    varchar(20) default NULL,
domains      text,
constraint modifications_Targets_Target_id_fk
foreign key (Target_id) references Targets
on delete cascade
);

DROP TABLE IF EXISTS isoform_modifications;
CREATE TABLE IF NOT EXISTS isoform_modifications
(
isoform_id varchar(30) default NULL,
mod_id     varchar(50) default NULL,
Date       timestamp   default CURRENT_TIMESTAMP,
constraint isoform_modifications_Isoforms_Isoform_id_fk
foreign key (isoform_id) references Isoforms
on delete cascade
);

create index idx_isoform_modifications_isoform_modifications_Isoforms_Isoform_id_fk
on isoform_modifications (isoform_id);

create index idx_isoform_modifications_isoform_modifications_modifications_Unique_modID_fk
on isoform_modifications (mod_id);

create index idx_modifications_modifications_Targets_Target_id_fk
on modifications (Target_id);

create index idx_modifications_modifications_Unique_modID_index
on modifications (Unique_modID);

create index idx_modifications_modifications_mod_id_index
on modifications (mod_id);

DROP TABLE IF EXISTS opentarget_association;
CREATE TABLE IF NOT EXISTS opentarget_association
(
target_id           VARCHAR(25),
disease_area        TEXT,
disease_name        VARCHAR(100),
overall_score       FLOAT,
genetic_association FLOAT,
known_drug          FLOAT,
litterature_mining  FLOAT,
animal_model        FLOAT,
affected_pathway    FLOAT,
rna_expression      FLOAT,
somatic_mutation    FLOAT,
CONSTRAINT opentarget_association_Targets_Target_id_fk 
FOREIGN KEY (target_id) REFERENCES Targets (Target_id) ON DELETE CASCADE
);

create index OT_Association_idx
on opentarget_association (target_id);

DROP TABLE IF EXISTS pathways;
CREATE TABLE IF NOT EXISTS pathways
(
Target_id       varchar(50)  default NULL,
pathway_dataset varchar(100) default NULL,
pathway_name    text,
date            timestamp    default CURRENT_TIMESTAMP not null,
constraint pathways_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_pathways_pathways_Targets_Target_id_fk
on pathways (Target_id);

DROP TABLE IF EXISTS pdb_bind;
CREATE TABLE IF NOT EXISTS pdb_bind
(
pdb_code         varchar(5) not null,
pub_year         integer     default NULL,
binding_type     varchar(20) default NULL,
binding_operator varchar(5)  default NULL,
binding_value    float       default NULL,
binding_units    varchar(10) default NULL,
lig_name         longtext,
type             varchar(50) default NULL,
version          varchar(30) default NULL,
primary key (pdb_code)
);

DROP TABLE IF EXISTS phenotype;
CREATE TABLE IF NOT EXISTS phenotype
(
Target_id      varchar(50) default NULL,
Allele_id      varchar(50) default NULL,
Phenotype      text,
Phenotype_desc text,
genotype       text,
organism       varchar(50) default NULL,
zygosity       varchar(50) default NULL,
date           timestamp   default CURRENT_TIMESTAMP not null,
Allele_type    text,
Allele_symbol  text,
constraint phenotype_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_phenotype_phenotype_Targets_Target_id_fk
on phenotype (Target_id);

DROP TABLE IF EXISTS protein_blast;
CREATE TABLE IF NOT EXISTS protein_blast
(
Query_target_id    varchar(10) not null,
Hit_gene_id        varchar(10) not null,
sequence_in_common float       not null,
similarity         float       not null,
last_updated       timestamp   default CURRENT_TIMESTAMP not null,
Hit_gene_name      varchar(50) default NULL,
Hit_gene_species   varchar(50) default NULL,
Unique_ID          varchar(50) default '' not null,
primary key (Unique_ID),
constraint protein_blast_Targets_Target_id_fk
foreign key (Query_target_id) references Targets
on delete cascade
);

create index idx_protein_blast_protein_blast_Hit_gene_id_index
on protein_blast (Hit_gene_id);

create index idx_protein_blast_protein_blast_Targets_Target_id_fk
on protein_blast (Query_target_id);

create index idx_protein_blast_protein_blast_Unique_ID_index
on protein_blast (Unique_ID);

DROP TABLE IF EXISTS protein_expression_levels;
CREATE TABLE IF NOT EXISTS protein_expression_levels
(
Target_id varchar(10) default NULL,
organ     varchar(30) default NULL,
tissue    varchar(20) default NULL,
cell      varchar(50) default NULL,
value     integer     default NULL,
date      timestamp   default CURRENT_TIMESTAMP,
Entry_id  varchar(100) not null,
primary key (Entry_id),
constraint protein_expression_levels_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_protein_expression_levels_protein_expression_levels_Targets_Target_id_fk
on protein_expression_levels (Target_id);

DROP TABLE IF EXISTS protein_expression_selectivity;
CREATE TABLE IF NOT EXISTS protein_expression_selectivity
(
Target_id           varchar(10) default '' not null,
Selectivity_entropy float       default NULL,
max_organ           varchar(30) default NULL,
date                timestamp   default CURRENT_TIMESTAMP,
primary key (Target_id),
constraint protein_expression_selectivity_Targets_Target_id_fk
foreign key (Target_id) references Targets
on update cascade
on delete cascade
);

create index idx_protein_expression_selectivity_protein_expression_selectivity_Target_id_index
on protein_expression_selectivity (Target_id);

DROP TABLE IF EXISTS purchasable_compounds;
CREATE TABLE IF NOT EXISTS purchasable_compounds
(
target_id      varchar(20) default NULL,
smiles         text,
affinity_type  varchar(50) default NULL,
affinity_value float       default NULL,
affinity_unit  varchar(10) default 'nM',
website        text,
price          text,
op             VARCHAR(10)
);

create index idx_purchasable_compounds_purchasable_compounds_target_id_index
on purchasable_compounds (target_id);

DROP TABLE IF EXISTS tcrd_id;
CREATE TABLE IF NOT EXISTS tcrd_id
(
    tcrd_id int,
    Target_id VARCHAR(25)
);
CREATE INDEX tcrd_id_tcrd_id_index ON tcrd_id (tcrd_id);
CREATE INDEX tcrd_id_Target_id_index ON tcrd_id (Target_id);


DROP TABLE IF EXISTS tcrd_target;
CREATE TABLE IF NOT EXISTS tcrd_target
(
    tcrd_id int,
    Pharos_class VARCHAR(10),
    protein_family text,
    protein_family_detail text
);
CREATE INDEX tcrd_target_tcrd_id_index ON tcrd_target (tcrd_id);

DROP TABLE IF EXISTS tcrd_info;
CREATE TABLE IF NOT EXISTS tcrd_info
(
    protein_id int,
    "Ab Count" int,
    "EBI Total Patent Count" int,
    "JensenLab PubMed Score" float,
    "MAb Count" int,
    "NCBI Gene PubMed Count" int,
    "PubTator Score" int
);
CREATE INDEX tcrd_info_protein_id_index ON tcrd_info (protein_id);

DROP TABLE IF EXISTS tcrd_patent;
CREATE TABLE IF NOT EXISTS tcrd_patent
(
    protein_id int,
    year int,
    count int
);
CREATE INDEX tcrd_patent_protein_id_index ON tcrd_patent (protein_id);

DROP TABLE IF EXISTS tcrd_disease;
CREATE TABLE IF NOT EXISTS tcrd_disease
(
    protein_id int,
    disease_id int,
    doid VARCHAR(20),
    score float,
    name text,
    parent VARCHAR(20)
);
CREATE INDEX tcrd_disease_protein_id_index ON tcrd_disease (protein_id);
CREATE INDEX tcrd_disease_disease_id_index ON tcrd_disease (disease_id);

DROP TABLE IF EXISTS tcrd_novelty;
CREATE TABLE IF NOT EXISTS tcrd_novelty
(
    score float,
    protein_id int
);
CREATE INDEX tcrd_novelty_protein_id_index ON tcrd_novelty (protein_id);

"""


def create_db():
	path_exist = False
	targetDB_path = Path()
	while not path_exist:
		targetDB_path = Path(askdirectory(title='Select where you want to save the targetDB file'))
		if targetDB_path.is_dir():
			path_exist = True
		else:
			print('[ERROR]: The directory you have entered does not exists')

	db_name = 'TargetDB_'+time.strftime("%d_%m_%y")+'.db'

	db_file_path = targetDB_path.joinpath(db_name)

	print('[DATABASE]: Creating targetDB database table structure')
	connector = sqlite3.connect(str(db_file_path))
	cursor = connector.cursor()
	cursor.executescript(creation_sql)
	connector.commit()
	print('[DATABASE]: targetDB empty database created')
	print('[FILE CREATED]: ',db_file_path)

	pck_path = Path(str(pkg_resources.resource_filename('targetDB.utils', ''))).parent
	data_pck_path = pck_path.joinpath('data')

	fill_db(data_pck_path,connector)
	print('[DATABASE]: targetDB table pre-filling completed')

	connector.close()
	return db_file_path


def fill_db(path_to_datafiles,connector):
	for f in path_to_datafiles.iterdir():
		df = pd.read_pickle(str(f))
		table_name = '_'.join(f.name.split('_')[:-3])
		print('[DATABASE]: Filling the table ',table_name)
		df.to_sql(table_name, con=connector, if_exists='append', index=False)
