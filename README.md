TargetDB
=========

TargetDB is a tool to quickly querry multiple publicly available databases and provide an integrated view of the information available about potential targets. A quick binding pocket search is also performed (using `fpocket`).

Installation
------------
### Python package installation
#### Pip installation

```
pip install targetDB
```

##### Python compatibility

python version **>= 3.4**

Preferred python distribution : [Anaconda 3](https://www.anaconda.com/download/)


#### Database files installation

This package relies on the use of sqlite database to properly function. Depending the mode you use you will be required to point at these files 
the first time you run the script.

This script can be used in two modes :
+ Report Only
+ Database Creation
>See later sections for more details

##### Report only

The only required databases is :
+ targetDB

You can download a copy of the database [HERE](https://github.com/sdecesco/targetDB/releases/download/v1.1.6/TargetDB_20_02_19.db.zip)

>This database contains all of the human genome genes that have a uniprot ID

##### Database creation

The list of required databases is :
+ ChEMBL v24

ChEMBL sqlite database can be directly downloaded [HERE](https://www.ebi.ac.uk/chembl/downloads)

>This mode will generate a targetDB database that can then be used in report mode

#### Other dependencies for the database creation mode

##### blast
This mode use blast locally to perform similarity search and sequence alignments 

information to download and install blast locally can be found [HERE](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

##### fpocket
In order to perform binding pocket searches and assess their druggability the program fpocket is used 

>Vincent Le Guilloux, Peter Schmidtke and Pierre Tuffery, "Fpocket: An open source platform for ligand pocket detection", BMC Bioinformatics, 2009, 10:168

>Peter Schmidtke, Xavier Barril "Understanding and Predicting Druggability. A High-Throughput Method for Detection of Drug Binding Sites", J. Med. Chem., 2010, 53 (15), pp 5858–5867

instructions to download and install fpocket can be found [HERE](https://github.com/Discngine/fpocket)

Current version of targetDB works with fpocket3

**targetDB will not be able to perform pocket search on windows as fpocket is not available on windows platform**

#### System compatibility

+ Linux
+ Windows*


*Binding pocket search is not possible on Windows (see comment above)

Usage
-----
targetDB package provides two command-line tools `target_DB` and `target_REPORT`

### `target_DB`

`target_DB` is the command used to generate/fill a targetDB database that can then be used in the report mode. If you don't want to use this mode you can download a pre-filled database and use it to generate reports
>Pre-filled databases will be generated when major updates of datasources are available

Here are the flag available:
```
One of the three following option is required:
        -g : enter a single gene name
        -i : Name of the input file with a list of genes (.txt - 1 gene per line)
        -l : Enter a list of gene name separated by a ","
        
        -update : Update record if already in database (default: No)
        -blastcore : Enter the value of processor core to use for the blast search (default=8)
        
        -v : Print detailed information in the terminal
        
        -update_config : use this if you want to update the config file
        
        -create_db : Use this option to generate an empty targetDB database (file name = TargetDB_Day_Month_Year.db)
```

When using `target_DB` for the first time you will be asked to enter information about:
+ ChEMBL sqlite database file
+ targetDB sqlite database file (not asked if you use the `-create_db` flag)
+ path to blast 
+ path to fpocket
+ folder to store database files (PDB files/Blast XML/Pockets) `/!\ This folder can occupy a lot of space /!\ `
+ email address (used for pubmed searches if none provided no pubmed search will be run)

Those informations will be stored in ~/.druggability/config.ini 

If you want to modify the informations in that file just use the `-update_config` flag

##### example 
input:
```
target_DB -g DYRK1A -v
```

output:
```
[BEGINNING OF GENE]: DYRK1A (Q13627)
[ISOFORMS SEQUENCE ALIGNMENT IN PROGRESS]
[ISOFORMS SEQUENCE ALIGNMENT DONE]
[FILE ALREADY THERE]: 2VX3.pdb
[...]
[FILE ALREADY THERE]: 5AIK.pdb
[PDB DOWNLOAD DONE]
[PDB INFO]: Extracting PDB information
[PDB INFO]: Done
[POCKET SEARCH]:  (2VX3.pdb)
[POCKET SEARCH DONE]: 2VX3
[...]
[POCKET SEARCH]:  (5AIK.pdb)
[POCKET SEARCH DONE]: 5AIK
[LIGANDS]: Extracting ligand informations
[3D BLAST]:DYRK1A(Q13627)
[3D BLAST FILE FOUND]:DYRK1A
[3D BLAST DONE]: Now parsing the data - DYRK1A(Q13627)
[PDB DOWNLOAD DONE]
[PROTEIN BLAST] DYRK1A(Q13627)
[PROTEIN BLAST FILE FOUND]:DYRK1A
[PROTEIN BLAST DONE]: Now parsing the data - DYRK1A(Q13627)
[DATABASE]: Start to write info into the database
[DATABASE]: Target table populated/updated
[DATABASE]: Domain table populated/updated
[DATABASE]: PDB Blast table populated/updated
[DATABASE]: Blast table populated/updated
[DATABASE]: Isoforms tables populated/updated
[DATABASE]: PDB tables populated/updated
[DATABASE]: fPockets tables populated/updated
[DATABASE]: alternate fPockets tables populated/updated
[DATABASE]: Bioactivities table populated
[DATABASE]: Protein expression levels tables populated/updated
[DATABASE]: Disease table populated/updated
[DATABASE]: Differential expression (tissues) table populated/updated
[DATABASE]: Differential expression (disease) table populated/updated
[DATABASE]: GWAS table populated/updated
[DATABASE]: phenotype table populated/updated
[DATABASE]: Pathway table populated/updated
[DATABASE]: Cross-references table populated/updated
[DATABASE]: Open-targets table populated/updated
[END OF GENE]: DYRK1A (in 54 sec.)
=======================================================================
```


### `target_REPORT`

`target_REPORT` is the command-line used to generate reports on one or many targets here are the available flags:

```
One of the three following option is required:
        -g : enter a single gene name
        -i : Name of the input file with a list of genes (.txt - 1 gene per line)
        -l : Enter a list of gene name separated by a ","
	
One of these two option is required
        -rl : produce a single output file in the form of a list with quantitative informations on the targets"
        -rs : produce an output file for each target listed with detailed informations (Caution : 1 File created per target
        
        -v : Print detailed information in the terminal
        
        -update_config : use this if you want to update the config file

```
When using `target_REPORT` for the first time you will be asked to enter information about:
+ targetDB sqlite database file (not asked if you already used `target_DB`)
+ path to save list output files 
+ path to save detailed target output files
+ email address (used for pubmed searches if none provided no pubmed search will be run) (not asked if you already used `target_DB`)

Those informations will be stored in ~/.druggability/config.ini 

If you want to modify the informations in that file just use the `-update_config` flag

##### example 
input:
```
target_REPORT -g DYRK1A -rs -v
```
output:
```
[FILE]: File for  DYRK1A  generated [ /path/to/file/DYRK1A_Q13627.xlsx ]
```
input:
```
target_REPORT -l DYRK1A,TREM2,APP,MAPT -rl -v
```
output:
```
[EXPORT]: Excel file: [Export_4_entries_10Aug2018_120530.xlsx] successfully generated
```

examples output files can be found in the folder `docs`:
+ `single_output.xlsx`
+ `list_output.xlsx`

more details on the information displayed and the datasources used can be found in those files
