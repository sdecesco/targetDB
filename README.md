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

This package relies on the use of sqlite database to properly function.

##### Required database

+ targetDB

You can download a copy of the database [HERE](https://github.com/sdecesco/targetDB/releases/download/v1.1.6/TargetDB_20_02_19.db.zip)

>This database contains all of the human genome genes that have a uniprot ID

#### System compatibility

+ Linux
+ Windows
+ Mac (not tested)


Usage
-----
targetDB package provides a user interface to use the tool: `targetDB`

### `targetDB`

When using `targetDB` for the first time you will be asked to enter information about:
+ targetDB sqlite database file
+ path to save list output files 
+ path to save detailed target output files
+ email address (used for pubmed searches if none provided no pubmed search will be run)

Those informations will be stored in ~/.targetdb/config.ini 

![Configuration panel](targetDB/resources/configuration.png)

Once created it will automatically start the main user interface (as seen below)

![Main interface](targetDB/resources/targetdb_gui.png)

Instructions to create a targetDB database from scratch
---

#### System compatibility
+ Linux

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
