TargetDB
=========

TargetDB is a tool to quickly querry multiple publicly available databases

Installation
------------
### Python package installation
#### Pip installation

```
pip install targetDB
```

##### Compatibility

python version >= 3.2

Preferred python package : [Anaconda 3](https://www.anaconda.com/download/)

****

### Database files installation
This package relies on the use of sqlite database to properly function. You will be required to point at these files 
the first time you run the script.

The list of required databases is :

+ ChEMBL version 24
+ targetDB
+ tcrd_v5.2.0

##### ChEMBL

ChEMBL sqlite database can be directly downloaded [HERE](ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24_1/chembl_24_1_sqlite.tar.gz)

##### targetDB

You can download a copy of the database [HERE](https://github.com/sdecesco/targetDB/releases/download/v0.1/TargetDB_v1.db.tar.gz)

##### TCRD
tcrd database can only be downloaded as a MySQL dump [HERE](http://juniper.health.unm.edu/tcrd/download/tcrd_v5.2.0.sql.gz)

to convert it into a sqlite database please follow the [instructions](https://github.com/webyrd/mediKanren/tree/master/pharos)

>#####Note
>
>dependencies on the TCRD database will be removed in future updates

Usage
-----
