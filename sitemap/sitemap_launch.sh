#!/bin/bash

# COMMAND TO LAUNCH PREPWIZARD
# REPLACE XX by the pdb name/code

$SCHRODINGER/utilities/prepwizard -j XXX_prep -LOCAL -HOST localhost:12 XXX.pdb XXX_prep.mae

# COMMAND FOR SITEMAP

$SCHORDINGER/sitemap -LOCAL -HOST localhost:12 -i sitemap_param.in -prot XXX_prep.mae


