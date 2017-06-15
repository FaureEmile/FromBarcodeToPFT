# FromBarcodeToPFT

Allow to jump from metabarcoding data to a matrix of PFT abundance


# About

This script was designed during my master 2 internship. It is designed to read a metabarcoding dataset, detecting phylogenetic lineage information, and then use a table of correspondance between phylogenetic lineages and Plankton Functional Types (PFT) in order to build a matrix of metabarcode abundance per PFT. It allows a jump between a classic metabarcoding dataset and a table that can serve as the input of a biogeochemical PFT type model, or any other type of study on PFT, for example statistical or biogeographical analysis.

In order to use this code as it is, your metabarcoding dataset should have one column giving the best guess of phylogenetic lineage for each metabarcode, and you should have a pre-built table linking phylogenetic lineages to plankton functional types. This code should be easy to run with any metabarcoding dataset by tweaking it a little bit to fit with the data construction.

In addition, you will find in this repository an R code showing some examples of which statistical and biogeographical analysis can be conducted using the PFT abundance matrix.

# Requirements

R 3.3.2

Ubuntu bash shell

# Usage

Datasets used in this project were different for eukaryotes and prokaryotes. 
Eukaryotes data came out of "real" metabarcoding, whereas for Prokaryotes miTAGs were used (see Sunagawa et al., 2015 in Science). My guess would be that most of metabarcoding datasets look like the one used for eukaryotes, but scripts working for both datasets are presented here.

### Eukaryotes

Order of run : 
1. Script_creationPFT.sh (which calls end_scriptPFT.sh)
2. Run_all_sp_scripts.sh
3. Suite_euka.sh

The last file created can then be used as input for the Rscript (read.table line 19)

### Prokaryotes

Order of run :
1. to 4. (no order needed) : Script_N2auto.sh, Script_N2hetero.sh, Script_Piconanohetero.sh, Script_Piconanophyto.sh
5. Suite_proka.sh

The last file created can then be used as input for the Rscript (read.table line 24)






