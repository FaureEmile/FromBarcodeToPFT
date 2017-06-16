# We use a table with 5 columns, one for the phylogenetic lineage and 4 for the PFT corresponding to the different size fractions used
# The next line creates the head of a script for each phylogenetic lineage read in the table

awk -F"\t" '{print "ESPECE=\""$1"\"\nSIZE1=\""$2"\"\nSIZE2=\""$3"\"\nSIZE3=\""$4"\"\nSIZE4=\""$5"\"\n" > "script_"$1".sh"}' pft_table.csv

# We obtained script heads that are specific to each phylogenetic lineage, and now we add the content of end_scriptPFT.sh to each of them
for i in script_*
do
cat end_scriptPFT.sh >> $i 
done
