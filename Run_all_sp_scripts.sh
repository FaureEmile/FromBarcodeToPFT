# You need to go to the directory where scripts for all your phylogenetic lineages are
cd /media/DATA/Emile/script_sp_euka

# Then just run all of them
for i in `find /media/DATA/Emile/script_sp_euka -name "script_*" -type f ! -size 0`
do
bash $i
done

