# We go where all the .OK files are (see end_scriptPFT.sh)
cd /media/DATA/Emile/OK_files

# For all taxogroup we want to get rid of doubled lines and values of 0 in abundance

for i in `find . -name "*.OK" -type f ! -size 0`
do
awk '!a[$0]++' $i | awk '$NF!=0' > /media/DATA/Emile/OK2_Files/$i"_"2
done


## Finally we pool the results into one file for all eukaryotes
cd /media/DATA/Emile/OK2_Files
rm -f alleuka.OK_2
cat *.OK_2 | sed '/^$/d' >> alleuka.OK2
mv alleuka.OK2 alleuka.OK_2


