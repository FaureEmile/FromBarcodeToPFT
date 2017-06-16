# We create a file for the focal lineage
mkdir /media/DATA/Emile/script_sp_euka/$ESPECE
cd /media/DATA/Emile/script_sp_euka/$ESPECE

# Save a temp file with all the lines of metaB corresponding to your lineage
zcat /media/DATA/Emile/globaldataset.barcode.v95ts.gz | awk -F"\t" -v sp=$ESPECE '{OFS="\t"} NR==1{print $0} NR>1 && $7==sp{print}' > temp.sp

# Then select only MetaB ID, taxogroup, and abundance at each station  
cut -f1,7,49- temp.sp > temp2.sp
rm temp.sp  # cleaning

# Get the headers
head -n1 temp2.sp | tr "\t" "\n" > temp.h

# For each line, we take all values and switch them into a column facing the headers
# Our two first lines are then ID and taxogroup
# For lines >2, we have station in column 1 and abundance in column 2
# We then use the info included in the script headers that are specific to each species (See Script_creationPFT.sh)
# to retrieve the PFT at each location for each size class :

n=0
tail -n+2 temp2.sp | sed '/^$/d' | while read i
do
        n=$(($n+1))
        rm -f temp.l.$n temp.h.l.$n
        echo -e "$i" | tr "\t" "\n" >> temp.l.$n
        paste temp.h temp.l.$n > temp.h.l.$n
        rm -f temp.l.$n

        rm -f temp.h.l.tot
        cat temp.h.l.$n | awk -F"\t" -v pft1=$SIZE1 -v pft2=$SIZE2 -v pft3=$SIZE3 -v pft4=$SIZE4 'BEGIN{OFS="\t"} 
        NR==1{idbc=$2}
        NR==2{sp=$2} 
        NR>2{ if(($1~"0.8-3") || ($1~"0.8-5")){print pft1,idbc, sp,$1, $2} 
        if(($1~"5-20")||($1~"3-20")) {print pft2,idbc, sp,$1, $2}
        if(($1~"20-180")||($1~"20-200")){print pft3,idbc, sp,$1, $2}
        if(($1~"180-2000")||($1~"200-2000")){print pft4,idbc, sp,$1, $2}  
        }' > temp.h.l.tot.$n
        rm -f temp.h.l.$n
done

# We concatenate all lines into a file with lineage name and .OK extention
find . -name "temp.h.l.tot.*" -exec cat {} \; > $ESPECE.OK

## Cleaning
find . -name "temp.h.l.tot.*" -exec rm {} \;
rm -f temp2.sp
rm -f temp.h
