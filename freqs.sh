cut -f 1 47_Genomas_snps_table | 
awk '
BEGIN{
ignore[32]=0; ignore[34]=0; ignore[44]=0;
# the next array contains the values to consider as missing "-" and "N"
miss["-"]=0; miss["N"]=0;
maxH=27; maxL=18;
# next, the maximum "-" values which we accept for HighP and LowP
maxMissH=20; maxMissL=14;
# HEADER
print "row_number\tSNPfiltered\tnts\tnts_HP\tnts_LP\tfreq_HP\tfreq_LP\tfreq_HP-freq_LP\tSNPoriginal";
}
{
# now we create the arrays needed for each line
delete bases; delete hp; delete lp; 
# next variables count how many "-" values we have read in current line
hmiss = 0; lmiss = 0;
if (NR==1) next;
# Therefore, if we have not skipped (NR>1) 
# then the current line IS NOT HEADER, so we process it:

# 1: split the string of nucleotides into single nucleotides
# v[1] = G
# v[2] = A
# ...
# v[48] = G
split($0,v,"");

# 2: now we are going to read the "v" array position by position
finalstring="";

for (i=1;i<=length(v);i++) {
        if (i in ignore) continue; 
        currvalue=v[i];
        finalstring=finalstring""currvalue;

        # in "bases" we count how many times we have found each nucleotide
        if (currvalue!="-"){ # ignore "-" values
        if (!(currvalue in miss)) {
                if (currvalue in bases) bases[currvalue]++; 
                else bases[currvalue]=1; 
        }

# now we process separately the HighP and LowP strains
        # if it is HighP:
        if (i<=maxH) {
                if (currvalue in miss) {
                        hmiss++; # another way hmiss=hmiss+1
                        if (hmiss>maxMissH) next; # if too many missings, skip this line
                } else {
                        # in "hp" we count how many times we have found each 
                        # nucleotide, for HighP
                        if (currvalue in hp) hp[currvalue]++;
                        else hp[currvalue]=1;
                }
        }

        # if it is LowP
        if (i>maxH) {
                # check if it is missing value
                if (currvalue in miss) {
                        lmiss++; # another way lmiss+=1
                        if (lmiss>maxMissL) next; # if too many missings, skip this line
                } else {
                        if (currvalue in lp) lp[currvalue]++; 
                        else lp[currvalue]=1;
                }
        }
        
print $0;
#       print i"-"v[i];
#       print length(bases)","length(hp)","length(lp);

} 

#print $0;
#print "AFTER LOOP "length(bases)","length(hp)","length(lp);

# After the "for" loop, we have already counted the number of occurences
# of each nucleotide for HighP and LowP

# Now we need to compute frequencies from those counts
# 3: for each nucleotide we have already read, we calculate its frequency
for (base in bases) {
        if (base in hp) hpval=hp[base]; else hpval=0;
        if (base in lp) lpval=lp[base]; else lpval=0;
        # frequency of nucleotide in HighP
        currhp=hpval/(maxH-hmiss); 
        # frequency of nucleotide in LowP
        currlp=lpval/(maxL-lmiss); 
        score=currhp-currlp; 
        if (score<0) score=score*(-1);
} 

#print length(lp);

print NR"\t"finalstring"\t"length(bases)"\t"length(hp)"\t"length(lp)"\t"currhp"\t"currlp"\t"score"\t"$0;

}'

#threshold=0.7; 
#if (score<-threshold || score>threshold
