#
# this awk script realizes following actions on a featureCounts.txt file:
# - round the count (some counting options provide not integer counts)
# - suppression of rRNA & tRNA genes (Clostridium difficiles 630 species, tmRNA suppressed too)
# - get the integer part of counts (when using the -F option of featureCounts)
#TODO - add counts of different fragments of a gene (based on the "_p" information into the gene name)
# run with:
# awk -f thisAwkScript.awk featureCounts_resultFile.txt
#
BEGIN{
  FS="\t";
  OFS="\t";
}
{
  if(($1!~/^#/)&&($1!~/CD630_r/)&&($1!~/CD630_t/)&&($1!~/CD630_s0600/)){ # _s0600: tmRNA
    print $1,$2,$3,$4,$5,$6,int($7)
  }
}
