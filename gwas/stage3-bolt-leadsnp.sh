chmod +x ${NEARESTGENE_BIN}

####
# Find Genome-wide significant loci / lead SNPs
#### 

>&2 echo "Processing ${JOINED_OUTPUT}"

# Should contain the name of our GWAS trait
HOLDING_FOLDER="$(basename $(dirname $JOINED_OUTPUT))"

####
# Find lead SNPs
####
cat < ${JOINED_OUTPUT} | awk -F $'\t' 'NR == 1 || $16 < '${GWAS_P_THRESHOLD}' {print $0}' > gwas.snps

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

Rscript \
    --vanilla ${BOLT_CLUMP} \
    --file "gwas.snps" \
    --leadsnp_file tmp.leadsnp \
    --locusid_file "${JOINED_OUTPUT_LOCUS_IDS}" \
    --pvalue=${GWAS_P_THRESHOLD} \
    --locusRadius ${GWAS_LOCUS_RADIUS}

# Add the info about which GWAS this came from directly to the leadsnp file
cat tmp.leadsnp | awk -F $'\t' -v folder="$HOLDING_FOLDER" 'NR==1 {print $0 FS "GWASTrait"}; NR > 1{print $0 FS folder}' > "${JOINED_OUTPUT_LEAD}"

# Annotate with the nearest gene
${NEARESTGENE_BIN} \
  --assembly ${ASSEMBLY} \
  --sites <(tail -n +2 ${JOINED_OUTPUT_LEAD} | awk '{print $2 ":" $3}') \
  > tmp.out

paste ${JOINED_OUTPUT_LEAD} tmp.out > ${ANNOTATED_JOINED_OUTPUT_LEAD}

####
# Find suggestive loci / lead SNPs
#### 

>&2 echo "Processing ${JOINED_OUTPUT}"

cat < ${JOINED_OUTPUT} | awk -F $'\t' 'NR == 1 || $16 < '${SUGGESTIVE_P_THRESHOLD}' {print $0}' > suggestive.snps

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

Rscript \
    --vanilla ${BOLT_CLUMP} \
    --file "suggestive.snps" \
    --leadsnp_file tmp.suggestive \
    --locusid_file discardable.file \
    --pvalue=${SUGGESTIVE_P_THRESHOLD} \
    --locusRadius ${GWAS_LOCUS_RADIUS}

# Add the info about which GWAS this came from directly to the leadsnp file
cat tmp.suggestive | awk -F $'\t' -v folder="$HOLDING_FOLDER" 'NR==1 {print $0 FS "GWASTrait"}; NR > 1{print $0 FS folder}' > "${JOINED_OUTPUT_SUGGESTIVE}"

# Annotate with the nearest gene
${NEARESTGENE_BIN} \
  --assembly ${ASSEMBLY} \
  --sites <(tail -n +2 ${JOINED_OUTPUT_SUGGESTIVE} | awk '{print $2 ":" $3}') \
  > tmp.out

paste ${JOINED_OUTPUT_SUGGESTIVE} tmp.out > ${ANNOTATED_JOINED_OUTPUT_SUGGESTIVE}