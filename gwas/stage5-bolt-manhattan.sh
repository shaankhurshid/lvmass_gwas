####
# Manhattan plot
#### 

# Fix extreme P values that R won't recognize. Since we're seeing P-values as
# extreme as 1E-3739, need to cap at a value that is permitted by LDSC. This
# code is from https://github.com/bulik/ldsc/issues/144#issue-411843650 and --
# treating the P-values like a string -- forces everything from 1E-9999 to 1E300
# down to 1E300.
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' ${INPUT_FILE}
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' ${INPUT_FILE}

# Require some trivial significance
TRIVIAL_P_THRESHOLD=1E-3
cat ${INPUT_FILE} | awk -F $'\t' 'NR == 1 || $16 < '${TRIVIAL_P_THRESHOLD}' {print $0}' > thinned.tsv

Rscript \
    --vanilla ${BOLT_MANHATTAN} \
    --file thinned.tsv \
    --trait ${PHENO_COL}

cp manhattan.pdf ${OUTPUT_FILE}