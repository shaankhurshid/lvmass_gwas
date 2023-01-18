# Assumes this is being run within docker / on dsub

set -e

DATAROOT="${LDSC_MOUNT}/projects/jamesp/data/ldsc"

FILE="${INPUT_FILE}"

echo ""
echo "Processing ${trait}"
echo "Finding files at ${FILE}"
echo "Will write output to ${OUTPUT_FILE}"

echo "Postprocessing"

# As of the time of this script, this was the BOLT-LMM header around which these
# column numbers were based:
#
# SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	CHISQ_LINREG	P_LINREG	BETA	SE	CHISQ_BOLT_LMM_INF	P_BOLT_LMM_INF	CHISQ_BOLT_LMM	P_BOLT_LMM
# Screen out erroneous SE and ChiSq values
zcat ${FILE} | awk -F $'\t' 'NR == 1{print $0}; NR > 1 && $12>0 && $15>0 {print $0}' > ${FILTERED}
# munge_sumstats chokes on similar column names even if you explicitly state the
# desired ones, so remove similarly named but unused columns
cat ${FILTERED} | awk -v samplesize=$(wc -l < ${PHENO_FILE}) -F $'\t' 'NR == 1{print $1 FS $2 FS $3 FS $5 FS $6 FS $7 FS $8 FS $11 FS $12 FS $16 FS "N"}; NR > 1 {print $1 FS $2 FS $3 FS $5 FS $6 FS $7 FS $8 FS $11 FS $12 FS $16 FS samplesize}' > filtered.tsv

# Since we're seeing P-values as extreme as 1E-3739, need to cap at a value that
# is permitted by LDSC. This code is from
# https://github.com/bulik/ldsc/issues/144#issue-411843650 and -- treating the
# P-values like a string -- forces everything from 1E-9999 to 1E300 down to
# 1E300.
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9][0-9]/1.0E-300/' filtered.tsv
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' filtered.tsv
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' filtered.tsv

echo "First 2 lines of the filtered file"
head -n 2 ${FILTERED}

echo "First 2 lines of ldsc input file"
head -n 2 filtered.tsv

echo "Munging"

# Note that ***A1 in ldsc is the the effect allele***:
# https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

munge_sumstats.py \
--sumstats filtered.tsv \
--p P_BOLT_LMM \
--a1 ALLELE1 \
--a2 ALLELE0 \
--merge-alleles ${DATAROOT}/w_hm3.snplist \
--out filtered.munged

cp filtered.munged.sumstats.gz ${MUNGED}

echo "Running LDSC"

ldsc.py \
--h2 ${MUNGED} \
--ref-ld-chr ${DATAROOT}/baselineLD_v2.2/baselineLD. \
--w-ld-chr ${DATAROOT}/weights_hm3_no_hla/weights. \
--out out.file

cp out.file.log ${OUTPUT_FILE}