#!/bin/bash

# Author: James Pirruccello <jamesp@broadinstitute.org>
# June 29, 2019.

# This script runs a UK Biobank GWAS using BOLT on GCP using dsub. First,
# please install dsub: https://github.com/DataBiosphere/dsub#install-dsub .
# Your dsub version must be at least 0.3.4

BILLING_PROJECT="broad-ml4cvd"

# This script assumes that you have access to the Broad-wide genetic data stored under
# gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183

# There is a manual step: please create a tab-delimited phenotype file at
# $SHARED_PROJECT_PATH/pheno.tsv . Output from this project will ultimately go to
# $PROJECT_PATH/output/* .

# VERSION="v23"
# PHENO_COL="lv_mass"
# SHARED_PROJECT_PATH="gs://ml4cvd/jamesp/lvmass/gwas/${VERSION}"
# PROJECT_PATH="${SHARED_PROJECT_PATH}"

VERSION="v19_seg"
PHENO_COL="lvmi_seg_adjusted"
SHARED_PROJECT_PATH="gs://ml4cvd/projects/skhurshid/ecg_lvh/gwas_actual/${VERSION}"
PROJECT_PATH=${SHARED_PROJECT_PATH}/${PHENO_COL}

# Choose one of these options for the creation of the GRM
# Larger filesets are probably more accurate but run more slowly.
# As a ballpark, ~16k SNPs take 90min for Step 1, while 
# 170k SNPs take 9.5 hours for Step 1.
#GENO_FILESET="pruned-ukbb-16k"
#GENO_FILESET="pruned-ukbb-37k"
# GENO_FILESET="pruned-ukbb-62k"
# GENO_FILESET="pruned-ukbb-141k"
GENO_FILESET="merged-ukbb"

# Allow this many preemptions. The script runs stage 2 twice: first preemptibly, and
# then again without preemption. It won't duplicate work, but this ensures that any 
# job that didn't finish with preemption eventually gets done before this script completes.
MAX_PREEMPTION=2

# # #
# Should not usually need to edit below here
# # #

# Project linker
SAMPLE_FILE="gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/ukb7089_imp_chr15_v3_s487395.sample"
SAMPLE_FILE_CHRX="gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/ukb7089_imp_chrX_v3_s486743.sample"

# Exclusion files
CONSENT_WITHDRAWN="gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/consent_withdrawn_20181016"
NOT_IMPUTED="gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/bolt.in_plink_but_not_imputed.FID_IID.968.txt"
NOT_IMPUTED_CHRX="gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/bolt.in_plink_but_not_imputed.FID_IID.13.txt"

# BOLT inputs
LDSCORES_FILE="gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz"
GENETIC_MAP_FILE="gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz"

# Which SAIGE docker image do you want to run?
SAIGE_DOCKER_IMAGE="gcr.io/ukbb-analyses/saige:0.35.8.8"

# Enable exit on error
set -o errexit

# Check for errors when we can
echo "Checking to make sure that ${SHARED_PROJECT_PATH}/pheno.tsv exists"
gsutil ls ${SHARED_PROJECT_PATH}/pheno.tsv

echo "Checking to make sure that ${SHARED_PROJECT_PATH}/exclusion.tsv exists"
gsutil ls ${SHARED_PROJECT_PATH}/exclusion.tsv

# Launch step 1 and block until completion
echo "Launching BOLT with up to ${MAX_PREEMPTION} preemptible attempts"
dsub \
    --project ${BILLING_PROJECT} \
    --provider google-v2 \
    --use-private-address \
    --regions us-central1 us-east1 us-west1 \
    --disk-type pd-standard \
    --disk-size 80 \
    --min-cores 64 \
    --min-ram 80 \
    --image ${SAIGE_DOCKER_IMAGE} \
    --skip \
    --wait \
    --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/bolt \
    --input PHENO_FILENAME=${SHARED_PROJECT_PATH}/pheno.tsv \
    --input EXCLUSION_FILENAME=${SHARED_PROJECT_PATH}/exclusion.tsv \
    --input CONSENT_WITHDRAWN=${CONSENT_WITHDRAWN} \
    --input NOT_IMPUTED=${NOT_IMPUTED} \
    --input LDSCORES_FILE=${LDSCORES_FILE} \
    --input GENETIC_MAP_FILE=${GENETIC_MAP_FILE} \
    --input SAMPLE_FILE=${SAMPLE_FILE} \
    --mount MOUNT_IMPUTED="gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183" \
    --mount MOUNT_PLINK_BUCKET="gs://ukbb_v2" \
    --env PLINK_MOUNTED_PATH="projects/jamesp/data/genotypes/${GENO_FILESET}" \
    --env PHENO_COL=${PHENO_COL} \
    --output GT_OUTPUT=${PROJECT_PATH}/bolt.gt.tsv.gz \
    --output IMPUTED_OUTPUT=${PROJECT_PATH}/bolt.imputed.tsv.gz \
    --name ${VERSION}-${PHENO_COL:0:30}-bolt \
    --label "pheno=${PHENO_COL,,}" \
    --label "version=${VERSION}" \
    --label 'stage=1' \
    --label 'script=bolt' \
    --script bolt.sh

# Launch step 1 and block until completion
echo "Launching BOLT for ChrX"
dsub \
    --project ${BILLING_PROJECT} \
    --provider google-v2 \
    --use-private-address \
    --regions us-central1 us-east1 us-west1 \
    --disk-type pd-standard \
    --disk-size 80 \
    --min-cores 64 \
    --min-ram 80 \
    --image ${SAIGE_DOCKER_IMAGE} \
    --skip \
    --wait \
    --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/boltX \
    --input PHENO_FILENAME=${SHARED_PROJECT_PATH}/pheno.tsv \
    --input EXCLUSION_FILENAME=${SHARED_PROJECT_PATH}/exclusion.tsv \
    --input CONSENT_WITHDRAWN=${CONSENT_WITHDRAWN} \
    --input NOT_IMPUTED=${NOT_IMPUTED} \
    --input NOT_IMPUTED_CHRX=${NOT_IMPUTED_CHRX} \
    --input LDSCORES_FILE=${LDSCORES_FILE} \
    --input GENETIC_MAP_FILE=${GENETIC_MAP_FILE} \
    --input SAMPLE_FILE=${SAMPLE_FILE_CHRX} \
    --input FAM_FILE="gs://ukbb_v2/projects/jamesp/data/ukb708_cal_chr1_v2_s488374.fixed.fam" \
    --mount MOUNT_IMPUTED="gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183" \
    --mount MOUNT_PLINK_BUCKET="gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b" \
    --env PHENO_COL=${PHENO_COL} \
    --output GT_OUTPUT=${PROJECT_PATH}/bolt-x.gt.tsv.gz \
    --output IMPUTED_OUTPUT=${PROJECT_PATH}/bolt-x.imputed.tsv.gz \
    --name ${VERSION}-${PHENO_COL:0:30}-boltX \
    --label "pheno=${PHENO_COL,,}" \
    --label "version=${VERSION}" \
    --label 'stage=1x' \
    --label 'script=boltx' \
    --script bolt-X.sh

   

# Since we can't split BOLT, don't allow preemption; too risky.
# --preemptible ${MAX_PREEMPTION} \
# --retries ${MAX_PREEMPTION} \
# Disable full GWAS

echo "Launching Step 2: Validity filtering and LDSC postprocessing"
dsub \
   --project ${BILLING_PROJECT} \
   --provider google-v2 \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --disk-type pd-standard \
   --disk-size 50 \
   --min-cores 8 \
   --min-ram 32 \
   --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/ldsc \
   --image gcr.io/ukbb-analyses/ldsc:1.0.1 \
   --preemptible ${MAX_PREEMPTION} \
   --retries ${MAX_PREEMPTION} \
   --skip \
   --wait \
   --input PHENO_FILE=${SHARED_PROJECT_PATH}/pheno.tsv \
   --input INPUT_FILE=${PROJECT_PATH}/bolt.imputed.tsv.gz \
   --mount LDSC_MOUNT=gs://ukbb_v2 \
   --output FILTERED=${PROJECT_PATH}/bolt.filtered.tsv \
   --output MUNGED=${PROJECT_PATH}/joined-ldsc-munged.tsv.gz \
   --output OUTPUT_FILE=${PROJECT_PATH}/ldsc.txt \
   --name ${VERSION}-${PHENO_COL:0:30}-ldsc \
   --label "pheno=${PHENO_COL,,}" \
   --label "version=${VERSION}" \
   --label 'stage=2' \
   --label 'script=ldsc' \
   --script stage2-bolt-ldsc.sh

echo "Launching Step 3: BOLT Lead SNPs"
dsub \
   --project ${BILLING_PROJECT} \
   --provider google-v2 \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --disk-type pd-standard \
   --disk-size 50 \
   --min-cores 8 \
   --min-ram 64 \
   --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/leadsnp \
   --image ${SAIGE_DOCKER_IMAGE} \
   --preemptible ${MAX_PREEMPTION} \
   --retries ${MAX_PREEMPTION} \
   --skip \
   --wait \
   --input NEARESTGENE_BIN="gs://ml4cvd/jamesp/bin/nearestgene" \
   --input BOLT_CLUMP=gs://ukbb_v2/projects/jamesp/bin/bolt-clump-by-p-and-distance.r \
   --input JOINED_OUTPUT=${PROJECT_PATH}/bolt.filtered.tsv \
   --env PROJECT_PATH=${PROJECT_PATH} \
   --env GWAS_P_THRESHOLD=5e-8 \
   --env GWAS_LOCUS_RADIUS=500000 \
   --env SUGGESTIVE_P_THRESHOLD=1e-6 \
   --env ASSEMBLY=37 \
   --output JOINED_OUTPUT_LEAD=${PROJECT_PATH}/bolt.filtered.tsv-leadsnps.tsv \
   --output JOINED_OUTPUT_LOCUS_IDS=${PROJECT_PATH}/bolt.filtered.tsv-locusids.tsv \
   --output JOINED_OUTPUT_SUGGESTIVE=${PROJECT_PATH}/bolt.filtered.tsv-suggestive-leadsnps.tsv \
   --output ANNOTATED_JOINED_OUTPUT_LEAD=${PROJECT_PATH}/annotated.bolt.filtered.tsv-leadsnps.tsv \
   --output ANNOTATED_JOINED_OUTPUT_SUGGESTIVE=${PROJECT_PATH}/annotated.bolt.filtered.tsv-suggestive-leadsnps.tsv \
   --name ${VERSION}-${PHENO_COL:0:30}-leadsnp \
   --label "pheno=${PHENO_COL,,}" \
   --label "version=${VERSION}" \
   --label 'stage=3' \
   --label 'script=leadsnp' \
   --script stage3-bolt-leadsnp.sh

echo "Launching Step 5: BOLT Manhattan Plot"
dsub \
   --project ${BILLING_PROJECT} \
   --provider google-v2 \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --disk-type pd-standard \
   --disk-size 20 \
   --min-cores 4 \
   --min-ram 64 \
   --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/manhattan \
   --image ${SAIGE_DOCKER_IMAGE} \
   --preemptible ${MAX_PREEMPTION} \
   --retries ${MAX_PREEMPTION} \
   --skip \
   --wait \
   --input BOLT_MANHATTAN=gs://ukbb_v2/projects/jamesp/bin/manhattan-bolt.r \
   --input INPUT_FILE=${PROJECT_PATH}/bolt.filtered.tsv \
   --env PHENO_COL=${PHENO_COL} \
   --output OUTPUT_FILE=${PROJECT_PATH}/manhattan.pdf \
   --name ${VERSION}-${PHENO_COL:0:30}-manhattan \
   --label "pheno=${PHENO_COL,,}" \
   --label "version=${VERSION}" \
   --label 'stage=5' \
   --label 'script=manhattan' \
   --script stage5-bolt-manhattan.sh

# echo "Launching Step 4: Lead-SNP PRS scoring" Save til the end so we still get
# Manhattan plot if this fails.
dsub \
   --project ${BILLING_PROJECT} \
   --provider google-v2 \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --disk-type pd-standard \
   --disk-size 30 \
   --min-cores 8 \
   --min-ram 16 \
   --logging ${SHARED_PROJECT_PATH}/${PHENO_COL}/dsub-logs/prs \
   --image ${SAIGE_DOCKER_IMAGE} \
   --preemptible ${MAX_PREEMPTION} \
   --retries ${MAX_PREEMPTION} \
   --skip \
   --wait \
   --input APPLYPRSBASIC=gs://ukbb_v2/projects/jamesp/bin/applyprsbasic \
   --input SAMPLE=gs://ukbb_v2/projects/jamesp/data/ukb7089_imp_chr15_v3_s487395.sample \
   --input INPUT_FILE=${PROJECT_PATH}/bolt.filtered.tsv-leadsnps.tsv \
   --input INPUT_FILE_SUGGESTIVE=${PROJECT_PATH}/bolt.filtered.tsv-suggestive-leadsnps.tsv \
   --mount MOUNTED_GENOME=gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183 \
   --env SOURCE=${PHENO_COL}-${VERSION} \
   --output OUTPUT_FILE=${PROJECT_PATH}/leadsnp.prs.tsv \
   --output OUTPUT_FILE_SUGGESTIVE=${PROJECT_PATH}/suggestive-leadsnp.prs.tsv \
   --name ${VERSION}-${PHENO_COL:0:30}-prs \
   --label "pheno=${PHENO_COL,,}" \
   --label "version=${VERSION}" \
   --label 'stage=4' \
   --label 'script=prs' \
   --script stage4-bolt-prs.sh
