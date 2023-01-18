#! /bin/bash

# Import BOLT
gsutil -m cp -n gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/bolt ~/
chmod +x ~/bolt
mkdir -p ~/lib/
gsutil -m cp -n gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/lib/libiomp5.so ~/lib/

~/bolt \
    --bed=${MOUNT_PLINK_BUCKET}/genotype/ukb_cal_chr{1:22}_v2.bed \
    --bed=${MOUNT_PLINK_BUCKET}/genotype/ukb_cal_chrX_v2.bed \
    --bim=${MOUNT_PLINK_BUCKET}/genotype/ukb_snp_chr{1:22}_v2.bim \
    --bim=${MOUNT_PLINK_BUCKET}/genotype/ukb_snp_chrX_v2.bim \
    --fam=${FAM_FILE} \
    --remove=${EXCLUSION_FILENAME} \
    --remove=${CONSENT_WITHDRAWN} \
    --remove=${NOT_IMPUTED} \
    --remove=${NOT_IMPUTED_CHRX} \
    --phenoFile=${PHENO_FILENAME} \
    --phenoCol=${PHENO_COL} \
    --covarFile=${PHENO_FILENAME} \
    --sampleFile=${SAMPLE_FILE} \
    --bgenFile=${MOUNT_IMPUTED}/imputed/ukb_imp_chrX_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --statsFileBgenSnps="${IMPUTED_OUTPUT}" \
    --covarMaxLevels=1000 \
    --covarCol=array_UKBB \
    --covarCol=male \
    --qCovarCol=age_at_mri \
    --qCovarCol=PC{1:5} \
    --LDscoresFile=${LDSCORES_FILE} \
    --geneticMapFile=${GENETIC_MAP_FILE} \
    --LDscoresMatchBp \
    --lmm \
    --numThreads=$(nproc) \
    --maxMissingPerSnp=0.05 \
    --maxMissingPerIndiv=0.5 \
    --statsFile="${GT_OUTPUT}" \
    --lmmForceNonInf \
    --verboseStats
