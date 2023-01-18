#! /bin/bash

# Import BOLT
gsutil -m cp -n gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/bolt ~/
chmod +x ~/bolt
mkdir -p ~/lib/
gsutil -m cp -n gs://ukbb_v2/projects/jamesp/lib/BOLT-LMM_v2.3.4/lib/libiomp5.so ~/lib/

~/bolt \
    --bfile=${MOUNT_PLINK_BUCKET}/${PLINK_MOUNTED_PATH} \
    --remove=${EXCLUSION_FILENAME} \
    --remove=${CONSENT_WITHDRAWN} \
    --remove=${NOT_IMPUTED} \
    --phenoFile=${PHENO_FILENAME} \
    --phenoCol=${PHENO_COL} \
    --covarFile=${PHENO_FILENAME} \
    --sampleFile=${SAMPLE_FILE} \
    --bgenFile=${MOUNT_IMPUTED}/imputed/ukb_imp_chr{1:22}_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --statsFileBgenSnps=${IMPUTED_OUTPUT} \
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
    --statsFile=${GT_OUTPUT} \
    --lmmForceNonInf \
    --verboseStats
