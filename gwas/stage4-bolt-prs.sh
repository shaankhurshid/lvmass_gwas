chmod +x ${APPLYPRSBASIC}

# SAIGE custom-layout: 5,4,5,0,1,10
# BOLT custom-layout: 4,4,5,1,2,10

${APPLYPRSBASIC} \
-bgen-template "${MOUNTED_GENOME}/imputed/ukb_imp_chr%s_v3.bgen" \
-input ${INPUT_FILE} \
-source ${SOURCE} \
-custom-layout 4,4,5,1,2,10 \
-sample ${SAMPLE} \
> ${OUTPUT_FILE}

${APPLYPRSBASIC} \
-bgen-template "${MOUNTED_GENOME}/imputed/ukb_imp_chr%s_v3.bgen" \
-input ${INPUT_FILE_SUGGESTIVE} \
-source ${SOURCE} \
-custom-layout 4,4,5,1,2,10 \
-sample ${SAMPLE} \
> ${OUTPUT_FILE_SUGGESTIVE}