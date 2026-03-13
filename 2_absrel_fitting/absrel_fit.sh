GENE="$1"

BASE_DIR="./${GENE}/hyphy_aBSREL_${GENE}"

hyphy CPU=1 absrel \
  --alignment "${BASE_DIR}/combined.data" \
  --multiple-hits Double+Triple \
  --output "${BASE_DIR}/${GENE}_absrel.json" > "${BASE_DIR}/output_${GENE}_absrel.txt"

