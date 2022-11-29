#!/usr/bin/env bash

################################################################################
# Script to collect GWAS summary statistics from IEU OpenGWAS
#
# Usage: bash ieu_fetch.sh <study_id>
#   <study_id> : a string with the study ID to download
#
# Author: David Piñeyro
# Date: 2022-02-23
################################################################################

STUDY="$1"

if [ -z "$STUDY" ]; then
    echo "[ERROR] No study ID detected."
    exit 1
fi

# Collecting all studies in JSON format using the API
curl -X GET "http://gwas-api.mrcieu.ac.uk/gwasinfo" -H  "accept: application/json" -H  "X-Api-Token: null" > all_meta.json

# Query all_studies.json
# Check for the study Id to directly collect it
if [[ $STUDY != "NA" ]] || [[ $STUDY != "" ]]; then
    echo $STUDY > selected_studies.txt
else
    echo "[ERROR] Please, include the study Id in the input file"
    exit 1
fi

# Use the obtained ids to download their GWAS summary statistics VCF files.
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    echo "Downloading study $gwas_id..."
    wget https://gwas.mrcieu.ac.uk/files/"$gwas_id"/"$gwas_id".vcf.gz
    echo "Study $gwas_id downloaded"
done < selected_studies.txt

# Collecting metadata for the downloaded studies
printf "note\tmr\tyear\tgroup_name\tauthor\tconsortium\tsex\tpriority\tpopulation\tunit\tnsnp\tsample_size\tbuild\ttrait\tid\tsubcategory\tcategory\tontology\n" > selected_studies_metadata.tsv
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    cat all_meta.json | jq -c ".[] | select(.id == \"$gwas_id\") | \"\(.note);;\(.mr);;\(.year);;\(.group_name);;\(.author);;\(.consortium);;\(.sex);;\(.priority);;\(.population);;\(.unit);;\(.nsnp);;\(.sample_size);;\(.build);;\(.trait);;\(.id);;\(.subcategory);;\(.category);;\(.ontology)\"" >> selected_studies_metadata.tsv
done < selected_studies.txt
# Remove quotes, brackets, backslash and change ";;" to \t
sed -e 's/"//g' -e 's/\[//g' -e 's/\]//g' -e 's/\\//g' -e 's/;;/\t/g' selected_studies_metadata.tsv > "$STUDY"_metadata_ieu.tsv

# Print sample_size as a separated file to be able to add it to the VCF GWAS column later
awk -F "\t" '{if (NR==2) print $12}' "$STUDY"_metadata_ieu.tsv > "$STUDY"_sample_size.txt
