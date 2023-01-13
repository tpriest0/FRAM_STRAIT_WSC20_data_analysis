NAME=$1

awk '{FS=OFS="\t"}{if ($6 ~ /ribosomal protein L11/ && $10 ~ /K02867/) print FILENAME,">"$2,">"$2"%%%RBP_L11"; \
else if ($6 ~ /ribosomal protein L14/ && $10 ~ /K02874/) print FILENAME,">"$2,">"$2"%%%RBP_L14"; \
else if ($6 ~ /ribosomal protein L15/ && $10 ~ /K02876/) print FILENAME,">"$2,">"$2"%%%RBP_L15"; \
else if ($6 ~ /ribosomal protein L16/ && $10 ~ /K02878/) print FILENAME,">"$2,">"$2"%%%RBP_L16"; \
else if ($6 ~ /ribosomal protein L18/ && $10 ~ /K02881/) print FILENAME,">"$2,">"$2"%%%RBP_L18"; \
else if ($6 ~ /ribosomal protein L22/ && $10 ~ /K02890/) print FILENAME,">"$2,">"$2"%%%RBP_L22"; \
else if ($6 ~ /ribosomal protein L24/ && $10 ~ /K02895/) print FILENAME,">"$2,">"$2"%%%RBP_L24"; \
else if ($6 ~ /ribosomal protein L2/ && $10 ~ /K02886/) print FILENAME,">"$2,">"$2"%%%RBP_L2"; \
else if ($6 ~ /ribosomal protein L3/ && $10 ~ /K02906/) print FILENAME,">"$2,">"$2"%%%RBP_L3"; \
else if ($6 ~ /ribosomal protein L4/ && $10 ~ /K02926/) print FILENAME,">"$2,">"$2"%%%RBP_L4"; \
else if ($6 ~ /ribosomal protein L5/ && $10 ~ /K02931/) print FILENAME,">"$2,">"$2"%%%RBP_L5"; \
else if ($6 ~ /ribosomal protein L6/ && $10 ~ /K02933/) print FILENAME,">"$2,">"$2"%%%RBP_L6"; \
else if ($6 ~ /ribosomal protein S3/ && $10 ~ /K02982/) print FILENAME,">"$2,">"$2"%%%RBP_S3"; \
else if ($6 ~ /ribosomal protein S8/ && $10 ~ /K02994/) print FILENAME,">"$2,">"$2"%%%RBP_S8"; \
else if ($6 ~ /ribosomal protein S17/ && $10 ~ /K02961/) print FILENAME,">"$2,">"$2"%%%RBP_S17"; \
else if ($6 ~ /ribosomal protein S19/ && $10 ~ /K02965/) print FILENAME,">"$2,">"$2"%%%RBP_S19"}' \
${NAME}/${NAME}_prokarya_gene_table_complete_with_taxonomy.txt | awk '{FS=OFS="\t"}{gsub("/.*","",$1)}1' > ${NAME}/${NAME}_SC_RBP.tmp

awk '{FS=OFS="\t"}FNR==NR{a[$2]=$3;next}{if ($1 in a) print a[$1]; else print $1}' ${NAME}/${NAME}_SC_RBP.tmp ${NAME}/${NAME}_prokarya_genes.fa | grep -A 1 '%%%RBP_' | sed 's/--//g' | sed '/^$/d' > ${NAME}/${NAME}_prokarya_SC_RBP_gene_nucleotide_seqs.fa

awk '{FS=OFS="\t"}FNR==NR{a[$2]=$3;next}{if ($1 in a) print a[$1]; else print $1}' ${NAME}/${NAME}_SC_RBP.tmp ${NAME}/${NAME}_prokarya_genes.faa | grep -A 1 '%%%RBP_' | sed 's/--//g' | sed '/^$/d' > ${NAME}/${NAME}_prokarya_SC_RBP_gene_amino_acid_seqs.faa

mkdir -p RBP_protein_analysis

if [[ -s ${NAME}/${NAME}_prokarya_SC_RBP_gene_nucleotide_seqs.fa && -s ${NAME}/${NAME}_prokarya_SC_RBP_gene_amino_acid_seqs.faa ]] ; then
    rm ${NAME}/${NAME}_SC_RBP.tmp &&
    cp ${NAME}/${NAME}_prokarya_SC_RBP_gene_amino_acid_seqs.faa RBP_protein_analysis/ &&
    cp ${NAME}/${NAME}_prokarya_SC_RBP_gene_nucleotide_seqs.fa RBP_protein_analysis/
else
    echo "Output files were empty, something went wrong :("
fi


