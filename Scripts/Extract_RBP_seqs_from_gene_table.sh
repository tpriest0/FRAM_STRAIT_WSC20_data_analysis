#!/bin/bash
# Specify variable from input
NAME=$1

# Using awk, extract gene names if the Pfam and KEGG annotation columns contain a matching ribosomal protein gene annotation
# and then place names into a new file along with reformatted version that include the ribosomal protein name along with the gene name combined
awk '{FS=OFS="\t"}{if ($6 ~ /ribosomal protein L11/ && $9 ~ /K02867/) print FILENAME,">"$2,">"$2"%%%RP_L11"; \
else if ($6 ~ /ribosomal protein L14/ && $9 ~ /K02874/) print FILENAME,">"$2,">"$2"%%%RP_L14"; \
else if ($6 ~ /ribosomal protein L15/ && $9 ~ /K02876/) print FILENAME,">"$2,">"$2"%%%RP_L15"; \
else if ($6 ~ /ribosomal protein L16/ && $9 ~ /K02878/) print FILENAME,">"$2,">"$2"%%%RP_L16"; \
else if ($6 ~ /ribosomal protein L18/ && $9 ~ /K02881/) print FILENAME,">"$2,">"$2"%%%RP_L18"; \
else if ($6 ~ /ribosomal protein L22/ && $9 ~ /K02890/) print FILENAME,">"$2,">"$2"%%%RP_L22"; \
else if ($6 ~ /ribosomal protein L24/ && $9 ~ /K02895/) print FILENAME,">"$2,">"$2"%%%RP_L24"; \
else if ($6 ~ /ribosomal protein L2/ && $9 ~ /K02886/) print FILENAME,">"$2,">"$2"%%%RP_L2"; \
else if ($6 ~ /ribosomal protein L3/ && $9 ~ /K02906/) print FILENAME,">"$2,">"$2"%%%RP_L3"; \
else if ($6 ~ /ribosomal protein L4/ && $9 ~ /K02926/) print FILENAME,">"$2,">"$2"%%%RP_L4"; \
else if ($6 ~ /ribosomal protein L5/ && $9 ~ /K02931/) print FILENAME,">"$2,">"$2"%%%RP_L5"; \
else if ($6 ~ /ribosomal protein L6/ && $9 ~ /K02933/) print FILENAME,">"$2,">"$2"%%%RP_L6"; \
else if ($6 ~ /ribosomal protein S3/ && $9 ~ /K02982/) print FILENAME,">"$2,">"$2"%%%RP_S3"; \
else if ($6 ~ /ribosomal protein S8/ && $9 ~ /K02994/) print FILENAME,">"$2,">"$2"%%%RP_S8"; \
else if ($6 ~ /ribosomal protein S17/ && $9 ~ /K02961/) print FILENAME,">"$2,">"$2"%%%RP_S17"; \
else if ($6 ~ /ribosomal protein S19/ && $9 ~ /K02965/) print FILENAME,">"$2,">"$2"%%%RP_S19"}' \
Metagenome_gene_annotation_tables/${NAME}_MG_gene_table.txt | awk '{FS=OFS="\t"}{gsub("/.*","",$1)}1' > Metagenome_gene_annotation_tables/${NAME}_RBP_gene_headers.txt

# Using the list of gene names generated above, use an awk array to extract the gene sequences, first from the amino acid
# gene files and then from the nucleotide gene files
awk '{FS=OFS="\t"}FNR==NR{a[$2]=$3;next}{if ($1 in a) print a[$1]; else print $1}' Metagenome_gene_annotation_tables/${NAME}_RBP_gene_headers.txt Metagenome_gene_sequences/${NAME}_MG_genes_amino_acid.faa | \
grep -A 1 '%%%RP_' | grep -v '^--' > Metagenome_gene_sequences/${NAME}_RBP_genes_amino_acid_seqs.faa

awk '{FS=OFS="\t"}FNR==NR{a[$2]=$3;next}{if ($1 in a) print a[$1]; else print $1}' Metagenome_gene_annotation_tables/${NAME}_RBP_gene_headers.txt Metagenome_gene_sequences/${NAME}_MG_genes_nucleotide.fa | \
grep -A 1 '%%%RP_' | grep -v '^--' > Metagenome_gene_sequences/${NAME}_RBP_genes_nucleotide.fa

# Check to see if the ribosomal protein gene sequence files exist and are not empty, and then delete the intermediate files
if [[ -s Metagenome_gene_sequences/${NAME}_RBP_genes_amino_acid_seqs.faa && -s Metagenome_gene_sequences/${NAME}_RBP_genes_nucleotide.fa ]] ; then
    rm Metagenome_gene_annotation_tables/${NAME}_RBP_gene_headers.txt
else
    echo "Output files were empty, something went wrong :("
fi
