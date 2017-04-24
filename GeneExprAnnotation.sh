#!/bin/bash

## I downloaded the data from the ensembl ftp to annotate the genes 

ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf  
# From the GTF file, I take the feature_types "gene" and print chr,start, end. There are a total of 63,677
awk 'BEGIN{FS="\t"; OFS="\t"} {if ($3 =="gene") print $1,$4,$5,$9} ' Homo_sapiens.GRCh37.75.gtf > allgenes.tsv
awk 'BEGIN{FS=";"; OFS="\t"} {print $1,$2} ' allgenes.tsv | sed 's/ gene_name "//g' | sed 's/"$//g' > allgenes.txt

# To consider only the "protein_coding"  genes # There are a total of 22,810 
grep "protein_coding" allgenes.tsv | awk 'BEGIN{FS=";"; OFS="\t"} {print $1,$2} ' | sed 's/ gene_name "//g' | sed 's/"$//g' > coding_genes.txt


##I downloaded the data from GENCODE (Relaese 19 GRCH37) that correspond to ensemble 74. 
http://www.gencodegenes.org/releases/
# From the GTF file, I take the feature_types "gene" and print chr,start, end. There are a total of 57,820
awk 'BEGIN{FS="\t"; OFS="\t"} {if ($3 =="gene") print $1,$4,$5,$9} ' gencode.v19.annotation.gtf > allgenes.tsv
awk 'BEGIN{FS=";"; OFS="\t"} {print $1,$3,$5} ' allgenes.tsv > allgenes.txt
awk 'BEGIN{FS=" "; OFS="\t"} {print $1,$2,$3,$5,$7,$9} ' allgenes.txt > allgenes2.txt

# To consider only the "protein_coding"  genes # There are a total of 20,345 
grep "protein_coding" allgenes2.txt | awk 'BEGIN{FS=";"; OFS="\t"} {print} '  > coding_genes.txt

##evidence-based annotation of the human genome (GRCh38), version 23 (Ensembl 81)
# From the GTF file, I take the feature_types "gene" and print chr,start, end. There are a total of 
awk 'BEGIN{FS="\t"; OFS="\t"} {if ($3 =="gene") print $1,$4,$5,$9} ' gencode.v23.annotation.gtf > allgenesCh38.tsv
awk 'BEGIN{FS=";"; OFS="\t"} {print $1,$2,$4} ' allgenesCh38.tsv > allgenesCh38.txt
awk 'BEGIN{FS=" "; OFS="\t"} {print $1,$2,$3,$5,$7,$9} ' allgenesCh38.txt > allgenesCh382.txt

# To consider only the "protein_coding"  genes # There are a total of 19,797
grep "protein_coding" allgenesCh382.txt | awk 'BEGIN{FS=";"; OFS="\t"} {print} '  > coding_genes.Ch38.txt

