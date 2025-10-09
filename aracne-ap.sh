# aracne-ap.sh
# Gene Regulatory Networks deconvolution using ARACNe-AP (Califano lab) for paper:
# "Two Cohorts, One Network: Consensus Master Regulators 
# Orchestrating Papillary Thyroid Carcinoma"
# Hugo Tovar, National Institute of Genomic Medicine, Mexico 
# hatovar@inmegen.gob.mx

# TCGA
java -Xmx50G -jar ARACNe-AP/dist/aracne.jar -e inmat_TCGA.txt -o tcga_network --tfs extra_data/TF_lambert.txt --pvalue 1E-8 --seed 1 --calculateThreshold

for i in {1..100}
do
java -Xmx200G -jar ARACNe-AP/dist/aracne.jar -e inmat_TCGA.txt -o tcga_network --tfs extra_data/TF_lambert.txt --pvalue 1E-8 --threads 50 --seed $i
done

java -Xmx50G -jar ARACNe-AP/dist/aracne.jar -o tcga_network --consolidate

mv tcga_network/network.txt ./tcga_tumor_network.txt


# GEO
# Using p-value threshold of 1E-4 for ARACNe-AP on the GEO dataset (GSE33630).
# This dataset includes only 49 samples, so a more relaxed threshold is appropriate to account
# for the lower statistical power and to ensure detection of biologically relevant interactions.
# The MI threshold difference between p=1E-8 and p=1E-4 is minimal (~1e-5), but can significantly
# impact the number of inferred edges due to the skewed MI distribution.

java -Xmx50G -jar ARACNe-AP/dist/aracne.jar -e inmat_GEO.txt -o geo_network --tfs extra_data/TF_lambert.txt --pvalue 1E-4 --seed 1 --calculateThreshold

for i in {1..100}
do
java -Xmx200G -jar ARACNe-AP/dist/aracne.jar -e inmat_GEO.txt -o geo_network --tfs extra_data/TF_lambert.txt --pvalue 1E-4 --threads 50 --seed $i
done

java -Xmx50G -jar ARACNe-AP/dist/aracne.jar -o geo_network --consolidate

mv geo_network/network.txt ./geo_tumor_network.txt

# The inferred GEO tumor network (n=49) yielded ~76k interactions, compared to ~234k from
# TCGA (n~350). This is consistent with the expected loss of statistical power at lower sample
# sizes. The conservative p=1E-4 threshold appears justified, as the resulting network size
# is in proportion.
