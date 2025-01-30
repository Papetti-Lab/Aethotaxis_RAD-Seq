### ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ###
### This file contains the command lines used in a Linux terminal to run the analysed presented in the paper "Sex differentiation and long-distance gene flow in the Antarctic fish Aethotaxis mitopteryx ###
### ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ###


		###########################
		### PROCESSING RAW DATA ###
		###########################

# A. mitopteryx samples were sequenced together with the samples used for the paper "Limited interspecific gene flow in the evolutionary history of the icefish Chionodraco spp."
# Demultiplexing and removal of PCR duplicates were the same as in that paper.

	### Demultiplexing and discarding low quality reads with Stacks ###

# Define input and output folders
IN=/path/to/raw/reads/folder
OUT=/path/to/output/the/processed/files/folder
# Command line for the first library
process_radtags -1 $IN/2791_S31_L002_R1_001.fastq.gz -2 $IN/2791_S31_L002_R2_001.fastq.gz -o $OUT -b $IN/barcodes_library_01.txt -e sbfI -r -c -q -t 143
# Command line for the second library
process_radtags -1 $IN/sample_S0_L006_R1_001.fastq.gz -2 $IN/sample_S0_L006_R2_001.fastq.gz -o $OUT -b $IN/barcodes_library_02.txt -e sbfI -r -c -q -t 143


	### Removal of PCR clones with Stacks ###

# Define input and output folders
IN=/path/to/cleaned/reads/generated/by/process_radtags/folder
OUT=/path/to/output/the/processed/files/folder

# Define an array with first input file names in a set of paired-end sequences
PAIR_1=($IN/*.1.fq.gz)
echo ${PAIR_1[@]} # check the file names
# Define an array with second input file names in a set of paired-end sequences
PAIR_2=($IN/*.2.fq.gz)
echo ${PAIR_2[@]} # check the file names

# Clone filter
for i in "${!PAIR_1[@]}" ; do
	clone_filter -P -1 ${PAIR_1[i]} -2 ${PAIR_2[i]} -i gzfastq -o $OUT > clone-filter-log.txt
 done





		#################################
		### ALIGN TO REFERENCE GENOME ###
		#################################

	### Scaffolding of the Aethotaxis genome using the D. mawsoni genome as a reference ###

# RagTag
REF=/path/to/Dissostichus_mawsoni/genome/file
QUERY=/path/to/Aethotaxis/genome/file
UNIMAP=/path/to/unimap/software
OUT=/path/to/output/folder

conda activate myenv # activate the python environment where ragtag is installed
ragtag.py scaffold $REF $QUERY -C -u -t 8 --aligner $UNIMAP -o $OUT


	### Alignment with BWA ###

	# Create genome index with BWA

# Define input/output folders
BWA=/path/to/bwa-0.7.17
GENOME=/path/to/genome/file
OUT=/path/to/output/folder
# Create BWA database
$BWA/bwa index -p $OUT/bwa_index/bwa $GENOME

	# Align with BWA mem

# Define input/output folders
BWA=/path/to/bwa-0.7.17
INDEX=/path/to/bwa_index/bwa
IN=/path/to/reads/to/be/aligned/folder
POPMAP=/path/to/popmap/file
OUT=/path/to/output/the/processed/files/folder

# Align
cat $POPMAP | cut -f 1 |
while read sample; do
	$BWA/bwa mem -t 8 -M $INDEX $IN/${sample}.1.1.fq.gz $IN/${sample}.2.2.fq.gz | samtools view -b -h | samtools sort --threads 8 -o $OUT/${sample}.bam
done




		###################
		### SNP CALLING ###
		###################

# Define input/output folders
IN=/path/to/bam/files/folder
POPMAP=/path/to/popmap/file
OUT=/path/to/output/the/processed/files/folder
	# Run gstacks
gstacks -I $IN -M $POPMAP -O $OUT --gt-alpha 0.01 --min-mapq 20 -t 8

# Define input/output folders
IN=/path/to/gstacks/folder
POPMAP=/path/to/popmap/file
OUT=/path/to/output/the/processed/files/folder
	# Run populations
populations -P $IN -M $POPMAP -O $OUT -r 0.75 --max-obs-het 0.75 --min-mac 2 --ordered-export --vcf




		#################
		### FILTERING ###
		#################

cd /path/to/folder/where/you/want/to/put/the/output
INPUT0=/path/to/stacks_populations/folder/populations.snps.vcf

# STEP 01: filter for minimum genotype quality and minimum depth:
# --minGQ 30
# --minDP 5
# --min-meanDP 10
# --mac 2
#
vcftools --vcf $INPUT0 --remove-indels --minGQ 30 --minDP 5 --min-meanDP 10 --mac 2 --max-missing 0.5 --recode --stdout --out step01 > step01.vcf
# Calculate per individual statistics
FILE=step01
vcftools --vcf $FILE.vcf --depth -c --out $FILE-depth_indv > $FILE-depth_indv.txt
vcftools --vcf $FILE.vcf --site-mean-depth -c --out $FILE-depth_site > $FILE-depth_site.txt
vcftools --vcf $FILE.vcf --missing-indv -c --out $FILE-miss_indv > $FILE-miss_indv.txt
vcftools --vcf $FILE.vcf --missing-site -c --out $FILE-miss_site > $FILE-miss_site.txt

# ---------------------- #
# Move to R for plotting #
# ---------------------- #

# ------------------------------------- #
# Back to Linux to filter for max depth #
# ------------------------------------- #
# STEP 02: filter for maximum depth
vcftools --vcf step01.vcf --maxDP 40 --max-meanDP 40 --mac 2 --max-missing 0.5 --recode --stdout --out step02 > step02.vcf

# --------------------------------------------- #
# Move to R for further filtering with SNPfiltR #
# --------------------------------------------- #
# STEP 03: filter for missing data
# STEP 04: thin for linkage



# -> From here on "step03" and "step04" indicate the names of the filtered SNP dataset used for the different analyses <- #



# Uncompress vcf.gz files produces by SNPfiltR
bgzip -d step04.vcf.gz



		#################
		### Admixture ###
		#################

# Convert from vcf to plink
cd /path/to/folder/where/you/want/to/put/the/output
IN=/path/to/step04.vcf # path to the input file
./plink --vcf $IN --const-fid --allow-extra-chr --make-bed --out step04
# Prepare the input for Admixture
FILE=step04 # name of the file without the file extension, not the file itself
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

# Run Admixture
ADM=/path/to/admixture/folder
for K in {1..5} # adjust the number of interation of the loop according to the number of clusters you want to test
do
	$ADM/admixture32 --cv $FILE.bed $K | tee log${K}.out
done

# Visualize cross-validation (CV) errors
grep -h CV log*.out
# Save cross-validation errors in a convenient format to plot them in R
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' | sort -k 1 -n > $FILE-cv_errors.txt




		########################
		### fineRADstructure ###
		########################

# Create a whitelist of loci with <= 10 SNPs based on the populations.sumstats.tsv file:
cd /path/to/folder/where/you/want/to/put/the/output
cat /path/to/stacks_populations/folder/populations.sumstats.tsv | grep -v "^#" | cut -f 1 | sort -n | uniq -c | awk '$1<=10 {print $2}' > whitelist10.txt

# Run populations to create the input file
IN=/path/to/gstacks/folder
POPMAP=/path/to/popmap/file
OUT=/path/to/folder/where/you/want/to/put/the/output
populations -P $IN -M $POPMAP -O $OUT -W whitelist10.txt -r 0.80 --max-obs-het 0.75 --min-mac 2 --filter-haplotype-wise --ordered-export --radpainter
# Rename the file with a better name for further processing
mv populations.haps.radpainter aetho.radpainter

# Calculate the co-ancestry matrix
FILE=aetho
FINERADSTR=/path/to/fineRADstructure
$FINERADSTR/RADpainter paint $FILE.radpainter
# Assign individuals to populations
$FINERADSTR/finestructure -x 100000 -y 100000 -z 1000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml
# Tree building
$FINERADSTR/finestructure -m T -x 100000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml ${FILE}_chunks.mcmcTree.xml

# ---------------------- #
# Move to R for plotting #
# ---------------------- #



		#############################
		### Divergence Statistics ###
		#############################

# Between sexes
IN=/path/to/gstacks/folder
POPMAP=/path/to/popmap-sex/file # popmap where individuals are labelled by sex
WHITELIST=/path/to/whitelist10.txt
OUT=/path/to/folder/where/you/want/to/put/the/output
populations -P $IN -M $POPMAP -O $OUT -W $WHITELIST -p 2 -r 0.80 --max-obs-het 0.75 --min-mac 2 --filter-haplotype-wise --fstats --smooth -t 8
# Between populations
IN=/path/to/gstacks/folder
POPMAP=/path/to/popmap-geo/file # popmap where individuals are labelled by geographical origin
WHITELIST=/path/to/whitelist10.txt
OUT=/path/to/folder/where/you/want/to/put/the/output
populations -P $IN -M $POPMAP -O $OUT -W $WHITELIST -p 2 -r 0.80 --max-obs-het 0.75 --min-mac 2 --filter-haplotype-wise --fstats --smooth -t 8

# Prepare file for plotting in R
sed '1,2d' populations.phistats_female-male.tsv > phi-plot.tsv
sed '1,2d' populations.phistats_AP-WS.tsv > phi-plot.tsv

# ---------------------- #
# Move to R for plotting #
# ---------------------- #




		#######################################
		### fineRADstructure without Scf 06 ###
		#######################################

# Create a whitelist of loci with <= 10 SNPs and without Scf 06 based on the populations.sumstats.tsv file:
cd /path/to/folder/where/you/want/to/put/the/output
cat /path/to/stacks_populations/folder/populations.sumstats.tsv | grep -v "^#" | awk '$2 != "HiC_scaffold_6_RagTag" ' | cut -f 1 | sort -n | uniq -c | awk '$1<=10 {print $2}' > whitelist10_No06.txt

# Run populations to create the input file
IN=/path/to/gstacks/folder
POPMAP=/path/to/popmap/file
OUT=/path/to/folder/where/you/want/to/put/the/output
populations -P $IN -M $POPMAP -O $OUT -W whitelist10_No06.txt -r 0.80 --max-obs-het 0.75 --min-mac 2 --filter-haplotype-wise --ordered-export --radpainter
# Rename the file with a better name for further processing
mv populations.haps.radpainter aetho_No06.radpainter

# Calculate the co-ancestry matrix
FILE=aetho_No06
FINERADSTR=/path/to/fineRADstructure
$FINERADSTR/RADpainter paint $FILE.radpainter
# Assign individuals to populations
$FINERADSTR/finestructure -x 100000 -y 100000 -z 1000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml
# Tree building
$FINERADSTR/finestructure -m T -x 100000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml ${FILE}_chunks.mcmcTree.xml

# ---------------------- #
# Move to R for plotting #
# ---------------------- #



		################################
		### Admixture without Scf 06 ###
		################################

cd /path/to/folder/where/you/want/to/put/the/output
IN=/path/to/step04.vcf
# Remove Scf 06
vcftools --vcf $IN --not-chr HiC_scaffold_6_RagTag --recode --stdout --out step04_No06 > step04_No06.vcf
# Convert from vcf to plink
IN=/path/to/step04_No06.vcf
./plink --vcf $IN --const-fid --allow-extra-chr --make-bed --out step04_No06
# Prepare the input for Admixture
FILE=step04_No06 # name of the file without the file extension, not the file itself
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

# Run Admixture
ADM=/path/to/admixture/folder
for K in {1..5} # adjust the number of interation of the loop according to the number of clusters you want to test
do
	$ADM/admixture32 --cv $FILE.bed $K | tee log${K}.out
done

# Visualize cross-validation (CV) errors
grep -h CV log*.out
# Save cross-validation errors in a convenient format to plot them in R
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' | sort -k 1 -n > $FILE-cv_errors.txt

# ---------------------- #
# Move to R for plotting #
# ---------------------- #



		####################
		### Homozigosity ###
		####################

	# Compute Homozigosity per sample per chromosome

IN=/path/to/folder/containing/input/file
FILE=step04 # name of the file without the file extension, not the file itself
# Define a variable containing the chromosome names
SCF="HiC_scaffold_1_RagTag HiC_scaffold_3_RagTag HiC_scaffold_4_RagTag HiC_scaffold_5_RagTag 
     HiC_scaffold_6_RagTag HiC_scaffold_8_RagTag HiC_scaffold_10_RagTag HiC_scaffold_11_RagTag 
     HiC_scaffold_13_RagTag HiC_scaffold_21_RagTag HiC_scaffold_22_RagTag HiC_scaffold_23_RagTag 
     HiC_scaffold_24_RagTag HiC_scaffold_25_RagTag HiC_scaffold_26_RagTag HiC_scaffold_28_RagTag 
     HiC_scaffold_30_RagTag HiC_scaffold_31_RagTag HiC_scaffold_32_RagTag HiC_scaffold_33_RagTag 
     HiC_scaffold_34_RagTag HiC_scaffold_35_RagTag HiC_scaffold_38_RagTag HiC_scaffold_41_RagTag"

# Compute homozigosity with VCFtools, put the results from different chromosomes into different folders, combine all the results into a single file
for scaffold in $SCF
do
	mkdir ${scaffold}
	cd ${scaffold}
	vcftools --vcf $IN/$FILE.vcf --chr ${scaffold} --het -c --out $FILE-${scaffold}-homozigosity > $FILE-${scaffold}-homozigosity.txt
	sed '1d' $FILE-${scaffold}-homozigosity.txt | awk -v value=${scaffold} '{print value, $0}'>> ../$FILE-homozigosity.txt # see explanation below
	cd ../
done
# Explanation of the command:
# sed... removes the first line from the file
# awk... add the scaffold name to the beginning of each line (in the end it creates a new first column with scaffold name)
# >>... append the modified files to ../$FILE-homozigosity.txt

# Add the header to step04-homozigosity.txt
echo -e "SCF\tINDV\tO_HOM\tE_HOM\tN_SITES\tF" | cat - step04-homozigosity.txt > headers && mv headers step04-homozigosity.txt
# Explanation: define the header, concatenate the header with the content of step04-homozigosity.txt into a temporary file called header, rename (move) headers back to step04-homozigosity.txt

# ---------------------- #
# Move to R for plotting #
# ---------------------- #



		#############
		### Depth ###
		#############

IN=/path/to/folder/containing/input/file
FILE=step04
# Select from the whole dataset only male individuals
# males.txt is a file where only male individuals are listed (one per row)
vcftools --vcf $IN/$FILE.vcf --keep males.txt --recode --stdout > $FILE-males.vcf
# Compute depth
vcftools --vcf $IN/$FILE-males.vcf --site-mean-depth --stdout --out males-depth_site
# Repeat for females
# males.txt is a file where only female individuals are listed (one per row)
vcftools --vcf $IN/$FILE.vcf --keep females.txt --recode --stdout > $FILE-females.vcf
vcftools --vcf $IN/$FILE-females.vcf --site-mean-depth --stdout --out females-depth_site

# ---------------------- #
# Move to R for plotting #
# ---------------------- #



		##############################
		### Linkage disequilibrium ###
		##############################

IN=/path/to/folder/containing/input/file
# All individuals
FILE=step03
# Convert from VCF to plink file format
./plink --vcf $IN/$FILE.vcf --const-fid --allow-extra-chr --make-bed --out $FILE
# Compute linkage decay
./plink --bfile $FILE --const-fid --allow-extra-chr --set-missing-var-ids @:# --r2 --ld-window-r2 0.0 --out $FILE

# Repeat for males (select from the whole dataset only male individuals, convert from vcf to plink and compute linkage decay)
# males.txt is a file where only male individuals are listed (one per row)
vcftools --vcf $IN/$FILE.vcf --keep males.txt --recode --stdout > $FILE-males.vcf
MALES=step03-males
./plink --vcf $MALES.vcf --const-fid --allow-extra-chr --make-bed --out $MALES
./plink --bfile $MALES --const-fid --allow-extra-chr --set-missing-var-ids @:# --r2 --ld-window-r2 0.0 --out $MALES

# Repeat for females  (select from the whole dataset only female individuals, convert from vcf to plink and compute linkage decay)
# females.txt is a file where only female individuals are listed (one per row)
vcftools --vcf $IN/$FILE.vcf --keep females.txt --recode --stdout > $FILE-females.vcf
FEMALES=step03-females
./plink --vcf $FEMALES.vcf --const-fid --allow-extra-chr --make-bed --out $FEMALES
./plink --bfile $FEMALES --const-fid --allow-extra-chr --set-missing-var-ids @:# --r2 --ld-window-r2 0.0 --out $FEMALES

# ---------------------- #
# Move to R for plotting #
# ---------------------- #




		##############
		### RADSEX ###
		##############

RADSEX=/path/to/RADsex/bin
IN=/path/to/clone_filtered_reads/
OUT=/path/to/folder/where/you/want/to/put/the/markers_table.tsv
cd /path/to/output/folder

# popmap-sex.txt is a popmap where individuals are labelled by sex

# Generating a table of marker depths for the entire dataset (paired-end reads are not supported, use only the first read)
$RADSEX/radsex process --input-dir $IN --output-file $OUT --threads 8

# Computes the distribution of marker depths in each individual
$RADSEX/radsex depth --markers-table markers_table.tsv --popmap popmap-sex.txt --output-file depth.txt

# Computing the distribution of markers between sexes (repeat with min-depth values 1, 2, 5 and 10)
for depth in {1,2,5,10}; do $RADSEX/radsex distrib --markers-table markers_table.tsv --output-file distrib_${depth}.tsv --popmap popmap-sex.txt --min-depth ${depth} --groups male,female; done

# Finding markers significantly associated with sex
for depth in {1,2,5,10}; do $RADSEX/radsex signif --markers-table markers_table.tsv --output-file signif_${depth}.tsv --popmap popmap-sex.txt --min-depth ${depth} --groups male,female; done

# Aligning markers to a genome
GENOME=/mnt/d/aetho/genome-aetho/unimapped/ragtag.scaffold.fasta
for depth in {1,2,5,10}; do $RADSEX/radsex map --markers-file markers_table.tsv --output-file map_${depth}.tsv --popmap popmap-sex.txt --genome-file $GENOME --min-depth ${depth} --groups male,female; done

# ---------------------- #
# Move to R for plotting #
# ---------------------- #
