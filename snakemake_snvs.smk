import pickle
import os

configfile: 'config.yaml'

###IMPORTENT###
# all the 
chrom = config["chromosomes"]
length = config["chromosome_length"]
datasets = config["datasets"]
genome_bed = config["genome_bed"]
blacklist = config["blacklist"]
two_bit = config["twobitgenome"]
exons = config["exon"]
window_sizes = config["window_size_kb"]
NumberWithDepth = config["NumberWithDepth"]
allelefrequency = config["allelefrequency"]
types = config["pattern_type"]
kmer = config["kmer_snvs"]
signatures = config["NumberOfSignatures"]
# methylation_data = config["methylation_data"]
# replication_time = config["replication_time"]
# size_partition = config["size_partition"]
# complex_structure = config["complex_structure"]
# recombination = config["recombination"]

#this makes the wildcards of the regions i want to investigate
def making_windows(chromosomes, length, window_sizes):
	regions = []
	for window_size in window_sizes:
		size = window_size*1000 # in snv it should already be at 10.000
		for chrom in chromosomes:
			for pos in range(size, int(length[chrom]), size):
				title = chrom+"_"+str(int((pos-size)/1000))+"kb_to_"+str(int(pos/1000))+"kb"
				regions.append(title)
	return regions

#regions = making_windows(chrom, length, window_sizes)


# Define the filename for your list data

# Check if the file exists
if os.path.exists("regions.pkl"):
    # Load the list from the file
    with open("regions.pkl", "rb") as file:
        regions = pickle.load(file)
else:
	regions = making_windows(chrom, length, window_sizes)
	with open("regions.pkl", "wb") as file:
		pickle.dump(regions, file)

def creating_radius(kmer):
	radius = (int(kmer)-1)/2
	return [radius]

# def making_windows2(chromosome, length, window_size):
# 	regions = []
# 	size = window_size*1000000
# 	for pos in range(size, int(length[chrom]), size):
# 		title = str(int((pos-size)/1000000))+"mb_to_"+str(int(pos/1000000))+"mb"
# 			regions.append(title)
# 	return regions
rule all:
	input:
		expand(["files/{datasets}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz",
		"files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", # This file contains a bedfile of all the regions that passes the restriction i have put on 80% of the individuals needs to have a coverage of 10x
		"files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed",
		"files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz",
		"files/{datasets}/derived_files/vcf_snvs/{chrom}_indel_{freq}.vcf.gz",
		"files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz"
		"{window_sizes}kb_windows/regions/{region}.bed",
		"{window_sizes}kb_windows/filtered_regions/{region}_{fraction}p.bed",
		"{window_sizes}kb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed",
		"{window_sizes}kb_windows/variants/snv_{region}_{freq}_{fraction}p.bed",
		"{window_sizes}kb_windows/snv_{kmer}mer/frequency_{freq}_at_{fraction}p/snv_counts_{region}_{kmer}mer.bed",
		"{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{kmer}mer.bed",
		"{window_sizes}kb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed",
		"{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_dataframe_{kmer}mer.rds",
		"{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
		], datasets = datasets, chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, types = types, kmer = kmer, signatures = signatures) #region = regions, window_sizes = window_sizes, kmer = kmer_indels,chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, size_partition = size_partition, complex_structure = complex_structure)

rule coverage_regions:
	input:
		seq_zipped = "files/{datasets}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz"
	params: 
		procent = lambda wildcards: float(int(wildcards.fraction)/100)
	output:
		bedfile = "files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", # This file contains a bedfile of all the regions that passes the restriction i have put on 100% of the individuals needs to have a coverage of 10x
	shell:"""
	temp_unzipped=$(mktemp -u)
    gunzip -c {input.seq_zipped} > $temp_unzipped
	python scripts/countingregions.py $temp_unzipped {params.procent} > {output.bedfile}
	gzip $temp_unzipped
	"""
rule vcf_snvs:
	input:
		raw_vcf = "files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz"
	conda: "envs/bcftools.yaml"
	output:
		filtered = "files/{datasets}/derived_files/vcf_snvs/{chrom}_indel_{freq}.vcf.gz",
	shell:"""
	tabix -f -p vcf {input.raw_vcf}
	bcftools filter -O z -o {output.filtered} -i 'AF<{wildcards.freq} && VRT=1' {input.raw_vcf}
	"""
rule aggregate_chromosomes:
	input:
		individual_coverage = expand("files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", datasets=datasets, chrom=chrom, fraction=NumberWithDepth),
		individual_vcf =  expand("files/{datasets}/derived_files/vcf_snvs/{chrom}_indel_{freq}.vcf.gz", datasets=datasets, chrom=chrom, freq = allelefrequency)
	output:
		summary_coverage = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
		summary_vcf= expand("files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
	shell:"""
	cat {input.individual_coverage} >> {output.summary_coverage}
	cat {input.individual_vcf} >> {output.summary_vcf}
	"""

# a rule which makes the MegaBases bedfile 
# Creating regions which are to be investigated

#might be faster in a for loop
if not os.path.exists("10kb_windows/regions/chr3_50_60kb.bed"):
	for region in regions:
		chromos = region.split("_")[0].split("k")[0]
		start = str(int(region.split("_")[1].split("k")[0])*1000)
		end = str(int(region.split("_")[3].split("k")[0])*1000)

		filename = f"10kb_windows/regions/{region}.bed"
		with open(filename, "w") as file:
			file.write(f"{chromos}	{start}	{end}\n")

# rule mega_bases:
# 	params: 
# 		chrom = lambda wildcards: wildcards.region.split("_")[0].split("k")[0],
# 		start = lambda wildcards: str(int(wildcards.region.split("_")[1].split("k")[0])*1000),
# 		end = lambda wildcards: str(int(wildcards.region.split("_")[3].split("k")[0])*1000)
# 	resources:
# 		threads=1,
# 		time=10,
# 		mem_mb=1000
# 	output:
# 		bedfiles = "{window_sizes}kb_windows/regions/{region}.bed"
# 	shell:"""
# 	printf '%s\t%s\t%s\n' {params.chrom} {params.start} {params.end} > {output.bedfiles}
# 	"""

rule filtering_regions:
	input:
		regions = "{window_sizes}kb_windows/regions/{region}.bed",
		coverage_accepted = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
		blacklist = blacklist,
		exons = exons
	conda: "envs/bedtools.yaml"
	resources:
		threads=1,
		time=120,
		mem_mb=5000
	output:
		tmp_cov = temporary("{window_sizes}kb_windows/tmp/tmp_coverage_{region}_{fraction}p.bed"),
		tmp_blacklist = temporary("{window_sizes}kb_windows/tmp/blacklist_{region}_{fraction}p.bed"),
		tmp_exons = temporary("{window_sizes}kb_windows/tmp/exons_{region}_{fraction}p.bed"),
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{region}_{fraction}p.bed"
	shell:"""
	bedtools intersect -a {input.regions} -b {input.coverage_accepted} > {output.tmp_cov}
	bedtools intersect -v -a {output.tmp_cov} -b {input.blacklist} > {output.tmp_blacklist}
	bedtools intersect -v -a {output.tmp_blacklist} -b {input.exons} > {output.tmp_exons}
	tmp=`bedtools intersect -wo -a {input.regions} -b {output.tmp_exons}| awk '{{s+=$7}} END {{print s}}'`
	num=$(expr {window_sizes} \* 10000 / 2)
	if [[ $tmp -ge $num ]]
	then 
		cp {output.tmp_exons} {output.filtered_regions}
	else
		touch {output.filtered_regions}
	fi
	"""

rule background_counter: #im not sure this works
	input:
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{region}_{fraction}p.bed",
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	params:
		radius  = lambda wildcards: int(creating_radius(wildcards.kmer)[0]),
	output:
		background = temporary("{window_sizes}kb_windows/tmp/background_{region}_{kmer}mer_{fraction}p.bed"),
		ss_background = "{window_sizes}kb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed"
	shell:"""
	check=`cat {input.filtered_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then 
		kmer_counter background --bed {input.filtered_regions} --radius {params.radius} {input.genome} > {output.background}
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.background} > {output.ss_background}
	else
		touch {output.background}
		touch {output.ss_background}
	fi
	"""

rule creating_variants:
	input:
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{region}_{fraction}p.bed",
		vcf_file = expand("files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
	conda: "envs/bedtools.yaml"
	output:
		variants =  "{window_sizes}kb_windows/variants/snv_{region}_{freq}_{fraction}p.bed"
	shell:"""
	check=`cat {input.filtered_regions} | wc -l`
	if [[ $check -gt 0 ]]
	then
		bedtools intersect -a {input.vcf_file} -b {input.filtered_regions} | awk -v OFS='\t' '{{print $1,$2,$4,$5}}' > {output.variants}
	else
		touch {output.variants}
	fi
	"""

rule snv_variant_counter:
	input:
		variants = "{window_sizes}kb_windows/variants/snv_{region}_{freq}_{fraction}p.bed",
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	params:
		radius  = lambda wildcards: int(int(wildcards.kmer)/2)
	output:
		kmer_count_snv = temporary("{window_sizes}kb_windows/tmp/snv_{kmer}mer/frequency_{freq}_at_{fraction}p/counts_{region}_{kmer}mer.bed"),
		ss_snv = "{window_sizes}kb_windows/snv_{kmer}mer/frequency_{freq}_at_{fraction}p/snv_counts_{region}_{kmer}mer.bed",
	shell:"""
	check=`cat {input.variants} | wc -l`
	if [[ $check -gt 0 ]]
	then
		kmer_counter snv {input.genome} {input.variants} > {output.kmer_count_snv}
		awk -v OFS='\t' '{{print "{wildcards.region}",$1,$2}}' {output.kmer_count_snv} > {output.ss_snv} 
	else
		touch {output.kmer_count_snv}
		touch {output.ss_snv} 
	fi
	"""
##Make a check for the directories
rule aggregate_regions:
	input:
		snvs = expand("{window_sizes}kb_windows/snv_{{kmer}}mer/frequency_{freq}_at_{fraction}p/snv_counts_{region}_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer),
		background = expand("{window_sizes}kb_windows/background_{{kmer}}mer/background_{region}_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer)
	output:
		summary_snv = expand("{window_sizes}kb_windows/snv_{{kmer}}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer),
		summary_background = expand("{window_sizes}kb_windows/background_{{kmer}}mer/combined/background_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer)
	params:
	shell:"""
	cat {input.snv} >> {output.summary_snv}
	cat {input.background} >> {output.summary_background}
	"""
###Now let do some nmf###

rule prepare_for_nmf:
	input:
		summary_snv = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{kmer}mer.bed",
		summary_background = "{window_sizes}kb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed"

	conda: "envs/callr.yaml"
	output:
		snv_dataframe = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_dataframe_{kmer}mer.rds",
	shell:"""
	Rscript scripts/creating_dataframes.R {input.summary_background} {input.summary_snv} {output.snv_dataframe}
	"""
rule modelselection:
	input:
		count_data = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/{types}_dataframe_{kmer}mer.rds"
	conda: "envs/nmf.yaml"
	resources:
		threads=4,
		time=800,
		mem_mb=10000
	output:
		model = "{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
	shell:"""
    Rscript scripts/opportunity_modelselection.R {wildcards.signatures} {input.count_data} {output.model}
    """