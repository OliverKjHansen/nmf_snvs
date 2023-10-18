# snakemake --profile slurm
#sort -k1,1 -k2,2n
configfile: 'config.yaml'

###IMPORTENT###
#file_path = config["file_path"]
# all the 
chrom = config["chromosomes"]
length = config["chromosome_length"]
datasets = config["datasets"]

#genome_bed = config["genome_bed"]
#blacklist = config["blacklist"]
#topmed = config["topmed_data"] # when i want to start filtering 
#coverage_files = config["coverage"] # use this to make bedfiles, maybe later
window_sizes = config["window_size_kb"]
# kmer_indels = config["kmer_indels"]
NumberWithDepth = config["NumberWithDepth"]
# two_bit = config["twobitgenome"]
allelefrequency = config["allelefrequency"]
# methylation_data = config["methylation_data"]
# replication_time = config["replication_time"]
# size_partition = config["size_partition"]
# complex_structure = config["complex_structure"]
# recombination = config["recombination"]

#this makes the wildcards of the regions i want to investigate
def making_windows(chromosomes, length, window_sizes):
	regions = []
	for window_size in window_sizes:
		size = window_size*1000000 # for an indel analysis only windows with the size of 10^6 is doable
		for chrom in chromosomes:
			for pos in range(size, int(length[chrom]), size):
				title = chrom+"_"+str(int((pos-size)/1000000))+"mb_to_"+str(int(pos/1000000))+"mb"
				regions.append(title)
	return regions

regions = making_windows(chrom, length, window_sizes)

def creating_breakpoints(kmer):
	before = int(kmer)/2
	after = before-1
	return [before, after]

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
		#"files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed",
		"files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz",
		"files/{datasets}/derived_files/vcf_snvs/{chrom}_indel_{freq}.vcf.gz",
		#"files/{datasets}/derived_files/vcf_indels/all_indels_{freq}.vcf.gz"
		], datasets = datasets, chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency) #region = regions, window_sizes = window_sizes, kmer = kmer_indels,chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, size_partition = size_partition, complex_structure = complex_structure)

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
	python scripts/countingregions.py {unzipped} {params.procent} > {output.bedfile}
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

# rule aggregate_chromosomes:
# 	input:
# 		individual_coverage = expand("files/{datasets}/accepted_coverage/{chrom}x10_{fraction}p.bed", datasets=datasets, chrom=chrom, fraction=NumberWithDepth),
# 		individual_vcf =  expand("files/{datasets}/vcf_indels/{chrom}_indel_{freq}.vcf.gz", datasets=datasets, chrom=chrom, freq = allelefrequency)
# 	output:
# 		summary_coverage = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
# 		summary_vcf= expand("files/{datasets}/derived_files/vcf_indels/all_indels_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
# 	shell:"""
# 	cat {output.individual_coverage} > {output.summary_coverage}
# 	cat {output.individual_vcf} > {output.summary_vcf}
# 	"""