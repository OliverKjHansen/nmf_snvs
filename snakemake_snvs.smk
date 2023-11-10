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
    chromo = {}
    chromo_splits = {}
    chr_split_length = 1000000 # change this value to control the size of the sub folders
    for window_size in window_sizes:
        size = window_size * 1000  # in snv it should already be at 10.000
        for chrom in chromosomes:
            chromo[chrom] = {}
            for split in range(chr_split_length, int(length[chrom]), chr_split_length):
                regions = []
                for pos in range((split-chr_split_length), split, size):
                    title = chrom+"_"+str(int((pos) / 1000)) + "kb_to_" + str(int((pos+size)/ 1000)) + "kb"
                    regions.append(title)
                chromo[chrom][f"{chrom}_{int((split-chr_split_length)/1000)}kb_{int(split/1000)}kb"] = regions
    return chromo

#regions_and_splits = making_windows(chrom, length, window_sizes)

# Define the filename for your list data

# Check if the file exists
if os.path.exists("regions_and_splits.pkl"):
    # Load the list from the file
    with open("regions_and_splits.pkl", "rb") as file:
        regions_and_splits = pickle.load(file)
else:
	regions_and_splits = making_windows(chrom, length, window_sizes)
	with open("regions_and_splits.pkl", "wb") as file:
		pickle.dump(regions_and_splits, file)

def creating_radius(kmer):
	radius = (int(kmer)-1)/2
	return [radius]

if os.path.exists("splits_list.pkl"):
    # Load the list from the file
    with open("splits_list.pkl", "rb") as file:
        splits_list = pickle.load(file)
else:
    splits_list = []
    for chromos in regions_and_splits:
        if not os.path.exists(f"10kb_windows/regions/{chromos}"):
            os.makedirs(f"10kb_windows/regions/{chromos}")
        splits = regions_and_splits[chromos]
        for split in splits:
            if not os.path.exists(f"10kb_windows/regions/{chromos}/{split}/"):
                os.makedirs(f"10kb_windows/regions/{chromos}/{split}/")
            splits_list.append(f"{chromos}/{split}")
            regions = regions_and_splits[chromos][split]
            for region in regions:
                #print(region)
                chromosome = region.split("_")[0].split("k")[0]
                start = str(int(region.split("_")[1].split("k")[0]) * 1000)
                end = str(int(region.split("_")[3].split("k")[0]) * 1000)
                filename = f"10kb_windows/regions/{chromos}/{split}/{region}_90p.bed"
                with open(filename, "w") as file:
                    file.write(f"{chromosome}\t{start}\t{end}\n")
    with open("splits_list.pkl", "wb") as file:
        pickle.dump(splits_list, file)

rule all:
	input:
		expand(["files/{datasets}/coverage_files/{chrom}.BRAVO_TOPMed_coverage_hg38.txt.gz",
		"files/{datasets}/derived_files/accepted_coverage/{chrom}x10_{fraction}p.bed", # This file contains a bedfile of all the regions that passes the restriction i have put on 80% of the individuals needs to have a coverage of 10x
		"files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed",
		"files/{datasets}/vcf_files/{chrom}.BRAVO_TOPMed_Freeze_8.vcf.gz",
		"files/{datasets}/derived_files/vcf_snvs/{chrom}_indel_{freq}.vcf.gz",
		"files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz",
		"{window_sizes}kb_windows/regions/{splits_list}/",
		"{window_sizes}kb_windows/filtered_regions/{splits_list}/dummyfile_{fraction}p.bed",
		"{window_sizes}kb_windows/background_{kmer}mer_{fraction}p/{splits_list}/dummyfile_background.bed",
		#"{window_sizes}kb_windows/variants_{freq}_{fraction}p/{splits_list}/dummy_snv.bed",
		"{window_sizes}kb_windows/snv_{kmer}mer_freq_{freq}_at_{fraction}p/{splits_list}/dummy_{kmer}mer.bed"
		# "{window_sizes}kb_windows/background_{kmer}mer/background_{region}_{kmer}mer_{fraction}p.bed",
		# "{window_sizes}kb_windows/variants/snv_{region}_{freq}_{fraction}p.bed",
		# "{window_sizes}kb_windows/snv_{kmer}mer/frequency_{freq}_at_{fraction}p/snv_counts_{region}_{kmer}mer.bed",
		# "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{kmer}mer.bed",
		# "{window_sizes}kb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed",
		# "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_dataframe_{kmer}mer.rds",
		# "{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
		], datasets = datasets, chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, types = types, kmer = kmer, signatures = signatures, splits_list = splits_list) #region = regions, window_sizes = window_sizes, kmer = kmer_indels,chrom = chrom, fraction = NumberWithDepth, freq = allelefrequency, size_partition = size_partition, complex_structure = complex_structure)

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
# 		bedfiles = "{window_sizes}kb_windows/regions/{chrom}/{region}.bed"
# 	shell:"""
# 	printf '%s\t%s\t%s\n' {params.chrom} {params.start} {params.end} > {output.bedfiles}
# 	"""

rule filtering_regions:
	input:
 		regions = "{window_sizes}kb_windows/regions/{splits_list}/",
 		coverage_accepted = expand("files/{datasets}/derived_files/accepted_coverage/all_coverage_x10_{fraction}p.bed", datasets=datasets, fraction=NumberWithDepth),
 		blacklist = blacklist,
		exons = exons
	conda: "envs/bedtools.yaml"
	params: 
	resources:
		threads=2,
		time=180,
		mem_mb=2500
	output:
		# tmp_cov = temporary("{window_sizes}kb_windows/tmp/tmp_coverage_{region}_{fraction}p.bed"),
		# tmp_blacklist = temporary("{window_sizes}kb_windows/tmp/blacklist_{region}_{fraction}p.bed"),
		# tmp_exons = temporary("{window_sizes}kb_windows/tmp/exons_{region}_{fraction}p.bed"),
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{splits_list}/dummyfile_{fraction}p.bed"
	shell:"""

	folder={input.regions}
	output_dir=$(dirname {output.filtered_regions})
	echo $output_dir

	for file in "$folder"/*
	do
		if [ -f "$file" ]
		then
			file_name=$(basename "$file")

			tmp_filter=$(mktemp)

			bedtools intersect -a "$file" -b {input.coverage_accepted} | \
    		bedtools intersect -v -a - -b {input.blacklist} | \
    		bedtools intersect -v -a - -b {input.exons} > "$tmp_filter"

			tmp=`bedtools intersect -wo -a "$file" -b "$tmp_filter" | awk '{{s+=$7}} END {{print s}}'`
			echo "$tmp"

			num=$(expr {window_sizes} \* 1000 / 2)
			echo "$num"

			output_file="$output_dir/$file_name"
			echo "$output_file"

			if [[ $tmp -ge $num ]]
			then 
				cp "$tmp_filer" "$output_file"
			else
				touch "$output_file"
			fi
			rm -f "$tmp_filter""
		fi
	done
	touch {output.filtered_regions}
	"""

# bedtools intersect -a "$file" -b {input.coverage_accepted} > "$tmp_cov"
# bedtools intersect -v -a "$tmp_cov" -b {input.blacklist} > "$tmp_blacklist"
# bedtools intersect -v -a "$tmp_blacklist" -b {input.exons} > "$tmp_exon"
# tmp=`bedtools intersect -wo -a "$file" -b "$tmp_exon"| awk '{{s+=$7}} END {{print s}}'`

rule background_counter: #im not sure this works tmp_bck=$(mktemp)
	input:
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{splits_list}/dummyfile_{fraction}p.bed",
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	params:
		radius  = lambda wildcards: int(creating_radius(wildcards.kmer)[0])
	resources:
		threads=2,
		time=240,
		mem_mb=3000
	output:
# 		background = temporary("{window_sizes}kb_windows/tmp/background_{region}_{kmer}mer_{fraction}p.bed"),
		ss_background = "{window_sizes}kb_windows/background_{kmer}mer_{fraction}p/{splits_list}.tsv",
		background_dummy "{window_sizes}kb_windows/background_{kmer}mer_{fraction}p/dummy_{splits_list}.txt",
	shell:"""
	folder=$(dirname {input.filtered_regions})
	mkdir -p $(dirname {ss_background})
	touch {ss_background}

	for file in "$folder"/*
	do
		if [ -f "$file" ]
		then
			file_name=$(basename "$file")
			check=`cat "$file" | wc -l`
			if [[ "$check" -gt 0 ]]
			then 
				kmer_counter background --bed "$file" --radius {params.radius} {input.genome} | \
				awk -v OFS='\t' -v file_name="$file_name" '{{print file_name,$1,$2}}' - >> {ss_background}
			fi
		fi
	done
	touch {output.background_dummy}
	"""

# rule creating_variants:
# 	input:
# 		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{splits_list}/dummyfile_{fraction}p.bed",
# 		vcf_file = expand("files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency)
# 	conda: "envs/bedtools.yaml"
# 	resources:
# 		threads=4,
# 		time=420,
# 		mem_mb=5000
# 	output:
# 		variants =  "{window_sizes}kb_windows/variants_{freq}_{fraction}p/{splits_list}/dummy_snv.bed"
# 	shell:"""
# 	folder=$(dirname {input.filtered_regions})
# 	output_snakemake={output.variants}
# 	output_dir=$(dirname "$output_snakemake")

# 	for file in "$folder"/*
# 	do
# 		if [ -f "$file" ]
# 		then
# 			file_name=$(basename "$file")
# 			output_file="$output_dir/$file_name"
# 			check=`cat "$file" | wc -l`
# 			if [[ "$check" -gt 0 ]]
# 			then 
# 				bedtools intersect -a {input.vcf_file} -b "$file" | awk -v OFS='\t' '{{print $1,$2,$4,$5}}' > "$output_file"
# 			else
# 				touch "$output_file"
# 			fi
# 		fi
# 	done
# 	touch {output.variants}
# 	"""

rule snv_variant_counter:
	input:
		filtered_regions = "{window_sizes}kb_windows/filtered_regions/{splits_list}/dummyfile_{fraction}p.bed",
		#variants = "{window_sizes}kb_windows/variants_{freq}_{fraction}p/{splits_list}/dummy_snv.bed",
		vcf_file = expand("files/{datasets}/derived_files/vcf_snvs/all_snvs_{freq}.vcf.gz", datasets=datasets, freq = allelefrequency),
		genome = two_bit
	conda: "envs/kmer_counter.yaml"
	resources:
		threads=4,
		time=420,
		mem_mb=5000
	params:
		radius  = lambda wildcards: int(creating_radius(wildcards.kmer)[0])
	output:
		#kmer_count_snv = temporary("{window_sizes}kb_windows/tmp/snv_{kmer}mer/frequency_{freq}_at_{fraction}p/counts_{region}_{kmer}mer.bed"),
		ss_snv = "{window_sizes}kb_windows/snv_{kmer}mer_freq_{freq}_at_{fraction}p/{splits_list}.tsv"
		dummy_snv = "{window_sizes}kb_windows/snv_{kmer}mer_freq_{freq}_at_{fraction}p/{splits_list}/dummy_{kmer}mer.txt"
	shell:"""
	folder=$(dirname {input.filtered_regions})
	mkdir -p $(dirname {output.ss_snv})
	touch {output.ss_snv}

	for file in "$folder"/*
	do
		if [ -f "$file" ]
		then
			file_name=$(basename "$file")
			check=`cat "$file" | wc -l`
			if [[ "$check" -gt 0 ]]
			then 
				bedtools intersect -a {input.vcf_file} -b "$file" | \
				awk -v OFS='\t' '{{print $1, $2, $4, $5}}' | \
				kmer_counter snv {input.genome} - | \
				awk -v OFS='\t' -v file_name="$file_name" '{{print file_name, $1, $2}}' >> {output.ss_snv}
			fi
		fi
	done
	touch {output.ss_snv}
	"""
# ##Make a check for the directories
# rule aggregate_regions:
# 	input:
# 		snvs = expand("{window_sizes}kb_windows/snv_{{kmer}}mer/frequency_{freq}_at_{fraction}p/snv_counts_{region}_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer),
# 		background = expand("{window_sizes}kb_windows/background_{{kmer}}mer/background_{region}_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, region = regions, window_sizes = window_sizes, kmer = kmer)
# 	output:
# 		summary_snv = expand("{window_sizes}kb_windows/snv_{{kmer}}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{{kmer}}mer.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer),
# 		summary_background = expand("{window_sizes}kb_windows/background_{{kmer}}mer/combined/background_{{kmer}}mer_{fraction}p.bed", fraction = NumberWithDepth, freq = allelefrequency, window_sizes = window_sizes, kmer = kmer)
# 	params:
# 	shell:"""
# 	cat {input.snv} >> {output.summary_snv}
# 	cat {input.background} >> {output.summary_background}
# 	"""
# ###Now let do some nmf###

# rule prepare_for_nmf:
# 	input:
# 		summary_snv = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_counts_{kmer}mer.bed",
# 		summary_background = "{window_sizes}kb_windows/background_{kmer}mer/combined/background_{kmer}mer_{fraction}p.bed"

# 	conda: "envs/callr.yaml"
# 	output:
# 		snv_dataframe = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/snv_dataframe_{kmer}mer.rds",
# 	shell:"""
# 	Rscript scripts/creating_dataframes.R {input.summary_background} {input.summary_snv} {output.snv_dataframe}
# 	"""
# rule modelselection:
# 	input:
# 		count_data = "{window_sizes}kb_windows/snv_{kmer}mer/combined/frequency_{freq}_at_{fraction}p/{types}_dataframe_{kmer}mer.rds"
# 	conda: "envs/nmf.yaml"
# 	resources:
# 		threads=4,
# 		time=800,
# 		mem_mb=10000
# 	output:
# 		model = "{window_sizes}mb_windows/models/frequency_{freq}_at_{fraction}p/{types}_{kmer}mer/{types}_{kmer}mer_{signatures}.rds"
# 	shell:"""
#     Rscript scripts/opportunity_modelselection.R {wildcards.signatures} {input.count_data} {output.model}
#     """