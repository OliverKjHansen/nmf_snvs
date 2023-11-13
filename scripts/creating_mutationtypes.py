import sys
# window_size = 3
# string = "chr1_10mb_to_5mb"

# #print(string.split("_")[1].split("m")[0])

# #print(string.split("_")[0].split("m")[0])
# # for idx in range(0,10, window_size):
# #     print("chr"+str(idx)+str((idx+window_size)))

# while True:
#     print("waiting")

base_table = {"A":"T", 
              "T":"A", 
              "G":"C", 
              "C":"G"}
  
def splitting_to_mutationstypes(file, ref, alt):
    with open(file, "r") as f:
        for line in f:
            #print(line)
            vcf_ref = line.split("\n")[0].split("\t")[2]
            vcf_alt = line.split("\n")[0].split("\t")[3]
            if vcf_ref.upper() == ref and vcf_alt.upper() == alt:
                print(line.split("\n")[0]) 
            if base_table[vcf_ref].upper() == ref and base_table[vcf_alt].upper() == alt:
                print(line.split("\n")[0])


# def counting_indels(variant_file):
#     insertion = 0
#     deletion = 0
#     with open(variant_file, "r") as f:
#         for line in f:
#             variant_indel = line.split("\n")[0].split("\t")
#             if len(variant_indel[2]) > len(variant_indel[3]):
#                 deletion += 1
#             if len(variant_indel[2]) < len(variant_indel[3]):
#                 insertion += 1
#         print("insetions", insertion)
#         print("deletions", deletion)
#         print("indel", deletion+insertion)


if __name__ == '__main__':
    variants = sys.argv[1]
    ref_mut = sys.argv[2]
    alt_mut = sys.argv[3]
    splitting_to_mutationstypes(variants,ref_mut.upper(),alt_mut.upper())


