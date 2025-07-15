''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  cuteSVTrio_Description.py
 * @author: Lixin
 * @date: Apr. 19th 2025
 * @version V0.1.0
'''
import argparse

VERSION = '0.1.0'

class cuteSVTriodp(object):
	'''
	Detailed descriptions of cuteSVTrio version and its parameters.
	'''

	USAGE="""\
		
	Current version: v%s
	Author: Lixin
	Contact: xinli01@stu.hit.edu.cn

	Suggestions:

	For PacBio CLR data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	200
		--diff_ratio_merging_DEL	0.5

	For PacBio CCS(HIFI) data:
		--max_cluster_bias_INS		1000
		--diff_ratio_merging_INS	0.9
		--max_cluster_bias_DEL	1000
		--diff_ratio_merging_DEL	0.5

	For ONT data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	100
		--diff_ratio_merging_DEL	0.3


	"""%(VERSION)

	# MinSizeDel = 'For current version of cuteSVTrio, it can detect deletions larger than this size.'

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="cuteSVTrio", 
		description=cuteSVTriodp.USAGE, 
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--version', '-v', 
		action = 'version', 
		version = '%(prog)s {version}'.format(version=VERSION))

	# **************Parameters of input******************
	parser.add_argument('-r', '--reference', 
		help = "The reference genome in fasta format.", 
		type = str)
	parser.add_argument('-o', '--output', 
		type = str, 
		help = "Output VCF format file.")
	parser.add_argument('-w', '--work_dir', 
		type = str, 
		help = "Work-directory for distributed jobs.")

	# ************** Other Parameters******************
	parser.add_argument('-e', '--execute_stage', 
		help = "The stage of this operation execution.[%(default)s]", 
		default = 0, 
		type = int)
	parser.add_argument('--performing_phasing',
		help = "The option of performing structural variant phasing.",
		action="store_true")
	parser.add_argument('--parents_phasing',
		help = "The option of performing parents' structural variant phasing.",
		action="store_true")
	parser.add_argument('--phase_all_ctgs',
		help = "Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}.",
		action="store_true")
	parser.add_argument('-f', '--family_mode', 
		help = "Mode of members in family.[%(default)s]", 
		default = "M1", 
		type = str)
	parser.add_argument('-io', "--input_offspring", 
		metavar="[BAM]", 
		type = str, 
		help ="Sorted .bam file of offspring in family from NGMLR or Minimap2.")
	parser.add_argument('-if', "--input_parent_1", 
		metavar="[BAM]", 
		type = str, 
		help ="Sorted .bam file of father or only parent in family from NGMLR or Minimap2.")
	parser.add_argument('-im', "--input_parent_2", 
		metavar="[BAM]", 
		type = str, 
		help ="Sorted .bam file of mother in family from NGMLR or Minimap2.")
	
	parser.add_argument('-t', '--threads', 
		help = "Number of threads to use.[%(default)s]", 
		default = 16, 
		type = int)
	parser.add_argument('-b', '--batches', 
		help = "Batch of genome segmentation interval.[%(default)s]", 
		default = 10000000, 
		type = int)
	# The description of batches needs to improve.
	parser.add_argument('-S', '--sample',
		help = "Sample name/id",
		default = "NULL",
		type = str)

	parser.add_argument('--retain_work_dir',
		help = "Enable to retain temporary folder and files.",
		action="store_true")
	
	parser.add_argument('--write_old_sigs',
		help = "Enable to write sigs file in temporary folder for legacy compatibilities.",
		action="store_true")

	parser.add_argument('--report_readid',
		help = "Enable to report supporting read ids for each SV.",
		action="store_true")
	
	parser.add_argument('--run_TRA',
		help = "Run calling of TRA/BND.",
		action="store_true")

	# **************Parameters in signatures collection******************
	GroupSignaturesCollect = parser.add_argument_group('Collection of SV signatures')
	GroupSignaturesCollect.add_argument('-sp', '--max_split_parts', 
		help = "Maximum number of split segments a read may be aligned before it is ignored. All split segments are considered when using -1. \
			(Recommand -1 when applying assembly-based alignment.)[%(default)s]", 
		default = 7, 
		type = int)
	GroupSignaturesCollect.add_argument('-mq', '--min_mapq', 
		help = "Minimum mapping quality value of alignment to be taken into account.[%(default)s]", 
		default = 20, 
		type = int)
	GroupSignaturesCollect.add_argument('-mr', '--min_read_len', 
		help = "Ignores reads that only report alignments with not longer than bp.[%(default)s]", 
		default = 500, 
		type = int)
	GroupSignaturesCollect.add_argument('-md', '--merge_del_threshold', 
		help = "Maximum distance of deletion signals to be merged. In our paper, I used -md 500 to process HG002 real human sample data.[%(default)s]", 
		default = 0, 
		type = int)
	GroupSignaturesCollect.add_argument('-mi', '--merge_ins_threshold', 
		help = "Maximum distance of insertion signals to be merged. In our paper, I used -mi 500 to process HG002 real human sample data.[%(default)s]", 
		default = 100, 
		type = int)
	GroupSignaturesCollect.add_argument('-include_bed', 
		help = "Optional given bed file. Only detect SVs in regions in the BED file. [NULL]",
		default = None,
        type = str)
	# The min_read_len in last version is 2000.
	# signatures with overlap need to be filtered

	# **************Parameters in clustering******************
	GroupSVCluster = parser.add_argument_group('Generation of SV clusters')
	GroupSVCluster.add_argument('-s', '--min_support_list', 
		help = "Minimum number of reads of each member of family that support a SV to be reported.It is recommended to divide the data coverage by 6.[%(default)s]",  
		type = str)
	GroupSVCluster.add_argument('-l', '--min_size', 
		help = "Minimum size of SV to be reported.[%(default)s]", 
		default = 30, 
		type = int)
	GroupSVCluster.add_argument('-L', '--max_size', 
		help = "Maximum size of SV to be reported. All SVs are reported when using -1. [%(default)s]", 
		default = 100000, 
		type = int)
	GroupSVCluster.add_argument('-sl', '--min_siglength', 
		help = "Minimum length of SV signal to be extracted.[%(default)s]", 
		default = 10, 
		type = int)
	GroupSVCluster.add_argument('--all_ins_singnature_reads',
		help = "Print all singnature reads of INS in output vcf.",
		action="store_true")

	# **************Parameters in genotyping******************
	GroupGenotype = parser.add_argument_group('Computing genotypes')
	GroupGenotype.add_argument('--gt_round', 
		help = "Maximum round of iteration for alignments searching if perform genotyping.[%(default)s]", 
		default = 500, 
		type = int)
	GroupGenotype.add_argument('--read_range', 
		help = "The interval range for counting reads distribution.[%(default)s]", 
		default = 1000, 
		type = int)
	# GroupGenotype.add_argument('--hom', 
	# 	help = "Threshold on allele frequency for homozygous.[%(default)s]", 
	# 	default = 0.8, 
	# 	type = float)
	# GroupGenotype.add_argument('--het', 
	# 	help = "Threshold on allele frequency for heterozygous.[%(default)s].", 
	# 	default = 0.2, 
	# 	type = float)

	# Just a parameter for debug.
	# Will be removed in future.
	# GroupSVCluster.add_argument('--preset',
	# 	help = "Parameter presets for different sequencing technologies (pbclr/pbccs/ont).[%(default)s]",
	# 	default = "pbccs",
	# 	type = str)

	# **************Parameters in force calling******************
	GroupGenotype = parser.add_argument_group('Force calling')
	GroupGenotype.add_argument('-Ivcf', #'--MERGED_VCF',
		help = "Optional given vcf file. Enable to perform force calling. [NULL]",
		default = None,
        type = str)

	# **************Parameters in force phasing******************
	GroupGenotype = parser.add_argument_group('Phasing')
	GroupGenotype.add_argument('-read_pos_interval',
		help = "Position division interval of the read position. [%(default)s]",
		default = 1000,
        type = int)

	# **************Advanced Parameters******************
	GroupAdvanced = parser.add_argument_group('Advanced')
	# ++++++Total++++++
	GroupAdvanced.add_argument('-p', '--sequencing_platform', 
		help = "The option of sequencing platform affects a series of parameters in the signature clustering.", 
		default = None, 
		type = str)
	# ++++++Total++++++
	GroupAdvanced.add_argument('-g', '--gold_standard_version', 
		help = "The gold standard version of variant detection fitting.", 
		default = "T2T", 
		type = str)
	# ++++++INS++++++
	GroupAdvanced.add_argument('--max_cluster_bias_INS', 
		help = "Maximum distance to cluster read together for insertion.[%(default)s]", 
		default = 50, 
		type = int)
	GroupAdvanced.add_argument('--diff_ratio_merging_INS', 
		help = "Do not merge breakpoints with basepair identity more than [%(default)s] for insertion.", 
		default = 0.05, 
		type = float)
	# GroupAdvanced.add_argument('--diff_ratio_filtering_INS', 
	# 	help = "Filter breakpoints with basepair identity less than [%(default)s] for insertion.", 
	# 	default = 0.6, 
	# 	type = float)

	# ++++++DEL++++++
	GroupAdvanced.add_argument('--max_cluster_bias_DEL', 
		help = "Maximum distance to cluster read together for deletion.[%(default)s]", 
		default = 100, 
		type = int)
	GroupAdvanced.add_argument('--diff_ratio_merging_DEL', 
		help = "Do not merge breakpoints with basepair identity more than [%(default)s] for deletion.", 
		default = 0.05, 
		type = float)
	# GroupAdvanced.add_argument('--diff_ratio_filtering_DEL', 
	# 	help = "Filter breakpoints with basepair identity less than [%(default)s] for deletion.", 
	# 	default = 0.7, 
	# 	type = float)

	# ++++++INV++++++
	GroupAdvanced.add_argument('--max_cluster_bias_INV', 
		help = "Maximum distance to cluster read together for inversion.[%(default)s]", 
		default = 500, 
		type = int)

	# ++++++DUP++++++
	GroupAdvanced.add_argument('--max_cluster_bias_DUP', 
		help = "Maximum distance to cluster read together for duplication.[%(default)s]", 
		default = 500, 
		type = int)

	# ++++++TRA++++++
	GroupAdvanced.add_argument('--max_cluster_bias_TRA', 
		help = "Maximum distance to cluster read together for translocation.[%(default)s]", 
		default = 50, 
		type = int)
	GroupAdvanced.add_argument('--diff_ratio_filtering_TRA', 
		help = "Filter breakpoints with basepair identity less than [%(default)s] for translocation.", 
		default = 0.6, 
		type = float)

	GroupAdvanced.add_argument('--remain_reads_ratio', 
		help = "The ratio of reads remained in cluster. Set lower when the alignment data have high quality but recommand over 0.5.[%(default)s]", 
		default = 1.0, 
		type = float)

	# parser.add_argument('-d', '--max_distance', 
	# 	help = "Maximum distance to group SV together..[%(default)s]", 
	# 	default = 1000, type = int)



	# These parameters are drawn lessons from pbsv v2.2.0
	# parser.add_argument('--min_del_size', 
	# 	help = "Minimum size of a deletion.[%(default)s]", 
	# 	default = 20, type = int)

	args = parser.parse_args(argv)
	return args

def Generation_VCF_header(file, contiginfo, sample, argv):
	# General header
	file.write("##fileformat=VCFv4.2\n")
	file.write("##source=cuteSVTrio-%s\n"%(VERSION))
	import time
	file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
	for i in contiginfo:
		file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

	# Specific header
	# ALT
	file.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
	file.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
	file.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
	file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
	file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

	# INFO
	file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
	file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
	file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
	file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
	file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
	file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
	file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
	file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
	# file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
	file.write("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
	file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
	file.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
	file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency.\">\n")
	file.write("##INFO=<ID=No_Mendel,Number=0,Type=Flag,Description=\"Not in accordance with Mendelian laws of inheritance.\">\n")
	file.write("##INFO=<ID=CorrectType,Number=1,Type=Integer,Description=\"The trio SV correct type.\">\n")
	file.write("##INFO=<ID=Denovo,Number=1,Type=Integer,Description=\"The high confidence de novo class.\">\n")
	file.write("##INFO=<ID=QUALLIST,Number=3,Type=Float,Description=\"Quality of all family members.\">\n")
	file.write("##INFO=<ID=FILTERLIST,Number=3,Type=String,Description=\"Filter flags for all family members.\">\n")
	file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
	# FORMAT
	# file.write("\n")
	file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	file.write("##FORMAT=<ID=HP_GT,Number=1,Type=String,Description=\"Haplotyped Genotype\">\n")
	file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
	file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
	file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
	file.write("##FORMAT=<ID=AP,Number=2,Type=Float,Description=\"ALT allele probability of first haplotype and second haplotype\">\n")
	file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")

	file.write("##CommandLine=\"cuteSVTrio %s\"\n"%(" ".join(argv)))
