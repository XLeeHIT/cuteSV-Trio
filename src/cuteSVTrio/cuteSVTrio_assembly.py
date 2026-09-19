


def msa_consensus_for_cluster(cluster_seqs, aligner=None):
    pass

def capture_reads_within_lcr(read, start, end) :
    pass

def redetect_nearby_sv(path, chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, hap_range, homo_range, family_mode, family_bams, minimum_support_reads_list, remap_merge_k, remap_minimizer_window, assembly_correction_threshold, assembly_correction_setting, family_assembly_reads, similarity_supplement_threshold, assembly_accelerate) :
    pass

# 超长的DEL合INS，可能会出现不同成员长度不匹配的情况
def make_extralong_sv_neat(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# minimap2错误的将邻接的两个sv的read比对到了一起，导致出现fp和fn
def correct_2svs_align_together(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 因为minimap2的比对可能会产生比对的断裂，导致原本的长变异消失，出现两个错误的短变异
def correct_2svs_align_crack(chr, tmp_svs, member_svs, family_member_ls) :
    pass

def verify_allele_svs_byassembly(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 按照ass的检测结果修正现有sv的pos
def refine_svs_pos(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 两个位置相近的杂合，minimap2的比对和cutesvTrio的聚类会是他们被错误的识别为纯和，需要按照ass的结果将两者分开
def split_error_homo_svs(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    pass

# 类等位基因中，如果两个sv的长度都在50附近，就是用assembly做筛选
def screening_allele_boundary_len(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 针对家庭纠错方法种后两类补充的sv，使用assembly重新验证其正确性
def revalidate_trio_correction(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    pass

# 使用assembly验证denovo，尤其是第一类denovo的正确性
def revalidate_denovo(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    pass

# 三重及以上等位基因中，如果assembly和原本都只支持的一个等位基因，但是不一致，那么修改为assembly支持的
# 三重及以上等位基因中，如果assembly不支持，但是原本支持一个，那去除该sv
# 三重及以上等位基因中，如果原本检测记录中没有变异，但是assembly支持了一个，很有可能是minimap2的比对错误导致的
def screening_threeup_allele(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 过滤前后紧密连接在一起的两个变异，使用assembly重新验证他们的正确性
# 效果不好，暂时不用
def screening_connected_svs(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 通过assembly中位置和长度几乎完全相同的记录来补充sv
def supple_samesv_bybyassembly(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 如果assembly支持的sv是一个左右完全没有原本检测变异支持的sv，那么反向证明了其可能的正确性
def reverse_supple_byunsupport_sv(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 如果原检测结果中两个变异接近，并且长度几乎一致，而assembly只支持一个，两个sv的支持信号互斥，那么只保留assembly支持的
def screening_two_near_svs(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    pass

# 使用lcr中确定的单倍型分支结果，重新验证纯和变异的正确性
def correct_homogt_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 使用lcr中确定的单倍型分支结果，重新验证杂合变异的正确性
def correct_hetegt_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 针对19号染色体进行针对CMRG的特异性优化
def correct_19_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 针对长并整齐的信号进行针对CMRG的特异性补充
def supple_long_neat_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    pass

# 过滤列表中的所有纯和变异，随机分配纯和变异的所有read，重复多次，如果中间出现了一定次数的杂合结果，那么将其修改为杂合
# 使用随机的方法将范围内的所有read随机的分配到两个hap中
def correct_homogt_allsv_random(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window) :
    pass

# 过滤列表中的所有纯和变异，随机分配纯和变异的所有read，重复多次，如果中间出现了一定次数的杂合结果，那么将其修改为杂合
# 将支持read分配给一个hap，范围内的其余read分配给另一个hap
def correct_homogt_allsv_support(path, chr, tmp_svs, fa_hap_ls, mo_hap_ls, homo_range, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list, family_assembly_reads, similarity_supplement_threshold, assembly_accelerate) :
    pass

# 将第二类denovo，或者无中生有的纯和结果提取出来，重新验证一遍
def correct_homogt_allsv_denovo(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list) :
    pass

# 将杂合结果提取出来，重新验证一遍
def correct_hetegt_allsv_denovo(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list) :
    pass

def sv_from_cigar(chr, cigar_str, ref_start, sv_min_size, origin_ref_start, sequence, reference):
    pass

# 按照同类sv的位置和长度的临近关系，去除相近的相似sv
def remove_redundant_sv(chr, sv_ls) :
    pass

def run_assembly(args) :
    return redetect_nearby_sv(*args)
