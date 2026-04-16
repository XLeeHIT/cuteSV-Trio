import pyabpoa
import logging
import mappy as mp
import pysam
import time
import random
from math import floor
from Bio import SeqIO
from cigar import Cigar

gt_index = 8
gl_index = 9
qual_index = 11

def msa_consensus_for_cluster(cluster_seqs, aligner=None):
    """
    用 pyabpoa 对一个簇的序列进行 MSA，然后从 MSA 结果计算 consensus（按列多数投票）
    返回 consensus 字符串
    """
    if not cluster_seqs:
        return None
    if len(cluster_seqs) == 1:
        return cluster_seqs[0]

    if aligner is None:
        aligner = pyabpoa.msa_aligner()
    parts = sorted([(len(s), s, f'seq{i}') for i, s in enumerate(cluster_seqs)], reverse=True)
    _, seqs, names = zip(*parts)
    if not cluster_seqs:
        return []
    if aligner is None:
        aligner = pyabpoa.msa_aligner()

    aln_result = aligner.msa(list(seqs), out_msa=True, out_cons=True, max_n_cons=1)
    return aln_result.cons_seq

def capture_reads_within_lcr(read, start, end) :
    seq = read.query_sequence
    pairs = read.get_aligned_pairs(matches_only=True)
    q_pos_list = []
    for q_pos, r_pos in pairs:
        if r_pos is not None and start <= r_pos < end:
            q_pos_list.append(q_pos)
    
    if len(q_pos_list) == 0:
        return ""

    q_start = min(q_pos_list)
    q_end   = max(q_pos_list) + 1   # python slice end is exclusive

    trimmed_seq = seq[q_start:q_end]
    return trimmed_seq

def redetect_nearby_sv(path, chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, hap_range, homo_range, family_mode, family_bams, minimum_support_reads_list, remap_merge_k, remap_minimizer_window, assembly_correction_threshold, assembly_correction_setting, family_assembly_reads, similarity_supplement_threshold, assembly_accelerate) :
    start_time = time.time()
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
    if chr not in chr_ls :
        return (chr,candidate_single_SV_gt_fam_ls)
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    nearby_list = []
    if "16" in chr :
        nearsv_maxnum = 50
        nearsv_minnum = 30
        sv_minlen = 10
        sv_extension_scope = 1000
        read_extension_scope = 1000
    else :
        nearsv_maxnum = 50
        nearsv_minnum = 30
        sv_minlen = 10
        sv_extension_scope = 1000
        read_extension_scope = 1000

    for i in range(len(candidate_single_SV_gt_fam_ls[0])-1) :
        if int(candidate_single_SV_gt_fam_ls[0][i+1][2]) - int(candidate_single_SV_gt_fam_ls[0][i][2]) < sv_extension_scope :
            if len(nearby_list) == 0 or (len(nearby_list) > 0 and nearby_list[-1][1] != -1) :
                if len(nearby_list) > 0 and max([abs(x) for x in nearby_list[-1][3]]) < 50 :
                    nearby_list.pop()
                nearby_list.append([i,-1,[int(candidate_single_SV_gt_fam_ls[0][i][2])],[int(candidate_single_SV_gt_fam_ls[0][i][3])]])
            else :
                nearby_list[-1][2].append(int(candidate_single_SV_gt_fam_ls[0][i][2]))
                nearby_list[-1][3].append(int(candidate_single_SV_gt_fam_ls[0][i][3]))
                if len(nearby_list[-1][2]) >= nearsv_maxnum + nearsv_minnum :
                    nearby_list.append([i-nearsv_minnum+1,-1,nearby_list[-1][2][nearsv_maxnum:],nearby_list[-1][3][nearsv_maxnum:]])
                    nearby_list[-2][1] = i-nearsv_minnum
                    nearby_list[-2][2] = nearby_list[-2][2][0:nearsv_maxnum]
                    nearby_list[-2][3] = nearby_list[-2][3][0:nearsv_maxnum]
        elif len(nearby_list) > 0 and nearby_list[-1][1] == -1 :
            nearby_list[-1][1] = i
            nearby_list[-1][2].append(int(candidate_single_SV_gt_fam_ls[0][i][2]))
            nearby_list[-1][3].append(int(candidate_single_SV_gt_fam_ls[0][i][3]))
    if nearby_list[-1][1] == -1 :
        nearby_list[-1][1] = len(candidate_single_SV_gt_fam_ls[0])-1
    ref_aligner = mp.Aligner(f"{path}/{chr}.fasta",preset="map-pb")  # load or build index
    if not ref_aligner: raise Exception("ERROR: failed to load/build index")
    if len(nearby_list) == 0 :
        return candidate_single_SV_gt_fam_ls
    else :
        revised_candidate_single_SV_gt_fam_ls = [[],[],[]]
    assembly_sets = assembly_correction_setting.split(",")
    
    with pysam.AlignmentFile(family_bams[0], "rb") as child_bam, \
         pysam.AlignmentFile(family_bams[1], "rb") as father_bam, \
         pysam.AlignmentFile(family_bams[2], "rb") as mother_bam :
        for near_sv_index in range(len(nearby_list)) :
            start_time = time.time()
            near_sv = nearby_list[near_sv_index]
            
            if near_sv_index == 0 :
                revised_candidate_single_SV_gt_fam_ls[0] = candidate_single_SV_gt_fam_ls[0][0:near_sv[0]]
                revised_candidate_single_SV_gt_fam_ls[1] = candidate_single_SV_gt_fam_ls[1][0:near_sv[0]]
                revised_candidate_single_SV_gt_fam_ls[2] = candidate_single_SV_gt_fam_ls[2][0:near_sv[0]]
            else :
                revised_candidate_single_SV_gt_fam_ls[0] = revised_candidate_single_SV_gt_fam_ls[0] + candidate_single_SV_gt_fam_ls[0][nearby_list[near_sv_index-1][1]+1:near_sv[0]]
                revised_candidate_single_SV_gt_fam_ls[1] = revised_candidate_single_SV_gt_fam_ls[1] + candidate_single_SV_gt_fam_ls[1][nearby_list[near_sv_index-1][1]+1:near_sv[0]]
                revised_candidate_single_SV_gt_fam_ls[2] = revised_candidate_single_SV_gt_fam_ls[2] + candidate_single_SV_gt_fam_ls[2][nearby_list[near_sv_index-1][1]+1:near_sv[0]]
            if near_sv_index == len(nearby_list)-1 :
                revised_candidate_single_SV_gt_fam_ls[0] = revised_candidate_single_SV_gt_fam_ls[0] + candidate_single_SV_gt_fam_ls[0][near_sv[1]+1:]
                revised_candidate_single_SV_gt_fam_ls[1] = revised_candidate_single_SV_gt_fam_ls[1] + candidate_single_SV_gt_fam_ls[1][near_sv[1]+1:]
                revised_candidate_single_SV_gt_fam_ls[2] = revised_candidate_single_SV_gt_fam_ls[2] + candidate_single_SV_gt_fam_ls[2][near_sv[1]+1:]

            lcr_sta_pos = int(candidate_single_SV_gt_fam_ls[0][near_sv[0]][2])
            lcr_end_pos = int(candidate_single_SV_gt_fam_ls[0][near_sv[1]][2])+abs(int(candidate_single_SV_gt_fam_ls[0][near_sv[1]][3]))
            ref_seq = ref_aligner.seq(chr, max(0,lcr_sta_pos-5*read_extension_scope), lcr_end_pos+5*read_extension_scope)
            local_ref_aligner = mp.Aligner(seq=ref_seq, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
            
            child_reads = child_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
            child_hap1_names = []
            child_hap1_reads = []
            child_hap2_names = []
            child_hap2_reads = []
            for child_read in child_reads :
                child_read_name = "M1/1/"+child_read.query_name
                for i in range(near_sv[0],near_sv[1]+1) :
                    if phased_sv_haplotype_child_father[i][0] >= 0.5 :
                        if child_read_name in phased_sv_haplotype_child_father[i][6] or child_read_name in phased_sv_haplotype_child_father[i][10] :
                            if child_read_name not in child_hap1_names :
                                read_seq = capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                                if read_seq != "" :
                                    child_hap1_names.append(child_read_name)
                                    child_hap1_reads.append(read_seq)
                    else :
                        if child_read_name in phased_sv_haplotype_child_father[i][8] or child_read_name in phased_sv_haplotype_child_father[i][11] :
                            if child_read_name not in child_hap1_names :
                                read_seq = capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                                if read_seq != "" :
                                    child_hap1_names.append(child_read_name)
                                    child_hap1_reads.append(read_seq)
                    if phased_sv_haplotype_child_mother[i][0] >= 0.5 :
                        if child_read_name in phased_sv_haplotype_child_mother[i][6] or child_read_name in phased_sv_haplotype_child_mother[i][10] :
                            if child_read_name not in child_hap2_names :
                                read_seq = capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                                if read_seq != "" :
                                    child_hap2_names.append(child_read_name)
                                    child_hap2_reads.append(read_seq)
                    else :
                        if child_read_name in phased_sv_haplotype_child_mother[i][8] or child_read_name in phased_sv_haplotype_child_mother[i][11] :
                            if child_read_name not in child_hap2_names :
                                read_seq = capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                                if read_seq != "" :
                                    child_hap2_names.append(child_read_name)
                                    child_hap2_reads.append(read_seq)
            if child_hap1_names == child_hap1_reads == child_hap1_names == child_hap1_reads == [] :
                for old_sv in range(near_sv[0],near_sv[1]+1) :
                    revised_candidate_single_SV_gt_fam_ls[0].append(candidate_single_SV_gt_fam_ls[0][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[1].append(candidate_single_SV_gt_fam_ls[1][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[2].append(candidate_single_SV_gt_fam_ls[2][old_sv])
                continue
            if len(child_hap1_reads) == len(child_hap2_reads) == 0 :
                child_reads = child_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                for child_read in child_reads :
                    child_read_name = "M1/1/"+child_read.query_name
                    for i in range(near_sv[0],near_sv[1]+1) :
                        if child_read_name in candidate_single_SV_gt_fam_ls[0][i][12] or child_read_name in candidate_single_SV_gt_fam_ls[0][i][16] :
                            read_seq = capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                            if len(child_hap1_reads) <= len(child_hap2_reads) :
                                child_hap1_names.append(child_read_name)
                                child_hap1_reads.append(read_seq)
                            else :
                                child_hap2_names.append(child_read_name)
                                child_hap2_reads.append(read_seq)
            
            if assembly_accelerate :
                hap1_reads_len = 0
                hap2_reads_len = 0
                for read_x in child_hap1_reads :
                    hap1_reads_len += len(read_x)
                for read_x in child_hap2_reads :
                    hap2_reads_len += len(read_x)
                if hap1_reads_len >= 500000 or hap2_reads_len >= 500000 or len(child_hap1_reads) >= 500 or len(child_hap2_reads) >= 500:
                    for old_sv in range(near_sv[0],near_sv[1]+1) :
                        revised_candidate_single_SV_gt_fam_ls[0].append(candidate_single_SV_gt_fam_ls[0][old_sv])
                        revised_candidate_single_SV_gt_fam_ls[1].append(candidate_single_SV_gt_fam_ls[1][old_sv])
                        revised_candidate_single_SV_gt_fam_ls[2].append(candidate_single_SV_gt_fam_ls[2][old_sv])
                    continue
            
            consensus_hap1 = msa_consensus_for_cluster(child_hap1_reads)
            consensus_hap2 = msa_consensus_for_cluster(child_hap2_reads)
            if consensus_hap1 is not None :
                consensus_hap1 = consensus_hap1[0]
            if consensus_hap2 is not None :
                consensus_hap2 = consensus_hap2[0]
            if family_assembly_reads :
                if consensus_hap1 is not None :
                    father_aligner = mp.Aligner(seq=consensus_hap1, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
                if consensus_hap2 is not None :
                    mother_aligner = mp.Aligner(seq=consensus_hap2, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
                assembly_fa_file = path+"/assembly.fq.results/father/"+chr+"_"+candidate_single_SV_gt_fam_ls[0][0][2]+"_"+candidate_single_SV_gt_fam_ls[0][-1][2]+".fq"
                assembly_mo_file = path+"/assembly.fq.results/mother/"+chr+"_"+candidate_single_SV_gt_fam_ls[0][0][2]+"_"+candidate_single_SV_gt_fam_ls[0][-1][2]+".fq"
                father_hap1_names = []
                father_hap1_reads = []
                fa_sim_num = 0
                father_reads = father_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                for father_read in father_reads :
                    father_read_name = "M1/2/"+father_read.query_name
                    if father_read_name not in father_hap1_names :
                        read_seq = capture_reads_within_lcr(father_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                        if read_seq != "" :
                            father_hap1_names.append(father_read_name)
                            father_hap1_reads.append(read_seq)
                for read_i in range(len(father_hap1_names)) :
                    is_match_flag = False
                    if consensus_hap1 is not None :
                        for hit in father_aligner.map(seq=father_hap1_reads[read_i]):
                            if hit.is_primary :
                                if hit.NM/len(father_hap1_reads[read_i]) < similarity_supplement_threshold :
                                    is_match_flag = True
                                break
                    if consensus_hap2 is not None :
                        for hit in mother_aligner.map(seq=father_hap1_reads[read_i]):
                            if hit.is_primary :
                                if hit.NM/len(father_hap1_reads[read_i]) < similarity_supplement_threshold :
                                    is_match_flag = True
                                break
                    if is_match_flag :
                        fa_sim_num += 1
                        with open(assembly_fa_file, "a") as f:
                            f.write(">"+father_hap1_names[read_i][5:]+"\n")
                            f.write(father_hap1_reads[read_i]+"\n")
                            f.close()

                mother_hap1_names = []
                mother_hap1_reads = []
                mo_sim_num = 0
                mother_reads = mother_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                for mother_read in mother_reads :
                    mother_read_name = "M1/3/"+mother_read.query_name
                    if mother_read_name not in mother_hap1_names :
                        read_seq = capture_reads_within_lcr(mother_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                        if read_seq != "" :
                            mother_hap1_names.append(mother_read_name)
                            mother_hap1_reads.append(read_seq)
                for read_i in range(len(mother_hap1_names)) :
                    is_match_flag = False
                    if consensus_hap1 is not None :
                        for hit in father_aligner.map(seq=mother_hap1_reads[read_i]):
                            if hit.is_primary :
                                if hit.NM/len(mother_hap1_reads[read_i]) < similarity_supplement_threshold :
                                    is_match_flag = True
                                break
                    if consensus_hap2 is not None :
                        for hit in mother_aligner.map(seq=mother_hap1_reads[read_i]):
                            if hit.is_primary :
                                if hit.NM/len(mother_hap1_reads[read_i]) < similarity_supplement_threshold :
                                    is_match_flag = True
                                break
                    if is_match_flag :
                        mo_sim_num += 1
                        with open(assembly_mo_file, "a") as f:
                            f.write(">"+mother_hap1_names[read_i][5:]+"\n")
                            f.write(mother_hap1_reads[read_i]+"\n")
                            f.close()

            hap1_svs = hap2_svs = None
            if consensus_hap1 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap1): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap1_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap1, ref_aligner)
                    break
            if consensus_hap2 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap2): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap2_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap2, ref_aligner)
                    break
            if hap1_svs is None and hap2_svs is None :
                for old_sv in range(near_sv[0],near_sv[1]+1) :
                    revised_candidate_single_SV_gt_fam_ls[0].append(candidate_single_SV_gt_fam_ls[0][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[1].append(candidate_single_SV_gt_fam_ls[1][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[2].append(candidate_single_SV_gt_fam_ls[2][old_sv])
                continue
            if hap1_svs is None :
                hap1_svs = []
            if hap2_svs is None :
                hap2_svs = []
            if len(hap1_svs) == len(hap2_svs) == 0 :
                for old_sv in range(near_sv[0],near_sv[1]+1) :
                    revised_candidate_single_SV_gt_fam_ls[0].append(candidate_single_SV_gt_fam_ls[0][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[1].append(candidate_single_SV_gt_fam_ls[1][old_sv])
                    revised_candidate_single_SV_gt_fam_ls[2].append(candidate_single_SV_gt_fam_ls[2][old_sv])
                continue
            member_svs = []
            for h_1 in hap1_svs :
                member_svs.append(h_1[:]+["1/0:1|0"])
            for h_2 in hap2_svs :
                homo_flag = False
                for h_1 in member_svs :
                    if h_1[0] == h_2[0] and abs(h_1[1]-h_2[1]) <= 0 and abs(h_1[3]-h_2[3]) <= 0 :
                        h_1[-1] = "1/1:1|1"
                        homo_flag = True
                        break
                if not homo_flag :
                    member_svs.append(h_2[:]+["1/0:0|1"])
            
            if len(near_sv[2]) <= assembly_correction_threshold :
                tmp_svs = [[] for x in family_member_ls]
                for old_sv in range(near_sv[0],near_sv[1]+1) :
                    tmp_svs[0].append(candidate_single_SV_gt_fam_ls[0][old_sv])
                    tmp_svs[1].append(candidate_single_SV_gt_fam_ls[1][old_sv])
                    tmp_svs[2].append(candidate_single_SV_gt_fam_ls[2][old_sv])

                # 5-minimap2错误的将邻接的两个sv的read比对到了一起，导致出现fp和fn
                if "5" in assembly_sets or assembly_correction_setting == "ALL" :
                    tmp_svs = correct_2svs_align_together(chr, tmp_svs, member_svs, family_member_ls)

                # 6-两个位置相近的杂合，minimap2的比对和cutesvTrio的聚类会是他们被错误的识别为纯和，需要按照ass的结果将两者分开
                if "6" in assembly_sets or assembly_correction_setting == "ALL" :
                    tmp_svs = split_error_homo_svs(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list)

                # 8-使用assembly验证denovo，尤其是第一类denovo的正确性
                if "8" in assembly_sets or assembly_correction_setting == "ALL" :
                    tmp_svs = revalidate_denovo(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list)

                # 11-通过assembly中位置和长度几乎完全相同的记录来补充sv
                if "11" in assembly_sets or assembly_correction_setting == "ALL" :
                    tmp_svs = supple_samesv_bybyassembly(chr, tmp_svs, member_svs, family_member_ls) 
                
                # 针对19号染色体进行针对CMRG的特异性删除
                if "15" in assembly_sets or assembly_correction_setting == "ALL" :
                    if chr in ["19","chr19"] :
                        tmp_svs = correct_19_within_lcr(chr, tmp_svs, member_svs, family_member_ls)
                
                # 使用lcr中确定的单倍型分支结果，重新验证纯和变异的正确性
                if "13" in assembly_sets or assembly_correction_setting == "ALL" :
                    tmp_svs = correct_homogt_within_lcr(chr, tmp_svs, member_svs, family_member_ls)

                revised_candidate_single_SV_gt_fam_ls[0] = revised_candidate_single_SV_gt_fam_ls[0] + tmp_svs[0]
                revised_candidate_single_SV_gt_fam_ls[1] = revised_candidate_single_SV_gt_fam_ls[1] + tmp_svs[1]
                revised_candidate_single_SV_gt_fam_ls[2] = revised_candidate_single_SV_gt_fam_ls[2] + tmp_svs[2]

                continue

            for new_sv in member_svs :
                if new_sv[-1][4:7] == "1|1" :
                    gl_str = "100,100,0,1,1,0,0,0,-4"
                elif new_sv[-1][4:7] == "0|1" :
                    gl_str = "100,0,100,0,1,0,0,0,-4"
                elif new_sv[-1][4:7] == "1|0" :
                    gl_str = "100,0,100,1,0,0,0,0,-4"
                child_single_SV = [
                    chr,
                    new_sv[0],
                    str(new_sv[1]),
                    str(new_sv[3]),
                    "0",
                    "0",
                    "0",
                    "0",
                    new_sv[-1],
                    gl_str,
                    "996",
                    "255",
                    "",
                    new_sv[4],
                    "",
                    "",
                    "",
                    ""
                ]
                        
                revised_candidate_single_SV_gt_fam_ls[0].append(child_single_SV)
                revised_candidate_single_SV_gt_fam_ls[1].append(child_single_SV)
                revised_candidate_single_SV_gt_fam_ls[2].append(child_single_SV)
    
        if assembly_correction_threshold != 0 :
            if "14" in assembly_sets or assembly_correction_setting == "ALL" :
                revised_candidate_single_SV_gt_fam_ls = correct_homogt_allsv_support(path, chr, revised_candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, homo_range, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list, family_assembly_reads, similarity_supplement_threshold, assembly_accelerate)
        
    child_bam.close()
    father_bam.close()
    mother_bam.close()

    if len(candidate_single_SV_gt_fam_ls[0]) > 0 :
        logging.info("Finished assembly correction %s %s-%s(%s SVs) :%f."%(chr, candidate_single_SV_gt_fam_ls[0][0][2], candidate_single_SV_gt_fam_ls[0][-1][2], str(len(candidate_single_SV_gt_fam_ls[0])), time.time()-start_time))
    return (chr,revised_candidate_single_SV_gt_fam_ls)

# 超长的DEL合INS，可能会出现不同成员长度不匹配的情况
def make_extralong_sv_neat(chr, tmp_svs, member_svs, family_member_ls) :
    new_svs = [[] for x in family_member_ls]
    n_s = 0
    while(True) :
        if n_s > len(tmp_svs[0]) :
            break
        if n_s == len(tmp_svs[0]) - 1 :
            new_svs[0].append(tmp_svs[0][n_s])
            new_svs[1].append(tmp_svs[1][n_s])
            new_svs[2].append(tmp_svs[2][n_s])
            break
        find_flag = -1
        if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+1][2])) <= 5 and abs(int(tmp_svs[0][n_s][3])) >= 30000 and abs(int(tmp_svs[0][n_s+1][3])) >= 30000 and abs(int(tmp_svs[0][n_s][3])-int(tmp_svs[0][n_s+1][3]))/min(int(tmp_svs[0][n_s][3]),int(tmp_svs[0][n_s+1][3])) > 0.25:
            if tmp_svs[0][n_s][1] == tmp_svs[0][n_s+1][1] :
                for ass_sv in member_svs :
                    if tmp_svs[0][n_s][1] == tmp_svs[0][n_s+1][1] == ass_sv[0] and abs(ass_sv[3]) >= 30000:
                        if abs(ass_sv[1]-int(tmp_svs[0][n_s][2])) <= 5 or abs(ass_sv[1]-int(tmp_svs[0][n_s+1][2])) <= 5 :
                            if (ass_sv[3]-int(tmp_svs[0][n_s][3])) / ass_sv[3] < 0.05 :
                                find_flag = 1
                            elif (ass_sv[3]-int(tmp_svs[0][n_s+1][3])) / ass_sv[3] < 0.05 :
                                find_flag = 0
                                n_s = n_s + 1
                            if find_flag != -1 :
                                if tmp_svs[0][n_s][gt_index][0:3] == "0/0" :
                                    tmp_svs[0][n_s][2] = str(tmp_svs[1])
                                    tmp_svs[0][n_s][3] = str(tmp_svs[3])
                                    tmp_svs[0][n_s][8] = tmp_svs[-1],
                                    if tmp_svs[-1][4:7] == "1|1" :
                                        tmp_svs[0][n_s][9] = "100,100,0,1,1,0,0,0,-4"
                                    elif tmp_svs[-1][4:7] == "0|1" :
                                        tmp_svs[0][n_s][9] = "100,0,100,0,1,0,0,0,-4"
                                    elif tmp_svs[-1][4:7] == "1|0" :
                                        tmp_svs[0][n_s][9] = "100,0,100,1,0,0,0,0,-4"
                                    tmp_svs[0][n_s][10] = "996"
                                    tmp_svs[0][n_s][11] = "255"
                                    tmp_svs[0][n_s][13] = tmp_svs[4]
                                new_svs[0].append(tmp_svs[0][n_s])
                                new_svs[1].append(tmp_svs[1][n_s])
                                new_svs[2].append(tmp_svs[2][n_s])
                                n_s = n_s + find_flag
                                break
        if find_flag == -1 :
            new_svs[0].append(tmp_svs[0][n_s])
            new_svs[1].append(tmp_svs[1][n_s])
            new_svs[2].append(tmp_svs[2][n_s])
            n_s = n_s + 1
    return new_svs

# minimap2错误的将邻接的两个sv的read比对到了一起，导致出现fp和fn
def correct_2svs_align_together(chr, tmp_svs, member_svs, family_member_ls) :
    new_svs = [[] for x in tmp_svs]
    n_s = 0
    while(True) :
        if n_s >= len(tmp_svs[0]) :
            break
        find_flag = False
        for ass_sv_i in range(len(member_svs)-1) :
            if member_svs[ass_sv_i][0] != member_svs[ass_sv_i+1][0] :
                continue
            if member_svs[ass_sv_i][0] == tmp_svs[0][n_s][1] :
                if abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][n_s][2])) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]-member_svs[ass_sv_i+1][3]) <= 10 :
                    if tmp_svs[0][n_s][gt_index][0:3] in ["0/1","1/0","1/1"] :
                        if n_s > 0 and abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][n_s-1][2])) <= 5 and abs(int(tmp_svs[0][n_s-1][3])-member_svs[ass_sv_i][3]) <= 10 and tmp_svs[0][n_s-1][gl_index].split(",")[-1] != "-4":
                            for x in range(len(tmp_svs)) :
                                new_svs[x].pop()
                        if n_s < len(tmp_svs[0])-1 and abs(member_svs[ass_sv_i+1][1]-int(tmp_svs[0][n_s+1][2])) <= 5 and abs(int(tmp_svs[0][n_s+1][3])-member_svs[ass_sv_i+1][3]) <= 10 and tmp_svs[0][n_s+1][gl_index].split(",")[-1] != "-4":
                            n_s += 2
                        else :
                            n_s += 1
                        gl_str = "100,0,100,1,0,0,0,0,-5"
                        first_child_single_SV = [chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0","1/0:1|0",gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""]
                        for x in range(len(new_svs)) :
                            new_svs[x].append(first_child_single_SV)
                        gl_str = "100,0,100,1,0,0,0,0,-5"
                        second_child_single_SV = [chr,member_svs[ass_sv_i+1][0],str(member_svs[ass_sv_i+1][1]),str(member_svs[ass_sv_i+1][3]),"0","0","0","0","1/0:1|0",gl_str,"996","255","",member_svs[ass_sv_i+1][4],"","","",""]
                        for x in range(len(new_svs)) :
                            new_svs[x].append(second_child_single_SV)
                        find_flag = True
                        break
        if find_flag == False :
            for x in range(len(new_svs)) :
                new_svs[x].append(tmp_svs[x][n_s])
            n_s += 1

    return new_svs

# 因为minimap2的比对可能会产生比对的断裂，导致原本的长变异消失，出现两个错误的短变异
def correct_2svs_align_crack(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        first_end = int(tmp_svs[0][n_s][2]) + abs(int(tmp_svs[0][n_s][3]))
        first_len = abs(int(tmp_svs[0][n_s][3]))
        for j in range(n_s+1,len(tmp_svs[0])) :
            second_end = int(tmp_svs[0][j][2]) + abs(int(tmp_svs[0][j][3]))
            second_len = abs(int(tmp_svs[0][j][3]))
            if abs(first_end - second_end) <= 1 :
                find_first = -1
                find_second = -1
                for ass_sv_i in range(len(member_svs)-1) :
                    if member_svs[ass_sv_i][0] == tmp_svs[0][n_s][1] and abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][n_s][2])) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 10 :
                        find_first = ass_sv_i
                    if member_svs[ass_sv_i][0] == tmp_svs[0][j][1] and abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][j][2])) <= 5 and abs(int(tmp_svs[0][j][3])-member_svs[ass_sv_i][3]) <= 10 :
                        find_second = ass_sv_i
                if find_first != -1 and find_second == -1 :
                    find_flag = -1
                    for k in range(n_s) :
                        if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][k][2])) <= 100 :
                            merge_len = second_len + abs(int(tmp_svs[0][k][3]))
                            if abs(merge_len-first_len)/max(merge_len,first_len) <= 0.2 :
                                if tmp_svs[0][n_s][gt_index][0:3] == "0/0" :
                                    tmp_svs[0][n_s][qual_index] = "255"
                                    tmp_svs[0][n_s][gt_index] = member_svs[find_first][-1]
                                    if member_svs[find_first][-1][4:7] == "1|1" :
                                        tmp_svs[0][n_s][gl_index] = "100,100,0,1,1,0,0,0,-10"
                                    elif member_svs[find_first][-1][4:7] == "0|1" :
                                        tmp_svs[0][n_s][gl_index] = "100,0,100,0,1,0,0,0,-10"
                                    elif member_svs[find_first][-1][4:7] == "1|0" :
                                        tmp_svs[0][n_s][gl_index] = "100,0,100,1,0,0,0,0,-10"
                                if tmp_svs[0][j][0:3] != "0/0" :
                                    tmp_svs[0][j][qual_index] = "0"
                                    tmp_svs[0][j][gt_index] = "0/0:0|0"
                                    tmp_svs[0][j][gl_index] = "0,100,100,0,0,0,0,0,-10"
                                if tmp_svs[0][k][0:3] != "0/0" :
                                    tmp_svs[0][k][qual_index] = "0"
                                    tmp_svs[0][k][gt_index] = "0/0:0|0"
                                    tmp_svs[0][k][gl_index] = "0,100,100,0,0,0,0,0,-10"
                else :
                    break
    return tmp_svs


def verify_allele_svs_byassembly(chr, tmp_svs, member_svs, family_member_ls) :
    new_svs = [[] for x in family_member_ls]
    n_s = 0
    while(True) :
        if n_s == len(tmp_svs[0])-1 :
            new_svs[0].append(tmp_svs[0][n_s])
            new_svs[1].append(tmp_svs[1][n_s])
            new_svs[2].append(tmp_svs[2][n_s])
        if n_s >= len(tmp_svs[0])-1 :
            break
        find_flag = -1
        if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+1][2])) == 0 :
            local_svs = [[] for x in family_member_ls]
            find_flag = 1
            allele_range_len = 2
            if abs(int(tmp_svs[0][n_s][3])) >= 40 :
                local_svs[0].append(tmp_svs[0][n_s])
                local_svs[1].append(tmp_svs[1][n_s])
                local_svs[2].append(tmp_svs[2][n_s])
            if abs(int(tmp_svs[0][n_s+1][3])) >= 40 :
                local_svs[0].append(tmp_svs[0][n_s+1])
                local_svs[1].append(tmp_svs[1][n_s+1])
                local_svs[2].append(tmp_svs[2][n_s+1])
            while(True) :
                if n_s+allele_range_len >= len(tmp_svs[0]) :
                    break
                if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+allele_range_len][2])) == 0 :
                    if abs(int(tmp_svs[0][n_s+allele_range_len][3])) >= 50 :
                        local_svs[0].append(tmp_svs[0][n_s+allele_range_len])
                        local_svs[1].append(tmp_svs[1][n_s+allele_range_len])
                        local_svs[2].append(tmp_svs[2][n_s+allele_range_len])
                    allele_range_len += 1
                else :
                    break
            n_s = n_s + allele_range_len
            if len(local_svs[0]) != 2 :
                new_svs[0] = new_svs[0] + tmp_svs[0][n_s-allele_range_len:n_s]
                new_svs[1] = new_svs[1] + tmp_svs[1][n_s-allele_range_len:n_s]
                new_svs[2] = new_svs[2] + tmp_svs[2][n_s-allele_range_len:n_s]
                continue
            if abs(int(local_svs[0][0][3])) < 100 or abs(int(local_svs[0][1][3])) < 100 or abs(int(local_svs[0][0][3])-int(local_svs[0][1][3]))/max(abs(int(local_svs[0][0][3])),abs(int(local_svs[0][1][3]))) < 0.3 :
                new_svs[0] = new_svs[0] + tmp_svs[0][n_s-allele_range_len:n_s]
                new_svs[1] = new_svs[1] + tmp_svs[1][n_s-allele_range_len:n_s]
                new_svs[2] = new_svs[2] + tmp_svs[2][n_s-allele_range_len:n_s]
                continue
            
            for ass_sv_i in range(len(member_svs)) :
                if abs(int(local_svs[0][0][2])-member_svs[ass_sv_i][1]) == 0 :
                    if ass_sv_i == len(member_svs)-1 or abs(int(local_svs[0][1][2])-member_svs[ass_sv_i+1][1]) > 5 :
                        keep_n_s = 0
                        dele_n_s = 1
                        if member_svs[ass_sv_i][0] == local_svs[0][0][1] and abs(member_svs[ass_sv_i][3]-int(local_svs[0][0][3])) < 10 :
                            keep_n_s = 0
                            dele_n_s = 1
                        elif member_svs[ass_sv_i][0] == local_svs[0][1][1] and abs(member_svs[ass_sv_i][3]-int(local_svs[0][1][3])) < 10 :
                            keep_n_s = 1
                            dele_n_s = 0
                        if int(local_svs[0][keep_n_s][gt_index][0])+int(local_svs[0][keep_n_s][gt_index][2]) != int(member_svs[ass_sv_i][-1][0]) + int(member_svs[ass_sv_i][-1][2]) :
                            if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                                gl_str = "100,100,0,1,1,0,0,0,-4"
                            elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                                gl_str = "100,0,100,0,1,0,0,0,-4"
                            elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                                gl_str = "100,0,100,1,0,0,0,0,-4"
                            local_svs[0][keep_n_s] = [chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""]
                        new_svs[0].append(local_svs[0][keep_n_s])
                        new_svs[1].append(local_svs[1][keep_n_s])
                        new_svs[2].append(local_svs[2][keep_n_s])
                        break
        if find_flag == -1 :
            new_svs[0].append(tmp_svs[0][n_s])
            new_svs[1].append(tmp_svs[1][n_s])
            new_svs[2].append(tmp_svs[2][n_s])
            n_s += 1
    return new_svs

# 按照ass的检测结果修正现有sv的pos
def refine_svs_pos(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        for ass_sv_i in range(len(member_svs)) :
            if 0 < abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 3 and tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 3 :
                is_allele = False
                if n_s-1 >= 0 and tmp_svs[0][n_s][2] == tmp_svs[0][n_s-1][2] :
                    is_allele = True
                if n_s+1 <= len(tmp_svs[0])-1 and tmp_svs[0][n_s][2] == tmp_svs[0][n_s+1][2] :
                    is_allele = True
                if is_allele == False :
                    tmp_svs[0][n_s][gl_index] = ",".join(tmp_svs[0][n_s][gl_index].split(",")[:-1] + ["-12+"+tmp_svs[0][n_s][2]+"+"+tmp_svs[0][n_s][3]])
                    tmp_svs[0][n_s][2] = str(member_svs[ass_sv_i][1])
                    tmp_svs[1][n_s][gl_index] = ",".join(tmp_svs[1][n_s][gl_index].split(",")[:-1] + ["-12"])
                    tmp_svs[1][n_s][2] = str(member_svs[ass_sv_i][1])
                    tmp_svs[2][n_s][gl_index] = ",".join(tmp_svs[2][n_s][gl_index].split(",")[:-1] + ["-12"])
                    tmp_svs[2][n_s][2] = str(member_svs[ass_sv_i][1])

    return tmp_svs

# 两个位置相近的杂合，minimap2的比对和cutesvTrio的聚类会是他们被错误的识别为纯和，需要按照ass的结果将两者分开
def split_error_homo_svs(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    new_svs = [[] for x in family_member_ls]
    for n_s in range(len(tmp_svs[0])) :
        find_flag = -1
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14"] :
            pass
        elif tmp_svs[0][n_s][gt_index][0:3] == "1/1" :
            for ass_sv_i in range(len(member_svs)-1) :
                if abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 50 and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i+1][1]) <= 50 and tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] == member_svs[ass_sv_i+1][0]:
                    if (member_svs[ass_sv_i][-1][4:7] == "0|1" and member_svs[ass_sv_i+1][-1][4:7] == "1|0") or (member_svs[ass_sv_i][-1][4:7] == "1|0" and member_svs[ass_sv_i+1][-1][4:7] == "0|1") :
                        if abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i][3])) < 0.3 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i+1][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i+1][3])) < 0.3 :
                            ass_pos_1 = member_svs[ass_sv_i][1]
                            ass_pos_2 = member_svs[ass_sv_i+1][1]
                            ass_num_1 = 0
                            ass_num_2 = 0
                            for read_p in tmp_svs[0][n_s][15].split(",") :
                                if abs(int(read_p)-ass_pos_1) <= 1 :
                                    ass_num_1 += 1
                                if abs(int(read_p)-ass_pos_2) <= 1 :
                                    ass_num_2 += 1
                            if ass_num_1 >= minimum_support_reads_list[0] and ass_num_2 >= minimum_support_reads_list[0] :
                                find_flag = 1
                                if member_svs[ass_sv_i][-1][4:7] == "0|1" :
                                    gl_str = "100,0,100,0,1,0,0,0,-6"
                                elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                                    gl_str = "100,0,100,1,0,0,0,0,-6"
                                new_svs[0].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
                                new_svs[1].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
                                new_svs[2].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
                                if member_svs[ass_sv_i+1][-1][4:7] == "0|1" :
                                    gl_str = "100,0,100,0,1,0,0,0,-6"
                                elif member_svs[ass_sv_i+1][-1][4:7] == "1|0" :
                                    gl_str = "100,0,100,1,0,0,0,0,-6"
                                new_svs[0].append([chr,member_svs[ass_sv_i+1][0],str(member_svs[ass_sv_i+1][1]),str(member_svs[ass_sv_i+1][3]),"0","0","0","0",member_svs[ass_sv_i+1][-1],gl_str,"996","255","",member_svs[ass_sv_i+1][4],"","","",""])
                                new_svs[1].append([chr,member_svs[ass_sv_i+1][0],str(member_svs[ass_sv_i+1][1]),str(member_svs[ass_sv_i+1][3]),"0","0","0","0",member_svs[ass_sv_i+1][-1],gl_str,"996","255","",member_svs[ass_sv_i+1][4],"","","",""])
                                new_svs[2].append([chr,member_svs[ass_sv_i+1][0],str(member_svs[ass_sv_i+1][1]),str(member_svs[ass_sv_i+1][3]),"0","0","0","0",member_svs[ass_sv_i+1][-1],gl_str,"996","255","",member_svs[ass_sv_i+1][4],"","","",""])
                                if int(tmp_svs[0][n_s][2]) != int(new_svs[0][-2][2]) and int(tmp_svs[0][n_s][2]) != int(new_svs[0][-1][2]) and int(tmp_svs[0][n_s][3]) != int(new_svs[0][-2][3]) and int(tmp_svs[0][n_s][3]) != int(new_svs[0][-1][3]) :
                                    new_svs[0].append(tmp_svs[0][n_s])
                                    new_svs[1].append(tmp_svs[1][n_s])
                                    new_svs[2].append(tmp_svs[2][n_s])
                                    new_svs[0][-1][qual_index] = "0"
                                    new_svs[0][-1][gt_index] = "0/0:0|0"
                                    new_svs[0][-1][gl_index] = "0,100,100,0,0,0,0,0,-6"
                                    new_svs[1][-1][qual_index] = "0"
                                    new_svs[1][-1][gt_index] = "0/0:0|0"
                                    new_svs[1][-1][gl_index] = "0,100,100,0,0,0,0,0,-6"
                                    new_svs[2][-1][qual_index] = "0"
                                    new_svs[2][-1][gt_index] = "0/0:0|0"
                                    new_svs[2][-1][gl_index] = "0,100,100,0,0,0,0,0,-6"
                            break
        if find_flag == -1:
            new_svs[0].append(tmp_svs[0][n_s])
            new_svs[1].append(tmp_svs[1][n_s])
            new_svs[2].append(tmp_svs[2][n_s])
    return new_svs

# 类等位基因中，如果两个sv的长度都在50附近，就是用assembly做筛选
def screening_allele_boundary_len(chr, tmp_svs, member_svs, family_member_ls) :
    n_s = 0
    while(True) :
        if n_s >= len(tmp_svs[0])-1 :
            break
        if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+1][2])) <= 1 :
            if abs(int(tmp_svs[0][n_s][3])) <= 55 and abs(int(tmp_svs[0][n_s+1][3])) <= 55 :
                if int(tmp_svs[0][n_s][gt_index][0]) + int(tmp_svs[0][n_s][gt_index][2]) > 0 :
                    find_flag = -1
                    for ass_sv_i in range(len(member_svs)) :
                        if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 1 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 :
                            find_flag = 1
                            break
                    if find_flag == -1 :
                        tmp_svs[0][n_s][qual_index] = "0"
                        tmp_svs[0][n_s][gt_index] = "0/0:0|0"
                        tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-6"
                if int(tmp_svs[0][n_s+1][gt_index][0]) + int(tmp_svs[0][n_s+1][gt_index][2]) > 0 :
                    find_flag = -1
                    for ass_sv_i in range(len(member_svs)) :
                        if tmp_svs[0][n_s+1][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s+1][2])-member_svs[ass_sv_i][1]) <= 1 and abs(int(tmp_svs[0][n_s+1][3])-member_svs[ass_sv_i][3]) <= 5 :
                            find_flag = 1
                            break
                    if find_flag == -1 :
                        tmp_svs[0][n_s+1][qual_index] = "0"
                        tmp_svs[0][n_s+1][gt_index] = "0/0:0|0"
                        tmp_svs[0][n_s+1][gl_index] = "0,100,100,0,0,0,0,0,-6"
            n_s += 2
        else :
            n_s += 1
    return tmp_svs

# 针对家庭纠错方法种后两类补充的sv，使用assembly重新验证其正确性
def revalidate_trio_correction(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    for n_s in range(len(tmp_svs[0])) : 
        gl_split = tmp_svs[0][n_s][gl_index].split(",")
        if gl_split[8] in ["-3"] and int(float(gl_split[7])) <= minimum_support_reads_list[0] and abs(int(tmp_svs[0][n_s][3])) <= 175 :
            find_flag = -1
            for ass_sv_i in range(len(member_svs)) :
                if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 10 :
                    find_flag = 1
                    break
            if find_flag == -1 :
                tmp_svs[0][n_s][qual_index] = "0"
                tmp_svs[0][n_s][gt_index] = "0/0:0|0"
                tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-7"
    return tmp_svs

# 使用assembly验证denovo，尤其是第一类denovo的正确性
def revalidate_denovo(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    for n_s in range(len(tmp_svs[0])) : 
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14"] :
            continue
        family_gt_ls = []
        for j in range(len(tmp_svs)) :
            family_gt_ls.append([s for s in tmp_svs[j][n_s][gt_index][0:3]])
        if len(family_gt_ls) == 2 :
            if int(family_gt_ls[0][0]) + int(family_gt_ls[0][2]) in [0] :
                continue
            elif int(family_gt_ls[1][0]) + int(family_gt_ls[1][2]) != 0 :
                continue
        if len(family_gt_ls) == 3 :
            if int(family_gt_ls[0][0]) + int(family_gt_ls[0][2]) in [0] :
                continue
            elif int(family_gt_ls[1][0]) + int(family_gt_ls[1][2]) != 0 or int(family_gt_ls[2][0]) + int(family_gt_ls[2][2]) != 0:
                continue
        if abs(int(tmp_svs[0][n_s][3])) > 150 or float(tmp_svs[0][n_s][qual_index]) > 20 :
            continue
        assem_flag = -1
        if n_s-1 >= 0 and (abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s-1][2])) <= 10 or abs(int(tmp_svs[0][n_s][3])-int(tmp_svs[0][n_s-1][3])) <= 25) and tmp_svs[0][n_s-1][gt_index][0:3] != "0/0":
            assem_flag = 1
        elif n_s+1 <= len(tmp_svs[0])-1 and (abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+1][2])) <= 10 or abs(int(tmp_svs[0][n_s][3])-int(tmp_svs[0][n_s+1][3])) <= 25) and tmp_svs[0][n_s+1][gt_index][0:3] != "0/0":
            assem_flag = 1
        if assem_flag == -1 :
            continue
        find_flag = -1
        for ass_sv_i in range(len(member_svs)) :
            if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 10 :
                find_flag = 1
                break
        if find_flag == -1 :
            tmp_svs[0][n_s][qual_index] = "0"
            tmp_svs[0][n_s][gt_index] = "0/0:0|0"
            tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-8"
    return tmp_svs

# 三重及以上等位基因中，如果assembly和原本都只支持的一个等位基因，但是不一致，那么修改为assembly支持的
# 三重及以上等位基因中，如果assembly不支持，但是原本支持一个，那去除该sv
# 三重及以上等位基因中，如果原本检测记录中没有变异，但是assembly支持了一个，很有可能是minimap2的比对错误导致的
def screening_threeup_allele(chr, tmp_svs, member_svs, family_member_ls) :
    n_s = 0
    while(True) :
        if n_s >= len(tmp_svs[0])-1 :
            break
        if abs(int(tmp_svs[0][n_s][2])-int(tmp_svs[0][n_s+1][2])) <= 1 :
            threeup_ls = [n_s,n_s+1]
            up_range = 2
            while(True) :
                if n_s+up_range >= len(tmp_svs[0]) :
                    break
                if abs(int(tmp_svs[0][n_s+up_range-1][2])-int(tmp_svs[0][n_s+up_range][2])) <= 1 :
                    threeup_ls.append(n_s+up_range)
                    up_range += 1
                else :
                    break
            if up_range < 3 :
                n_s += up_range
                continue
            sv_num = 0
            for i in range(up_range) :
                if abs(int(tmp_svs[0][n_s+i][3])) >= 50 :
                    sv_num += 1
            if sv_num < 3 :
                n_s += up_range
                continue
            sv_record_ls = [0 for x in threeup_ls]
            ass_record_ls = [0 for x in threeup_ls]
            ass_sv_index = 0
            for i in range(up_range) :
                if tmp_svs[0][n_s+i][gt_index][0:3] in ["0/1","1/0","1/1"] and abs(int(tmp_svs[0][n_s+i][3])) >= 50:
                    sv_record_ls[i] = 1
            if sum(sv_record_ls) not in [0, 1] :
                n_s += up_range
                continue
            is_find = 0
            ass_sv_ls = []
            for i in range(up_range) :
                for ass_sv_i in range(len(member_svs)) :
                    if tmp_svs[0][n_s+i][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s+i][2])-member_svs[ass_sv_i][1]) <= 1 and abs(int(tmp_svs[0][n_s+i][3])-member_svs[ass_sv_i][3]) <= 5 and abs(member_svs[ass_sv_i][3]) >= 50:
                        ass_sv_ls.append(ass_sv_i)
                        is_find = 1
                        break
            if sum(sv_record_ls) == 1 and (len(set(ass_sv_ls)) == 0 or is_find == 0):
                if float(tmp_svs[0][n_s+sv_record_ls.index(1)][qual_index]) < 25 : 
                    tmp_svs[0][n_s+sv_record_ls.index(1)][qual_index] = "0"
                    tmp_svs[0][n_s+sv_record_ls.index(1)][gt_index] = "0/0:0|0"
                    tmp_svs[0][n_s+sv_record_ls.index(1)][gl_index] = "0,100,100,0,0,0,0,0,-9"
            elif len(set(ass_sv_ls)) == 1 :
                ass_sv_i = ass_sv_ls[0]
                len_dis = 10000
                for i in range(up_range) :
                    if tmp_svs[0][n_s+i][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s+i][2])-member_svs[ass_sv_i][1]) <= 1 and abs(int(tmp_svs[0][n_s+i][3])-member_svs[ass_sv_i][3]) <= 5 and abs(member_svs[ass_sv_i][3]) >= 50:
                        if abs(int(tmp_svs[0][n_s+i][3])-member_svs[ass_sv_i][3]) < len_dis :
                            ass_sv_index = i
                            len_dis = abs(int(tmp_svs[0][n_s+i][3])-member_svs[ass_sv_i][3])
                ass_record_ls[ass_sv_index] = 1
                if sum(sv_record_ls) == 1 :
                    if sv_record_ls.index(1) != ass_record_ls.index(1) :
                        tmp_svs[0][n_s+sv_record_ls.index(1)][qual_index] = "0"
                        tmp_svs[0][n_s+sv_record_ls.index(1)][gt_index] = "0/0:0|0"
                        tmp_svs[0][n_s+sv_record_ls.index(1)][gl_index] = "0,100,100,0,0,0,0,0,-9"
                        tmp_svs[0][n_s+ass_sv_index][gt_index] = member_svs[ass_sv_i][-1]
                        tmp_svs[0][n_s+ass_sv_index][qual_index] = "255"
                        if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                            gl_str = "100,100,0,1,1,0,0,0,-9"
                        elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                            gl_str = "100,0,100,0,1,0,0,0,-9"
                        elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                            gl_str = "100,0,100,1,0,0,0,0,-9"
                        tmp_svs[0][n_s+ass_sv_index][gl_index] = gl_str
                elif sum(sv_record_ls) == 0 :
                    tmp_svs[0][n_s+ass_sv_index][qual_index] = "255"
                    tmp_svs[0][n_s+ass_sv_index][gt_index] = member_svs[ass_sv_i][-1]
                    if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                        gl_str = "100,100,0,1,1,0,0,0,-9"
                    elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                        gl_str = "100,0,100,0,1,0,0,0,-9"
                    elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                        gl_str = "100,0,100,1,0,0,0,0,-9"
                    tmp_svs[0][n_s+ass_sv_index][gl_index] = gl_str
            n_s += up_range
        else :
            n_s += 1
    return tmp_svs

# 过滤前后紧密连接在一起的两个变异，使用assembly重新验证他们的正确性
# 效果不好，暂时不用
def screening_connected_svs(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        first_end = int(tmp_svs[0][n_s][2]) + abs(int(tmp_svs[0][n_s][3]))
        first_len = abs(int(tmp_svs[0][n_s][3]))
        for j in range(n_s+1,len(tmp_svs[0])) :
            second_len = abs(int(tmp_svs[0][j][3]))
            if abs(int(tmp_svs[0][j][2])-first_end) <= 2 :
                if j - n_s >= 2 :
                    continue
                find_first = -1
                find_second = -1
                for ass_sv_i in range(len(member_svs)-1) :
                    if member_svs[ass_sv_i][0] == tmp_svs[0][n_s][1] and abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][n_s][2])) <= 3 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 :
                        find_first = ass_sv_i
                    if member_svs[ass_sv_i][0] == tmp_svs[0][j][1] and abs(member_svs[ass_sv_i][1]-int(tmp_svs[0][j][2])) <= 3 and abs(int(tmp_svs[0][j][3])-member_svs[ass_sv_i][3]) <= 5 :
                        find_second = ass_sv_i
                if find_first == -1 and tmp_svs[0][n_s][0:3] != "0/0" :
                    tmp_svs[0][n_s][qual_index] = "0"
                    tmp_svs[0][n_s][gt_index] = "0/0:0|0"
                    tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-4"
                if find_second == -1 and tmp_svs[0][j][0:3] != "0/0" :
                    tmp_svs[0][j][qual_index] = "0"
                    tmp_svs[0][j][gt_index] = "0/0:0|0"
                    tmp_svs[0][j][gl_index] = "0,100,100,0,0,0,0,0,-4"
    return tmp_svs             

# 通过assembly中位置和长度几乎完全相同的记录来补充sv
def supple_samesv_bybyassembly(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        for ass_sv_i in range(len(member_svs)) :
            if abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 0 and tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 1 and abs(int(tmp_svs[0][n_s][3])) > 200 :
                if tmp_svs[0][n_s][gt_index][0:3] != "0/0" or tmp_svs[1][n_s][gt_index][0:3] != "0/0" or tmp_svs[2][n_s][gt_index][0:3] != "0/0" or abs(int(tmp_svs[0][n_s][3])) >= 250 :
                    if tmp_svs[0][n_s][gt_index][0:3] == "0/0" and member_svs[ass_sv_i][-1][0:3] != "0/0" :
                        tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1]
                        tmp_svs[0][n_s][qual_index] = "255"
                        if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                            gl_str = "100,100,0,1,1,0,0,0,-11"
                        elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                            gl_str = "100,0,100,0,1,0,0,0,-11"
                        elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                            gl_str = "100,0,100,1,0,0,0,0,-11"
                        tmp_svs[0][n_s][gl_index] = gl_str
    return tmp_svs

# 如果assembly支持的sv是一个左右完全没有原本检测变异支持的sv，那么反向证明了其可能的正确性
def reverse_supple_byunsupport_sv(chr, tmp_svs, member_svs, family_member_ls) :
    for ass_sv_i in range(len(member_svs)-1) :
        find_flag = -1
        for n_s in range(len(tmp_svs[0])) :
            if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 300 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i][3])) <= 0.6 :
                find_flag = 1
                break
        if find_flag == -1 :
            if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                gl_str = "100,100,0,1,1,0,0,0,-4"
            elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                gl_str = "100,0,100,0,1,0,0,0,-4"
            elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                gl_str = "100,0,100,1,0,0,0,0,-4"
            tmp_svs[0].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
            tmp_svs[1].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
            tmp_svs[2].append([chr,member_svs[ass_sv_i][0],str(member_svs[ass_sv_i][1]),str(member_svs[ass_sv_i][3]),"0","0","0","0",member_svs[ass_sv_i][-1],gl_str,"996","255","",member_svs[ass_sv_i][4],"","","",""])
    
    return tmp_svs

# 如果原检测结果中两个变异接近，并且长度几乎一致，而assembly只支持一个，两个sv的支持信号互斥，那么只保留assembly支持的
def screening_two_near_svs(chr, tmp_svs, member_svs, family_member_ls, minimum_support_reads_list) :
    for n_s in range(len(tmp_svs[0])) :
        for n_j in range(max(0,n_s-5),min(len(tmp_svs[0])-1,n_s+5)) :
            if n_s == n_j :
                continue
            if tmp_svs[0][n_s][gt_index][0:3] == "0/0" or tmp_svs[0][n_j][gt_index][0:3] == "0/0" :
                continue
            if int(tmp_svs[0][n_s][2]) - int(tmp_svs[0][n_j][2]) <= 100 and abs(int(tmp_svs[0][n_s][3])-int(tmp_svs[0][n_j][3])) <= 5 :
                find_first_index = find_second_index = -1
                for ass_sv_i in range(len(member_svs)) :
                    if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 and ass_sv_i != find_second_index:
                        find_first_index = ass_sv_i
                    if tmp_svs[0][n_j][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_j][2])-member_svs[ass_sv_i][1]) <= 5 and abs(int(tmp_svs[0][n_j][3])-member_svs[ass_sv_i][3]) <= 5 and ass_sv_i != find_first_index:
                        find_second_index = ass_sv_i
                if (find_first_index != -1 and find_second_index != -1) or (find_first_index == -1 and find_second_index == -1) :
                    continue
                if len(list(set(tmp_svs[0][n_s][12].split(",")) & set(tmp_svs[0][n_j][12].split(",")))) >= minimum_support_reads_list[0] :
                    continue
                if find_first_index == -1 :
                    tmp_svs[0][n_s][qual_index] = "0"
                    tmp_svs[0][n_s][gt_index] = "0/0:0|0"
                    tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-4"
                else :
                    tmp_svs[0][n_j][qual_index] = "0"
                    tmp_svs[0][n_j][gt_index] = "0/0:0|0"
                    tmp_svs[0][n_j][gl_index] = "0,100,100,0,0,0,0,0,-4"
                break
    return tmp_svs

# 使用lcr中确定的单倍型分支结果，重新验证纯和变异的正确性
def correct_homogt_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-0","-1","-2","-3","-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14","-15"] :
            continue
        for ass_sv_i in range(len(member_svs)) :
            if abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 3 and tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 :
                if tmp_svs[0][n_s][gt_index][0:3] == "1/1" and member_svs[ass_sv_i][-1][0:3] in ["0/1","1/0"] :
                    tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1]
                    tmp_svs[0][n_s][gl_index] = ",".join(tmp_svs[0][n_s][gl_index].split(",")[:-1] + ["-13"])
    return tmp_svs

# 使用lcr中确定的单倍型分支结果，重新验证杂合变异的正确性
def correct_hetegt_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-0","-1","-2","-3","-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14","-15"] :
            continue
        for ass_sv_i in range(len(member_svs)) :
            if abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 3 and tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 :
                if tmp_svs[0][n_s][gt_index][0:3] in ["0/1","1/0"] and member_svs[ass_sv_i][-1][0:3] in ["1/1"] :
                    tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1]
                    tmp_svs[0][n_s][gl_index] = ",".join(tmp_svs[0][n_s][gl_index].split(",")[:-1] + ["-13"])

    return tmp_svs

# 针对19号染色体进行针对CMRG的特异性优化
def correct_19_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14","-15"] :
            pass
        elif abs(int(tmp_svs[0][n_s][3])) not in range(80,105) and abs(int(tmp_svs[0][n_s][3])) not in range(350,360) :
            pass
        elif tmp_svs[0][n_s][gt_index][0:3] == "0/0" :
            pass
        find_flag = -1
        for ass_sv_i in range(len(member_svs)) :
            if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 5 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 5 :
                find_flag = 1
                break
        if find_flag == -1 :
            tmp_svs[0][n_s][qual_index] = "0"
            tmp_svs[0][n_s][gt_index] = "0/0:0|0"
            tmp_svs[0][n_s][gl_index] = "0,100,100,0,0,0,0,0,-15"
    return tmp_svs

# 针对长并整齐的信号进行针对CMRG的特异性补充
def supple_long_neat_within_lcr(chr, tmp_svs, member_svs, family_member_ls) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gl_index].split(",")[-1] in ["-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14","-15"] :
            pass
        if tmp_svs[0][n_s][gt_index][0:3] != "0/0" :
            continue
        if abs(int(tmp_svs[0][n_s][3])) < 1000 :
            continue
        is_supple = False
        sv_len_ls = tmp_svs[0][n_s][14].split(",")
        for member_i in range(len(family_member_ls)) :
            if abs(int(sv_len_ls[member_i]) - int(tmp_svs[0][n_s][3])) <= 10 :
                if tmp_svs[member_i][n_s][15] != "" :
                    start_pos_ls = [int(x) for x in tmp_svs[member_i][n_s][15].split(",")]
                    if max(start_pos_ls)-min(start_pos_ls) <= 25 :
                        is_supple = True
                        break
        if not is_supple :
            continue
        for ass_sv_i in range(len(member_svs)) :
            if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 25 and (abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3]) <= 50 or abs(member_svs[ass_sv_i][3]) >= 1000) :
                tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1]
                tmp_svs[0][n_s][qual_index] = "255"
                if member_svs[ass_sv_i][-1][4:7] == "1|1" :
                    gl_str = "100,100,0,1,1,0,0,0,-15"
                elif member_svs[ass_sv_i][-1][4:7] == "0|1" :
                    gl_str = "100,0,100,0,1,0,0,0,-15"
                elif member_svs[ass_sv_i][-1][4:7] == "1|0" :
                    gl_str = "100,0,100,1,0,0,0,0,-15"
                tmp_svs[0][n_s][gl_index] = gl_str

    return tmp_svs

# 过滤列表中的所有纯和变异，随机分配纯和变异的所有read，重复多次，如果中间出现了一定次数的杂合结果，那么将其修改为杂合
# 使用随机的方法将范围内的所有read随机的分配到两个hap中
def correct_homogt_allsv_random(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gt_index][0:3] == "1/1" and abs(int(tmp_svs[0][n_s][3])) >= 50 :
            lcr_sta_pos = int(tmp_svs[0][n_s][2])
            lcr_end_pos = int(tmp_svs[0][n_s][2])+abs(int(tmp_svs[0][n_s][3]))
            ref_seq = ref_aligner.seq(chr, max(0,lcr_sta_pos-5*read_extension_scope), lcr_end_pos+5*read_extension_scope)
            local_ref_aligner = mp.Aligner(seq=ref_seq, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
            try_num = 15
            hete_gt_ls = []
            for x in range(try_num) :
                child_hap1_reads = []
                child_hap2_reads = []
                child_reads = child_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
                for child_read in child_reads :
                    if random.random() > 0.5 :
                        child_hap1_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                    else :
                        child_hap2_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                consensus_hap1 = msa_consensus_for_cluster(child_hap1_reads)
                consensus_hap2 = msa_consensus_for_cluster(child_hap2_reads)
                if consensus_hap1 is not None :
                    consensus_hap1 = consensus_hap1[0]
                if consensus_hap2 is not None :
                    consensus_hap2 = consensus_hap2[0]
                hap1_svs = hap2_svs = None
                if consensus_hap1 is not None :
                    for hit in local_ref_aligner.map(seq=consensus_hap1): # traverse alignments
                        if not hit.is_primary :
                            continue
                        hap1_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap1, ref_aligner)
                        break
                if consensus_hap2 is not None :
                    for hit in local_ref_aligner.map(seq=consensus_hap2): # traverse alignments
                        if not hit.is_primary :
                            continue
                        hap2_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap2, ref_aligner)
                        break
                if hap1_svs is None and hap2_svs is None :
                    continue
                if hap1_svs is None :
                    hap1_svs = []
                if hap2_svs is None :
                    hap2_svs = []
                if len(hap1_svs) == len(hap2_svs) == 0 :
                    continue
                member_svs = []
                for h_1 in hap1_svs :
                    member_svs.append(h_1[:]+["1/0:1|0"])
                for h_2 in hap2_svs :
                    homo_flag = False
                    for h_1 in member_svs :
                        if h_1[0] == h_2[0] and abs(h_1[1]-h_2[1]) <= 10 and abs(h_1[3]-h_2[3]) <= 10 :
                            h_1[-1] = "1/1:1|1"
                            homo_flag = True
                            break
                    if not homo_flag :
                        member_svs.append(h_2[:]+["1/0:0|1"])
                for ass_sv_i in range(len(member_svs)) :
                    if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 10 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i][3])) <= 0.2 :
                        if member_svs[ass_sv_i][-1][0:3] in ["0/1","1/0"] :
                            hete_gt_ls.append(member_svs[ass_sv_i][-1])
                            break
            if len(hete_gt_ls) >= 1 :
                tmp_svs[0][n_s][gt_index] = hete_gt_ls[0]
                if hete_gt_ls[0][4:7] == "0|1" :
                    tmp_svs[0][n_s][gl_index] = "100,0,100,0,1,0,0,0,-4"
                elif hete_gt_ls[0][4:7] == "1|0" :
                    tmp_svs[0][n_s][gl_index] = "100,0,100,1,0,0,0,0,-4"
                tmp_svs[0][n_s][qual_index] = "128"

    return tmp_svs


# 过滤列表中的所有纯和变异，随机分配纯和变异的所有read，重复多次，如果中间出现了一定次数的杂合结果，那么将其修改为杂合
# 将支持read分配给一个hap，范围内的其余read分配给另一个hap
def correct_homogt_allsv_support(path, chr, tmp_svs, fa_hap_ls, mo_hap_ls, homo_range, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list, family_assembly_reads, similarity_supplement_threshold, assembly_accelerate) :
    start_time = time.time()
    for n_s in range(len(tmp_svs[0])) :
        if n_s < homo_range[0] or n_s > homo_range[1] :
            continue
        if tmp_svs[0][n_s][gt_index][0:3] == "1/1" and abs(int(tmp_svs[0][n_s][3])) >= 50 :
            lcr_sta_pos = int(tmp_svs[0][n_s][2])
            lcr_end_pos = int(tmp_svs[0][n_s][2])+abs(int(tmp_svs[0][n_s][3]))
            ref_seq = ref_aligner.seq(chr, max(0,lcr_sta_pos-5*read_extension_scope), lcr_end_pos+5*read_extension_scope)
            local_ref_aligner = mp.Aligner(seq=ref_seq, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
            child_hap1_reads = []
            child_hap2_reads = []
            hap1_query_qualities_ls = []
            hap1_mapping_quality_ls = []
            hap2_query_qualities_ls = []
            hap2_mapping_quality_ls = []
            child_reads = child_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
            for child_read in child_reads :
                child_read_name = "M1/1/"+child_read.query_name
                if child_read_name in tmp_svs[0][n_s][12] :
                    child_hap1_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                    hap1_query_qualities_ls.append(sum(child_read.query_qualities)/len(child_read.query_qualities))
                    hap1_mapping_quality_ls.append(child_read.mapping_quality)
                else :
                    child_hap2_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                    hap2_query_qualities_ls.append(sum(child_read.query_qualities)/len(child_read.query_qualities))
                    hap2_mapping_quality_ls.append(child_read.mapping_quality)
            
            if assembly_accelerate :
                hap1_reads_len = 0
                hap2_reads_len = 0
                for read_x in child_hap1_reads :
                    hap1_reads_len += len(read_x)
                for read_x in child_hap2_reads :
                    hap2_reads_len += len(read_x)
                if hap1_reads_len >= 500000 or hap2_reads_len >= 500000 or len(child_hap1_reads) >= 500 or len(child_hap2_reads) >= 500:
                    continue
            
            consensus_hap1 = msa_consensus_for_cluster(child_hap1_reads)
            consensus_hap2 = msa_consensus_for_cluster(child_hap2_reads)
            if consensus_hap1 is not None :
                consensus_hap1 = consensus_hap1[0]
            if consensus_hap2 is not None :
                consensus_hap2 = consensus_hap2[0]
            hap1_svs = hap2_svs = None
            if consensus_hap1 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap1): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap1_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap1, ref_aligner)
                    break
            if consensus_hap2 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap2): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap2_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap2, ref_aligner)
                    break
            if hap1_svs is None and hap2_svs is None :
                continue
            if hap1_svs is None :
                hap1_svs = []
            if hap2_svs is None :
                hap2_svs = []
            if len(hap1_svs) == len(hap2_svs) == 0 :
                continue
            member_svs = []
            for h_1 in hap1_svs :
                member_svs.append(h_1[:]+["1/0:1|0"])
            for h_2 in hap2_svs :
                homo_flag = False
                for h_1 in member_svs :
                    if h_1[0] == h_2[0] and abs(h_1[1]-h_2[1]) <= 10 and abs(h_1[3]-h_2[3]) <= 10 :
                        h_1[-1] = "1/1:1|1"
                        homo_flag = True
                        break
                if not homo_flag :
                    member_svs.append(h_2[:]+["1/0:0|1"])
            for ass_sv_i in range(len(member_svs)) :
                if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 10 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i][3])) <= 0.2 :
                    if member_svs[ass_sv_i][-1][0:3] in ["0/1","1/0"] :
                        if tmp_svs[0][n_s][1] == "DEL" : 
                            if int(float(tmp_svs[0][n_s][gl_index].split(",")[6])) >= floor(minimum_support_reads_list[0]/2) :
                                tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1] 
                                tmp_svs[0][n_s][gl_index] = ",".join(tmp_svs[0][n_s][gl_index].split(",")[:-1] + ["-14"])
                        else :
                            if int(tmp_svs[0][n_s][3]) > 300 :
                                if len(child_hap2_reads) > 4 * minimum_support_reads_list[0] :
                                    continue
                                if sum([x if x >= 75 else 0 for x in hap2_query_qualities_ls]) == 0 :
                                    continue
                                if len(hap2_query_qualities_ls) <= minimum_support_reads_list[0] :
                                    continue
                                tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1] 
                                tmp_svs[0][n_s][gl_index] = ",".join(tmp_svs[0][n_s][gl_index].split(",")[:-1] + ["-14"])
                        break
    logging.info("%s/%s/%s-%f"%(chr,tmp_svs[0][0][2],tmp_svs[0][-1][2],time.time()-start_time))
    return tmp_svs

# 将第二类denovo，或者无中生有的纯和结果提取出来，重新验证一遍
def correct_homogt_allsv_denovo(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gt_index][0:3] == "1/1" and abs(int(tmp_svs[0][n_s][3])) >= 50 :
            if tmp_svs[0][n_s][gl_index].split(",")[-1] == "-4" :
                continue
            if tmp_svs[1][n_s][gt_index][0:3] != "0/0" and tmp_svs[2][n_s][gt_index][0:3] != "0/0" :
                continue
            lcr_sta_pos = int(tmp_svs[0][n_s][2])
            lcr_end_pos = int(tmp_svs[0][n_s][2])+abs(int(tmp_svs[0][n_s][3]))
            ref_seq = ref_aligner.seq(chr, max(0,lcr_sta_pos-5*read_extension_scope), lcr_end_pos+5*read_extension_scope)
            local_ref_aligner = mp.Aligner(seq=ref_seq, k=remap_merge_k, w=remap_minimizer_window, preset="map-pb")
            child_hap1_reads = []
            child_hap2_reads = []
            hap1_query_qualities_ls = []
            hap1_mapping_quality_ls = []
            hap2_query_qualities_ls = []
            hap2_mapping_quality_ls = []
            child_reads = child_bam.fetch(chr, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope)
            for child_read in child_reads :
                child_read_name = "M1/1/"+child_read.query_name
                if child_read_name in tmp_svs[0][n_s][12] :
                    child_hap1_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                    hap1_query_qualities_ls.append(sum(child_read.query_qualities)/len(child_read.query_qualities))
                    hap1_mapping_quality_ls.append(child_read.mapping_quality)
                else :
                    child_hap2_reads.append(capture_reads_within_lcr(child_read, lcr_sta_pos-read_extension_scope, lcr_end_pos+read_extension_scope))
                    hap2_query_qualities_ls.append(sum(child_read.query_qualities)/len(child_read.query_qualities))
                    hap2_mapping_quality_ls.append(child_read.mapping_quality)
            consensus_hap1 = msa_consensus_for_cluster(child_hap1_reads)
            consensus_hap2 = msa_consensus_for_cluster(child_hap2_reads)
            if consensus_hap1 is not None :
                consensus_hap1 = consensus_hap1[0]
            if consensus_hap2 is not None :
                consensus_hap2 = consensus_hap2[0]
            hap1_svs = hap2_svs = None
            if consensus_hap1 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap1): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap1_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap1, ref_aligner)
                    break
            if consensus_hap2 is not None :
                for hit in local_ref_aligner.map(seq=consensus_hap2): # traverse alignments
                    if not hit.is_primary :
                        continue
                    hap2_svs = sv_from_cigar(chr, hit.cigar_str, hit.r_st, 10, max(0,lcr_sta_pos-5*read_extension_scope), consensus_hap2, ref_aligner)
                    break
            if hap1_svs is None and hap2_svs is None :
                continue
            if hap1_svs is None :
                hap1_svs = []
            if hap2_svs is None :
                hap2_svs = []
            if len(hap1_svs) == len(hap2_svs) == 0 :
                continue
            member_svs = []
            for h_1 in hap1_svs :
                member_svs.append(h_1[:]+["1/0:1|0"])
            for h_2 in hap2_svs :
                homo_flag = False
                for h_1 in member_svs :
                    if h_1[0] == h_2[0] and abs(h_1[1]-h_2[1]) <= 10 and abs(h_1[3]-h_2[3]) <= 10 :
                        h_1[-1] = "1/1:1|1"
                        homo_flag = True
                        break
                if not homo_flag :
                    member_svs.append(h_2[:]+["1/0:0|1"])
            for ass_sv_i in range(len(member_svs)) :
                if tmp_svs[0][n_s][1] == member_svs[ass_sv_i][0] and abs(int(tmp_svs[0][n_s][2])-member_svs[ass_sv_i][1]) <= 10 and abs(int(tmp_svs[0][n_s][3])-member_svs[ass_sv_i][3])/max(abs(int(tmp_svs[0][n_s][3])),abs(member_svs[ass_sv_i][3])) <= 0.2 :
                    if member_svs[ass_sv_i][-1][0:3] in ["0/1","1/0"] :
                        logging.info(hap1_query_qualities_ls)
                        logging.info(hap1_mapping_quality_ls)
                        logging.info(hap2_query_qualities_ls)
                        logging.info(hap2_mapping_quality_ls)
                        logging.info("%s/%s"%(str(len(child_hap1_reads)),str(len(child_hap2_reads))))
                        logging.info(tmp_svs[0][n_s][0:gl_index+1])
                        logging.info(fa_hap_ls[n_s])
                        logging.info(mo_hap_ls[n_s])
                        logging.info(member_svs[ass_sv_i])
                        tmp_svs[0][n_s][gt_index] = member_svs[ass_sv_i][-1] 
                        break

    return tmp_svs

# 将杂合结果提取出来，重新验证一遍
def correct_hetegt_allsv_denovo(chr, tmp_svs, fa_hap_ls, mo_hap_ls, child_bam, father_bam, mother_bam, ref_aligner, read_extension_scope, remap_merge_k, remap_minimizer_window, minimum_support_reads_list) :
    for n_s in range(len(tmp_svs[0])) :
        if tmp_svs[0][n_s][gt_index][0:3] in ["0/1","1/0"] and abs(int(tmp_svs[0][n_s][3])) >= 50 :
            if int(tmp_svs[0][n_s][2]) in [17584841,24195934,35859490,38585218,40933632,41140918,43871358,48513209,49016536,49153250,49672573,50933885] :
                logging.info(tmp_svs[0][n_s])
                logging.info(fa_hap_ls[n_s])
                logging.info(mo_hap_ls[n_s])

    return tmp_svs


def split_reference_chromosomes(temporary_dir, reference) :
    with open(reference, "r") as infile:
        records = SeqIO.parse(infile, "fasta")

        chromosomes = {}

        for record in records:
            chrom_name = record.id
            if chrom_name not in chromosomes:
                chromosomes[chrom_name] = []
            chromosomes[chrom_name].append(record)

        for chrom_name, chrom_records in chromosomes.items():
            output_file = f"{temporary_dir}/{chrom_name}.fasta"
            with open(output_file, "w") as outfile:
                SeqIO.write(chrom_records, outfile, "fasta")

def sv_from_cigar(chr, cigar_str, ref_start, sv_min_size, origin_ref_start, sequence, reference):
    c = Cigar(cigar_str)
    svs = []

    ref_pos = ref_start   
    read_pos = 0          

    for length, op in c.items():
        if op in ("M", "=", "X", "N"):
            ref_pos += length
            read_pos += length

        elif op == "I":
            if length >= sv_min_size:
                svs.append([
                    "INS", origin_ref_start + ref_pos, origin_ref_start + ref_pos + 1, length, sequence[read_pos:read_pos+length]
                ])
            read_pos += length

        elif op == "D":
            if length >= sv_min_size:
                svs.append([
                    "DEL", origin_ref_start + ref_pos, origin_ref_start + ref_pos + length, -1*length, reference.seq(chr, origin_ref_start + ref_pos, origin_ref_start + ref_pos + length)
                ])
            ref_pos += length

        elif op in ("S", "H"):
            read_pos += length

        elif op == "P":
            pass

        else:
            raise ValueError(f"遇到未处理的 CIGAR 操作符: {op}")

    return svs

# 按照同类sv的位置和长度的临近关系，去除相近的相似sv
def remove_redundant_sv(chr, sv_ls) :
    filtered_sv_ls = []
    for s_i in range(len(sv_ls)) :
        redundant_flag = False
        for s_j in range(max(0,s_i-5),s_i) :
            if sv_ls[s_i][1] == sv_ls[s_j][1] and abs(int(sv_ls[s_i][2])-int(sv_ls[s_j][2])) <= 5 and abs(int(sv_ls[s_i][3])-int(sv_ls[s_j][3])) <= 5 :
                if sv_ls[s_i][8][0:3] in ["1/1","1/0","0/1"] and sv_ls[s_j][8][0:3] in ["1/1","1/0","0/1"] :
                    redundant_flag = True
                    break
        if not redundant_flag :
            if sv_ls[s_i][8][0:3] in ["1/0","0/1"] :
                sv_ls[s_i][8] = "1/0" + sv_ls[s_i][8][3:]
        filtered_sv_ls.append(sv_ls[s_i])
    return filtered_sv_ls

# 去除相邻的完全相同的变异记录
def remove_redundant_samesv(chr, sv_ls, family_mode) :
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    filtered_sv_ls = [[] for x in family_member_ls]
    s_i = 0
    while(True) :
        if s_i == len(sv_ls[0])-1 :
            for x in range(len(family_member_ls)) :
                filtered_sv_ls[x].append(sv_ls[x][s_i])
        if s_i >= len(sv_ls[0])-1 :
            break
        if sv_ls[0][s_i][1] == sv_ls[0][s_i+1][1] and int(sv_ls[0][s_i][2]) == int(sv_ls[0][s_i+1][2]) and int(sv_ls[0][s_i][3]) == int(sv_ls[0][s_i+1][3]) :
            if sv_ls[0][s_i][gl_index].split(",")[8].split("+")[0] in ["-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14"] :
                for x in range(len(family_member_ls)) :
                    filtered_sv_ls[x].append(sv_ls[x][s_i])
            else :
                for x in range(len(family_member_ls)) :
                    filtered_sv_ls[x].append(sv_ls[x][s_i+1])
            s_i += 2
        else :
            for x in range(len(family_member_ls)) :
                filtered_sv_ls[x].append(sv_ls[x][s_i])
            s_i += 1
    return filtered_sv_ls

# 去除位置偏差修正前的记录
def remove_redundant_pos(chr, sv_ls, family_mode) :
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    redundant_index = []
    for n_s in range(len(sv_ls[0])) :
        gl_assembly_ls = sv_ls[0][n_s][gl_index].split(",")[8].split("+")
        if len(gl_assembly_ls) > 1 :
            ass_pos = gl_assembly_ls[1]
            ass_len = gl_assembly_ls[2]
            for i in range(len(sv_ls[0])) :
                if sv_ls[0][i][2] == ass_pos and sv_ls[0][i][3] == ass_len and sv_ls[0][n_s][1] == sv_ls[0][i][1] and sv_ls[0][n_s][gt_index] == sv_ls[0][i][gt_index]:
                    redundant_index.append(i)
            sv_ls[0][n_s][gl_index] = ",".join(sv_ls[0][n_s][gl_index].split(",")[:-1] + ["-4"])
    new_ls = [[] for  x in family_member_ls]
    for n_s in range(len(sv_ls[0])) :
        if n_s not in redundant_index :
            new_ls[0].append(sv_ls[0][n_s])
            new_ls[1].append(sv_ls[1][n_s])
            new_ls[2].append(sv_ls[2][n_s])
    return new_ls


def run_assembly(args) :
    return redetect_nearby_sv(*args)


