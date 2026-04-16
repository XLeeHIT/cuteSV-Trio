from cuteSVTrio.cuteSVTrio_genotype import cal_CIPOS, overlap_cover, assign_gt, unsolvable_correction, modify_genotype, increase_sigs_through_pedigree, check_gt_consistent, inconformity_mendel_modify
from cuteSVTrio.cuteSVTrio_mendel import resolution_mendel,correction_mendel
import numpy as np
from math import ceil
import logging
import pickle
import copy
import time

'''
*******************************************
                TO DO LIST
*******************************************
    1. Identify DP with samfile pointer;
    2. Add CIPOS, CILEN and/or CIEND;
    3. Determine (IM)PRECISE type.
*******************************************

'''
def resolution_DEL(path, chr, read_count, threshold_gloab, max_cluster_bias,
                   minimum_support_reads_list, gt_round, remain_reads_ratio, merge_del_threshold, 
                   read_pos_interval, family_mode, performing_phasing):

    '''
    cluster DEL
    ********************************************************************************************
    path:	DEL.sigs
    chr:	chromosome id
    svtype:	<DEL>
    '''
    start_time = time.time()
    further_threshold_gloab = threshold_gloab*1.5
    max_cluster_len_bias = 2
    if remain_reads_ratio > 1:
        remain_reads_ratio = 1
    semi_del_cluster = list()
    semi_del_cluster.append([0,0,''])
    candidate_single_SV = list()
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    
    seqs = []
    for family_member in family_member_ls :
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"sigindex"), 'rb') as f:
            sigs_index=pickle.load(f)
            f.close()
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"DEL"), 'rb') as f:
            if chr not in sigs_index["DEL"].keys() :
                f.close()
                continue
            f.seek(sigs_index["DEL"][chr])
            f_seqs = pickle.load(f)
            seqs += f_seqs
            f.close()

    if len(seqs) == 0 :
        return
    seqs = sorted(seqs, key=lambda x: x[0])
    for seq in seqs :
        pos = int(seq[0])
        indel_len = int(seq[1])
        read_id = seq[2]
        if pos - semi_del_cluster[-1][0] > max_cluster_bias:
            if len(semi_del_cluster) >= read_count:
                if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
                    pass
                else:
                    del_tmp = sorted(semi_del_cluster, key=lambda x: x[1])
                    global_len = [i[1] for i in del_tmp]
                    global_pos = [i[0] for i in del_tmp]
                    global_first = np.mean(global_pos[0:ceil(len(global_pos)/2)])
                    global_last = np.mean(global_pos[ceil(len(global_pos)/2)-1:])
                    if abs(global_last-global_first) <= merge_del_threshold * 2:
                        semi_del_cluster_ls = [semi_del_cluster]
                    else :
                        DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = max_cluster_len_bias * np.mean(global_len)
                        semi_del_cluster_ls = [[del_tmp[0]]]
                        last_len = del_tmp[0][1]
                        for i in del_tmp[1:]:
                            if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
                                semi_del_cluster_ls.append([])
                            semi_del_cluster_ls[-1].append(i)
                            last_len = i[1]
                    
                    semi_del_cluster_ls = [sorted(i, key=lambda x: x[0]) for i in semi_del_cluster_ls]
                    del_tmp_ls = semi_del_cluster_ls
                    semi_del_cluster_ls = list()
                    for del_tmp in del_tmp_ls :
                        semi_del_cluster_ls.append([del_tmp[0]])
                        for semi_read in del_tmp[1:] :
                            if semi_read[0] - semi_del_cluster_ls[-1][-1][0] > max_cluster_bias :
                                if len(semi_del_cluster_ls[-1]) < read_count :
                                    semi_del_cluster_ls.remove(semi_del_cluster_ls[-1])
                                semi_del_cluster_ls.append([])
                            semi_del_cluster_ls[-1].append(semi_read)
                    
                    for semi_del_cluster_new in semi_del_cluster_ls :
                        generate_del_cluster(semi_del_cluster_new, 
                                            chr, 
                                            read_count, 
                                            threshold_gloab, 
                                            further_threshold_gloab,
                                            # threshold_local, 
                                            min(minimum_support_reads_list), 
                                            candidate_single_SV,
                                            gt_round,
                                            remain_reads_ratio,
                                            merge_del_threshold,
                                            family_mode,
                                            performing_phasing)
            semi_del_cluster = []
            semi_del_cluster.append([pos, indel_len, read_id])
        else:
            if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
                semi_del_cluster = []
                semi_del_cluster.append([pos, indel_len, read_id])
            else:
                semi_del_cluster.append([pos, indel_len, read_id])
    
    if len(semi_del_cluster) >= read_count:
        if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
            pass
        else:
            del_tmp = sorted(semi_del_cluster, key=lambda x: x[1])
            global_len = [i[1] for i in del_tmp]
            global_pos = [i[0] for i in del_tmp]
            global_first = np.mean(global_pos[0:ceil(len(global_pos)/2)])
            global_last = np.mean(global_pos[ceil(len(global_pos)/2)-1:])
            if abs(global_last-global_first) <= merge_del_threshold :
                semi_del_cluster_ls = [semi_del_cluster]
            else :
                DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = max_cluster_len_bias * np.mean(global_len)
                semi_del_cluster_ls = [[del_tmp[0]]]
                last_len = del_tmp[0][1]
                for i in del_tmp[1:]:
                    if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
                        semi_del_cluster_ls.append([])
                    semi_del_cluster_ls[-1].append(i)
                    last_len = i[1]
            
            semi_del_cluster_ls = [sorted(i, key=lambda x: x[0]) for i in semi_del_cluster_ls]
            del_tmp_ls = semi_del_cluster_ls
            semi_del_cluster_ls = list()
            for del_tmp in del_tmp_ls :
                semi_del_cluster_ls.append([del_tmp[0]])
                for semi_read in del_tmp[1:] :
                    if semi_read[0] - semi_del_cluster_ls[-1][-1][0] > max_cluster_bias :
                        if len(semi_del_cluster_ls[-1]) < read_count :
                            semi_del_cluster_ls.remove(semi_del_cluster_ls[-1])
                        semi_del_cluster_ls.append([])
                    semi_del_cluster_ls[-1].append(semi_read)

            for semi_del_cluster_new in semi_del_cluster_ls :
                generate_del_cluster(semi_del_cluster_new, 
                                    chr, 
                                    read_count, 
                                    threshold_gloab, 
                                    further_threshold_gloab,
                                    # threshold_local, 
                                    min(minimum_support_reads_list), 
                                    candidate_single_SV,
                                    gt_round,
                                    remain_reads_ratio,
                                    merge_del_threshold,
                                    family_mode,
                                    performing_phasing)

    candidate_single_SV_fam_ls = []
    candidate_single_SV_fam_ls.append(candidate_single_SV)
    for i in range(1,len(family_member_ls)) :
        candidate_single_SV_fam_ls.append([])
    for i in range(len(candidate_single_SV)) :
        fam_candidate_other = [[]]
        for j in range(len(family_member_ls)-1) :
            fam_candidate_other.append(['' for x in candidate_single_SV[i]])
        fam_reads = [[] for j in family_member_ls]
        reads_pos = [[] for j in family_member_ls]
        for j in range(len(candidate_single_SV[i][8])) :
            fam_reads[int(candidate_single_SV[i][8][j][3])-1].append(candidate_single_SV[i][8][j])
            if performing_phasing :
                reads_pos[int(candidate_single_SV[i][8][j][3])-1].append(candidate_single_SV[i][10][j])
        candidate_single_SV_fam_ls[0][i][8] = fam_reads[0]
        candidate_single_SV_fam_ls[0][i][10] = reads_pos[0]
        candidate_single_SV_fam_ls[0][i][4] = len(fam_reads[0])
        for j in range(1,len(family_member_ls)) :
            fam_candidate_other[j][8] = fam_reads[j]
            fam_candidate_other[j][10] = reads_pos[j]
            fam_candidate_other[j][4] = len(fam_reads[j])
            candidate_single_SV_fam_ls[j].append(fam_candidate_other[j])
    candidate_single_SV_gt_fam_ls = []
    for i in range(len(family_member_ls)) :
        family_member = family_member_ls[i]
        candidate_single_SV_gt_fam_ls.append(call_gt(path, chr, candidate_single_SV_fam_ls[i], candidate_single_SV_fam_ls[0], max_cluster_bias, 'DEL', family_mode, family_member, minimum_support_reads_list, performing_phasing))
    standard_list = 0
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) != 0 :
            standard_list = i
            break
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) == 0 :
            candidate_single_SV_gt_fam_ls[i] = copy.deepcopy(candidate_single_SV_gt_fam_ls[standard_list])
            for s_v in candidate_single_SV_gt_fam_ls[i] :
                s_v[8] = "0/0"
                s_v[9] = "0,100,100,0,0,0,0,0,0"
                s_v[10] = 996
                s_v[11] = 0
                s_v[12] = ""
                s_v[15] = ""
                s_v[16] = ""
                s_v[17] = ""
    
    increase_sigs_through_pedigree(candidate_single_SV_gt_fam_ls, 'DEL', minimum_support_reads_list, family_mode)    # 使用ESS算法
    
    unsolvable_correction(candidate_single_SV_gt_fam_ls, 'DEL', family_mode)
    correction_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)

    if not performing_phasing and family_mode == "M1" :
        resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)

    logging.info("Finished calling %s:%s:%f."%(chr, "DEL", time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)


def generate_del_cluster(semi_del_cluster, chr, read_count, 
                         threshold_gloab, further_threshold_gloab, minimum_support_reads, candidate_single_SV, 
                         gt_round, remain_reads_ratio, merge_del_threshold, family_mode, performing_phasing):

    is_print_flag = False 
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    read_tag = dict()
    for element in semi_del_cluster:
        if element[2] not in read_tag:
            read_tag[element[2]] = element
        else:
            if element[1] > read_tag[element[2]][1]:
                read_tag[element[2]] = element

    if len(read_tag) < read_count:
        return

    read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
    global_len = [i[1] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = threshold_gloab * np.mean(global_len)

    last_len = read_tag2SortedList[0][1]

    allele_collect = list()
    allele_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
        [read_tag2SortedList[0][2]]])

    for i in read_tag2SortedList[1:]:
        if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[]])

        allele_collect[-1][0].append(i[0])
        allele_collect[-1][1].append(i[1])
        allele_collect[-1][3].append(i[2])
        last_len = i[1]
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    if is_print_flag :
        logging.info(allele_collect)
    is_segment_flag = True
    while(is_segment_flag) :
        is_segment_flag = False
        further_allele_collect = list()
        for allele in allele_collect :
            DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = max(10,further_threshold_gloab * np.mean(allele[1]))
            further_allele_collect.append([[allele[0][0]],[allele[1][0]],[],[allele[3][0]]])
            first_len = allele[1][0]
            for i in range(1,len(allele[1])) :
                if allele[1][i] - first_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP :
                    is_segment_flag = True
                    further_allele_collect[-1][2].append(len(further_allele_collect[-1][0]))
                    further_allele_collect.append([[],[],[],[]])
                    first_len = allele[1][i]
                further_allele_collect[-1][0].append(allele[0][i])
                further_allele_collect[-1][1].append(allele[1][i])
                further_allele_collect[-1][3].append(allele[3][i])
            further_allele_collect[-1][2].append(len(further_allele_collect[-1][0]))
        allele_collect = further_allele_collect
    if is_print_flag :
        logging.info(allele_collect)
    allele_sort = sorted(allele_collect, key = lambda x:x[2])
    len_coefficient_var_threshold = 0.5
    for allele in allele_sort:
        if allele[2][0] >= minimum_support_reads:
            if np.std(allele[1])/np.mean(allele[1]) > len_coefficient_var_threshold :
                continue
            allele_list = list()
            var_list = list()
            remain_allele_num = max(int(remain_reads_ratio * allele[2][0]), 1)
            pos_mean = np.mean(allele[0])
            for i in range(len(allele[0])):
                var_list.append((abs(allele[0][i] - pos_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[0][var_list[i][1]])
            breakpointStart = np.mean(allele_list)
            search_threshold = allele_list[0]
            allele_list = list()
            var_list = list()
            len_mean = np.mean(allele[1])
            for i in range(len(allele[1])):
                var_list.append((abs(allele[1][i] - len_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[1][var_list[i][1]])
            signalLen = np.median(allele_list)
            CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
            signalLen_STD = np.std(allele[1])
            CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))

            fam_allele_len_ls = [[] for x in family_member_ls]
            for seq_i in range(len(allele[0])) :
                fam_allele_len_ls[int(allele[3][seq_i][3])-1].append(allele[1][seq_i])
            for fam_i in range(len(fam_allele_len_ls)) :
                if len(fam_allele_len_ls[fam_i]) == 0 :
                    fam_allele_len_ls[fam_i].append(0)
            sv_len_str = ""
            for fam_i in range(len(fam_allele_len_ls)) :
                sv_len_str = sv_len_str + str(int(np.mean(fam_allele_len_ls[fam_i]))) + ","

            if performing_phasing :
                candidate_single_SV.append([chr, 
                                            "DEL", 
                                            max(1,int(breakpointStart)), 
                                            int(-signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(search_threshold),
                                            allele[3],
                                            sv_len_str[:-1],
                                            allele[0]])
            else :
                candidate_single_SV.append([chr, 
                                            "DEL", 
                                            max(1,int(breakpointStart)), 
                                            int(-signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(search_threshold),
                                            allele[3],
                                            sv_len_str[:-1],
                                            ''])
    

def resolution_INS(path, chr, read_count, threshold_gloab, 
                   max_cluster_bias, minimum_support_reads_list, gt_round, remain_reads_ratio, merge_INS_threshold, 
                   read_pos_interval, family_mode, performing_phasing, all_ins_singnature_reads):
    
    '''
    cluster INS
    ********************************************************************************************
    path:	INS.sigs
    chr:	chromosome id
    svtype:	<INS>
    '''
    start_time = time.time()
    further_threshold_gloab = threshold_gloab*7.5
    max_cluster_len_bias = 1.2
    if remain_reads_ratio > 1:
        remain_reads_ratio = 1
    semi_ins_cluster = list()
    semi_ins_cluster.append([0,0,'',''])
    candidate_single_SV = list()
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

    seqs = []
    for family_member in family_member_ls :
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"sigindex"), 'rb') as f:
            sigs_index=pickle.load(f)
            f.close()
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"INS"), 'rb') as f:
            if chr not in sigs_index["INS"].keys() :
                f.close()
                continue
            f.seek(sigs_index["INS"][chr])
            f_seqs = pickle.load(f)
            seqs += f_seqs
            f.close()
    if len(seqs) == 0 :
        return
    seqs = sorted(seqs, key=lambda x: x[0])
    for seq in seqs:
        pos = int(seq[0])
        indel_len = int(seq[1])
        read_id = seq[2]
        try:
            ins_seq = seq[3]
        except:
            ins_seq = ''
        if pos - semi_ins_cluster[-1][0] > max_cluster_bias:
            if len(semi_ins_cluster) >= read_count:
                if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
                    pass
                else:
                    ins_tmp = sorted(semi_ins_cluster, key=lambda x: x[1])
                    global_len = [i[1] for i in ins_tmp]
                    global_pos = [i[0] for i in ins_tmp]
                    global_first = np.mean(global_pos[0:ceil(len(global_pos)/2)])
                    global_last = np.mean(global_pos[ceil(len(global_pos)/2)-1:])
                    if abs(global_last-global_first) <= merge_INS_threshold * 2 :
                        semi_ins_cluster_ls = [semi_ins_cluster]
                    else :
                        DISCRETE_THRESHOLD_LEN_CLUSTER_ins_TEMP = max_cluster_len_bias * np.mean(global_len)
                        semi_ins_cluster_ls = [[ins_tmp[0]]]
                        last_len = ins_tmp[0][1]
                        for i in ins_tmp[1:]:
                            if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_ins_TEMP:
                                semi_ins_cluster_ls.append([])
                            semi_ins_cluster_ls[-1].append(i)
                            last_len = i[1]
                    
                    semi_ins_cluster_ls = [sorted(i, key=lambda x: x[0]) for i in semi_ins_cluster_ls]
                    ins_tmp_ls = semi_ins_cluster_ls
                    semi_ins_cluster_ls = list()
                    for ins_tmp in ins_tmp_ls :
                        semi_ins_cluster_ls.append([ins_tmp[0]])
                        for semi_read in ins_tmp[1:] :
                            if semi_read[0] - semi_ins_cluster_ls[-1][-1][0] > max_cluster_bias :
                                if len(semi_ins_cluster_ls[-1]) < read_count :
                                    semi_ins_cluster_ls.remove(semi_ins_cluster_ls[-1])
                                semi_ins_cluster_ls.append([])
                            semi_ins_cluster_ls[-1].append(semi_read)
                    
                    for semi_ins_cluster_new in semi_ins_cluster_ls :
                        generate_ins_cluster(semi_ins_cluster_new, 
                                            chr, 
                                            read_count, 
                                            threshold_gloab, 
                                            further_threshold_gloab,
                                            # threshold_local, 
                                            min(minimum_support_reads_list), 
                                            candidate_single_SV,
                                            gt_round,
                                            remain_reads_ratio,
                                            merge_INS_threshold,
                                            family_mode,
                                            performing_phasing,
                                            all_ins_singnature_reads)
            semi_ins_cluster = []
            semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])
        else:
            if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
                semi_ins_cluster = []
                semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])
            else:
                semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])

    if len(semi_ins_cluster) >= read_count:
        if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
            pass
        else:
            ins_tmp = sorted(semi_ins_cluster, key=lambda x: x[1])
            global_len = [i[1] for i in ins_tmp]
            global_pos = [i[0] for i in ins_tmp]
            global_first = np.mean(global_pos[0:ceil(len(global_pos)/2)])
            global_last = np.mean(global_pos[ceil(len(global_pos)/2)-1:])
            if abs(global_last-global_first) <= merge_INS_threshold :
                semi_ins_cluster_ls = [semi_ins_cluster]
            else :
                DISCRETE_THRESHOLD_LEN_CLUSTER_ins_TEMP = max_cluster_len_bias * np.mean(global_len)
                semi_ins_cluster_ls = [[ins_tmp[0]]]
                last_len = ins_tmp[0][1]
                for i in ins_tmp[1:]:
                    if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_ins_TEMP:
                        semi_ins_cluster_ls.append([])
                    semi_ins_cluster_ls[-1].append(i)
                    last_len = i[1]
            
            semi_ins_cluster_ls = [sorted(i, key=lambda x: x[0]) for i in semi_ins_cluster_ls]
            ins_tmp_ls = semi_ins_cluster_ls
            semi_ins_cluster_ls = list()
            for ins_tmp in ins_tmp_ls :
                semi_ins_cluster_ls.append([ins_tmp[0]])
                for semi_read in ins_tmp[1:] :
                    if semi_read[0] - semi_ins_cluster_ls[-1][-1][0] > max_cluster_bias :
                        if len(semi_ins_cluster_ls[-1]) < read_count :
                            semi_ins_cluster_ls.remove(semi_ins_cluster_ls[-1])
                        semi_ins_cluster_ls.append([])
                    semi_ins_cluster_ls[-1].append(semi_read)
            
            for semi_ins_cluster_new in semi_ins_cluster_ls :
                generate_ins_cluster(semi_ins_cluster_new, 
                                    chr, 
                                    read_count, 
                                    threshold_gloab, 
                                    further_threshold_gloab,
                                    # threshold_local, 
                                    min(minimum_support_reads_list), 
                                    candidate_single_SV,
                                    gt_round,
                                    remain_reads_ratio,
                                    merge_INS_threshold,
                                    family_mode,
                                    performing_phasing,
                                    all_ins_singnature_reads)
    candidate_single_SV_fam_ls = []
    candidate_single_SV_fam_ls.append(candidate_single_SV)
    for i in range(1,len(family_member_ls)) :
        candidate_single_SV_fam_ls.append([])
    for i in range(len(candidate_single_SV)) :
        fam_candidate_other = [[]]
        for j in range(len(family_member_ls)-1) :
            fam_candidate_other.append(['' for x in candidate_single_SV[i]])
        fam_reads = [[] for j in family_member_ls]
        reads_pos = [[] for j in family_member_ls]
        for j in range(len(candidate_single_SV[i][8])) :
            fam_reads[int(candidate_single_SV[i][8][j][3])-1].append(candidate_single_SV[i][8][j])
            if performing_phasing :
                reads_pos[int(candidate_single_SV[i][8][j][3])-1].append(candidate_single_SV[i][11][j])
        candidate_single_SV_fam_ls[0][i][8] = fam_reads[0]
        candidate_single_SV_fam_ls[0][i][11] = reads_pos[0]
        candidate_single_SV_fam_ls[0][i][4] = len(fam_reads[0])
        for j in range(1,len(family_member_ls)) :
            fam_candidate_other[j][8] = fam_reads[j]
            fam_candidate_other[j][11] = reads_pos[j]
            fam_candidate_other[j][4] = len(fam_reads[j])
            candidate_single_SV_fam_ls[j].append(fam_candidate_other[j])
    
    candidate_single_SV_gt_fam_ls = []
    for i in range(len(family_member_ls)) :
        family_member = family_member_ls[i]
        candidate_single_SV_gt_fam_ls.append(call_gt(path, chr, candidate_single_SV_fam_ls[i], candidate_single_SV_fam_ls[0], max_cluster_bias, 'INS', family_mode, family_member, minimum_support_reads_list, performing_phasing))
    standard_list = 0
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) != 0 :
            standard_list = i
            break
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) == 0 :
            candidate_single_SV_gt_fam_ls[i] = copy.deepcopy(candidate_single_SV_gt_fam_ls[standard_list])
            for s_v in candidate_single_SV_gt_fam_ls[i] :
                s_v[8] = "0/0"
                s_v[9] = "0,100,100,0,0,0,0,0,0"
                s_v[10] = 996
                s_v[11] = 0
                s_v[12] = ""
                s_v[15] = ""
                s_v[16] = ""
                s_v[17] = ""

    increase_sigs_through_pedigree(candidate_single_SV_gt_fam_ls, 'INS', minimum_support_reads_list, family_mode)
    
    unsolvable_correction(candidate_single_SV_gt_fam_ls, 'INS', family_mode)
    
    correction_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
    if not performing_phasing and family_mode == "M1" :
        resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)

    logging.info("Finished calling %s:%s:%f."%(chr, "INS", time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)
    

def generate_ins_cluster(semi_ins_cluster, chr, read_count, 
                         threshold_gloab, further_threshold_gloab, minimum_support_reads, candidate_single_SV, 
                         gt_round, remain_reads_ratio, merge_INS_threshold, family_mode, performing_phasing, all_ins_singnature_reads):
        
    is_print_flag = False 
    for i in semi_ins_cluster :
        if chr == '7' and i[0] > 158995120 and i[0] < 158996120 :
            is_print_flag = True
    # Remove duplicates
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]

    read_name_dict = dict()
    for i in range(len(semi_ins_cluster)) :
        element = semi_ins_cluster[i]
        if element[2] not in read_name_dict:
            read_name_dict[element[2]] = [element]
        else :
            read_name_dict[element[2]].append(element)
    for read_name in read_name_dict :
        read_name_dict[read_name] = sorted(read_name_dict[read_name], key = lambda x:x[0])
    len_max = 0
    for read_name in read_name_dict :
        pre_merge_ls = read_name_dict[read_name]
        read_merge_list = [pre_merge_ls[0]+[0]]
        pre_pos = pre_merge_ls[0][0]
        merge_dis_no_flag = False
        if pre_merge_ls[0][1] > len_max :
            len_max = pre_merge_ls[0][1] 
        for i in range(1,len(pre_merge_ls)) :
            if pre_merge_ls[i][1] > len_max :
                len_max = pre_merge_ls[i][1]
            if pre_merge_ls[i][0] - pre_pos > merge_INS_threshold :
                if merge_dis_no_flag :
                    read_merge_list[-1][4] = 0
                read_merge_list.append(pre_merge_ls[i]+[0])
                merge_dis_no_flag = False
            else :
                read_merge_list[-1][1] += pre_merge_ls[i][1]
                read_merge_list[-1][3] += pre_merge_ls[i][3]
                read_merge_list[-1][4] = read_merge_list[-1][4] + pre_merge_ls[i][0] - pre_pos
                if pre_merge_ls[i][1] >= 50 :
                    merge_dis_no_flag = True
            pre_pos = pre_merge_ls[i][0]
        if merge_dis_no_flag :
            read_merge_list[-1][4] = 0
        read_name_dict[read_name] = read_merge_list
    if len_max > 500 :
        for read_name in read_name_dict :
            pre_merge_ls = read_name_dict[read_name]
            if abs(sum([x[1] for x in pre_merge_ls])-len_max)/len_max < 0.1 : 
                read_name_dict[read_name] = [[pre_merge_ls[0][0],sum([x[1] for x in pre_merge_ls]),pre_merge_ls[0][2],''.join([x[3] for x in pre_merge_ls]),sum([x[4] for x in pre_merge_ls])]]

    # kepp longest in read
    read_tag = dict()
    for read_name in read_name_dict :
        pre_merge_ls = read_name_dict[read_name]
        read_tag[read_name] = pre_merge_ls[0]
        for i in range(1,len(pre_merge_ls)) :
            if pre_merge_ls[i][1] > read_tag[read_name][1] :
                read_tag[read_name] = pre_merge_ls[i]
    if len(read_tag) < read_count:
        return
    read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
    # start&end breakpoint
    global_len = [i[1] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    last_len = read_tag2SortedList[0][1]

    allele_collect = list()
    allele_collect.append([[read_tag2SortedList[0][0]],
                            [read_tag2SortedList[0][1]],
                            [], 
                            [read_tag2SortedList[0][2]],
                            [read_tag2SortedList[0][3]],
                            [read_tag2SortedList[0][4]]])

    for i in read_tag2SortedList[1:]:
        if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[],[],[]])

        allele_collect[-1][0].append(i[0])
        allele_collect[-1][1].append(i[1])
        allele_collect[-1][3].append(i[2])
        allele_collect[-1][4].append(i[3])
        allele_collect[-1][5].append(i[4])
        last_len = i[1]
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    is_segment_flag = True
    while(is_segment_flag) :
        is_segment_flag = False
        further_allele_collect = list()
        for allele in allele_collect :
            DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = further_threshold_gloab * np.mean(allele[1])
            further_allele_collect.append([[allele[0][0]],[allele[1][0]],[],[allele[3][0]],[allele[4][0]],[allele[5][0]]])
            first_len = allele[1][0]
            for i in range(1,len(allele[1])) :
                if allele[1][i] - first_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP :
                    is_segment_flag = True
                    further_allele_collect[-1][2].append(len(further_allele_collect[-1][0]))
                    further_allele_collect.append([[],[],[],[],[],[]])
                    first_len = allele[1][i]
                further_allele_collect[-1][0].append(allele[0][i])
                further_allele_collect[-1][1].append(allele[1][i])
                further_allele_collect[-1][3].append(allele[3][i])
                further_allele_collect[-1][4].append(allele[4][i])
                further_allele_collect[-1][5].append(allele[5][i])
            further_allele_collect[-1][2].append(len(further_allele_collect[-1][0]))
        allele_collect = further_allele_collect
        break
    allele_sort = sorted(allele_collect, key = lambda x:x[2])
    dis_coefficient_var_threshold = 0.01
    len_coefficient_var_threshold = 0.3
    for allele in allele_sort:
        if allele[2][0] >= minimum_support_reads:
            if sum(allele[5]) / len(allele[5]) >= merge_INS_threshold/2 and np.mean(allele[2]) <= 150:
                fam_alleles = [[] for i in family_member_ls]
                for i in range(len(allele[0])) :
                    fam_index = int(allele[3][i][3])-1
                    fam_alleles[fam_index].append(allele[5][i])
                is_merge_indel = []
                is_merge_indel_num = 0
                for i in range(len(fam_alleles)) :
                    if len(fam_alleles[i]) > 0 and np.mean(fam_alleles[i]) > 0 and np.std(fam_alleles[i])/np.mean(fam_alleles[i]) <= dis_coefficient_var_threshold :
                        is_merge_indel.append(True)
                        is_merge_indel_num += 1
                    else :
                        is_merge_indel.append(False)
                fam_merge_ls = []
                if sum(is_merge_indel) > 0 and is_merge_indel_num >= 2:
                    for i in range(len(fam_alleles)) :
                        if is_merge_indel[i] :
                            fam_merge_ls += fam_alleles[i]
                    if len(fam_merge_ls) > 0 and np.mean(fam_merge_ls) > 0 and np.std(fam_merge_ls)/np.mean(fam_merge_ls) <= dis_coefficient_var_threshold :
                        filtered_allele = [[] for i in allele]
                        for i in range(len(allele[0])) :
                            if not is_merge_indel[int(allele[3][i][3])-1] :
                                filtered_allele[0].append(allele[0][i])
                                filtered_allele[1].append(allele[1][i])
                                filtered_allele[3].append(allele[3][i])
                                filtered_allele[4].append(allele[4][i])
                        filtered_allele[2].append(len(filtered_allele[0]))
                        allele = filtered_allele
            if allele[2][0] < minimum_support_reads:
                continue
            if np.std(allele[1])/np.mean(allele[1]) > len_coefficient_var_threshold :
                continue
            allele_list = list()
            var_list = list()
            remain_allele_num = max(int(remain_reads_ratio * allele[2][0]), 1)
            pos_mean = np.mean(allele[0])
            for i in range(len(allele[0])):
                var_list.append((abs(allele[0][i] - pos_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[0][var_list[i][1]])
            breakpointStart = np.mean(allele_list)

            allele_list = list()
            var_list = list()
            len_mean = np.mean(allele[1])
            for i in range(len(allele[1])):
                var_list.append((abs(allele[1][i] - len_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[1][var_list[i][1]])
            signalLen = np.median(allele_list)

            CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
            signalLen_STD = np.std(allele[1])
            CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))
            ideal_ins_seq = '<INS>'
            for pos,i in zip(allele[0],allele[4]):
                if len(i) >= int(signalLen):
                    breakpointStart = pos
                    ideal_ins_seq = i[0:int(signalLen)]
                    break
            if ideal_ins_seq == '<INS>':
                continue
            
            # keep sew of INS for te
            if all_ins_singnature_reads :
                if family_mode == "M1" :
                    ideal_ins_seq_family_ls = [[],[],[]]
                else :
                    ideal_ins_seq_family_ls = [[],[]]
                for i in range(len(allele[4])) :
                    ideal_ins_seq_family_ls[int(allele[3][i][3])-1].append(allele[4][i])
                ideal_ins_seq = ""
                for i in range(len(ideal_ins_seq_family_ls)) :
                    ideal_ins_seq = ideal_ins_seq + ";".join(ideal_ins_seq_family_ls[i]) + ";;"
                ideal_ins_seq = ideal_ins_seq[:-2]

            fam_allele_len_ls = [[] for x in family_member_ls]
            for seq_i in range(len(allele[0])) :
                fam_allele_len_ls[int(allele[3][seq_i][3])-1].append(allele[1][seq_i])
            for fam_i in range(len(fam_allele_len_ls)) :
                if len(fam_allele_len_ls[fam_i]) == 0 :
                    fam_allele_len_ls[fam_i].append(0)
            sv_len_str = ""
            for fam_i in range(len(fam_allele_len_ls)) :
                sv_len_str = sv_len_str + str(int(np.mean(fam_allele_len_ls[fam_i]))) + ","
            if performing_phasing :
                candidate_single_SV.append([chr, 
                                            "INS", 
                                            max(1,int(breakpointStart)), 
                                            int(signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(breakpointStart), 
                                            allele[3],
                                            ideal_ins_seq,
                                            sv_len_str[:-1],
                                            allele[0]])
            else :
                candidate_single_SV.append([chr, 
                                            "INS", 
                                            max(1,int(breakpointStart)), 
                                            int(signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(breakpointStart), 
                                            allele[3],
                                            ideal_ins_seq,
                                            sv_len_str[:-1],
                                            ''])


def run_del(args):
    return resolution_DEL(*args)

def run_ins(args):
    return resolution_INS(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, candidate_info_SV, max_cluster_bias, svtype, family_mode, family_member, minimum_support_reads_list, performing_phasing):
    with open("%s%s.%s.%s.pickle"%(temporary_dir,family_mode,family_member,"sigindex"), 'rb') as f:
        sigs_index=pickle.load(f)
        f.close()    
    if chr not in sigs_index["reads"].keys():
        return []
    readsfile = open("%s%s.%s.reads.pickle"%(temporary_dir,family_mode,family_member), 'rb')
    readsfile.seek(sigs_index["reads"][chr])
    reads_list=pickle.load(readsfile)
    readsfile.close()

    svs_list = list()
    for item in candidate_info_SV :
        pos_bias = 1000 if item[3] < 1000 else max_cluster_bias
        svs_list.append((max(item[7] - pos_bias, 0), item[7] + abs(item[3]) + pos_bias, item[3]))
    iteration_dict, primary_num_dict, cover_dict, overlap_dict, cover_pos_dict = overlap_cover(svs_list, reads_list, performing_phasing) # both key(sv idx), value(set(read id))
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"

    read_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][8]
    
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, cover_pos_dict, svtype, family_member, minimum_support_reads_list[int(family_member)-1], performing_phasing)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        assert len(assign_list[i]) == 8, "assign genotype error"
        #0-chr|1-sv_type|2-pos|3-len|4-support read num|5-CIPOS|6-CILEN|7-DR|8-GT|9-GL|10-GQ|11-QUAL|12-support read name|13-read squence|14-sv len str|15-support read pos|16-not support read name|17-not support read pos
        candidate_single_SV_gt.append([candidate_info_SV[i][0], 
                                    candidate_info_SV[i][1], 
                                    str(candidate_info_SV[i][2]), 
                                    str(candidate_info_SV[i][3]), 
                                    str(candidate_single_SV[i][4]), 
                                    candidate_info_SV[i][5],
                                    candidate_info_SV[i][6],
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][8])])
        if svtype == 'INS':
            candidate_single_SV_gt[i].append(candidate_info_SV[i][9])
        else :
            candidate_single_SV_gt[i].append([])
        if svtype == 'INS':
            candidate_single_SV_gt[i].append(candidate_info_SV[i][10])
            if performing_phasing :
                candidate_single_SV_gt[i].append(','.join([str(x) for x in candidate_single_SV[i][11]]))
                candidate_single_SV_gt[i].append(','.join(assign_list[i][6]))
                candidate_single_SV_gt[i].append(','.join([str(x) for x in assign_list[i][7]]))
            else :
                candidate_single_SV_gt[i].append('')
                candidate_single_SV_gt[i].append('')
                candidate_single_SV_gt[i].append('')
        else :
            candidate_single_SV_gt[i].append(candidate_info_SV[i][9])
            if performing_phasing :
                candidate_single_SV_gt[i].append(','.join([str(x) for x in candidate_single_SV[i][10]]))
                candidate_single_SV_gt[i].append(','.join(assign_list[i][6]))
                candidate_single_SV_gt[i].append(','.join([str(x) for x in assign_list[i][7]]))
            else :
                candidate_single_SV_gt[i].append('')
                candidate_single_SV_gt[i].append('')
                candidate_single_SV_gt[i].append('')
    return candidate_single_SV_gt