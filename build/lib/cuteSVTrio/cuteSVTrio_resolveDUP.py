from cuteSVTrio.cuteSVTrio_genotype import cal_CIPOS, overlap_cover, assign_gt, unsolvable_correction, modify_genotype, increase_sigs_through_pedigree, check_gt_consistent, inconformity_mendel_modify
from cuteSVTrio.cuteSVTrio_mendel import resolution_mendel
import numpy as np
from math import ceil
import logging
import pickle
import copy
import time

def resolution_DUP(path, chr, read_count, max_cluster_bias, minimum_support_reads_list, sv_size, MaxSize, gt_round, read_pos_interval, family_mode, performing_phasing):
    start_time = time.time()
    semi_dup_cluster = list()
    semi_dup_cluster.append([0, 0, ''])
    candidate_single_SV = list()
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

    seqs = []
    for family_member in family_member_ls :
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"sigindex"), 'rb') as f:
            sigs_index=pickle.load(f)
            #if family_member == "1" :
            #    logging.info(sigs_index["INS"])
            f.close()
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"DUP"), 'rb') as f:
            if chr not in sigs_index["DUP"].keys() :
                f.close()
                continue
            f.seek(sigs_index["DUP"][chr])
            f_seqs = pickle.load(f)
            seqs += f_seqs
            #if family_member == "1" and chr == "GL000243.1":
            #    logging.info(f_seqs)
            f.close()

    if len(seqs) == 0 :
        return
    seqs = sorted(seqs, key=lambda x: x[0])
    #if chr == "chr1" :
    #    logging.info(seqs)
    for seq in seqs:
        pos_1 = int(seq[0])
        pos_2 = int(seq[1])
        read_id = seq[2]
        # if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias or pos_2 - semi_dup_cluster[-1][1] > max_cluster_bias:
        if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias:
            if len(semi_dup_cluster) >= read_count:
                if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                    pass
                else:
                    generate_dup_cluster(semi_dup_cluster, 
                                         chr, 
                                         read_count, 
                                         max_cluster_bias, 
                                         sv_size, 
                                         candidate_single_SV,
                                         MaxSize,
                                         gt_round,
                                         family_mode,
                                         performing_phasing)
            semi_dup_cluster = []
            semi_dup_cluster.append([pos_1, pos_2, read_id])
        else:
            if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                semi_dup_cluster = []
                semi_dup_cluster.append([pos_1, pos_2, read_id])
            else:
                semi_dup_cluster.append([pos_1, pos_2, read_id])

    if len(semi_dup_cluster) >= read_count:
        if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
            pass
        else:
            generate_dup_cluster(semi_dup_cluster, 
                                 chr, 
                                 read_count, 
                                 max_cluster_bias, 
                                 sv_size, 
                                 candidate_single_SV,
                                 MaxSize,
                                 gt_round,
                                 family_mode,
                                 performing_phasing)
    candidate_single_SV_fam_ls = []
    candidate_single_SV_fam_ls.append(candidate_single_SV)
    #logging.info(candidate_single_SV[0])
    #logging.info(candidate_single_SV[1])
    #logging.info(candidate_single_SV[2])
    #logging.info(candidate_single_SV[3])
    #logging.info(candidate_single_SV[4])
    for i in range(1,len(family_member_ls)) :
        candidate_single_SV_fam_ls.append([])
    for i in range(len(candidate_single_SV)) :
        fam_candidate_other = [[]]
        for j in range(len(family_member_ls)-1) :
            fam_candidate_other.append(['' for x in candidate_single_SV[i]])
        fam_reads = [[] for j in family_member_ls]
        reads_pos = [[] for j in family_member_ls]
        for j in range(len(candidate_single_SV[i][4])) :
            fam_reads[int(candidate_single_SV[i][4][j][3])-1].append(candidate_single_SV[i][4][j])
            if performing_phasing :
                reads_pos[int(candidate_single_SV[i][4][j][3])-1].append(candidate_single_SV[i][5][j])
        #for j in range(len(family_member_ls)) :
        #    if len(fam_reads[j]) != len(reads_pos[j]) :
        #        logging.info("%s/%s"%(str(fam_reads[j]),str(reads_pos[j])))
        candidate_single_SV_fam_ls[0][i][4] = fam_reads[0]
        candidate_single_SV_fam_ls[0][i][5] = reads_pos[0]
        for j in range(1,len(family_member_ls)) :
            fam_candidate_other[j][4] = fam_reads[j]
            fam_candidate_other[j][5] = reads_pos[j]
            candidate_single_SV_fam_ls[j].append(fam_candidate_other[j])
    
    #if chr == "chr1" :
    #    logging.info(candidate_single_SV_fam_ls)
    candidate_single_SV_gt_fam_ls = []
    for i in range(len(family_member_ls)) :
        family_member = family_member_ls[i]
        #gt_candidate_sv = call_gt(path, chr, candidate_single_SV_fam_ls[i], 100, 'INS', family_mode, family_member)
        candidate_single_SV_gt_fam_ls.append(call_gt(path, chr, candidate_single_SV_fam_ls[i], candidate_single_SV_fam_ls[0], 1000, family_mode, family_member, minimum_support_reads_list, performing_phasing))
    #logging.info("DUP/%s/%d/%s"%(chr,len(candidate_single_SV_gt_fam_ls[0]),str(candidate_single_SV_gt_fam_ls[0][0:10])))
    #logging.info(candidate_single_SV_gt_fam_ls)
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
                #logging.info(s_v)
    #logging.info("%s/%s"%(chr,str(candidate_single_SV_gt_fam_ls)))
    #if chr == "1" :
    #    logging.info(candidate_single_SV_gt_fam_ls)
    increase_sigs_through_pedigree(candidate_single_SV_gt_fam_ls, 'DUP', minimum_support_reads_list, family_mode)
    unsolvable_correction(candidate_single_SV_gt_fam_ls, 'DUP', family_mode)
    if not performing_phasing and family_mode == "M1" :
        resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
    #如果非INDEL的sv参与phasing，那么可能会在phasing部分进行denovo分析
    #if not performing_phasing :
    #    resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, performing_phasing)
    #unsolvable_correction(candidate_single_SV_gt_fam_ls, 'DUP', family_mode, family_member)
    #not_accord_mendel_num,run_mcmc_num,run_mcmc_ls=resolution_mendel(candidate_single_SV_gt_fam_ls,'DUP',family_mode, gap_thres)
    #resolution_denovo(candidate_single_SV_gt_fam_ls, 'DUP', family_mode)
    #modify_genotype(candidate_single_SV_gt_fam_ls, 'DUP')
    logging.info("Finished calling %s:%s:%f."%(chr, "DUP", time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)

def generate_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, 
    sv_size, candidate_single_SV, MaxSize, gt_round, family_mode, performing_phasing):
    # calculate support reads
    #logging.info(semi_dup_cluster)
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    support_read = list(set([i[2] for i in semi_dup_cluster]))
    if len(support_read) < read_count:
        return

    semi_dup_cluster.sort(key = lambda x:x[1])
    allele_collect = []
    allele_collect.append([semi_dup_cluster[0]])
    last_len = semi_dup_cluster[0][1]
    for i in semi_dup_cluster[1:]:
        if i[1] - last_len > max_cluster_bias:
            allele_collect.append([])
        allele_collect[-1].append(i)
        last_len = i[1]
    for i in allele_collect:
        support_read = []
        pos_ls = []
        fam_allele_len_ls = [[] for x in family_member_ls]
        for j in i :
            if j[2] not in support_read :
                #logging.info(j)
                support_read.append(j[2])
                pos_ls.append(j[0])
                fam_allele_len_ls[int(j[2][3])-1].append(abs(int(j[1])-int(j[0])))
        for fam_i in range(len(fam_allele_len_ls)) :
            if len(fam_allele_len_ls[fam_i]) == 0 :
                fam_allele_len_ls[fam_i].append(0)
        sv_len_str = ""
        for fam_i in range(len(fam_allele_len_ls)) :
            sv_len_str = sv_len_str + str(int(np.mean(fam_allele_len_ls[fam_i]))) + ","
        #logging.info(i)
        if len(support_read) < read_count:
            continue
        low_b = int(len(i)*0.4)
        up_b = int(len(i)*0.6)

        if low_b == up_b:
            breakpoint_1 = i[low_b][0]
            breakpoint_2 = i[low_b][1]
        else:
            breakpoint_1 = [i[0] for i in i[low_b:up_b]]
            breakpoint_2 = [i[1] for i in i[low_b:up_b]]
            breakpoint_1 = int(sum(breakpoint_1)/len(i[low_b:up_b]))
            breakpoint_2 = int(sum(breakpoint_2)/len(i[low_b:up_b]))


        if sv_size <= breakpoint_2 - breakpoint_1 <= MaxSize or (sv_size <= breakpoint_2 - breakpoint_1 and MaxSize == -1):
            if performing_phasing :
                candidate_single_SV.append([chr,
                                            'DUP', 
                                            max(1,breakpoint_1), 
                                            max(1,breakpoint_2),
                                            support_read,
                                            pos_ls,
                                            sv_len_str[:-1]])
            else :
                candidate_single_SV.append([chr,
                                            'DUP', 
                                            max(1,breakpoint_1), 
                                            max(1,breakpoint_2),
                                            support_read,
                                            '',
                                            sv_len_str[:-1]])

def run_dup(args):
    return resolution_DUP(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, candidate_info_SV, max_cluster_bias, family_mode, family_member, minimum_support_reads_list, performing_phasing):
    # reads_list = list() # [(10000, 10468, 0, 'm54238_180901_011437/52298335/ccs'), ...]
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
    for item in candidate_info_SV:
        new_cluster_bias = max(max_cluster_bias, item[3] - item[2])
        svs_list.append((max(item[2] - new_cluster_bias, 0), item[3] + new_cluster_bias))
    for item in candidate_info_SV:
        new_cluster_bias = max(max_cluster_bias, item[3] - item[2])
        svs_list.append((max(item[3] - new_cluster_bias, 0), item[3] + item[3] - item[2] + new_cluster_bias))
    
    iteration_dict, primary_num_dict, cover_dict, overlap_dict, cover_pos_dict = overlap_cover(svs_list, reads_list, performing_phasing) # both key(sv idx), value(set(read id))
    assert len(cover_dict) == 2 * len(candidate_single_SV), "overlap length error"
    candidate_single_SV_length = len(candidate_single_SV)
    for idx in range(candidate_single_SV_length):
        for i in range(len(cover_dict[idx + candidate_single_SV_length])) :
            if cover_dict[idx + candidate_single_SV_length][i] not in cover_dict[idx] :
                cover_dict[idx].append(cover_dict[idx + candidate_single_SV_length][i])
                if performing_phasing :
                    cover_pos_dict[idx].append(cover_pos_dict[idx + candidate_single_SV_length][i])
    for idx in range(candidate_single_SV_length, candidate_single_SV_length * 2, 1):
        cover_dict.pop(idx)
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"

    read_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][4]
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, cover_pos_dict, "DUP", family_member, minimum_support_reads_list[int(family_member)-1], performing_phasing)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    #0-chr|1-sv_type|2-pos|3-len|4-support read num|5-''|6-''|7-DR|8-GT|9-GL|10-GQ|11-QUAL|12-support read name|13-[]|14-sv len str|15-support read pos|16-not support read name|17-not support read pos
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_info_SV[i][0], 
                                    candidate_info_SV[i][1], 
                                    str(candidate_info_SV[i][2]), 
                                    str(candidate_info_SV[i][3] - candidate_info_SV[i][2]), 
                                    str(len(candidate_single_SV[i][4])),
                                    '',
                                    '',
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][4])])
        candidate_single_SV_gt[i].append([])
        candidate_single_SV_gt[i].append(candidate_info_SV[i][6])

        if performing_phasing :
            candidate_single_SV_gt[i].append(','.join([str(x) for x in candidate_single_SV[i][5]]))
            candidate_single_SV_gt[i].append(','.join(assign_list[i][6]))
            candidate_single_SV_gt[i].append(','.join([str(x) for x in assign_list[i][7]]))
        else :
            candidate_single_SV_gt[i].append('')
            candidate_single_SV_gt[i].append('')
            candidate_single_SV_gt[i].append('')
    
    
    
    
    return candidate_single_SV_gt	