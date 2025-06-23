from cuteSVTrio.cuteSVTrio_genotype import cal_GL_2,normalize_log10_probs
from cuteSVTrio.cuteSVTrio_mendel import resolution_mendel
import numpy as np
from math import log10,ceil
import logging
import random
import copy
import time

#consistency_threshold用来限制两个基因型概率的差值
consistency_threshold = 0.15
#genotype_threshold用来确定某个基因型概率是否接近0，1
genotype_threshold = 0.1

# 以下代码功能为安装pos分区域记录read和sv之间的关系，read为read name，sv为编译列表中的索引位置
def make_read_to_sv_list(chr, candidate_single_SV_gt, family_mode, read_pos_interval) :
    read_pos_max = 0
    for i in candidate_single_SV_gt :
        if int(i[2]) > read_pos_max :
            read_pos_max = int(i[2])
    read_pos_max =int(read_pos_max * 1.01)
    #logging.info(candidate_single_SV_gt_fam_ls[0][:10])
    read_to_var_list = [[[],[]] for j in range(ceil(read_pos_max/read_pos_interval))]
    read_to_non_list = [[[],[]] for j in range(ceil(read_pos_max/read_pos_interval))]
    for j in range(len(candidate_single_SV_gt)) :
        if candidate_single_SV_gt[j][12] == '' :
            continue
        #if chr == "1" :
        #    logging.info("%s/%s"%(candidate_single_SV_gt[j][12],candidate_single_SV_gt[j][15]))
        read_names = candidate_single_SV_gt[j][12].split(",")
        read_poss = [int(x) for x in candidate_single_SV_gt[j][15].split(",")]
        for k in range(len(read_names)) :
            #logging.info("%s/%s"%(str(read_pos_max),str(read_poss[k])))
            read_interval_index = ceil(read_poss[k]/read_pos_interval) - 1
            if read_names[k] in read_to_var_list[read_interval_index][0] :
                read_to_var_list[read_interval_index][1][read_to_var_list[read_interval_index][0].index(read_names[k])].append(j)
            else :
                read_to_var_list[read_interval_index][0].append(read_names[k])
                read_to_var_list[read_interval_index][1].append([j])
        if candidate_single_SV_gt[j][16] == "" :
            continue
        read_names = candidate_single_SV_gt[j][16].split(",")
        read_poss = [int(x) for x in candidate_single_SV_gt[j][17].split(",")]
        for k in range(len(read_names)) :
            #logging.info("%s/%s"%(str(read_pos_max),str(read_poss[k])))
            read_interval_index = ceil(read_poss[k]/read_pos_interval) - 1
            if read_names[k] in read_to_non_list[read_interval_index][0] :
                read_to_non_list[read_interval_index][1][read_to_non_list[read_interval_index][0].index(read_names[k])].append(j)
            else :
                read_to_non_list[read_interval_index][0].append(read_names[k])
                read_to_non_list[read_interval_index][1].append([j])
    return read_to_var_list,read_to_non_list

def confirm_haplotype_source(fam_genotype_ls, family_mode, family_member) :
    if family_mode == "M1" :
        ch_gt = fam_genotype_ls[0]
        fa_gt = fam_genotype_ls[1]
        mo_gt = fam_genotype_ls[2]
        if family_member == "1" :
            dis_1 = min(abs(ch_gt[0]-fa_gt[0]),abs(ch_gt[0]-fa_gt[1])) + min(abs(ch_gt[1]-mo_gt[0]),abs(ch_gt[1]-mo_gt[1]))
            dis_2 = min(abs(ch_gt[1]-fa_gt[0]),abs(ch_gt[1]-fa_gt[1])) + min(abs(ch_gt[0]-mo_gt[0]),abs(ch_gt[0]-mo_gt[1]))
            # 0代表孩子的第一个基因来自父本，第二个来自母本；1代表孩子的第一个基因来自母本，第二个来自父本
            if dis_1 < consistency_threshold and dis_2 > 1-consistency_threshold :
                return True,0
            elif dis_2 < consistency_threshold and dis_1 > 1-consistency_threshold :
                return True,1
            else :
                return False,-1
        elif family_member == "2" :
            dis_1 = min(abs(fa_gt[0]-ch_gt[0])+min(abs(mo_gt[0]-ch_gt[1]),abs(mo_gt[1]-ch_gt[1])),abs(fa_gt[0]-ch_gt[1])+min(abs(mo_gt[0]-ch_gt[0]),abs(mo_gt[1]-ch_gt[0])))
            dis_2 = min(abs(fa_gt[1]-ch_gt[0])+min(abs(mo_gt[0]-ch_gt[1]),abs(mo_gt[1]-ch_gt[1])),abs(fa_gt[1]-ch_gt[1])+min(abs(mo_gt[0]-ch_gt[0]),abs(mo_gt[1]-ch_gt[0])))
            # 0代表父本的第一个基因遗传向孩子，1代表父本的第二个基因遗传向孩子
            if dis_1 < consistency_threshold and dis_2 > 1-consistency_threshold :
                return True,0
            elif dis_2 < consistency_threshold and dis_1 > 1-consistency_threshold :
                return True,1
            else :
                return False,-1
        else :
            dis_1 = min(abs(mo_gt[0]-ch_gt[0])+min(abs(fa_gt[0]-ch_gt[1]),abs(fa_gt[1]-ch_gt[1])),abs(mo_gt[0]-ch_gt[1])+min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[1]-ch_gt[0])))
            dis_2 = min(abs(mo_gt[1]-ch_gt[0])+min(abs(fa_gt[0]-ch_gt[1]),abs(fa_gt[1]-ch_gt[1])),abs(mo_gt[1]-ch_gt[1])+min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[1]-ch_gt[0])))
            # 0代表母本的第一个基因遗传向孩子，1代表母本的第二个基因遗传向孩子
            if dis_1 < consistency_threshold and dis_2 > 1-consistency_threshold :
                return True,0
            elif dis_2 < consistency_threshold and dis_1 > 1-consistency_threshold :
                return True,1
            else :
                return False,-1
    elif family_mode == "M2" :
        ch_gt = fam_genotype_ls[0]
        fa_gt = fam_genotype_ls[1]
        if family_member == "1" :
            dis_1 = min(abs(ch_gt[0]-fa_gt[0]),abs(ch_gt[0]-fa_gt[1]))
            dis_2 = max(abs(ch_gt[0]-fa_gt[0]),abs(ch_gt[0]-fa_gt[1]))
            dis_3 = min(abs(ch_gt[1]-fa_gt[0]),abs(ch_gt[1]-fa_gt[1]))
            dis_4 = max(abs(ch_gt[1]-fa_gt[0]),abs(ch_gt[1]-fa_gt[1]))
            # 0代表孩子的第一个基因来自父本，第二个来自母本；1代表孩子的第一个基因来自母本，第二个来自父本
            if dis_2 < consistency_threshold and dis_3 > 1-consistency_threshold :
                return True,0
            elif dis_4 < consistency_threshold and dis_1 > 1-consistency_threshold :
                return True,1
            else :
                return False,-1
        elif family_member == "2" :
            dis_1 = min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[0]-ch_gt[1]))
            dis_2 = max(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[0]-ch_gt[1]))
            dis_3 = min(abs(fa_gt[1]-ch_gt[0]),abs(fa_gt[1]-ch_gt[1]))
            dis_4 = max(abs(fa_gt[1]-ch_gt[0]),abs(fa_gt[1]-ch_gt[1]))
            # 0代表父本的第一个基因遗传向孩子，1代表父本的第二个基因遗传向孩子
            if dis_2 < consistency_threshold and dis_3 > 1-consistency_threshold :
                return True,0
            elif dis_4 < consistency_threshold and dis_1 > 1-consistency_threshold :
                return True,1
            else :
                return False,-1

# [gt1,gt2] 如果是纯和，则返回True
def gt_homozygous(genotype_ls) :
    return max(abs(genotype_ls[0]-1),abs(genotype_ls[1]-1)) < genotype_threshold or max(abs(genotype_ls[0]-0),abs(genotype_ls[1]-0)) < genotype_threshold

def genetic_phasing_family(chr, candidate_single_SV_gt_fam_ls, family_mode, read_pos_interval, minimum_support_reads_list, phase_all_ctgs, parents_phasing) :
    start_time = time.time()
    #logging.info("Phasing starting of %s:%f."%(chr, time.time()-start_time))
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    #chr_ls = []
    if phase_all_ctgs == True or (phase_all_ctgs == False and chr in chr_ls) :
        #logging.info("phasing %s"%(chr))
        #return (chr,candidate_single_SV_gt_fam_ls)
        read_to_var_list, read_to_non_list = make_read_to_sv_list(chr, candidate_single_SV_gt_fam_ls[0], family_mode, read_pos_interval)
        phased_sv_haplotype_child_father = []
        phased_sv_haplotype_child_mother = []
        phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother = genetic_phasing_member(chr, 
                                                                                                   candidate_single_SV_gt_fam_ls, 
                                                                                                   phased_sv_haplotype_child_father, 
                                                                                                   phased_sv_haplotype_child_mother, 
                                                                                                   family_mode, 
                                                                                                   "1", 
                                                                                                   read_pos_interval, 
                                                                                                   read_to_var_list, 
                                                                                                   read_to_non_list, 
                                                                                                   1)
        phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother = genetic_phasing_member(chr, 
                                                                                                   candidate_single_SV_gt_fam_ls, 
                                                                                                   phased_sv_haplotype_child_father, 
                                                                                                   phased_sv_haplotype_child_mother, 
                                                                                                   family_mode, 
                                                                                                   "1", 
                                                                                                   read_pos_interval, 
                                                                                                   read_to_var_list, 
                                                                                                   read_to_non_list, 
                                                                                                   2)
        if parents_phasing :
            read_to_var_list, read_to_non_list = make_read_to_sv_list(chr, candidate_single_SV_gt_fam_ls[1], family_mode, read_pos_interval)
            phased_sv_haplotype_father_inher = []
            phased_sv_haplotype_father_forgo = []
            phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo = genetic_phasing_member(chr, 
                                                                                                       candidate_single_SV_gt_fam_ls, 
                                                                                                       phased_sv_haplotype_father_inher, 
                                                                                                       phased_sv_haplotype_father_forgo, 
                                                                                                       family_mode, 
                                                                                                       "2", 
                                                                                                       read_pos_interval, 
                                                                                                       read_to_var_list, 
                                                                                                       read_to_non_list, 
                                                                                                       1)
            phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo = genetic_phasing_member(chr, 
                                                                                                       candidate_single_SV_gt_fam_ls, 
                                                                                                       phased_sv_haplotype_father_inher, 
                                                                                                       phased_sv_haplotype_father_forgo, 
                                                                                                       family_mode, 
                                                                                                       "2", 
                                                                                                       read_pos_interval, 
                                                                                                       read_to_var_list, 
                                                                                                       read_to_non_list, 
                                                                                                       2)
            phased_sv_haplotype_mother_inher = []
            phased_sv_haplotype_mother_forgo = []
            if family_mode == "M1" :
                read_to_var_list, read_to_non_list = make_read_to_sv_list(chr, candidate_single_SV_gt_fam_ls[2], family_mode, read_pos_interval)
                phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo = genetic_phasing_member(chr, 
                                                                                                           candidate_single_SV_gt_fam_ls, 
                                                                                                           phased_sv_haplotype_mother_inher, 
                                                                                                           phased_sv_haplotype_mother_forgo, 
                                                                                                           family_mode, 
                                                                                                           "3", 
                                                                                                           read_pos_interval, 
                                                                                                           read_to_var_list, 
                                                                                                           read_to_non_list, 
                                                                                                           1)
                phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo = genetic_phasing_member(chr, 
                                                                                                           candidate_single_SV_gt_fam_ls, 
                                                                                                           phased_sv_haplotype_mother_inher, 
                                                                                                           phased_sv_haplotype_mother_forgo, 
                                                                                                           family_mode, 
                                                                                                           "3", 
                                                                                                           read_pos_interval, 
                                                                                                           read_to_var_list, 
                                                                                                           read_to_non_list, 
                                                                                                           2)
        else :
            phased_sv_haplotype_father_inher = []
            phased_sv_haplotype_father_forgo = []
            phased_sv_haplotype_mother_inher = []
            phased_sv_haplotype_mother_forgo = []
        
        phasing_candidate_fam_SV(chr,
                                 candidate_single_SV_gt_fam_ls, 
                                 phased_sv_haplotype_child_father, 
                                 phased_sv_haplotype_child_mother, 
                                 phased_sv_haplotype_father_inher, 
                                 phased_sv_haplotype_father_forgo, 
                                 phased_sv_haplotype_mother_inher, 
                                 phased_sv_haplotype_mother_forgo, 
                                 family_mode,
                                 True,
                                 parents_phasing)
        if family_mode == "M1" :
            if parents_phasing :
                resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, False, minimum_support_reads_list)
            else :
                resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
    else :
        phasing_candidate_fam_SV(chr,
                                 candidate_single_SV_gt_fam_ls, 
                                 [], 
                                 [], 
                                 [], 
                                 [], 
                                 [], 
                                 [], 
                                 family_mode,
                                 False,
                                 parents_phasing)
        if family_mode == "M1" :
            resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
        logging.info("no phasing %s"%(chr))
        return (chr,candidate_single_SV_gt_fam_ls)
    
    logging.info("Finished phasing %s:%f."%(chr, time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)

def genetic_no_phasing_family(chr, candidate_single_SV_gt_fam_ls, family_mode, read_pos_interval, minimum_support_reads_list, phase_all_ctgs) :
    start_time = time.time()
    phasing_candidate_fam_SV(chr,
                             candidate_single_SV_gt_fam_ls, 
                             [], 
                             [], 
                             [], 
                             [], 
                             [], 
                             [], 
                             family_mode,
                             False,
                             True)
    if family_mode == "M1" :
        resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
    logging.info("Finished no phasing %s:%f."%(chr, time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)


def family_nearest_match(fam_genotype_ls, family_mode, family_member) :
    if family_mode == "M1" :
        ch_gt = fam_genotype_ls[0]
        fa_gt = fam_genotype_ls[1]
        mo_gt = fam_genotype_ls[2]
        if family_member == "1" :
            dis_1 = min(abs(ch_gt[0]-fa_gt[0]),abs(ch_gt[0]-fa_gt[1])) + min(abs(ch_gt[1]-mo_gt[0]),abs(ch_gt[1]-mo_gt[1]))
            dis_2 = min(abs(ch_gt[1]-fa_gt[0]),abs(ch_gt[1]-fa_gt[1])) + min(abs(ch_gt[0]-mo_gt[0]),abs(ch_gt[0]-mo_gt[1]))
        elif family_member == "2" :
            dis_1 = min(abs(fa_gt[0]-ch_gt[0])+min(abs(mo_gt[0]-ch_gt[1]),abs(mo_gt[1]-ch_gt[1])),abs(fa_gt[0]-ch_gt[1])+min(abs(mo_gt[0]-ch_gt[0]),abs(mo_gt[1]-ch_gt[0])))
            dis_2 = min(abs(fa_gt[1]-ch_gt[0])+min(abs(mo_gt[0]-ch_gt[1]),abs(mo_gt[1]-ch_gt[1])),abs(fa_gt[1]-ch_gt[1])+min(abs(mo_gt[0]-ch_gt[0]),abs(mo_gt[1]-ch_gt[0])))
        else :
            dis_1 = min(abs(mo_gt[0]-ch_gt[0])+min(abs(fa_gt[0]-ch_gt[1]),abs(fa_gt[1]-ch_gt[1])),abs(mo_gt[0]-ch_gt[1])+min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[1]-ch_gt[0])))
            dis_2 = min(abs(mo_gt[1]-ch_gt[0])+min(abs(fa_gt[0]-ch_gt[1]),abs(fa_gt[1]-ch_gt[1])),abs(mo_gt[1]-ch_gt[1])+min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[1]-ch_gt[0])))
        #0代表子代的第一个基因遗传自父本，或父本第一个基因遗传向子代；1代表子代的第二个基因遗传自父本，或父本第二个基因遗传向子代
        if dis_1 < dis_2 :
            return True,0
        else :
            return True,1
    elif family_mode == "M2" :
        ch_gt = fam_genotype_ls[0]
        fa_gt = fam_genotype_ls[1]
        if family_member == "1" :
            dis_1 = min(abs(ch_gt[0]-fa_gt[0]),abs(ch_gt[0]-fa_gt[1]))
            dis_2 = min(abs(ch_gt[1]-fa_gt[0]),abs(ch_gt[1]-fa_gt[1]))
        elif family_member == "2" :
            dis_1 = min(abs(fa_gt[0]-ch_gt[0]),abs(fa_gt[0]-ch_gt[1]))
            dis_2 = min(abs(fa_gt[1]-ch_gt[0]),abs(fa_gt[1]-ch_gt[1]))
        if dis_1 < dis_2 :
            return True,0
        else :
            return True,1
        
def recalculate_block_no(phased_sv_haplotype_1, phased_sv_haplotype_2) :
    # 以下处理只考虑单独的杂合
    new_block_no = 1
    for sv_i in range(len(phased_sv_haplotype_1)) :
        if phased_sv_haplotype_1[sv_i][1] in [1,2,3] :
            phased_sv_haplotype_1[sv_i][3] = 0
            phased_sv_haplotype_2[sv_i][3] = 0
            continue
        elif phased_sv_haplotype_1[sv_i][1] == 0 :
            phased_sv_haplotype_1[sv_i][3] = -1
            phased_sv_haplotype_2[sv_i][3] = -1
            continue
        if sv_i == 0 :
            if phased_sv_haplotype_1[sv_i][1] in [4,5,6] :
                phased_sv_haplotype_1[sv_i][3] = new_block_no
                phased_sv_haplotype_2[sv_i][3] = new_block_no
        else :
            if phased_sv_haplotype_1[sv_i][1] in [4,5,6] :
                if phased_sv_haplotype_1[sv_i-1][1] not in [4,5,6] :
                    new_block_no += 1
                phased_sv_haplotype_1[sv_i][3] = new_block_no
                phased_sv_haplotype_2[sv_i][3] = new_block_no
    return new_block_no

def genetic_phasing_member(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_father, phased_sv_haplotype_mother, family_mode, family_member, read_pos_interval, read_to_var_list, read_to_non_list, phase_type) :
    if family_mode == "M1" :
        if family_member not in ["1","2","3"] :
            logging.info("Wrong family member!")
            return
    elif family_mode == "M1" :
        if family_member not in ["1","2"] :
            logging.info("Wrong family member!")
            return
    family_index = int(family_member) - 1
    #logging.info(read_to_sv_dict)
    read_to_vars = read_to_var_list
    read_to_nons = read_to_non_list
    gl_index = 9
    var_read_nameds_index = 12
    var_read_poss_index = 15
    non_read_nameds_index = 16
    non_read_poss_index = 17
    block_no = 1
    if phase_type == 1 :
        anchor_flag = 1
        phased_flag = 2
        break_flag = 3
    else :
        anchor_flag = 4
        phased_flag = 5
        break_flag = 6
    # 处理read phasing+genetic phasing
    # 使用遗传phasing确定基础锚点
    for j in range(len(candidate_single_SV_gt_fam_ls[family_index])) :
        #if j % 5000 == 0 :
        #    logging.info("%d/%d"%(j,len(candidate_single_SV_gt_fam_ls[family_index])))
        if phase_type == 2 and phased_sv_haplotype_father[j][1] in [1,2,3] :
            continue
        fam_genotype_ls = []
        if family_mode == "M1" :
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][j][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][j][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[2][j][gl_index].split(",")[3:5]])
        elif family_mode == "M2" :
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][j][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][j][gl_index].split(",")[3:5]])
        candidate_single_SV_ls = candidate_single_SV_gt_fam_ls[family_index][j]
        # 纯和
        if gt_homozygous(fam_genotype_ls[family_index]) :
            if phase_type == 1 :
                # 如果phase孩子
                # 0-基因型概率,1-phase状态,2-父母本来源,3-block序号,4-pos,5-type,6-变异reads列表,7-变异pos列表,8-无变异reads列表,9-无变异pos列表,10-phase得到的变异reads列表，11-phase得到的无变异reads列表
                # 父母本来源：0-前一个基因来自父本，后一个基因来自母本；1-前一个基因来自母本，后一个基因来自父本
                # 如果phase父母
                # 0-基因型概率,1-phase状态,2-遗传单倍体去向,3-block序号,4-pos,5-type,6-变异reads列表,7-变异pos列表,8-无变异reads列表,9-无变异pos列表,10-phase得到的变异reads列表，11-phase得到的无变异reads列表
                # 遗传单倍体去向：0-前一个基因遗传向孩子，1-后一个基因遗传向孩子
                # phase状态：0-未phase,1-遗传phase ,2-关联phase,3-断点,4-read phase anchor,5-read phase sv,6-read phase断点,7-就近匹配的纯和
                phased_sv_haplotype_father.append([0,0,-1,0,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]])
                phased_sv_haplotype_mother.append([0,0,-1,0,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]])
            else :
                continue
        else :
            if phase_type == 1 :
                is_con, gen_source = confirm_haplotype_source(fam_genotype_ls,family_mode,family_member)
            else :
                is_con, gen_source = family_nearest_match(fam_genotype_ls,family_mode,family_member)
            if not is_con :
                phased_sv_haplotype_father.append([0,0,-1,0,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]])
                phased_sv_haplotype_mother.append([0,0,-1,0,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]])
            else :
                if fam_genotype_ls[family_index][gen_source] < genotype_threshold :
                    if phase_type == 1 :
                        hap_t = [fam_genotype_ls[family_index][gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        hap_t[8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                        hap_t[9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                        phased_sv_haplotype_father.append(hap_t) 
                        hap_t = [fam_genotype_ls[family_index][1-gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        hap_t[6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                        hap_t[7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                        phased_sv_haplotype_mother.append(hap_t)
                    else :
                        phased_sv_haplotype_father[j] = [fam_genotype_ls[family_index][gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        phased_sv_haplotype_mother[j] = [fam_genotype_ls[family_index][1-gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                    #phased_sv_haplotype_father.append([fam_genotype_ls[0][gen_source],1,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],candidate_single_SV_ls[non_read_nameds_index].split(","),candidate_single_SV_ls[non_read_poss_index].split(",")])
                    #phased_sv_haplotype_mother.append([fam_genotype_ls[0][1-gen_source],1,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],candidate_single_SV_ls[var_read_nameds_index].split(","),candidate_single_SV_ls[var_read_poss_index].split(","),[],[]])
                else :
                    if phase_type == 1 :
                        hap_t = [fam_genotype_ls[family_index][gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        hap_t[6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                        hap_t[7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                        phased_sv_haplotype_father.append(hap_t)
                        hap_t = [fam_genotype_ls[family_index][1-gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        hap_t[8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                        hap_t[9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                        phased_sv_haplotype_mother.append(hap_t)
                    else :
                        phased_sv_haplotype_father[j] = [fam_genotype_ls[family_index][gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                        phased_sv_haplotype_mother[j] = [fam_genotype_ls[family_index][1-gen_source],anchor_flag,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],[],[],[],[]]
                    #phased_sv_haplotype_father.append([fam_genotype_ls[0][gen_source],1,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],candidate_single_SV_ls[var_read_nameds_index].split(","),candidate_single_SV_ls[var_read_poss_index].split(","),[],[]])
                    #phased_sv_haplotype_mother.append([fam_genotype_ls[0][1-gen_source],1,gen_source,block_no,candidate_single_SV_ls[1],candidate_single_SV_ls[2],[],[],candidate_single_SV_ls[non_read_nameds_index].split(","),candidate_single_SV_ls[non_read_poss_index].split(",")])
                if phase_type == 1 :
                    block_no += 1      
    #import sys
    #logging.info("%s/%f/%f"%(chr,float(sys.getsizeof(phased_sv_haplotype_father)),float(sys.getsizeof(candidate_single_SV_gt_fam_ls))))
    phased_sv_haplotype_father_forward = copy.deepcopy(phased_sv_haplotype_father)
    phased_sv_haplotype_mother_forward = copy.deepcopy(phased_sv_haplotype_mother)
    phased_sv_haplotype_father_reverse = copy.deepcopy(phased_sv_haplotype_father)
    phased_sv_haplotype_mother_reverse = copy.deepcopy(phased_sv_haplotype_mother)
    homozygous_phased_num = 0
    # 正向
    for sv_i in range(len(phased_sv_haplotype_father_forward)) :
        #if sv_i % 5000 == 0 :
        #    logging.info("%d/%d"%(sv_i,len(phased_sv_haplotype_father_forward)))
        if phase_type == 2 and phased_sv_haplotype_father_forward[sv_i][1] in [1,2,3] :
            continue
        genotype_ls = [float(x) for x in candidate_single_SV_gt_fam_ls[family_index][sv_i][gl_index].split(",")[3:5]]
        if gt_homozygous(genotype_ls) :
            if sv_i == 0:
                #phased_sv_haplotype_father_forward[sv_i][1] = 3
                pass
            elif ((len(phased_sv_haplotype_father_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_father_forward[sv_i][8]) == 0) and
                  (len(phased_sv_haplotype_mother_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_forward[sv_i][8]) == 0)) :
                    if phased_sv_haplotype_father_forward[sv_i-1][1] == 2 and phased_sv_haplotype_mother_forward[sv_i-1][1] == 2 and phase_type == 1:
                        phased_sv_haplotype_father_forward[sv_i-1][1] = 3
                        phased_sv_haplotype_mother_forward[sv_i-1][1] = 3
                    elif phased_sv_haplotype_father_forward[sv_i-1][1] == 5 and phased_sv_haplotype_mother_forward[sv_i-1][1] == 5 and phase_type == 2:
                        phased_sv_haplotype_father_forward[sv_i-1][1] = 6
                        phased_sv_haplotype_mother_forward[sv_i-1][1] = 6
            else :
                #homozygous_phased_num += 1
                phased_sv_haplotype_father_forward[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_forward[sv_i][1] = phased_flag
            continue
        if phased_sv_haplotype_father_forward[sv_i][1] in [0,4] :
            # 确定基因型父母本来源
            # 如果父母本的基因型确定结果有冲突，可能需要纠错
            if phased_sv_haplotype_father_forward[sv_i][1] == 0 and (len(phased_sv_haplotype_father_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_father_forward[sv_i][8]) == 0) and (len(phased_sv_haplotype_mother_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_forward[sv_i][8]) == 0) :
                if sv_i == 0:
                    #phased_sv_haplotype_father_forward[sv_i][1] = 3
                    pass
                elif ((len(phased_sv_haplotype_father_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_father_forward[sv_i][8]) == 0 and phased_sv_haplotype_father_forward[sv_i-1][1] in [2,5]) and
                      (len(phased_sv_haplotype_mother_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_forward[sv_i][8]) == 0 and phased_sv_haplotype_mother_forward[sv_i-1][1] in [2,5])) :
                    if phased_sv_haplotype_father_forward[sv_i-1][1] == 2 and phased_sv_haplotype_mother_forward[sv_i-1][1] == 2 and phase_type == 1:
                        phased_sv_haplotype_father_forward[sv_i-1][1] = 3
                        phased_sv_haplotype_mother_forward[sv_i-1][1] = 3
                    elif phased_sv_haplotype_father_forward[sv_i-1][1] == 5 and phased_sv_haplotype_mother_forward[sv_i-1][1] == 5 and phase_type == 2:
                        phased_sv_haplotype_father_forward[sv_i-1][1] = 6
                        phased_sv_haplotype_mother_forward[sv_i-1][1] = 6
                continue
            candidate_single_SV_ls = candidate_single_SV_gt_fam_ls[family_index][sv_i]
            if phased_sv_haplotype_father_forward[sv_i][1] == 4 and len(phased_sv_haplotype_father_forward[sv_i][6]) == 0 and len(phased_sv_haplotype_father_forward[sv_i][8]) == 0 :
                if phased_sv_haplotype_father_forward[sv_i][0] < genotype_threshold :
                    phased_sv_haplotype_father_forward[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                    phased_sv_haplotype_father_forward[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                    phased_sv_haplotype_mother_forward[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                    phased_sv_haplotype_mother_forward[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                else :
                    phased_sv_haplotype_father_forward[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                    phased_sv_haplotype_father_forward[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                    phased_sv_haplotype_mother_forward[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                    phased_sv_haplotype_mother_forward[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                block_no += 1
            # 有可能杂合的位点处同时变异或同时非变异，这种情况可能需要处理
            elif len(phased_sv_haplotype_father_forward[sv_i][8]) == 0 or len(phased_sv_haplotype_father_forward[sv_i][6]) / len(phased_sv_haplotype_father_forward[sv_i][8]) >= 1 :
                #if len(phased_sv_haplotype_father_forward[sv_i][6]) / len(phased_sv_haplotype_father_forward[sv_i][8]) <= 1.5 :
                #    logging.info(phased_sv_haplotype_father_forward[sv_i])
                phased_sv_haplotype_father_forward[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_forward[sv_i][1] = phased_flag
                phased_sv_haplotype_father_forward[sv_i][3] = block_no
                phased_sv_haplotype_mother_forward[sv_i][3] = block_no
                phased_sv_haplotype_father_forward[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                phased_sv_haplotype_father_forward[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                phased_sv_haplotype_father_forward[sv_i][8] = []
                phased_sv_haplotype_father_forward[sv_i][9] = []
                phased_sv_haplotype_mother_forward[sv_i][6] = []
                phased_sv_haplotype_mother_forward[sv_i][7] = []
                phased_sv_haplotype_mother_forward[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                phased_sv_haplotype_mother_forward[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                if genotype_ls[0] > genotype_ls[1] :
                    phased_sv_haplotype_father_forward[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_father_forward[sv_i][2] = 0
                    phased_sv_haplotype_mother_forward[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_mother_forward[sv_i][2] = 0
                else :
                    phased_sv_haplotype_father_forward[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_father_forward[sv_i][2] = 1
                    phased_sv_haplotype_mother_forward[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_mother_forward[sv_i][2] = 1
            else :
                #if len(phased_sv_haplotype_father_forward[sv_i][8]) != 0 and len(phased_sv_haplotype_father_forward[sv_i][6]) / len(phased_sv_haplotype_father_forward[sv_i][8]) > 0.7 :
                #    logging.info(phased_sv_haplotype_father_forward[sv_i])
                phased_sv_haplotype_father_forward[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_forward[sv_i][1] = phased_flag
                phased_sv_haplotype_father_forward[sv_i][3] = block_no
                phased_sv_haplotype_mother_forward[sv_i][3] = block_no
                phased_sv_haplotype_father_forward[sv_i][6] = []
                phased_sv_haplotype_father_forward[sv_i][7] = []
                phased_sv_haplotype_father_forward[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                phased_sv_haplotype_father_forward[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                phased_sv_haplotype_mother_forward[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                phased_sv_haplotype_mother_forward[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                phased_sv_haplotype_mother_forward[sv_i][8] = []
                phased_sv_haplotype_mother_forward[sv_i][9] = []
                if genotype_ls[0] < genotype_ls[1] :
                    phased_sv_haplotype_father_forward[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_father_forward[sv_i][2] = 1
                    phased_sv_haplotype_mother_forward[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_mother_forward[sv_i][2] = 1
                else :
                    phased_sv_haplotype_father_forward[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_father_forward[sv_i][2] = 0
                    phased_sv_haplotype_mother_forward[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_mother_forward[sv_i][2] = 0
        else :
            if phase_type == 1 :
                block_no = phased_sv_haplotype_father_forward[sv_i][3]
        if phased_sv_haplotype_father_forward[sv_i][0] < genotype_threshold :
            fa_sv_read_names = phased_sv_haplotype_father_forward[sv_i][8]
            fa_sv_read_poss = phased_sv_haplotype_father_forward[sv_i][9]
            mo_sv_read_names = phased_sv_haplotype_mother_forward[sv_i][6]
            mo_sv_read_poss = phased_sv_haplotype_mother_forward[sv_i][7]
        else :
            fa_sv_read_names = phased_sv_haplotype_father_forward[sv_i][6]
            fa_sv_read_poss = phased_sv_haplotype_father_forward[sv_i][7]
            mo_sv_read_names = phased_sv_haplotype_mother_forward[sv_i][8]
            mo_sv_read_poss = phased_sv_haplotype_mother_forward[sv_i][9]
        #next_sv_read_names = set()
        #next_sv_read_poss = set()
        #if phase_type == 2 :
        #    logging.info(len(fa_sv_read_names))
        for read_j in range(len(fa_sv_read_names)) :
            read_name = fa_sv_read_names[read_j]
            read_pos = int(fa_sv_read_poss[read_j])
            if read_name in read_to_vars[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_vars[ceil(read_pos/read_pos_interval)-1][1][read_to_vars[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv < sv_i :
                        continue
                    if phased_sv_haplotype_father_forward[sv][1] in [0,4] and read_name not in phased_sv_haplotype_father_forward[sv][6] :
                        if phased_sv_haplotype_father_forward[sv][1] != 1:
                            phased_sv_haplotype_father_forward[sv][6].append(read_name)
                            phased_sv_haplotype_father_forward[sv][3] = block_no
                        phased_sv_haplotype_father_forward[sv][10].append(read_name)
            elif read_name in read_to_nons[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_nons[ceil(read_pos/read_pos_interval)-1][1][read_to_nons[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv < sv_i :    
                        continue
                    if phased_sv_haplotype_father_forward[sv][1] in [0,4] and read_name not in phased_sv_haplotype_father_forward[sv][8] :
                        if phased_sv_haplotype_father_forward[sv][1] != 1:
                            phased_sv_haplotype_father_forward[sv][8].append(read_name)
                            phased_sv_haplotype_father_forward[sv][3] = block_no
                        phased_sv_haplotype_father_forward[sv][11].append(read_name)

        for read_j in range(len(mo_sv_read_names)) :
            read_name = mo_sv_read_names[read_j]
            if len(mo_sv_read_names) != len(mo_sv_read_poss) :
                logging.info(len(mo_sv_read_names))
                logging.info(len(mo_sv_read_poss))
                logging.info(phased_sv_haplotype_mother_forward[sv_i])
                logging.info(candidate_single_SV_gt_fam_ls[family_index][sv_i])
            #logging.info("%d/%d"%(len(mo_sv_read_names),len(mo_sv_read_poss)))
            read_pos = int(mo_sv_read_poss[read_j])
            if read_name in read_to_vars[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_vars[ceil(read_pos/read_pos_interval)-1][1][read_to_vars[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv < sv_i :
                        continue
                    if phased_sv_haplotype_mother_forward[sv][1] in [0,4] and read_name not in phased_sv_haplotype_mother_forward[sv][6] :
                        if phased_sv_haplotype_mother_forward[sv][1] != 1:
                            phased_sv_haplotype_mother_forward[sv][6].append(read_name)
                            phased_sv_haplotype_mother_forward[sv][3] = block_no
                        phased_sv_haplotype_mother_forward[sv][10].append(read_name)
            elif read_name in read_to_nons[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_nons[ceil(read_pos/read_pos_interval)-1][1][read_to_nons[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv < sv_i :
                        continue
                    if phased_sv_haplotype_mother_forward[sv][1] in [0,4] and read_name not in phased_sv_haplotype_mother_forward[sv][8] :
                        if phased_sv_haplotype_mother_forward[sv][1] != 1:
                            phased_sv_haplotype_mother_forward[sv][8].append(read_name)
                            phased_sv_haplotype_mother_forward[sv][3] = block_no  
                        phased_sv_haplotype_mother_forward[sv][11].append(read_name)
    #logging.info("%d/%s/%d"%(phase_type,chr,homozygous_phased_num))
    # 反向
    for sv_i in range(len(phased_sv_haplotype_father_reverse)-1,-1,-1) :
        #if sv_i % 5000 == 0 :
        #    logging.info("%d/%d"%(sv_i,len(phased_sv_haplotype_father_reverse)))
        if phase_type == 2 and phased_sv_haplotype_father_reverse[sv_i][1] in [1,2,3] :
            continue
        genotype_ls = [float(x) for x in candidate_single_SV_gt_fam_ls[family_index][sv_i][gl_index].split(",")[3:5]]
        if gt_homozygous(genotype_ls) :
            if sv_i == len(phased_sv_haplotype_father_reverse)-1:
                #phased_sv_haplotype_father_reverse[sv_i][1] = 3
                pass
            elif ((len(phased_sv_haplotype_father_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_father_reverse[sv_i][8]) == 0) and
                  (len(phased_sv_haplotype_mother_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_reverse[sv_i][8]) == 0)) :
                if phased_sv_haplotype_father_reverse[sv_i+1][1] == 2 and phased_sv_haplotype_mother_reverse[sv_i+1][1] == 2 and phase_type == 1:
                    phased_sv_haplotype_father_reverse[sv_i+1][1] = 3
                    phased_sv_haplotype_mother_reverse[sv_i+1][1] = 3
                elif phased_sv_haplotype_father_reverse[sv_i+1][1] == 5 and phased_sv_haplotype_mother_reverse[sv_i+1][1] == 5 and phase_type == 2:
                    phased_sv_haplotype_father_reverse[sv_i+1][1] = 6
                    phased_sv_haplotype_mother_reverse[sv_i+1][1] = 6
            else :
                phased_sv_haplotype_father_reverse[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_reverse[sv_i][1] = phased_flag
            continue
        if phased_sv_haplotype_father_reverse[sv_i][1] in [0,4] :
            # 确定基因型父母本来源
            # 如果父母本的基因型确定结果有冲突，可能需要纠错
            if phased_sv_haplotype_father_reverse[sv_i][1] == 0 and (len(phased_sv_haplotype_father_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_father_reverse[sv_i][8]) == 0) and (len(phased_sv_haplotype_mother_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_reverse[sv_i][8]) == 0) :
                if sv_i == len(phased_sv_haplotype_father_reverse)-1:
                    #phased_sv_haplotype_father_reverse[sv_i][1] = 3
                    pass
                elif ((len(phased_sv_haplotype_father_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_father_reverse[sv_i][8]) == 0 and phased_sv_haplotype_father_reverse[sv_i+1][1] in [2,5]) and
                      (len(phased_sv_haplotype_mother_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_mother_reverse[sv_i][8]) == 0 and phased_sv_haplotype_mother_reverse[sv_i+1][1] in [2,5])) :
                    if phased_sv_haplotype_father_reverse[sv_i+1][1] == 2 and phased_sv_haplotype_mother_reverse[sv_i+1][1] == 2 and phase_type == 1:
                        phased_sv_haplotype_father_reverse[sv_i+1][1] = 3
                        phased_sv_haplotype_mother_reverse[sv_i+1][1] = 3
                    elif phased_sv_haplotype_father_reverse[sv_i+1][1] == 5 and phased_sv_haplotype_mother_reverse[sv_i+1][1] == 5 and phase_type == 2:
                        phased_sv_haplotype_father_reverse[sv_i+1][1] = 6
                        phased_sv_haplotype_mother_reverse[sv_i+1][1] = 6
                continue
            candidate_single_SV_ls = candidate_single_SV_gt_fam_ls[family_index][sv_i]
            if phased_sv_haplotype_father_reverse[sv_i][1] == 4 and len(phased_sv_haplotype_father_reverse[sv_i][6]) == 0 and len(phased_sv_haplotype_father_forward[sv_i][8]) == 0 :
                if phased_sv_haplotype_father_reverse[sv_i][0] < genotype_threshold :
                    phased_sv_haplotype_father_reverse[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                    phased_sv_haplotype_father_reverse[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                    phased_sv_haplotype_mother_reverse[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                    phased_sv_haplotype_mother_reverse[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                else :
                    phased_sv_haplotype_father_reverse[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                    phased_sv_haplotype_father_reverse[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                    phased_sv_haplotype_mother_reverse[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                    phased_sv_haplotype_mother_reverse[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                block_no += 1
            # 有可能杂合的位点处同时变异或同时非变异，这种情况可能需要处理
            if len(phased_sv_haplotype_father_reverse[sv_i][8]) == 0 or len(phased_sv_haplotype_father_reverse[sv_i][6]) / len(phased_sv_haplotype_father_reverse[sv_i][8]) >= 1 :
                #if len(phased_sv_haplotype_father_reverse[sv_i][6]) / len(phased_sv_haplotype_father_reverse[sv_i][8]) <= 1.5 :
                #    logging.info(phased_sv_haplotype_father_reverse[sv_i])
                phased_sv_haplotype_father_reverse[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_reverse[sv_i][1] = phased_flag
                phased_sv_haplotype_father_reverse[sv_i][3] = block_no
                phased_sv_haplotype_mother_reverse[sv_i][3] = block_no
                phased_sv_haplotype_father_reverse[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                phased_sv_haplotype_father_reverse[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                phased_sv_haplotype_father_reverse[sv_i][8] = []
                phased_sv_haplotype_father_reverse[sv_i][9] = []
                phased_sv_haplotype_mother_reverse[sv_i][6] = []
                phased_sv_haplotype_mother_reverse[sv_i][7] = []
                phased_sv_haplotype_mother_reverse[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                phased_sv_haplotype_mother_reverse[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                if genotype_ls[0] > genotype_ls[1] :
                    phased_sv_haplotype_father_reverse[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_father_reverse[sv_i][2] = 0
                    phased_sv_haplotype_mother_reverse[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_mother_reverse[sv_i][2] = 0
                else :
                    phased_sv_haplotype_father_reverse[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_father_reverse[sv_i][2] = 1
                    phased_sv_haplotype_mother_reverse[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_mother_reverse[sv_i][2] = 1
            else :
                #if len(phased_sv_haplotype_father_reverse[sv_i][8]) != 0 and len(phased_sv_haplotype_father_reverse[sv_i][6]) / len(phased_sv_haplotype_father_reverse[sv_i][8]) > 0.7 :
                #    logging.info(phased_sv_haplotype_father_reverse[sv_i])
                phased_sv_haplotype_father_reverse[sv_i][1] = phased_flag
                phased_sv_haplotype_mother_reverse[sv_i][1] = phased_flag
                phased_sv_haplotype_father_reverse[sv_i][3] = block_no
                phased_sv_haplotype_mother_reverse[sv_i][3] = block_no
                phased_sv_haplotype_father_reverse[sv_i][6] = []
                phased_sv_haplotype_father_reverse[sv_i][7] = []
                phased_sv_haplotype_father_reverse[sv_i][8] = candidate_single_SV_ls[non_read_nameds_index].split(",") if candidate_single_SV_ls[non_read_nameds_index] != '' else []
                phased_sv_haplotype_father_reverse[sv_i][9] = candidate_single_SV_ls[non_read_poss_index].split(",") if candidate_single_SV_ls[non_read_poss_index] != '' else []
                phased_sv_haplotype_mother_reverse[sv_i][6] = candidate_single_SV_ls[var_read_nameds_index].split(",") if candidate_single_SV_ls[var_read_nameds_index] != '' else []
                phased_sv_haplotype_mother_reverse[sv_i][7] = candidate_single_SV_ls[var_read_poss_index].split(",") if candidate_single_SV_ls[var_read_poss_index] != '' else []
                phased_sv_haplotype_mother_reverse[sv_i][8] = []
                phased_sv_haplotype_mother_reverse[sv_i][9] = []
                if genotype_ls[0] < genotype_ls[1] :
                    phased_sv_haplotype_father_reverse[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_father_reverse[sv_i][2] = 1
                    phased_sv_haplotype_mother_reverse[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_mother_reverse[sv_i][2] = 1
                else :
                    phased_sv_haplotype_father_reverse[sv_i][0] = genotype_ls[1]
                    phased_sv_haplotype_father_reverse[sv_i][2] = 0
                    phased_sv_haplotype_mother_reverse[sv_i][0] = genotype_ls[0]
                    phased_sv_haplotype_mother_reverse[sv_i][2] = 0
        else :
            if phase_type == 1 :
                block_no = phased_sv_haplotype_father_reverse[sv_i][3]
        if phased_sv_haplotype_father_reverse[sv_i][0] < genotype_threshold :
            fa_sv_read_names = phased_sv_haplotype_father_reverse[sv_i][8]
            fa_sv_read_poss = phased_sv_haplotype_father_reverse[sv_i][9]
            mo_sv_read_names = phased_sv_haplotype_mother_reverse[sv_i][6]
            mo_sv_read_poss = phased_sv_haplotype_mother_reverse[sv_i][7]
        else :
            fa_sv_read_names = phased_sv_haplotype_father_reverse[sv_i][6]
            fa_sv_read_poss = phased_sv_haplotype_father_reverse[sv_i][7]
            mo_sv_read_names = phased_sv_haplotype_mother_reverse[sv_i][8]
            mo_sv_read_poss = phased_sv_haplotype_mother_reverse[sv_i][9]
        #next_sv_read_names = set()
        #next_sv_read_poss = set()
        # key为sv索引值，value为二元列表，其中第一元为支持read name列表，第二元为不支持read name列表
        #support_reads_ls = {}
        for read_j in range(len(fa_sv_read_names)) :
            read_name = fa_sv_read_names[read_j]
            read_pos = int(fa_sv_read_poss[read_j])
            if read_name in read_to_vars[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_vars[ceil(read_pos/read_pos_interval)-1][1][read_to_vars[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv > sv_i :
                        continue
                    if phased_sv_haplotype_father_reverse[sv][1] in [0,4] and read_name not in phased_sv_haplotype_father_reverse[sv][6] :
                        if phased_sv_haplotype_father_reverse[sv][1] != 1:
                            phased_sv_haplotype_father_reverse[sv][6].append(read_name)
                            phased_sv_haplotype_father_reverse[sv][3] = block_no
                        phased_sv_haplotype_father_reverse[sv][10].append(read_name)
                        #if sv not in support_reads_ls :
                        #    support_reads_ls[sv] = [[read_name],[]]
                        #else :
                        #    support_reads_ls[sv][0].append(read_name)
            elif read_name in read_to_nons[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_nons[ceil(read_pos/read_pos_interval)-1][1][read_to_nons[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv > sv_i :
                        continue
                    if phased_sv_haplotype_father_reverse[sv][1] in [0,4] and read_name not in phased_sv_haplotype_father_reverse[sv][8] :
                        if phased_sv_haplotype_father_reverse[sv][1] != 1:
                            phased_sv_haplotype_father_reverse[sv][8].append(read_name)
                            phased_sv_haplotype_father_reverse[sv][3] = block_no
                        phased_sv_haplotype_father_reverse[sv][11].append(read_name)
                        #if sv not in support_reads_ls :
                        #    support_reads_ls[sv] = [[],[read_name]]
                        #else :
                        #    support_reads_ls[sv][1].append(read_name)
        #for sv in support_reads_ls.keys() :
        #    if abs(len(support_reads_ls[sv][0])-len(support_reads_ls[sv][1]))/max(len(support_reads_ls[sv][0]),len(support_reads_ls[sv][1])) <= 0.5 and len(support_reads_ls[sv][0]) != 0 and len(support_reads_ls[sv][1]) != 0 :
        #        #logging.info("%d/%d"%(len(support_reads_ls[sv][0]),len(support_reads_ls[sv][1])))
        #        if len(support_reads_ls[sv][0]) <= len(support_reads_ls[sv][1]) :
        #            for read_name in support_reads_ls[sv][0] :
        #                phased_sv_haplotype_father_reverse[sv][6].remove(read_name)
        #                phased_sv_haplotype_father_reverse[sv][10].remove(read_name)
        #                phased_sv_haplotype_father_reverse[sv][8].append(read_name)
        #                phased_sv_haplotype_father_reverse[sv][11].append(read_name)
        #        else :
        #            for read_name in support_reads_ls[sv][1] :
        #                phased_sv_haplotype_father_reverse[sv][8].remove(read_name)
        #                phased_sv_haplotype_father_reverse[sv][11].remove(read_name)
        #                phased_sv_haplotype_father_reverse[sv][6].append(read_name)
        #                phased_sv_haplotype_father_reverse[sv][10].append(read_name)
        #support_reads_ls = {}
        for read_j in range(len(mo_sv_read_names)) :
            read_name = mo_sv_read_names[read_j]
            read_pos = int(mo_sv_read_poss[read_j])
            if read_name in read_to_vars[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_vars[ceil(read_pos/read_pos_interval)-1][1][read_to_vars[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv > sv_i :
                        continue
                    if phased_sv_haplotype_mother_reverse[sv][1] in [0,4] and read_name not in phased_sv_haplotype_mother_reverse[sv][6] :
                        if phased_sv_haplotype_mother_reverse[sv][1] != 1:
                            phased_sv_haplotype_mother_reverse[sv][6].append(read_name)
                            phased_sv_haplotype_mother_reverse[sv][3] = block_no
                        phased_sv_haplotype_mother_reverse[sv][10].append(read_name)
                        #if sv not in support_reads_ls :
                        #    support_reads_ls[sv] = [[read_name],[]]
                        #else :
                        #    support_reads_ls[sv][0].append(read_name)
            elif read_name in read_to_nons[ceil(read_pos/read_pos_interval)-1][0] :
                support_svs = read_to_nons[ceil(read_pos/read_pos_interval)-1][1][read_to_nons[ceil(read_pos/read_pos_interval)-1][0].index(read_name)]
                for sv in support_svs :
                    if sv > sv_i :
                        continue
                    if phased_sv_haplotype_mother_reverse[sv][1] in [0,4] and read_name not in phased_sv_haplotype_mother_reverse[sv][8] :
                        if phased_sv_haplotype_mother_reverse[sv][1] != 1:
                            phased_sv_haplotype_mother_reverse[sv][8].append(read_name)
                            phased_sv_haplotype_mother_reverse[sv][3] = block_no
                        phased_sv_haplotype_mother_reverse[sv][11].append(read_name)
                        #if sv not in support_reads_ls :
                        #    support_reads_ls[sv] = [[],[read_name]]
                        #else :
                        #    support_reads_ls[sv][1].append(read_name)
        #for sv in support_reads_ls.keys() :
        #    if abs(len(support_reads_ls[sv][0])-len(support_reads_ls[sv][1]))/max(len(support_reads_ls[sv][0]),len(support_reads_ls[sv][1])) <= 0.5 and len(support_reads_ls[sv][0]) != 0 and len(support_reads_ls[sv][1]) != 0 :
        #        #logging.info("%d/%d"%(len(support_reads_ls[sv][0]),len(support_reads_ls[sv][1])))
        #        if len(support_reads_ls[sv][0]) <= len(support_reads_ls[sv][1]) :
        #            for read_name in support_reads_ls[sv][0] :
        #                phased_sv_haplotype_mother_reverse[sv][6].remove(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][10].remove(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][8].append(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][11].append(read_name)
        #        else :
        #            for read_name in support_reads_ls[sv][1] :
        #                phased_sv_haplotype_mother_reverse[sv][8].remove(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][11].remove(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][6].append(read_name)
        #                phased_sv_haplotype_mother_reverse[sv][10].append(read_name)

    # 整合父本/整合遗传单倍体
    sv_i = 0
    inconsistent_num = 0
    while(sv_i <= len(phased_sv_haplotype_father_forward)-1) :
        genotype_ls = [float(x) for x in candidate_single_SV_gt_fam_ls[family_index][sv_i][gl_index].split(",")[3:5]]
        if (phased_sv_haplotype_father_forward[sv_i][1] == 0 and phased_sv_haplotype_father_reverse[sv_i][1] == 0) or (phased_sv_haplotype_father_forward[sv_i][1] == 1 and phased_sv_haplotype_father_reverse[sv_i][1] == 1) :
            sv_i += 1
            continue
        if not gt_homozygous(genotype_ls) :
            if phased_sv_haplotype_father_forward[sv_i][1] in [2,3,4,5,6] and phased_sv_haplotype_father_reverse[sv_i][1] not in [2,3,4,5,6] :
                phased_sv_haplotype_father[sv_i] = copy.deepcopy(phased_sv_haplotype_father_forward[sv_i])
            elif phased_sv_haplotype_father_forward[sv_i][1] not in [2,3,4,5,6] and phased_sv_haplotype_father_reverse[sv_i][1] in [2,3,4,5,6] :
                phased_sv_haplotype_father[sv_i] = copy.deepcopy(phased_sv_haplotype_father_reverse[sv_i])
            elif phased_sv_haplotype_father_forward[sv_i][1] in [2,3,4,5,6] and phased_sv_haplotype_father_reverse[sv_i][1] in [2,3,4,5,6] :
                if abs(phased_sv_haplotype_father_forward[sv_i][0]-phased_sv_haplotype_father_reverse[sv_i][0]) > consistency_threshold :
                    # 理论上这里需要添加正反向结果不一致的处理过程，但是根据实际测试，几乎没有这种情况，所以不再处理
                    inconsistent_num += 1
                phased_sv_haplotype_father[sv_i] = copy.deepcopy(phased_sv_haplotype_father_forward[sv_i])
        else :
            if phased_sv_haplotype_father_forward[sv_i][1] in [2,3,4,5,6] :
                phased_sv_haplotype_father[sv_i] = copy.deepcopy(phased_sv_haplotype_father_forward[sv_i])
                if phased_sv_haplotype_father_reverse[sv_i][1] in [2,3,4,5,6] :
                    for read_i in range(len(phased_sv_haplotype_father_reverse[sv_i][6])) :
                        if phased_sv_haplotype_father_reverse[sv_i][6][read_i] not in phased_sv_haplotype_father[sv_i][6] :
                            phased_sv_haplotype_father[sv_i][6].append(phased_sv_haplotype_father_reverse[sv_i][6][read_i])
                    for read_i in range(len(phased_sv_haplotype_father_reverse[sv_i][10])) :
                        if phased_sv_haplotype_father_reverse[sv_i][10][read_i] not in phased_sv_haplotype_father[sv_i][10] :
                            phased_sv_haplotype_father[sv_i][10].append(phased_sv_haplotype_father_reverse[sv_i][10][read_i])
                    for read_i in range(len(phased_sv_haplotype_father_reverse[sv_i][8])) :
                        if phased_sv_haplotype_father_reverse[sv_i][8][read_i] not in phased_sv_haplotype_father[sv_i][8] :
                            phased_sv_haplotype_father[sv_i][8].append(phased_sv_haplotype_father_reverse[sv_i][8][read_i])
                    for read_i in range(len(phased_sv_haplotype_father_reverse[sv_i][11])) :
                        if phased_sv_haplotype_father_reverse[sv_i][11][read_i] not in phased_sv_haplotype_father[sv_i][11] :
                            phased_sv_haplotype_father[sv_i][11].append(phased_sv_haplotype_father_reverse[sv_i][11][read_i])
            else :
                phased_sv_haplotype_father[sv_i] = copy.deepcopy(phased_sv_haplotype_father_reverse[sv_i])
        if phased_sv_haplotype_father_forward[sv_i][1] in [4,5] and phased_sv_haplotype_father_reverse[sv_i][1] in [4,5] and phased_sv_haplotype_father_forward[sv_i][1] != phased_sv_haplotype_father_reverse[sv_i][1] :
            phased_sv_haplotype_father[sv_i][1] = 5
        sv_i += 1
    # 整合母本/整合弃绝单倍体
    sv_i = 0
    inconsistent_num = 0
    while(sv_i <= len(phased_sv_haplotype_mother_forward)-1) :
        genotype_ls = [float(x) for x in candidate_single_SV_gt_fam_ls[family_index][sv_i][gl_index].split(",")[3:5]]
        if (phased_sv_haplotype_mother_forward[sv_i][1] == 0 and phased_sv_haplotype_mother_reverse[sv_i][1] == 0) or (phased_sv_haplotype_mother_forward[sv_i][1] == 1 and phased_sv_haplotype_mother_reverse[sv_i][1] == 1) :
            sv_i += 1
            continue
        if not gt_homozygous(genotype_ls) :
            if phased_sv_haplotype_mother_forward[sv_i][1] in [2,3,4,5,6] and phased_sv_haplotype_mother_reverse[sv_i][1] not in [2,3,4,5,6] :
                phased_sv_haplotype_mother[sv_i] = copy.deepcopy(phased_sv_haplotype_mother_forward[sv_i])
            elif phased_sv_haplotype_mother_forward[sv_i][1] not in [2,3,4,5,6] and phased_sv_haplotype_mother_reverse[sv_i][1] in [2,3,4,5,6] :
                phased_sv_haplotype_mother[sv_i] = copy.deepcopy(phased_sv_haplotype_mother_reverse[sv_i])
            elif phased_sv_haplotype_mother_forward[sv_i][1] in [2,3,4,5,6] and phased_sv_haplotype_mother_reverse[sv_i][1] in [2,3,4,5,6] :
                if abs(phased_sv_haplotype_mother_forward[sv_i][0]-phased_sv_haplotype_mother_reverse[sv_i][0]) > consistency_threshold :
                    inconsistent_num += 1
                phased_sv_haplotype_mother[sv_i] = copy.deepcopy(phased_sv_haplotype_mother_forward[sv_i])
        else :
            if phased_sv_haplotype_mother_forward[sv_i][1] in [2,3,4,5,6] :
                phased_sv_haplotype_mother[sv_i] = copy.deepcopy(phased_sv_haplotype_mother_forward[sv_i])
                if phased_sv_haplotype_mother_reverse[sv_i][1] in [2,3,4,5,6] :
                    for read_i in range(len(phased_sv_haplotype_mother_reverse[sv_i][6])) :
                        if phased_sv_haplotype_mother_reverse[sv_i][6][read_i] not in phased_sv_haplotype_mother[sv_i][6] :
                            phased_sv_haplotype_mother[sv_i][6].append(phased_sv_haplotype_mother_reverse[sv_i][6][read_i])
                    for read_i in range(len(phased_sv_haplotype_mother_reverse[sv_i][10])) :
                        if phased_sv_haplotype_mother_reverse[sv_i][10][read_i] not in phased_sv_haplotype_mother[sv_i][10] :
                            phased_sv_haplotype_mother[sv_i][10].append(phased_sv_haplotype_mother_reverse[sv_i][10][read_i])
                    for read_i in range(len(phased_sv_haplotype_mother_reverse[sv_i][8])) :
                        if phased_sv_haplotype_mother_reverse[sv_i][8][read_i] not in phased_sv_haplotype_mother[sv_i][8] :
                            phased_sv_haplotype_mother[sv_i][8].append(phased_sv_haplotype_mother_reverse[sv_i][8][read_i])
                    for read_i in range(len(phased_sv_haplotype_mother_reverse[sv_i][11])) :
                        if phased_sv_haplotype_mother_reverse[sv_i][11][read_i] not in phased_sv_haplotype_mother[sv_i][11] :
                            phased_sv_haplotype_mother[sv_i][11].append(phased_sv_haplotype_mother_reverse[sv_i][11][read_i])
            else :
                phased_sv_haplotype_mother[sv_i] = copy.deepcopy(phased_sv_haplotype_mother_reverse[sv_i])
        if phased_sv_haplotype_mother_forward[sv_i][1] in [4,5] and phased_sv_haplotype_mother_reverse[sv_i][1] in [4,5] and phased_sv_haplotype_mother_forward[sv_i][1] != phased_sv_haplotype_mother_reverse[sv_i][1] :
            phased_sv_haplotype_mother[sv_i][1] = 5
        sv_i += 1
    # 确定纯和
    for sv_i in range(len(phased_sv_haplotype_father)) :
        genotype_ls = [float(x) for x in candidate_single_SV_gt_fam_ls[family_index][sv_i][gl_index].split(",")[3:5]]
        if phased_sv_haplotype_father[sv_i][1] == 0 :
            continue
        if genotype_ls[0] < genotype_threshold and genotype_ls[1] < genotype_threshold :
            if len(phased_sv_haplotype_father[sv_i][8]) > len(phased_sv_haplotype_mother[sv_i][8]) :
                phased_sv_haplotype_father[sv_i][0] = min(genotype_ls)
                phased_sv_haplotype_mother[sv_i][0] = max(genotype_ls)
            else :
                phased_sv_haplotype_father[sv_i][0] = max(genotype_ls)
                phased_sv_haplotype_mother[sv_i][0] = min(genotype_ls)
        if genotype_ls[0] > 1-genotype_threshold and genotype_ls[1] > 1-genotype_threshold :
            if len(phased_sv_haplotype_father[sv_i][6]) > len(phased_sv_haplotype_mother[sv_i][6]) :
                phased_sv_haplotype_father[sv_i][0] = max(genotype_ls)
                phased_sv_haplotype_mother[sv_i][0] = min(genotype_ls)
            else :
                phased_sv_haplotype_father[sv_i][0] = min(genotype_ls)
                phased_sv_haplotype_mother[sv_i][0] = max(genotype_ls)

    # 以下处理只考虑单独的杂合
    if phase_type == 2 :
        _ = recalculate_block_no(phased_sv_haplotype_father, phased_sv_haplotype_mother)
    
    # 对没有phase的变异，基于就近匹配的方式确定
    for sv_i in range(len(phased_sv_haplotype_father)) :
        if phased_sv_haplotype_father[sv_i][1] == 0 and phased_sv_haplotype_mother[sv_i][1] == 0 :
            fam_genotype_ls = []
            if family_mode == "M1" :
                fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][j][gl_index].split(",")[3:5]])
                fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][j][gl_index].split(",")[3:5]])
                fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[2][j][gl_index].split(",")[3:5]])
            elif family_mode == "M2" :
                fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][j][gl_index].split(",")[3:5]])
                fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][j][gl_index].split(",")[3:5]])
            nearest_match_res = family_nearest_match(fam_genotype_ls,family_mode,family_member)
            phased_sv_haplotype_father[sv_i][0] = nearest_match_res[0]
            phased_sv_haplotype_mother[sv_i][0] = nearest_match_res[1]
            phased_sv_haplotype_father[sv_i][1] = 7
            phased_sv_haplotype_mother[sv_i][1] = 7
            phased_sv_haplotype_father[sv_i][3] = -1
            phased_sv_haplotype_mother[sv_i][3] = -1
    
    return phased_sv_haplotype_father,phased_sv_haplotype_mother

def make_family_block_record(phased_sv_haplotype_child_father, phased_sv_haplotype_father_inher, phased_sv_haplotype_mother_inher) :
    block_record_ls = [[],[],[]]
    # 记录子代block信息
    start_pos = -1
    pre_block_no = 0
    for sv_i in range(len(phased_sv_haplotype_child_father)) :
        if sv_i == 0 :
            if phased_sv_haplotype_child_father[sv_i][3] not in [0,-1] :
                start_pos = sv_i
                pre_block_no = phased_sv_haplotype_child_father[sv_i][3]
            continue
        if phased_sv_haplotype_child_father[sv_i][3] == 0 :
            if start_pos != -1 :
                block_record_ls[0].append([pre_block_no,start_pos,sv_i-1])
                #logging.info("%d/%d/%d/%d"%(start_pos,sv_i-1,phased_sv_haplotype_child_father[start_pos][3],phased_sv_haplotype_child_father[sv_i-1][3]))
                start_pos = -1
                pre_block_no = 0
            continue
        elif phased_sv_haplotype_child_father[sv_i][3] == -1 :
            if start_pos != -1 :
                block_record_ls[0].append([pre_block_no,start_pos,sv_i-1])
                #logging.info("%d/%d/%d/%d"%(start_pos,sv_i-1,phased_sv_haplotype_child_father[start_pos][3],phased_sv_haplotype_child_father[sv_i-1][3]))
                start_pos = -1
                pre_block_no = 0
        elif phased_sv_haplotype_child_father[sv_i][3] > pre_block_no :
            if start_pos != -1 :
                block_record_ls[0].append([pre_block_no,start_pos,sv_i-1])
                #logging.info("%d/%d/%d/%d"%(start_pos,sv_i-1,phased_sv_haplotype_child_father[start_pos][3],phased_sv_haplotype_child_father[sv_i-1][3]))
            start_pos = sv_i
            pre_block_no = phased_sv_haplotype_child_father[sv_i][3]
        #logging.info("%d/%d/%d/%d"%(start_pos,sv_i,phased_sv_haplotype_child_father[start_pos][3],phased_sv_haplotype_child_father[sv_i-1][3]))
    # 记录父本block信息
    start_pos = -1
    pre_block_no = 0
    for sv_i in range(len(phased_sv_haplotype_father_inher)) :
        if sv_i == 0 :
            if phased_sv_haplotype_father_inher[sv_i][3] not in [0,-1] :
                start_pos = sv_i
                pre_block_no = phased_sv_haplotype_father_inher[sv_i][3]
            continue
        if phased_sv_haplotype_father_inher[sv_i][3] == 0 :
            if start_pos != -1 :
                block_record_ls[1].append([pre_block_no,start_pos,sv_i-1])
                start_pos = -1
                pre_block_no = 0
            continue
        elif phased_sv_haplotype_father_inher[sv_i][3] == -1 :
            if start_pos != -1 :
                block_record_ls[1].append([pre_block_no,start_pos,sv_i-1])
                start_pos = -1
                pre_block_no = 0
        elif phased_sv_haplotype_father_inher[sv_i][3] > pre_block_no :
            if start_pos != -1 :
                block_record_ls[1].append([pre_block_no,start_pos,sv_i-1])
            start_pos = sv_i
            pre_block_no = phased_sv_haplotype_father_inher[sv_i][3]
    # 记录母本block信息
    start_pos = -1
    pre_block_no = 0
    for sv_i in range(len(phased_sv_haplotype_mother_inher)) :
        if sv_i == 0 :
            if phased_sv_haplotype_mother_inher[sv_i][3] not in [0,-1] :
                start_pos = sv_i
                pre_block_no = phased_sv_haplotype_mother_inher[sv_i][3]
            continue
        if phased_sv_haplotype_mother_inher[sv_i][3] == 0 :
            if start_pos != -1 :
                block_record_ls[2].append([pre_block_no,start_pos,sv_i-1])
                start_pos = -1
                pre_block_no = 0
            continue
        elif phased_sv_haplotype_mother_inher[sv_i][3] == -1 :
            if start_pos != -1 :
                block_record_ls[2].append([pre_block_no,start_pos,sv_i-1])
                start_pos = -1
                pre_block_no = 0
        elif phased_sv_haplotype_mother_inher[sv_i][3] > pre_block_no :
            if start_pos != -1 :
                block_record_ls[2].append([pre_block_no,start_pos,sv_i-1])
            start_pos = sv_i
            pre_block_no = phased_sv_haplotype_mother_inher[sv_i][3]
    
    return block_record_ls

def change_haplotype_res(phased_sv_haplotype_1, phased_sv_haplotype_2, is_exchange, block_no, new_phase_state, new_block_no) :
    for sv_i in range(len(phased_sv_haplotype_1)) :
        if phased_sv_haplotype_1[sv_i][3] == block_no :
            if is_exchange :
                t = phased_sv_haplotype_1[sv_i]
                phased_sv_haplotype_1[sv_i] = phased_sv_haplotype_2[sv_i]
                phased_sv_haplotype_2[sv_i] = t
            phased_sv_haplotype_1[sv_i][1] = new_phase_state
            phased_sv_haplotype_2[sv_i][1] = new_phase_state
            phased_sv_haplotype_1[sv_i][3] = new_block_no
            phased_sv_haplotype_2[sv_i][3] = new_block_no

def genetic_phasing_optimize(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, family_mode) :
    if family_mode != "M1" :
        return
    gl_index = 9
    # 使用父母本pt=1的phase信息优化子代，使用子代pt=1的phase信息优化父母本
    change_num_ls = [0,0,0,0]
    for sv_i in range(len(candidate_single_SV_gt_fam_ls[0])) :
        fam_genotype_ls = []
        if family_mode == "M1" :
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][sv_i][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][sv_i][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[2][sv_i][gl_index].split(",")[3:5]])
        if phased_sv_haplotype_father_inher[sv_i][1] in [1,2,3] and not gt_homozygous(fam_genotype_ls[1]) and phased_sv_haplotype_child_father[sv_i][1] not in [1,2,3] and not gt_homozygous(fam_genotype_ls[0]) :
            change_num_ls[0] += 1
            block_no = phased_sv_haplotype_child_father[sv_i][3]
            if abs(phased_sv_haplotype_father_inher[sv_i][0]-phased_sv_haplotype_child_father[sv_i][0]) < consistency_threshold :
                change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_no,2,0)
            else :
                change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_no,2,0)
        if phased_sv_haplotype_mother_inher[sv_i][1] in [1,2,3] and not gt_homozygous(fam_genotype_ls[2]) and phased_sv_haplotype_child_mother[sv_i][1] not in [1,2,3] and not gt_homozygous(fam_genotype_ls[0]) :
            change_num_ls[1] += 1
            block_no = phased_sv_haplotype_child_mother[sv_i][3]
            if abs(phased_sv_haplotype_mother_inher[sv_i][0]-phased_sv_haplotype_child_mother[sv_i][0]) < consistency_threshold :
                change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_no,2,0)
            else :
                change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_no,2,0)
        if phased_sv_haplotype_child_father[sv_i][1] in [1,2,3] and not gt_homozygous(fam_genotype_ls[0]) and phased_sv_haplotype_father_inher[sv_i][1] not in [1,2,3] and not gt_homozygous(fam_genotype_ls[1]) :
            change_num_ls[2] += 1
            block_no = phased_sv_haplotype_father_inher[sv_i][3]
            if abs(phased_sv_haplotype_child_father[sv_i][0]-phased_sv_haplotype_father_inher[sv_i][0]) < consistency_threshold :
                change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,False,block_no,2,0)
            else :
                change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,True,block_no,2,0)
        if phased_sv_haplotype_child_mother[sv_i][1] in [1,2,3] and not gt_homozygous(fam_genotype_ls[0]) and phased_sv_haplotype_mother_inher[sv_i][1] not in [1,2,3] and not gt_homozygous(fam_genotype_ls[2]) :
            change_num_ls[3] += 1
            block_no = phased_sv_haplotype_mother_inher[sv_i][3]
            if abs(phased_sv_haplotype_child_mother[sv_i][0]-phased_sv_haplotype_mother_inher[sv_i][0]) < consistency_threshold :
                change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,False,block_no,2,0)
            else :
                change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,True,block_no,2,0)
    # 使用pt=2的phase信息优化家系的phase效果
    exchange_happend_flag = True
    #exchange_num = 0
    while(exchange_happend_flag) :
        exchange_num = 0
        exchange_happend_flag = False
        block_record_ls = make_family_block_record(phased_sv_haplotype_child_father, phased_sv_haplotype_father_inher, phased_sv_haplotype_mother_inher)
        # 使用父母本pt=2的phase信息优化子代
        for b_i in range(len(block_record_ls[0])-1) :
            if phased_sv_haplotype_father_inher[block_record_ls[0][b_i][2]][3] == phased_sv_haplotype_father_inher[block_record_ls[0][b_i+1][1]][3] and phased_sv_haplotype_father_inher[block_record_ls[0][b_i][2]][3] not in [-1,0]:
                #logging.info("%d/%d/%d/%d"%(block_record_ls[0][b_i][2],block_record_ls[0][b_i+1][1],phased_sv_haplotype_child_father[block_record_ls[0][b_i][2]][3],phased_sv_haplotype_child_father[block_record_ls[0][b_i+1][1]][3]))
                exchange_happend_flag = True
                exchange_num += 1
                new_block_no = block_record_ls[0][b_i][0]
                if abs(phased_sv_haplotype_child_father[block_record_ls[0][b_i][2]][0]-phased_sv_haplotype_father_inher[block_record_ls[0][b_i][2]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_record_ls[0][b_i][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_record_ls[0][b_i][0],5,new_block_no)
                if abs(phased_sv_haplotype_child_father[block_record_ls[0][b_i+1][1]][0]-phased_sv_haplotype_father_inher[block_record_ls[0][b_i+1][1]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_record_ls[0][b_i+1][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_record_ls[0][b_i+1][0],5,new_block_no)
                for j in range(block_record_ls[0][b_i][2]+1,block_record_ls[0][b_i+1][1]) :
                    phased_sv_haplotype_child_father[j][1] = 5
                    phased_sv_haplotype_child_mother[j][1] = 5
                    phased_sv_haplotype_child_father[j][3] = new_block_no
                    phased_sv_haplotype_child_mother[j][3] = new_block_no
                
            elif phased_sv_haplotype_mother_inher[block_record_ls[0][b_i][2]][3] == phased_sv_haplotype_mother_inher[block_record_ls[0][b_i+1][1]][3] and phased_sv_haplotype_mother_inher[block_record_ls[0][b_i][2]][3] not in [-1,0]:
                exchange_happend_flag = True
                exchange_num += 1
                new_block_no = block_record_ls[0][b_i][0]
                if abs(phased_sv_haplotype_child_mother[block_record_ls[0][b_i][2]][0]-phased_sv_haplotype_mother_inher[block_record_ls[0][b_i][2]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_record_ls[0][b_i][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_record_ls[0][b_i][0],5,new_block_no)
                if abs(phased_sv_haplotype_child_mother[block_record_ls[0][b_i+1][1]][0]-phased_sv_haplotype_mother_inher[block_record_ls[0][b_i+1][1]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,False,block_record_ls[0][b_i+1][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_child_father,phased_sv_haplotype_child_mother,True,block_record_ls[0][b_i+1][0],5,new_block_no)
                for j in range(block_record_ls[0][b_i][2]+1,block_record_ls[0][b_i+1][1]) :
                    phased_sv_haplotype_child_father[j][1] = 5
                    phased_sv_haplotype_child_mother[j][1] = 5
                    phased_sv_haplotype_child_father[j][3] = new_block_no
                    phased_sv_haplotype_child_mother[j][3] = new_block_no
        _ = recalculate_block_no(phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother)
        _ = recalculate_block_no(phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo)
        _ = recalculate_block_no(phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo)
        block_record_ls = make_family_block_record(phased_sv_haplotype_child_father, phased_sv_haplotype_father_inher, phased_sv_haplotype_mother_inher)
        # 使用子代pt=2的phase信息优化父本
        for b_i in range(len(block_record_ls[1])-1) :
            if phased_sv_haplotype_child_father[block_record_ls[1][b_i][2]][3] == phased_sv_haplotype_child_father[block_record_ls[1][b_i+1][1]][3] and phased_sv_haplotype_child_father[block_record_ls[1][b_i][2]][3] not in [-1,0]:
                exchange_happend_flag = True
                exchange_num += 1
                new_block_no = block_record_ls[1][b_i][0]
                if abs(phased_sv_haplotype_father_inher[block_record_ls[1][b_i][2]][0]-phased_sv_haplotype_child_father[block_record_ls[1][b_i][2]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,False,block_record_ls[1][b_i][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,True,block_record_ls[1][b_i][0],5,new_block_no)
                if abs(phased_sv_haplotype_father_inher[block_record_ls[1][b_i+1][1]][0]-phased_sv_haplotype_child_father[block_record_ls[1][b_i+1][1]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,False,block_record_ls[1][b_i+1][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_father_inher,phased_sv_haplotype_father_forgo,True,block_record_ls[1][b_i+1][0],5,new_block_no)
                for j in range(block_record_ls[1][b_i][2]+1,block_record_ls[1][b_i+1][1]) :
                    phased_sv_haplotype_father_inher[j][1] = 5
                    phased_sv_haplotype_father_forgo[j][1] = 5
                    phased_sv_haplotype_father_inher[j][3] = new_block_no
                    phased_sv_haplotype_father_forgo[j][3] = new_block_no
        # 使用子代pt=2的phase信息优化母本
        for b_i in range(len(block_record_ls[2])-1) :
            if phased_sv_haplotype_child_mother[block_record_ls[2][b_i][2]][3] == phased_sv_haplotype_child_mother[block_record_ls[2][b_i+1][1]][3] and phased_sv_haplotype_child_mother[block_record_ls[2][b_i][2]][3] not in [-1,0]:
                exchange_happend_flag = True
                exchange_num += 1
                new_block_no = block_record_ls[2][b_i][0]
                if abs(phased_sv_haplotype_mother_inher[block_record_ls[2][b_i][2]][0]-phased_sv_haplotype_child_mother[block_record_ls[2][b_i][2]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,False,block_record_ls[2][b_i][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,True,block_record_ls[2][b_i][0],5,new_block_no)
                if abs(phased_sv_haplotype_mother_inher[block_record_ls[2][b_i+1][1]][0]-phased_sv_haplotype_child_mother[block_record_ls[2][b_i+1][1]][0]) < consistency_threshold :
                    change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,False,block_record_ls[2][b_i+1][0],5,new_block_no)
                else :
                    change_haplotype_res(phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_forgo,True,block_record_ls[2][b_i+1][0],5,new_block_no)
                for j in range(block_record_ls[2][b_i][2]+1,block_record_ls[2][b_i+1][1]) :
                    phased_sv_haplotype_mother_inher[j][1] = 5
                    phased_sv_haplotype_mother_forgo[j][1] = 5
                    phased_sv_haplotype_mother_inher[j][3] = new_block_no
                    phased_sv_haplotype_mother_forgo[j][3] = new_block_no
        _ = recalculate_block_no(phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother)
        _ = recalculate_block_no(phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo)
        _ = recalculate_block_no(phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo)
        #logging.info("%s/%d"%(chr,exchange_num))
        #break

def evaluate_phasing(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_father, phased_sv_haplotype_mother, family_mode, family_member) :
    family_index = int(family_member) - 1
    gl_index = 9
    heterozygous_num = 0
    confirm_num = 0
    heterozygous_phase_1 = 0
    heterozygous_phase_2 = 0
    max_block = 0
    for j in range(len(candidate_single_SV_gt_fam_ls[family_index])) :
        fam_genotype_ls = []
        if family_mode == "M1" :
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[0][j][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[1][j][gl_index].split(",")[3:5]])
            fam_genotype_ls.append([float(x) for x in candidate_single_SV_gt_fam_ls[2][j][gl_index].split(",")[3:5]])
        if phased_sv_haplotype_father[j][3] > max_block :
            max_block = phased_sv_haplotype_father[j][3]
        if not gt_homozygous(fam_genotype_ls[family_index]) :
            heterozygous_num += 1
            if phased_sv_haplotype_father[j][1] in [1,2,3] :
                heterozygous_phase_1 += 1
            if phased_sv_haplotype_father[j][1] in [4,5,6] :
                heterozygous_phase_2 += 1
            is_con, gen_source = confirm_haplotype_source(fam_genotype_ls,family_mode,family_member)
            if is_con :
                confirm_num += 1
    if heterozygous_num != 0 and chr in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]:
        logging.info("%s/%s/%d/%f/%f/%f/%d"%(family_member,chr,heterozygous_num,confirm_num/heterozygous_num,heterozygous_phase_1/heterozygous_num,heterozygous_phase_2/heterozygous_num,max_block))
    
def verify_haplotype_consistency(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, family_mode, family_member) :
    if chr != '1' :
        return
    if family_member == "1" :
        for sv_i in range(0,100) :
            if abs(phased_sv_haplotype_child_father[sv_i][0]-phased_sv_haplotype_father_inher[sv_i][0]) > consistency_threshold or abs(phased_sv_haplotype_child_mother[sv_i][0]-phased_sv_haplotype_mother_inher[sv_i][0]) > consistency_threshold :
                #logging.info("%s/%f/%f/%f/%f/%d"%(chr,phased_sv_haplotype_child_father[sv_i][0],phased_sv_haplotype_father_inher[sv_i][0],phased_sv_haplotype_child_mother[sv_i][0],phased_sv_haplotype_mother_inher[sv_i][0],sv_i))
                logging.info("%f/%f/%f/%f/%f/%f"%(phased_sv_haplotype_child_father[sv_i][0],phased_sv_haplotype_child_mother[sv_i][0],phased_sv_haplotype_father_inher[sv_i][0],phased_sv_haplotype_father_forgo[sv_i][0],phased_sv_haplotype_mother_inher[sv_i][0],phased_sv_haplotype_mother_forgo[sv_i][0]))
                logging.info("%s/%s/%s/%s"%(str(phased_sv_haplotype_child_father[sv_i]),str(phased_sv_haplotype_child_mother[sv_i]),str(phased_sv_haplotype_father_inher[sv_i]),str(phased_sv_haplotype_mother_inher[sv_i])))
                logging.info("%s/%s/%s"%(str(candidate_single_SV_gt_fam_ls[0][sv_i]),str(candidate_single_SV_gt_fam_ls[1][sv_i]),str(candidate_single_SV_gt_fam_ls[2][sv_i])))

# 将phasing结果转移到变异结果列表中
def phasing_candidate_fam_SV(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, family_mode, phase_flag, parents_phasing) :
    gt_index = 8
    gl_index = 9
    if family_mode == "M1" :
        family_member_ls = ["1","2","3"]
        for j in range(len(candidate_single_SV_gt_fam_ls[0])) :
            if phase_flag :
                single_SV_gt = candidate_single_SV_gt_fam_ls[0][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":" + str(round(phased_sv_haplotype_child_father[j][0])) + "|" + str(round(phased_sv_haplotype_child_mother[j][0]))
                single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + "," + str(phased_sv_haplotype_child_father[j][0]) + "," + str(round(phased_sv_haplotype_child_mother[j][0])) + "," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                if parents_phasing :
                    single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":" + str(round(phased_sv_haplotype_father_inher[j][0])) + "|" + str(round(phased_sv_haplotype_father_forgo[j][0]))
                    single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + "," + str(phased_sv_haplotype_father_inher[j][0]) + "," + str(round(phased_sv_haplotype_father_forgo[j][0])) + "," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                    single_SV_gt = candidate_single_SV_gt_fam_ls[2][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":" + str(round(phased_sv_haplotype_mother_inher[j][0])) + "|" + str(round(phased_sv_haplotype_mother_forgo[j][0]))
                    single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + "," + str(phased_sv_haplotype_mother_inher[j][0]) + "," + str(round(phased_sv_haplotype_mother_forgo[j][0])) + "," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                else :
                    single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":."
                    #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + ",.,.," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                    single_SV_gt = candidate_single_SV_gt_fam_ls[2][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":."
                    #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + ",.,.," + ",".join(single_SV_gt[gl_index].split(",")[5:])
            else :
                single_SV_gt = candidate_single_SV_gt_fam_ls[0][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":./."
                #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:5]) + ",".join(single_SV_gt[gl_index].split(",")[5:]) + ",0"
                single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":./."
                #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:5]) + ",".join(single_SV_gt[gl_index].split(",")[5:]) + ",0"
                single_SV_gt = candidate_single_SV_gt_fam_ls[2][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":./."
                #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:5]) + ",".join(single_SV_gt[gl_index].split(",")[5:]) + ",0"
    elif family_mode == "M2" :
        family_member_ls = ["1","2"]
        for j in range(len(candidate_single_SV_gt_fam_ls[0])) :
            if phase_flag :
                single_SV_gt = candidate_single_SV_gt_fam_ls[0][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":" + str(round(phased_sv_haplotype_child_father[j][0])) + "|" + str(round(phased_sv_haplotype_child_mother[j][0]))
                single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + "," + str(phased_sv_haplotype_child_father[j][0]) + "," + str(round(phased_sv_haplotype_child_mother[j][0])) + "," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                if parents_phasing :
                    single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":" + str(round(phased_sv_haplotype_father_inher[j][0])) + "|" + str(round(phased_sv_haplotype_father_forgo[j][0]))
                    single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:3]) + "," + str(phased_sv_haplotype_father_inher[j][0]) + "," + str(round(phased_sv_haplotype_father_forgo[j][0])) + "," + ",".join(single_SV_gt[gl_index].split(",")[5:])
                else :
                    single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                    single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":."
            else :
                single_SV_gt = candidate_single_SV_gt_fam_ls[0][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":./."
                #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:5]) + ",".join(single_SV_gt[gl_index].split(",")[5:]) + ",0"
                single_SV_gt = candidate_single_SV_gt_fam_ls[1][j]
                single_SV_gt[gt_index] = single_SV_gt[gt_index] + ":./."
                #single_SV_gt[gl_index] = ",".join(single_SV_gt[gl_index].split(",")[0:5]) + ",".join(single_SV_gt[gl_index].split(",")[5:]) + ",0"

def correct_gt_phasing_2_calling(candidate_single_SV,haplotype_father_single,haplotype_mother_single,f_prob,m_prob) :
    len_index = 3
    gt_index = 8
    gl_index = 9
    gq_index = 10
    qual_index = 11
    prob = list(normalize_log10_probs([log10(max((1-f_prob)*(1-m_prob),pow(0.1,100))), log10(max((1-f_prob)*m_prob+f_prob*(1-m_prob),pow(0.1,100))), log10(max(f_prob*m_prob,pow(0.1,100)))]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
    candidate_single_SV[gt_index] = str(max(round(f_prob),round(m_prob)))+"/"+str(min(round(f_prob),round(m_prob)))+":"+str(round(f_prob))+"|"+str(round(m_prob))
    candidate_single_SV[gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(PL[0], PL[1], PL[2], f_prob, m_prob, int(candidate_single_SV[gl_index].split(",")[5]), len(haplotype_father_single[11])+len(haplotype_mother_single[11]), len(haplotype_father_single[10])+len(haplotype_mother_single[10]))
    candidate_single_SV[gq_index] = str(max(GQ))
    candidate_single_SV[qual_index] = str(QUAL)

def phasing_sv_error_detection_correction(chr, candidate_single_SV_gt_fam_ls, phased_sv_haplotype_child_father, phased_sv_haplotype_child_mother, phased_sv_haplotype_father_inher, phased_sv_haplotype_father_forgo, phased_sv_haplotype_mother_inher, phased_sv_haplotype_mother_forgo, family_mode, minimum_support_reads_list) :
    change_threshold = 0.5
    bias_threshold = 0.15
    len_index = 3
    gt_index = 8
    gl_index = 9
    gq_index = 10
    qual_index = 11
    if family_mode == "M1" :
        haplotype_father_ls = [phased_sv_haplotype_child_father,phased_sv_haplotype_father_inher,phased_sv_haplotype_father_inher]
        haplotype_mother_ls = [phased_sv_haplotype_child_mother,phased_sv_haplotype_mother_inher,phased_sv_haplotype_mother_inher]
    elif family_mode == "M2" :
        haplotype_father_ls = [phased_sv_haplotype_child_father,phased_sv_haplotype_father_inher]
        haplotype_mother_ls = [phased_sv_haplotype_child_mother,phased_sv_haplotype_mother_inher]

    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        haplotype_father = haplotype_father_ls[i]
        haplotype_mother = haplotype_mother_ls[i]
        for j in range(len(candidate_single_SV_gt_fam_ls[i])) :
            if (len(haplotype_father[j][10]) != 0 or len(haplotype_father[j][11]) != 0) and max([int(x) for x in candidate_single_SV_gt_fam_ls[i][j][gl_index].split(",")[0:3]]) != 100 :
                is_change_f_flag = False
                is_change_m_flag = False
                _, GL_P_F_1, _, _ = cal_GL_2(len(haplotype_father[j][11]), len(haplotype_father[j][10]))
                _, GL_P_M_1, _, _ = cal_GL_2(len(haplotype_mother[j][11]), len(haplotype_mother[j][10]))
                if len(haplotype_father[j][10]) != len(haplotype_father[j][11]) and max(len(haplotype_father[j][10]),len(haplotype_father[j][11])) > minimum_support_reads_list[i]*1 :
                    if round(GL_P_F_1) == 0 and round(haplotype_father[j][0]) == 1 and abs(GL_P_F_1-haplotype_father[j][0]) >= change_threshold:
                        is_change_f_flag = True
                if len(haplotype_mother[j][10]) != len(haplotype_mother[j][11]) and max(len(haplotype_mother[j][10]),len(haplotype_mother[j][11])) > minimum_support_reads_list[i]*1 :
                    if round(GL_P_M_1) == 0 and round(haplotype_mother[j][0]) == 1 and abs(GL_P_M_1-haplotype_mother[j][0]) >= change_threshold:
                        is_change_m_flag = True
                if (is_change_f_flag or is_change_m_flag) and abs(int(candidate_single_SV_gt_fam_ls[i][j][len_index])) > 100:
                    #logging.info(candidate_single_SV_gt_fam_ls[i][j])
                    if is_change_f_flag and not is_change_m_flag :
                        f_prob = GL_P_F_1
                        m_prob = haplotype_mother[j][0]
                    elif not is_change_f_flag and is_change_m_flag :
                        f_prob = haplotype_father[j][0]
                        m_prob = GL_P_M_1
                    else :
                        f_prob = GL_P_F_1
                        m_prob = GL_P_M_1
                    #logging.info("%f/%f/%f"%(max((1-f_prob)*(1-m_prob),pow(0.1,100)), max((1-f_prob)*m_prob+f_prob*(1-m_prob),pow(0.1,100)), max(f_prob*m_prob,pow(0.1,100))))
                    correct_gt_phasing_2_calling(candidate_single_SV_gt_fam_ls[i][j],haplotype_father[j],haplotype_mother[j],f_prob,m_prob)
                    continue
                #if haplotype_father[j][0]<0.5 and abs(GL_P_F_1-haplotype_father[j][0])<=bias_threshold and round(GL_P_F_1)==round(haplotype_father[j][0]) and len(haplotype_father[j][11])>minimum_support_reads_list[i] :
                #    if round(GL_P_M_1) == 0 and round(haplotype_mother[j][0]) == 1 and abs(GL_P_F_1-GL_P_M_1)<=bias_threshold and len(haplotype_mother[j][11])>minimum_support_reads_list[i] :
                #        f_prob,m_prob = haplotype_father[j][0],GL_P_M_1
                #        correct_gt_phasing_2_calling(candidate_single_SV_gt_fam_ls[i][j],haplotype_father[j],haplotype_mother[j],f_prob,m_prob)
                #        continue
                #if haplotype_mother[j][0]<0.5 and abs(GL_P_M_1-haplotype_mother[j][0])<=bias_threshold and round(GL_P_M_1)==round(haplotype_mother[j][0]) and len(haplotype_mother[j][11])>minimum_support_reads_list[i] :
                #    if round(GL_P_F_1) == 0 and round(haplotype_father[j][0]) == 1 and abs(GL_P_F_1-GL_P_M_1)<=bias_threshold and len(haplotype_father[j][11])>minimum_support_reads_list[i] :
                #        f_prob,m_prob = GL_P_F_1,haplotype_mother[j][0]
                #        correct_gt_phasing_2_calling(candidate_single_SV_gt_fam_ls[i][j],haplotype_father[j],haplotype_mother[j],f_prob,m_prob)
                #        continue
                
                





def run_phasing(args) :
    #logging.info("1")
    return genetic_phasing_family(*args)