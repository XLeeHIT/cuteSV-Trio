from cuteSVTrio.cuteSVTrio_genotype import cal_GL_3, threshold_ref_count, count_coverage
import numpy as np
import logging
import pickle
import copy
import random
from math import log10

gap_thres = 0.2 # 这个值是确定两个单支基因是否一致的阈值
mutation_rate = 0.00000001
err = 0.1
chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

# 确定该SV是否为denovo变异，如果是，返回True
def resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, nearby_matching, minimum_support_reads_list) :
    #pass
    gl_index = 9
    gt_index = 8
    for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
        family_gt_ls = []
        family_gl_ls = []
        for j in range(len(candidate_single_SV_gt_fam_ls)) :
            family_gt_ls.append([float(s) for s in candidate_single_SV_gt_fam_ls[j][i][gl_index].split(",")[3:5]])
            family_gl_ls.append([float(s) for s in candidate_single_SV_gt_fam_ls[j][i][gl_index].split(",")[5:8]])
        # 要加上对孟德尔一致性的检验，因为不是所有的突变都引起孟德尔失序性，但是本研究中只讨论引起孟德尔失序性的突变，这一点值得讨论
        # 如果有phasing过程，则可以严格复现孟德尔遗传过程检测突变
        # 如果没有phasing过程，则只能通过孟德尔失序性检测突变，即就近匹配
        # 无论是否phasing，通过该函数的都认为是突变
        if mendel_assert(family_gt_ls, family_mode, nearby_matching) :
            #pass
            continue
        if resolution_denovo(family_gt_ls, family_gl_ls[0], family_mode) :
            gl_split_ls = candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")
            fa_split_ls = candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")
            mo_split_ls = candidate_single_SV_gt_fam_ls[2][i][gl_index].split(",")
            family_genotype_ls = []
            for j in range(len(candidate_single_SV_gt_fam_ls)) :
                family_genotype_ls.append(candidate_single_SV_gt_fam_ls[j][i][gt_index][0:3])
            # denovo=1:0/0->0/1
            # denovo=2:0/1->1/1
            # denovo=3:1/1->0/1
            # denovo=4:0/1->0/0
            if gl_split_ls[8] in ["-1","-2","-3"] or fa_split_ls[8] in ["-1","-2","-3"] or mo_split_ls[8] in ["-1","-2","-3"] :
                continue
            if family_genotype_ls[0] in ["0/1","1/0"] and "0/0" == family_genotype_ls[1] and "0/0" == family_genotype_ls[2] :
                if int(float(gl_split_ls[7])) >= minimum_support_reads_list[0] and int(float(fa_split_ls[7])) <= 0 and int(float(mo_split_ls[7])) <= 0 :
                    gl_split_ls[8] = "1"
            elif "1/1" == family_genotype_ls[0] and ("0/0" == family_genotype_ls[1] or "0/0" == family_genotype_ls[2]):
                if 0 <= int(float(gl_split_ls[6])) <= 2 and (0 <= int(float(fa_split_ls[7])) <= 3 or 0 <= int(float(mo_split_ls[7])) <= 3) :
                    gl_split_ls[8] = "2"
            elif family_genotype_ls[0] in ["0/1","1/0"] and "1/1" == family_genotype_ls[1] and "1/1" == family_genotype_ls[2] :
                if int(float(gl_split_ls[7])) >= minimum_support_reads_list[0] and 0 <= int(float(fa_split_ls[6])) <= 4 and int(float(fa_split_ls[7])) >= minimum_support_reads_list[1] and 0 <= int(float(mo_split_ls[6])) <= 4 and int(float(mo_split_ls[7])) >= minimum_support_reads_list[2] :
                    gl_split_ls[8] = "3"
            elif "0/0" == family_genotype_ls[0] and ("1/1" == family_genotype_ls[1] or "1/1" == family_genotype_ls[2]) :
                if int(float(gl_split_ls[7])) == 0 :
                    gl_split_ls[8] = "4"
            else :
                #gl_split_ls[8] = "5"
                pass
            candidate_single_SV_gt_fam_ls[0][i][gl_index] = ",".join(gl_split_ls)

# 判断一家三口是否符mendel遗传定律，符合mendel遗传定律返回True
# 输入：family_gt_ls：一家三口的基因型概率，如[[1.0, 1.0], [1.0, 0.0], [1.0, 0.0]]
# 如果没有phasing过程复现mendel遗传过程，则只能使用就近匹配确定基因来源 ； 如果已经通过phasing过程复现mendel遗传过程，则基因来源已经确定
def mendel_assert(family_gt_ls, family_mode, nearby_matching) :
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    if family_mode == "M1" :
        if nearby_matching :
            min_f_m = (min([abs(i-family_gt_ls[0][0]) for i in family_gt_ls[1]]),min([abs(i-family_gt_ls[0][1]) for i in family_gt_ls[2]]))
            min_m_f = (min([abs(i-family_gt_ls[0][1]) for i in family_gt_ls[1]]),min([abs(i-family_gt_ls[0][0]) for i in family_gt_ls[2]]))
            (min_gen,min_index) = (min_f_m,1) if sum(min_f_m) < sum(min_m_f) else (min_m_f,2)

            # 修改三口人的基因顺序，方便处理
            if min_index == 2 : # 区分基因来源，保证子女的基因前一个来自父亲，后一个来自母亲
                family_gt_ls[0] = [family_gt_ls[0][1],family_gt_ls[0][0]]
            if abs(family_gt_ls[0][0]-family_gt_ls[1][0]) > abs(family_gt_ls[0][0]-family_gt_ls[1][1]) :
                family_gt_ls[1] = [family_gt_ls[1][1],family_gt_ls[1][0]]
            if abs(family_gt_ls[0][1]-family_gt_ls[2][0]) > abs(family_gt_ls[0][1]-family_gt_ls[2][1]) :
                family_gt_ls[2] = [family_gt_ls[2][1],family_gt_ls[2][0]]

        return abs(family_gt_ls[0][0]-family_gt_ls[1][0]) <= gap_thres and abs(family_gt_ls[0][1]-family_gt_ls[2][0]) <= gap_thres 
    else :
        if nearby_matching :
            min_f_1 = min([abs(i-family_gt_ls[0][0]) for i in family_gt_ls[1]])
            min_f_2 = min([abs(i-family_gt_ls[0][1]) for i in family_gt_ls[1]])
            (min_gen,min_index) = (min_f_1,1) if min_f_1 < min_f_2 else (min_f_2,2)

            # 修改三口人的基因顺序，方便处理
            if min_index == 2 :
                family_gt_ls[0] = [family_gt_ls[0][1],family_gt_ls[0][0]]
            if abs(family_gt_ls[0][0]-family_gt_ls[1][0]) > abs(family_gt_ls[0][0]-family_gt_ls[1][1]) :
                family_gt_ls[1] = [family_gt_ls[1][1],family_gt_ls[1][0]]
        
        return abs(family_gt_ls[0][0]-family_gt_ls[1][0]) <= gap_thres

def resolution_denovo(family_gt_ls, family_gl_offspring, family_mode) :
    return True
    mutation_rate_ls = [0.25*(1-mutation_rate)*(1-mutation_rate),0.25*mutation_rate*(1-mutation_rate),0.25*mutation_rate*(1-mutation_rate),0.25*mutation_rate*mutation_rate]
    [p_data,c0,c1] = family_gl_offspring 
    mutation_posterior_ls = []
    for j in range(4) :
        if j == 0 :
            (f,m) = (family_gt_ls[1][0],family_gt_ls[2][0])
        elif j == 1 :
            (f,m) = (family_gt_ls[1][0],family_gt_ls[2][1])
        elif j == 2 :
            (f,m) = (family_gt_ls[1][1],family_gt_ls[2][0])
        else :
            (f,m) = (family_gt_ls[1][1],family_gt_ls[2][1])
        for k in range(4) :
            if k == 0 :
                f_value = f
                m_value = m
            elif k == 1 :
                f_value = 1-f
                m_value = m
            elif k == 2 :
                f_value = f
                m_value = 1-m
            else :
                f_value = 1-f
                m_value = 1-m
            p_data_f_m = pow(0.5*f_value*err+0.5*(1-f_value)*(1-err)+0.5*m_value*err+0.5*(1-m_value)*(1-err),int(c1))*pow(0.5*f_value*(1-err)+0.5*(1-f_value)*err+0.5*m_value*(1-err)+0.5*(1-m_value)*err,int(c0))
            #logging.info(candidate_single_SV_ls[0][i][gl_index])
            p_f_m_data = p_data_f_m * mutation_rate_ls[k]
            mutation_posterior_ls.append(p_f_m_data)
    not_denovo_prob = mutation_posterior_ls[0] + mutation_posterior_ls[4] + mutation_posterior_ls[8] + mutation_posterior_ls[12]
    #denovo_prob = sum(mutation_posterior_ls) - not_denovo_prob
    #return denovo_prob > not_denovo_prob
    if sum(mutation_posterior_ls) == 0 :
        return False
    not_denovo_prob = min(1,not_denovo_prob / sum(mutation_posterior_ls))
    denovo_prob = 1 - not_denovo_prob
    return not_denovo_prob - denovo_prob < 0.3