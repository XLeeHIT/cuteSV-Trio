from cuteSVTrio.cuteSVTrio_genotype import cal_GL_3, threshold_ref_count, count_coverage, unsolvable_correction, modify_genotype
from cuteSVTrio.cuteSVTrio_mendel import resolution_mendel
import numpy as np
import logging
import pickle
import copy
import pysam
import time

'''
            *********Description*********
            *    TYPE A:        N[chr:pos[    *
            *    TYPE B:        N]chr:pos]    *
            *    TYPE C:        [chr:pos[N    *
            *    TYPE D:        ]chr:pos]N    *
            *****************************
            '''

def resolution_TRA(path, chr_1, read_count, overlap_size, max_cluster_bias, minimum_support_reads_list, bam_path_list, gt_round, read_pos_interval, family_mode, performing_phasing, close_NSS, close_ESS):
    start_time = time.time()
    semi_tra_cluster = list()
    semi_tra_cluster.append([0,0,'','N'])
    #candidate_single_SV = list()
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    candidate_single_SV = []

    seqs = []
    for family_member in family_member_ls :
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"sigindex"), 'rb') as f:
            sigs_index=pickle.load(f)
            #if family_member == "1" :
            #    logging.info(sigs_index["INS"])
            f.close()
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"TRA"), 'rb') as f:
            if chr_1 not in sigs_index["TRA"].keys() :
                f.close()
                continue
            f.seek(sigs_index["TRA"][chr_1])
            f_seqs = pickle.load(f)
            seqs += f_seqs
            #if family_member == "1" and chr == "GL000243.1":
            #    logging.info(f_seqs)
            f.close()
    #logging.info(len(seqs))
    if len(seqs) == 0 :
        return
    seqs = sorted(seqs, key=lambda x: x[0])
    chr_2=seqs[0][2]

    for seq in seqs:
        if seq[2]!=chr_2:
            #logging.info("Finished %s-%s:%s:%f."%(chr_1, chr_2, "TRA/BND", time.time()-start_time))
            if len(semi_tra_cluster) >= read_count:
                if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_tra_cluster(semi_tra_cluster, 
                                              chr_1, 
                                              chr_2, 
                                              read_count, 
                                              overlap_size, 
                                              max_cluster_bias, 
                                              candidate_single_SV,
                                              gt_round)
            semi_tra_cluster = list()
            semi_tra_cluster.append([0,0,'','N'])
            chr_2=seq[2]
        pos_1 = int(seq[1])
        pos_2 = int(seq[3])
        read_id = seq[4]
        BND_type = seq[0]
    
        if pos_1 - semi_tra_cluster[-1][0] > max_cluster_bias or BND_type != semi_tra_cluster[-1][3]:
            if len(semi_tra_cluster) >= read_count:
                if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_tra_cluster(semi_tra_cluster, 
                                              chr_1, 
                                              chr_2, 
                                              read_count, 
                                              overlap_size, 
                                              max_cluster_bias, 
                                              candidate_single_SV,
                                              gt_round)
            semi_tra_cluster = []
            semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type])
        else:
            if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                semi_tra_cluster = []
                semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type])
            else:
                semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type])

    if len(semi_tra_cluster) >= read_count:
        if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
            pass
        else:
            generate_semi_tra_cluster(semi_tra_cluster, 
                                      chr_1, 
                                      chr_2, 
                                      read_count, 
                                      overlap_size, 
                                      max_cluster_bias, 
                                      candidate_single_SV,
                                      gt_round)
    #logging.info(candidate_single_SV)
    candidate_single_SV_fam_ls = [[] for x in family_member_ls]
    for i in range(len(candidate_single_SV)) :
        candidate_SV_fam = [copy.deepcopy(candidate_single_SV[i]) for x in family_member_ls]
        for j in candidate_SV_fam :
            j[12] = []
        #logging.info(candidate_single_SV[i])
        for j in range(len(candidate_single_SV[i][12])) :
            #logging.info(candidate_single_SV[i][12][j])
            #logging.info(candidate_single_SV[i][12][j][3])
            candidate_SV_fam[int(candidate_single_SV[i][12][j][3])-1][12].append(candidate_single_SV[i][12][j])
        retain_flag = False
        #logging.info(candidate_SV_fam)
        for j in range(len(candidate_SV_fam)) :
            if len(candidate_SV_fam[j][12]) >= minimum_support_reads_list[j] :
                retain_flag = True
        if retain_flag :
            for j in range(len(candidate_SV_fam)) :
                if len(candidate_SV_fam[j][12]) >= minimum_support_reads_list[j] :
                    candidate_SV_fam[j][12] = ",".join(candidate_SV_fam[j][12])
                    candidate_single_SV_fam_ls[j].append(candidate_SV_fam[j])
                else :
                    candidate_SV_fam[j][12] = ""
                    candidate_single_SV_fam_ls[j].append(candidate_SV_fam[j])
    # 形成gt
    for i in range(len(candidate_single_SV_fam_ls)) :
        for j in range(len(candidate_single_SV_fam_ls[i])) :
            #logging.info("%s/%d/%d"%(chr_1,i,j))
            if candidate_single_SV_fam_ls[i][j][12] == "" :
                DV, DR, GT, GL, GQ, QUAL = "0","0","0/0","0,100,100,0,0,0,0,0,0","0","0"
            else :
                DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path_list[i], 
                                                   int(candidate_single_SV_fam_ls[i][j][2]), 
                                                   int(candidate_single_SV_fam_ls[i][j][4]), 
                                                   candidate_single_SV_fam_ls[i][j][0], 
                                                   candidate_single_SV_fam_ls[i][j][3], 
                                                   set(candidate_single_SV_fam_ls[i][j][12].split(",")), 
                                                   max_cluster_bias, 
                                                   gt_round, 
                                                   close_NSS)
            candidate_single_SV_fam_ls[i][j][7:12] = [str(DR),str(GT),str(GL),str(GQ),str(QUAL)]
    # Todo:如果所有个体的sv质量控制都是q5，那么删除该记录，否则会有大量的多余BND 
    #logging.info(candidate_single_SV_fam_ls)
    resolution_mendel(candidate_single_SV_fam_ls, family_mode, True, minimum_support_reads_list)
    #if not performing_phasing :
    #    resolution_mendel(candidate_single_SV_fam_ls, family_mode, performing_phasing)



    logging.info("Finished calling %s:%s:%f."%(chr_1, "TRA/BND", time.time()-start_time))
    return (chr_1,candidate_single_SV_fam_ls)

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, candidate_single_SV, gt_round):
    BND_type = semi_tra_cluster[0][3]
    semi_tra_cluster = sorted(semi_tra_cluster, key = lambda x:x[1])
    read_tag = dict()
    temp = list()
    # p1, p2, count
    last_len = semi_tra_cluster[0][1]
    temp.append([semi_tra_cluster[0][0], semi_tra_cluster[0][1], [semi_tra_cluster[0][2]]])
    read_tag[semi_tra_cluster[0][2]] = 0
    for element in semi_tra_cluster:
        if element[1] - last_len > max_cluster_bias:
            temp.append([element[0],element[1],[element[2]]])
            last_len = element[1]
        else:
            temp[-1][0] += element[0]
            temp[-1][1] += element[1]
            temp[-1][2].append(element[2])
            last_len = element[1]

        if element[2] not in read_tag:
            read_tag[element[2]] = 0
    if len(read_tag) < read_count:
        return

    temp = sorted(temp, key = lambda x:-len(set(x[2])))

    if len(temp) > 1 and len(set(temp[1][2])) >= 0.5*read_count:
        if len(set(temp[0][2]))+len(set(temp[1][2])) >= len(semi_tra_cluster)*overlap_size:
            # candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
            # candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
            BND_pos_1 = "%s:%s"%(chr_2, int(temp[0][1]/len(temp[0][2])))
            BND_pos_2 = "%s:%s"%(chr_2, int(temp[1][1]/len(temp[1][2])))
            if BND_type == 'A':
                TRA_1 = "N[%s["%(BND_pos_1)
                TRA_2 = "N[%s["%(BND_pos_2)
            elif BND_type == 'B':
                TRA_1 = "N]%s]"%(BND_pos_1)
                TRA_2 = "N]%s]"%(BND_pos_2)
            elif BND_type == 'C':
                TRA_1 = "[%s[N"%(BND_pos_1)
                TRA_2 = "[%s[N"%(BND_pos_2)
            elif BND_type == 'D':
                TRA_1 = "]%s]N"%(BND_pos_1)
                TRA_2 = "]%s]N"%(BND_pos_2)
            else:
                return

            #if action:
            #    import time
            #    # time_start = time.time()
            #    DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[0][0]/len(temp[0][2])), 
            #                                int(temp[0][1]/len(temp[0][2])), chr_1, chr_2, set(temp[0][2]), 
            #                                max_cluster_bias, gt_round)
            #    # cost_time = time.time() - time_start
            #    # print("BND", chr_1, chr_2, int(temp[0][0]/len(temp[0][2])), int(temp[0][1]/len(temp[0][2])), DR, DV, QUAL, "%.4f"%cost_time)
            #else:
            #    DR = '.'
            #    GT = './.'
            #    GL = '.,.,.'
            #    GQ = "."
            #    QUAL = "."
            #0-chr1|1-TRA_type|2-ref pos|3-chr2|4-alt pos|5-support read num|6-none|7-DR|8-GT|9-GL|10-GQ|11-QUAL|12-support read name|13-none|14-none(todo)|15-none|16-none|17-none
            candidate_single_SV.append([chr_1, 
                                        TRA_1, 
                                        str(int(temp[0][0]/len(temp[0][2]))), 
                                        chr_2, 
                                        str(int(temp[0][1]/len(temp[0][2]))), 
                                        str(len(set(temp[0][2]))),
                                        '',
                                        '.',
                                        './.',
                                        '.,.,.,.,.,.,.,.',
                                        ".",
                                        ".",
                                        list(set(temp[0][2])),
                                        '',
                                        '',
                                        '',
                                        ''])

            #if action:
            #    import time
            #    # time_start = time.time()
            #    DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[1][0]/len(temp[1][2])), 
            #                                int(temp[1][1]/len(temp[1][2])), chr_1, chr_2, set(temp[1][2]), 
            #                                max_cluster_bias, gt_round)
            #    # cost_time = time.time() - time_start
            #    # print("BND", chr_1, chr_2, int(temp[1][0]/len(temp[1][2])), int(temp[1][1]/len(temp[1][2])), DR, DV, QUAL, "%.4f"%cost_time)
            #else:
            #    DR = '.'
            #    GT = './.'
            #    GL = '.,.,.'
            #    GQ = "."
            #    QUAL = "."
            candidate_single_SV.append([chr_1, 
                                        TRA_2, 
                                        str(int(temp[1][0]/len(temp[1][2]))), 
                                        chr_2, 
                                        str(int(temp[1][1]/len(temp[1][2]))), 
                                        str(len(set(temp[1][2]))),
                                        '',
                                        '.',
                                        './.',
                                        '.,.,.,.,.,.,.,.',
                                        ".",
                                        ".",
                                        list(set(temp[1][2])),
                                        '',
                                        '',
                                        '',
                                        ''])
    else:
        if len(set(temp[0][2])) >= len(semi_tra_cluster)*overlap_size:
            # print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
            # candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
            BND_pos = "%s:%s"%(chr_2, int(temp[0][1]/len(temp[0][2])))
            if BND_type == 'A':
                TRA = "N[%s["%(BND_pos)
            elif BND_type == 'B':
                TRA = "N]%s]"%(BND_pos)
            elif BND_type == 'C':
                TRA = "[%s[N"%(BND_pos)
            elif BND_type == 'D':
                TRA = "]%s]N"%(BND_pos)
            else:
                return

            #if action:
            #    import time
            #    # time_start = time.time()
            #    DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[0][0]/len(temp[0][2])), 
            #                                int(temp[0][1]/len(temp[0][2])), chr_1, chr_2, set(temp[0][2]), 
            #                                max_cluster_bias, gt_round)
            #    # cost_time = time.time() - time_start
            #    # print("BND", chr_1, chr_2, int(temp[0][0]/len(temp[0][2])), int(temp[0][1]/len(temp[0][2])), DR, DV, QUAL, "%.4f"%cost_time)
            #else:
            #    DR = '.'
            #    GT = './.'
            #    GL = '.,.,.'
            #    GQ = "."
            #    QUAL = "."
            candidate_single_SV.append([chr_1, 
                                        TRA, 
                                        str(int(temp[0][0]/len(temp[0][2]))), 
                                        chr_2, 
                                        str(int(temp[0][1]/len(temp[0][2]))), 
                                        str(len(set(temp[0][2]))),
                                        '',
                                        '.',
                                        './.',
                                        '.,.,.,.,.,.,.,.',
                                        ".",
                                        ".",
                                        list(set(temp[0][2])),
                                        '',
                                        '',
                                        '',
                                        ''])


def run_tra(args):
    return resolution_TRA(*args)

def call_gt(bam_path, pos_1, pos_2, chr_1, chr_2, read_id_list, max_cluster_bias, gt_round, close_NSS):
    import pysam
    bamfile = pysam.AlignmentFile(bam_path)
    querydata = set()
    search_start = max(int(pos_1) - max_cluster_bias, 0)
    search_end = min(int(pos_1) + max_cluster_bias, bamfile.get_reference_length(chr_1))

    up_bound = threshold_ref_count(len(read_id_list))

    status = count_coverage(chr_1, 
                            search_start, 
                            search_end, 
                            bamfile, 
                            querydata, 
                            up_bound, 
                            gt_round)

    if status == -1:
        DR = '.'
        GT = "./."
        GL = '.,.,.,.,.,.,.,.'
        GQ = "."
        QUAL = "."

    elif status == 1:
        DR = 0
        for query in querydata:
            if query not in read_id_list:
                DR += 1
        GT, GL, GQ, QUAL = cal_GL_3(DR, len(read_id_list), "TRA", 0, None, close_NSS)

    else:
        search_start = max(int(pos_2) - max_cluster_bias, 0)
        search_end = min(int(pos_2) + max_cluster_bias, bamfile.get_reference_length(chr_2))
        status_2 = count_coverage(chr_2, 
                                    search_start, 
                                    search_end, 
                                    bamfile, 
                                    querydata, 
                                    up_bound, 
                                    gt_round)
        # status_2 judgement
        DR = 0
        for query in querydata:
            if query not in read_id_list:
                DR += 1
        GT, GL, GQ, QUAL = cal_GL_3(DR, len(read_id_list), "TRA", 0, None, close_NSS)

    bamfile.close()
    return len(read_id_list), DR, GT, GL, GQ, QUAL
