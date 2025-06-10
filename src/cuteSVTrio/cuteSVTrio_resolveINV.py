from cuteSVTrio.cuteSVTrio_genotype import cal_CIPOS, overlap_cover, assign_gt, unsolvable_correction, modify_genotype, increase_sigs_through_pedigree, check_gt_consistent, inconformity_mendel_modify
from cuteSVTrio.cuteSVTrio_mendel import resolution_mendel
import numpy as np
from math import ceil
import logging
import pickle
import copy
import time

def resolution_INV(path, chr, read_count, max_cluster_bias, minimum_support_reads_list, sv_size, MaxSize, gt_round, read_pos_interval, family_mode, performing_phasing):
    '''
    cluster INV
    ************************************************************************
    path:	INV.sigs
    chr:	chromosome id
    svtype:	<INV>
    '''
    # INV检测还有一些非常相近的长度接近的fp可以解决，进一步提高presicion
    # Initialization of some temporary variables
    start_time = time.time()
    semi_inv_cluster = list()
    semi_inv_cluster.append([0,0,'',''])
    candidate_single_SV = list()
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    # Load inputs & cluster breakpoint from each signature read 
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

    seqs = []
    for family_member in family_member_ls :
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"sigindex"), 'rb') as f:
            sigs_index=pickle.load(f)
            #if family_member == "1" :
            #    logging.info(sigs_index["INS"])
            f.close()
        with open("%s%s.%s.%s.pickle"%(path,family_mode,family_member,"INV"), 'rb') as f:
            if chr not in sigs_index["INV"].keys() :
                f.close()
                continue
            f.seek(sigs_index["INV"][chr])
            f_seqs = pickle.load(f)
            seqs += f_seqs
            f.close()
        
    if len(seqs) == 0 :
        return
    #for i in range(len(seqs)) :
    #    seqs[i][0] = "++"
    seqs = sorted(seqs, key=lambda x: (int(x[1]),int(x[2])))
    for seq in seqs:
        #logging.info(seq)
        strand = "++"
        breakpoint_1_in_read = int(seq[1])
        breakpoint_2_in_read = int(seq[2])
        #if chr == "2" and breakpoint_1_in_read > 30473576 and breakpoint_1_in_read < 30474576 :
        #    logging.info("%s/%s/%s"%(str(strand),str(breakpoint_1_in_read),str(breakpoint_2_in_read)))
        read_id = seq[3]

        # print("new")
        # print(seq[1], seq[2], seq[3], seq[4], seq[5])
        # print(semi_inv_cluster)

        if breakpoint_1_in_read - semi_inv_cluster[-1][0] > max_cluster_bias or breakpoint_2_in_read - semi_inv_cluster[-1][1] > max_cluster_bias or strand != semi_inv_cluster[-1][-1]:
            if len(semi_inv_cluster) >= read_count:
                if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_inv_cluster(semi_inv_cluster, 
                                              chr, 
                                              read_count, 
                                              max_cluster_bias,
                                              sv_size, 
                                              candidate_single_SV, 
                                              MaxSize,
                                              gt_round,
                                              family_mode,
                                              performing_phasing)
            semi_inv_cluster = []
            semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])
        else:
            if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                semi_inv_cluster = []
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])
            else:
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])

    if len(semi_inv_cluster) >= read_count:
        if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
            pass
        else:
            generate_semi_inv_cluster(semi_inv_cluster, 
                                      chr, 
                                      read_count, 
                                      max_cluster_bias,
                                      sv_size, 
                                      candidate_single_SV, 
                                      MaxSize,
                                      gt_round,
                                      family_mode,
                                      performing_phasing)
    #logging.info(candidate_single_SV)
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
        for j in range(len(candidate_single_SV[i][6])) :
            fam_reads[int(candidate_single_SV[i][6][j][3])-1].append(candidate_single_SV[i][6][j])
            if performing_phasing :
                reads_pos[int(candidate_single_SV[i][6][j][3])-1].append(candidate_single_SV[i][8][j])
        #for j in range(len(family_member_ls)) :
        #    if len(fam_reads[j]) != len(reads_pos[j]) :
        #        logging.info("%s/%s"%(str(fam_reads[j]),str(reads_pos[j])))
        candidate_single_SV_fam_ls[0][i][6] = fam_reads[0]
        candidate_single_SV_fam_ls[0][i][8] = reads_pos[0]
        for j in range(1,len(family_member_ls)) :
            fam_candidate_other[j][6] = fam_reads[j]
            fam_candidate_other[j][8] = reads_pos[j]
            candidate_single_SV_fam_ls[j].append(fam_candidate_other[j])

    #if chr == "2" :
    #    for i in range(len(candidate_single_SV_fam_ls[0])) :
    #        if candidate_single_SV_fam_ls[0][i][2] > 30473576 and candidate_single_SV_fam_ls[0][i][2] < 30474576 :
    #            logging.info(candidate_single_SV_fam_ls[0][i])
    #            logging.info(candidate_single_SV_fam_ls[1][i])
    #            logging.info(candidate_single_SV_fam_ls[2][i])
    candidate_single_SV_gt_fam_ls = []
    for i in range(len(family_member_ls)) :
        family_member = family_member_ls[i]
        #gt_candidate_sv = call_gt(path, chr, candidate_single_SV_fam_ls[i], 1000, 'INS', family_mode, family_member)
        candidate_single_SV_gt_fam_ls.append(call_gt(path, chr, candidate_single_SV_fam_ls[i], candidate_single_SV_fam_ls[0], 1000, family_mode, family_member, minimum_support_reads_list, performing_phasing))
    #logging.info("INV/%s/%d/%s"%(chr,len(candidate_single_SV_gt_fam_ls[0]),str(candidate_single_SV_gt_fam_ls[0][0:10])))
    standard_list = 0
    #logging.info(candidate_single_SV_gt_fam_ls[0][0:10])
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) != 0 :
            standard_list = i
            break
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        if len(candidate_single_SV_gt_fam_ls[i]) == 0 :
            candidate_single_SV_gt_fam_ls[i] = copy.deepcopy(candidate_single_SV_gt_fam_ls[standard_list])
            for s_v in candidate_single_SV_gt_fam_ls[i] :
                #logging.info(s_v)
                s_v[8] = "0/0"
                s_v[9] = "0,100,100,0,0,0,0,0,0"
                s_v[10] = 996
                s_v[11] = 0
                s_v[12] = ""
                s_v[15] = ""
                s_v[16] = ""
                s_v[17] = ""
    #if chr == "2" :
    #    for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
    #        if int(candidate_single_SV_gt_fam_ls[0][i][2]) > 30473576 and int(candidate_single_SV_gt_fam_ls[0][i][2]) < 30474576 :
    #            logging.info(candidate_single_SV_gt_fam_ls[0][i])
    #            logging.info(candidate_single_SV_gt_fam_ls[1][i])
    #            logging.info(candidate_single_SV_gt_fam_ls[2][i])
    increase_sigs_through_pedigree(candidate_single_SV_gt_fam_ls, 'INV', minimum_support_reads_list, family_mode)
    #if chr == "2" :
    #    for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
    #        if int(candidate_single_SV_gt_fam_ls[0][i][2]) > 30473576 and int(candidate_single_SV_gt_fam_ls[0][i][2]) < 30474576 :
    #            logging.info(candidate_single_SV_gt_fam_ls[0][i])
    #            logging.info(candidate_single_SV_gt_fam_ls[1][i])
    #            logging.info(candidate_single_SV_gt_fam_ls[2][i])
    unsolvable_correction(candidate_single_SV_gt_fam_ls, 'INV', family_mode)
    if not performing_phasing and family_mode == "M1" :
        resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, True, minimum_support_reads_list)
    #if not performing_phasing :
    #    resolution_mendel(candidate_single_SV_gt_fam_ls, family_mode, performing_phasing)
    #unsolvable_correction(candidate_single_SV_gt_fam_ls, 'INV', family_mode, family_member)
    #not_accord_mendel_num,run_mcmc_num,run_mcmc_ls=resolution_mendel(candidate_single_SV_gt_fam_ls,'INV',family_mode, gap_thres)
    #resolution_denovo(candidate_single_SV_gt_fam_ls, 'INV', family_mode)
    #modify_genotype(candidate_single_SV_gt_fam_ls, 'INV')
    
    logging.info("Finished calling %s:%s:%f."%(chr, "INV", time.time()-start_time))
    return (chr,candidate_single_SV_gt_fam_ls)

def generate_semi_inv_cluster(semi_inv_cluster, chr, read_count, max_cluster_bias, 
    sv_size, candidate_single_SV, MaxSize, gt_round, family_mode, performing_phasing):
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    strand = semi_inv_cluster[0][-1]

    #logging.info(semi_inv_cluster)
    read_id = [i[2] for i in semi_inv_cluster]
    pos_list = [i[0] for i in semi_inv_cluster]
    support_read = len(list(set(read_id)))
    if support_read < read_count:
        return

    inv_cluster_b2 = sorted(semi_inv_cluster, key = lambda x:x[1])

    # breakpoint_1 = np.mean(breakpoint_1_candidate)
    last_bp = inv_cluster_b2[0][1]
    temp_count = 1
    # max_count = 0
    temp_sum_b1 = inv_cluster_b2[0][0]
    temp_sum_b2 = last_bp

    # max_sum = 0
    temp_id = dict()
    #logging.info(inv_cluster_b2[0])
    temp_id[inv_cluster_b2[0][2]] = [0,inv_cluster_b2[0][0],inv_cluster_b2[0][1],inv_cluster_b2[0][2]]
    #logging.info(inv_cluster_b2[0])

    for i in inv_cluster_b2[1:]:
        if i[1] - last_bp > max_cluster_bias:
            if temp_count >= read_count:
                max_count_id = len(temp_id)

                breakpoint_1 = round(temp_sum_b1 / temp_count)
                breakpoint_2 = round(temp_sum_b2 / temp_count)
                inv_len = breakpoint_2 - breakpoint_1
                if inv_len >= sv_size and max_count_id >= read_count:
                    # candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
                    if inv_len <= MaxSize or MaxSize == -1:
                        #logging.info(temp_id)
                        #logging.info([chr, svtype, breakpoint_1, inv_len, max_count_id,strand,list(temp_id.keys()),breakpoint_2])
                        #tmp = [temp_id[key] for key in temp_id.keys()]
                        #logging.info(tmp)
                        fam_allele_len_ls = [[] for x in family_member_ls]
                        for key in temp_id.keys() :
                            fam_allele_len_ls[int(temp_id[key][3][3])-1].append(abs(int(temp_id[key][2])-int(temp_id[key][1])))
                        for fam_i in range(len(fam_allele_len_ls)) :
                            if len(fam_allele_len_ls[fam_i]) == 0 :
                                fam_allele_len_ls[fam_i].append(0)
                        sv_len_str = ""
                        for fam_i in range(len(fam_allele_len_ls)) :
                            sv_len_str = sv_len_str + str(int(np.mean(fam_allele_len_ls[fam_i]))) + ","
                        candidate_single_SV.append([chr, 
                                                    "INV", 
                                                    breakpoint_1, 
                                                    inv_len, 
                                                    max_count_id,
                                                    strand,
                                                    list(temp_id.keys()),
                                                    breakpoint_2,
                                                    [temp_id[key][1] for key in temp_id.keys()],
                                                    sv_len_str[:-1]])
                        # print(chr, svtype, str(int(breakpoint_1)), str(int(inv_len)), str(max_count_id), str(DR), str(GT), strand)

            temp_id = dict()
            temp_count = 1
            temp_sum_b1 = i[0]
            temp_sum_b2 = i[1]
            temp_id[i[2]] = [0,i[0],i[1],i[2]]
        else:
            if i[2] not in temp_id:
                temp_id[i[2]] = [0,i[0],i[1],i[2]]
            else:
                temp_id[i[2]][0] += 1
            temp_count += 1
            temp_sum_b1 += i[0]
            temp_sum_b2 += i[1]
        last_bp = i[1]
    if temp_count >= read_count:
        max_count_id = len(temp_id)
        breakpoint_1 = round(temp_sum_b1 / temp_count)
        breakpoint_2 = round(temp_sum_b2 / temp_count)
        inv_len = breakpoint_2 - breakpoint_1
        if inv_len >= sv_size and max_count_id >= read_count:
            # candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
            if inv_len <= MaxSize or MaxSize == -1:
                #tmp = [temp_id[key] for key in temp_id.keys()]
                #logging.info(tmp)
                fam_allele_len_ls = [[] for x in family_member_ls]
                for key in temp_id.keys() :
                    fam_allele_len_ls[int(temp_id[key][3][3])-1].append(abs(int(temp_id[key][2])-int(temp_id[key][1])))
                for fam_i in range(len(fam_allele_len_ls)) :
                    if len(fam_allele_len_ls[fam_i]) == 0 :
                        fam_allele_len_ls[fam_i].append(0)
                sv_len_str = ""
                for fam_i in range(len(fam_allele_len_ls)) :
                    sv_len_str = sv_len_str + str(int(np.mean(fam_allele_len_ls[fam_i]))) + ","
                candidate_single_SV.append([chr, 
                                            "INV", 
                                            breakpoint_1, 
                                            inv_len, 
                                            max_count_id,
                                            strand,
                                            list(temp_id.keys()),
                                            breakpoint_2,
                                            [temp_id[key][1] for key in temp_id.keys()],
                                            sv_len_str[:-1]])
                # print(chr, svtype, str(int(breakpoint_1)), str(int(inv_len)), str(max_count_id), str(DR), str(GT), strand)

def run_inv(args):
    return resolution_INV(*args)

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
    #for item in candidate_info_SV:
    #    svs_list.append((max(item[2] - max_cluster_bias, 0), item[7] + max_cluster_bias))
    for item in candidate_info_SV:
        svs_list.append((max(item[2] - max_cluster_bias, 0), item[2] + max_cluster_bias))
    for item in candidate_info_SV:
        svs_list.append((max(item[7] - max_cluster_bias, 0), item[7] + max_cluster_bias))
    iteration_dict, primary_num_dict, cover_dict, overlap_dict, cover_pos_dict = overlap_cover(svs_list, reads_list, performing_phasing) # both key(sv idx), value(set(read id))
    assert len(cover_dict) == 2 * len(candidate_info_SV), "overlap length error"
    candidate_single_SV_length = len(candidate_info_SV)
    for idx in range(candidate_single_SV_length):
        for i in range(len(cover_dict[idx + candidate_single_SV_length])) :
            if cover_dict[idx + candidate_single_SV_length][i] not in cover_dict[idx] :
                cover_dict[idx].append(cover_dict[idx + candidate_single_SV_length][i])
                if performing_phasing :
                    cover_pos_dict[idx].append(cover_pos_dict[idx + candidate_single_SV_length][i])
    for idx in range(candidate_single_SV_length, candidate_single_SV_length * 2, 1):
        cover_dict.pop(idx)
    assert len(cover_dict) == len(candidate_info_SV), "overlap length error"

    read_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][6]
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, cover_pos_dict, "INV", family_member, minimum_support_reads_list[int(family_member)-1], performing_phasing)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    #0-chr|1-sv_type|2-pos|3-len|4-support read num|5-strand|6-''|7-DR|8-GT|9-GL|10-GQ|11-QUAL|12-support read name|13-[]|14-sv len|15-support read pos|16-not support read name|17-not support read pos
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_info_SV[i][0], 
                                    candidate_info_SV[i][1], 
                                    str(int(candidate_info_SV[i][2])), 
                                    str(int(candidate_info_SV[i][3])), 
                                    str(candidate_single_SV[i][4]), 
                                    candidate_single_SV[i][5],
                                    '',
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][6])])
        candidate_single_SV_gt[i].append([])
        candidate_single_SV_gt[i].append(candidate_single_SV[i][9])
        if performing_phasing :
            candidate_single_SV_gt[i].append(','.join([str(x) for x in candidate_single_SV[i][8]]))
            candidate_single_SV_gt[i].append(','.join(assign_list[i][6]))
            candidate_single_SV_gt[i].append(','.join([str(x) for x in assign_list[i][7]]))
        else :
            candidate_single_SV_gt[i].append('')
            candidate_single_SV_gt[i].append('')
            candidate_single_SV_gt[i].append('')

    return candidate_single_SV_gt