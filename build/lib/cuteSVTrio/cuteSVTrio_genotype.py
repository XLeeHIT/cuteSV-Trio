
from cuteSVTrio.cuteSVTrio_Description import Generation_VCF_header
from math import log10,ceil
import numpy as np
import pysam
import pickle
import logging
import random

err = 0.1
prior = float(1/3)
Genotype = ["0/0", "0/1", "1/1"]
deta_limit = -0.001
chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]


def log10sumexp(log10_probs):
    # Normalization of Genotype likelihoods
    m = max(log10_probs)
    return m + log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    # Adjust the Genotype likelihoods
    log10_probs = np.array(log10_probs)
    lse = log10sumexp(log10_probs)
    return np.minimum(log10_probs - lse, 0.0)

def rescale_read_counts(c0, c1, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    Total = c0 + c1
    if Total > max_allowed_reads:
        c0 = int(max_allowed_reads * float(c0/Total))
        c1 = max_allowed_reads - c0
    return c0, c1

def cal_GL_3(c0, c1, svtype, minimum_support_reads, sv_len):
    t0, t1 = c0, c1
    if c0 + c1 != 0 and c1 < minimum_support_reads and sv_len is None and svtype != "TRA":
        return '0/0',"0,100,100,0,0,0,"+str(t0)+","+str(t1)+",0",996,0
    if svtype not in ["DEL","INS"] :
        if c0 == 0 and c1 == 0 :
            return '0/0',"0,100,100,0,0,0,0,0,-1",996,0
    else :
        if c0 == 0 and c1 == 0 : 
            if sv_len is not None and int(sv_len) >= 1500:
                return '1/1',"100,100,0,1,1,0,0,0,-1",996,255    # NSS
            else :
                return '0/0',"0,100,100,0,0,0,0,0,-1",996,0    # no NSS
            
    c0, c1 = rescale_read_counts(c0, c1) # DR, DV
    
    ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
    ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
    ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
    p_data = ori_GL00+ori_GL01+ori_GL11
    p_data = int(np.around(-10*log10(p_data)))
    a = 1
    b = GL_P[0] - GL_P[2] - 1
    c = GL_P[2]
    roots = np.roots([a, b, c])
    return Genotype[prob.index(max(prob))], "%d,%d,%d,%f,%f,%d,%f,%f,0"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, p_data, t0, t1), max(GQ), QUAL

def cal_Gl_3_sim(c0, c1) :
    c0, c1 = rescale_read_counts(c0, c1) # DR, DV
    
    ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
    ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
    ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
    p_data = ori_GL00+ori_GL01+ori_GL11
    p_data = int(np.around(-10*log10(p_data)))
    a = 1
    b = GL_P[0] - GL_P[2] - 1
    c = GL_P[2]
    roots = np.roots([a, b, c])
    return "%f,%f,%f"%(log10(GL_P[0]), log10(GL_P[1]), log10(GL_P[2]))

#cal_GL_2得到单支基因型的编译可能性
def cal_GL_2(c0, c1) :
    c0, c1 = rescale_read_counts(c0, c1) # DR, DV
    ori_GL0 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior))
    ori_GL1 = np.float64(pow(err, c0)*pow((1-err), c1)*prior)
    prob = list(normalize_log10_probs([log10(ori_GL0), log10(ori_GL1)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
    return GL_P[0], GL_P[1], PL[0], PL[1], max(GQ), QUAL

def cal_CIPOS(std, num):
    pos = int(1.96 * std / num ** 0.5)
    return "-%d,%d"%(pos,pos)

def threshold_ref_count(num):
    if num <= 2:
        return 20*num
    elif 3 <= num <= 5:
        return 9*num 
    elif 6 <= num <= 15:
        return 7*num
    else:
        return 5*num

def count_coverage(chr, s, e, f, read_count, up_bound, itround):
    status = 0
    iteration = 0
    primary_num = 0
    for i in f.fetch(chr, s, e):
        iteration += 1
        if i.flag not in [0,16]:
            continue
        primary_num += 1
        if i.reference_start < s and i.reference_end > e:
            read_count.add(i.query_name)
            if len(read_count) >= up_bound:
                status = 1
                break
        if iteration >= itround:
            if float(primary_num/iteration) <= 0.2:
                status = 1
            else:
                status = -1
            break

    return status

def overlap_cover(svs_list, reads_list, performing_phasing):
    sort_list = list()
    idx = 0
    for i in reads_list:
        sort_list.append([i[0], 1, idx, i[2], i[3]])
        sort_list.append([i[1], 2, idx, i[2], i[3]])
        idx += 1
    idx = 0
    for i in svs_list:
        sort_list.append([i[0], 3, idx])
        sort_list.append([i[1], 0, idx])
        idx += 1
    sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    read_set = set()
    overlap_dict = dict()
    cover_dict = dict()
    for node in sort_list:
        if node[1] == 1: # set2(read) left
            read_set.add(node[2])
            for x in svs_set:
                if svs_list[x][1] == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set2(read) right
            read_set.remove(node[2])
        elif node[1] == 3: # set1(sv) left
            svs_set.add(node[2])
            overlap_dict[node[2]] = set()
            for x in read_set:
                overlap_dict[node[2]].add(x)
            cover_dict[node[2]] = set()
            for x in read_set:
                cover_dict[node[2]].add(x)
        elif node[1] == 0: # set1(sv) right
            svs_set.remove(node[2])
            temp_set = set()
            for x in read_set:
                temp_set.add(x)
            cover_dict[node[2]] = cover_dict[node[2]] & temp_set
    overlap2_dict = dict()
    cover2_dict = dict()
    cover_pos_dict = dict()
    iteration_dict = dict()
    primary_num_dict = dict()
    for idx in cover_dict:
        iteration_dict[idx] = len(overlap_dict[idx])
        primary_num_dict[idx] = 0
        for x in overlap_dict[idx]:
            if reads_list[x][2] == 1:
                primary_num_dict[idx] += 1
        cover2_dict[idx] = list()
        if performing_phasing :
            cover_pos_dict[idx] = list()
        for x in cover_dict[idx]:
            if reads_list[x][2] == 1:
                if reads_list[x][3] not in cover2_dict[idx] :
                    cover2_dict[idx].append(reads_list[x][3])
                    if performing_phasing :
                        cover_pos_dict[idx].append(reads_list[x][0])
        overlap2_dict[idx] = set()
        for x in overlap_dict[idx]:
            if reads_list[x][2] == 1:
                overlap2_dict[idx].add(reads_list[x][3])
    return iteration_dict, primary_num_dict, cover2_dict, overlap2_dict, cover_pos_dict

def assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, cover_pos_dict, svtype, family_member, minimum_support_reads, performing_phasing):
    assign_list = list()
    for idx in read_id_dict:
        cover_not_sv_read = []
        cover_not_sv_pos = []
        iteration = iteration_dict[idx]
        primary_num = primary_num_dict[idx]
        read_count = cover_dict[idx]
        if performing_phasing :
            pos_count = cover_pos_dict[idx]
        DR = 0
        for j in range(len(read_count)):
            query = read_count[j]
            if query not in read_id_dict[idx]:
                DR += 1
                if query[3] == family_member :
                    if performing_phasing :
                        cover_not_sv_read.append(query)
                        cover_not_sv_pos.append(pos_count[j])
        GT, GL, GQ, QUAL = cal_GL_3(DR, len(read_id_dict[idx]), svtype, minimum_support_reads, None)
        if performing_phasing :
            if len(cover_not_sv_read) != len(cover_not_sv_pos) :
                logging.info("%d/%d"%(len(cover_not_sv_read),len(cover_not_sv_pos)))
            assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL, cover_not_sv_read, cover_not_sv_pos])
        else :
            assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL, [], []])
    return assign_list

def duipai(svs_list, reads_list, iteration_dict, primary_num_dict, cover2_dict, overlap2_dict):
    print('start duipai')
    idx = 0
    correct_cover = 0
    correct_overlap = 0
    bb = set()
    for i in svs_list:
        cover = set()
        overlap = set()
        primary_num = 0
        iteration = 0
        for j in reads_list:
            if (j[0] <= i[0] and j[1] > i[0]) or (i[0] <= j[0] < i[1]):
                iteration += 1
                if j[2] == 1:
                    overlap.add(j[3])
                    primary_num += 1
                    if i[0] >= j[0] and i[1] <= j[1]:
                        cover.add(j[3])
        flag = 0
        if iteration != iteration_dict[idx]:
            print('Iteration error %d:%d(now) %d(ans)'%(idx, iteration_dict[idx], iteration))
        if primary_num != primary_num_dict[idx]:
            print('Primary_num error %d:%d(now) %d(ans)'%(idx, primary_num_dict[idx], primary_num))
        if len(cover) == len(cover2_dict[idx]): flag += 1
        if len(cover - cover2_dict[idx]) == 0: flag += 1
        if len(cover2_dict[idx] - cover) == 0: flag += 1
        if flag != 3:
            print(idx)
            print(cover)
            print(cover2_dict[idx])
            print(cover - cover2_dict[idx])
        else:
            correct_cover += 1
        flag = 0
        if len(overlap) == len(overlap2_dict[idx]): flag += 1
        if len(overlap - overlap2_dict[idx]) == 0: flag += 1
        if len(overlap2_dict[idx] - overlap) == 0: flag += 1
        if flag != 3:
            print(idx)
            print(overlap)
            print(overlap2_dict[idx])
            print(overlap - overlap2_dict[idx])
        else:
            correct_overlap += 1
        idx += 1
    print('Correct iteration cover %d; overlap %d'%(correct_cover, correct_overlap))

def generate_output(args, semi_result_ls, chrom, temporary_dir):
    
    '''
    Generation of VCF format file.
    VCF version: 4.2
    '''

    # genotype_trigger = TriggerGT[args.genotype]
    f=open("%s%s.results/%s.pickle"%(temporary_dir,args.family_mode,chrom), "wb")
    for i in range(len(semi_result_ls)) :
        semi_result_ls[i].sort(key = lambda x:int(x[2]))
    
    fa_file = pysam.FastaFile(args.reference)
    try:
        ref_chrom=fa_file.fetch(chrom)
    except:
        raise Exception("No corresponding contig in reference with %s."%(chrom))
    fa_file.close()
    lines=[]
    BATCH_SIZE=1000
    for i in range(len(semi_result_ls[0])) :
        res_i = semi_result_ls[0][i]
        if res_i[1] in ["DEL", "INS"]:
            if abs(int(float(res_i[3]))) > args.max_size and args.max_size != -1:
                continue
            if abs(int(float(res_i[3]))) < args.min_size:
                continue
            if res_i[1] == "INS":
                cal_end = int(res_i[2]) + 1
            else:
                cal_end = int(res_i[2]) + abs(int(float(res_i[3])))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if res_i[8] == "0/0" else "PRECISE", 
                SVTYPE = res_i[1], 
                SVLEN = res_i[3], 
                END = str(cal_end), 
                CIPOS = res_i[5], 
                CILEN = res_i[6], 
                RE = res_i[4],
                RNAMES = res_i[12] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(int(res_i[4]) / (int(res_i[4]) + int(res_i[7])), 4))
            except:
                info_list += ";AF=."
            if res_i[1] =="DEL":
                info_list += ";STRAND=+-"
            if res_i[9].split(",")[8] in ["-0","-1","-2","-3","-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14"] :
                info_list = info_list + ";CorrectType=" + res_i[9].split(",")[8][1:]
            if res_i[9].split(",")[8] in ["1","2","3"] :
                info_list = info_list + ";Denovo=" + res_i[9].split(",")[8]
            info_list = info_list + ";QUALLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                info_list = info_list + str(res_j[11]) + ","
            info_list = info_list[:-1] + ";FILTERLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                if res_j[11] == "." or res_j[11] == None:
                    filter_lable = "PASS"
                else:
                    filter_lable = "PASS" if float(res_j[11]) >= 5.0 else "q5"
                info_list = info_list + filter_lable + ","
            info_list = info_list[:-1]
            if res_i[11] == "." or res_i[11] == None:
                filter_lable = "PASS"
            else:
                filter_lable = "PASS" if float(res_i[11]) >= 5.0 else "q5"
            if args.performing_phasing :
                GT_FORMAT = "GT:HP_GT:DR:DV:PL:AP:GQ"
            else :
                GT_FORMAT = "GT:DR:DV:PL:AP:GQ"
            vcf_line = "{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}".format(
                CHR = res_i[0], 
                POS = str(int(res_i[2])), 
                ID = "cuteSVTrio.%s.<SVID>"%(res_i[1]),
                REF = ref_chrom[max(int(res_i[2])-1, 0)] if res_i[1] == 'INS' else ref_chrom[max(int(res_i[2])-1, 0):int(res_i[2])-int(res_i[3])],
                ALT = "%s"%(ref_chrom[max(int(res_i[2])-1, 0)]+res_i[13] if res_i[1] == 'INS' else ref_chrom[max(int(res_i[2])-1, 0)]), 
                INFO = info_list, 
                FORMAT = GT_FORMAT, 
                QUAL = res_i[11],
                PASS = filter_lable)
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                PL = "%s,%s,%s"%(res_j[9].split(",")[0],res_j[9].split(",")[1],res_j[9].split(",")[2])
                AP = "%s,%s"%(res_j[9].split(",")[3],res_j[9].split(",")[4])
                vcf_line = vcf_line + "\t{GT}:{DR}:{RE}:{PL}:{AP}:{GQ}".format(
                                        GT = res_j[8],
                                        DR = res_j[7],
                                        RE = res_j[4],
                                        PL = PL,
                                        AP = AP,
                                        GQ = res_j[10])
            lines.append((res_i[1],vcf_line+"\n"))
        elif res_i[1] == "DUP":
            if abs(int(float(res_i[3]))) > args.max_size and args.max_size != -1:
                continue
            cal_end = int(res_i[2]) + 1 + abs(int(float(res_i[3])))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND=-+;RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if res_i[8] == "0/0" else "PRECISE", 
                SVTYPE = res_i[1], 
                SVLEN = res_i[3], 
                END = str(cal_end), 
                RE = res_i[4],
                RNAMES = res_i[12] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(int(res_i[4]) / (int(res_i[4]) + int(res_i[7])), 4))
            except:
                info_list += ";AF=."
            if res_i[9].split(",")[8] in ["-0","-1","-2","-3","-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14"] :
                info_list = info_list + ";CorrectType=" + res_i[9].split(",")[8][1:]
            if res_i[9].split(",")[8] in ["1","2","3"] :
                info_list = info_list + ";Denovo=" + res_i[9].split(",")[8]
            info_list = info_list + ";QUALLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                info_list = info_list + str(res_j[11]) + ","
            info_list = info_list[:-1] + ";FILTERLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                if res_j[11] == "." or res_j[11] == None:
                    filter_lable = "PASS"
                else:
                    filter_lable = "PASS" if float(res_j[11]) >= 5.0 else "q5"
                info_list = info_list + filter_lable + ","
            info_list = info_list[:-1]
            if res_i[11] == "." or res_i[11] == None:
                filter_lable = "PASS"
            else:
                filter_lable = "PASS" if float(res_i[11]) >= 5.0 else "q5"
            if args.performing_phasing :
                GT_FORMAT = "GT:HP_GT:DR:DV:PL:AP:GQ"
            else :
                GT_FORMAT = "GT:DR:DV:PL:AP:GQ"
            vcf_line = "{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}".format(
                CHR = res_i[0], 
                POS = str(int(res_i[2]) + 1), 
                ID = "cuteSVTrio.%s.<SVID>"%(res_i[1]),
                REF = ref_chrom[int(res_i[2])],
                ALT = "<%s>"%(res_i[1]), 
                INFO = info_list, 
                FORMAT = GT_FORMAT, 
                QUAL = res_i[11],
                PASS = filter_lable)
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                PL = "%s,%s,%s"%(res_j[9].split(",")[0],res_j[9].split(",")[1],res_j[9].split(",")[2])
                AP = "%s,%s"%(res_j[9].split(",")[3],res_j[9].split(",")[4])
                vcf_line = vcf_line + "\t{GT}:{DR}:{RE}:{PL}:{AP}:{GQ}".format(
                                        GT = res_j[8],
                                        DR = res_j[7],
                                        RE = res_j[4],
                                        PL = PL,
                                        AP = AP,
                                        GQ = res_j[10])
            lines.append((res_i[1],vcf_line+"\n"))
        elif res_i[1] == "INV":
            if abs(int(float(res_i[3]))) > args.max_size and args.max_size != -1:
                continue
            cal_end = int(res_i[2]) + 1 + abs(int(float(res_i[3])))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND={STRAND};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if res_i[8] == "0/0" else "PRECISE", 
                SVTYPE = res_i[1], 
                SVLEN = res_i[3], 
                END = str(cal_end), 
                RE = res_i[4],
                STRAND = res_i[5],
                RNAMES = res_i[12] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(int(res_i[4]) / (int(res_i[4]) + int(res_i[7])), 4))
            except:
                info_list += ";AF=."
            if res_i[9].split(",")[8] in ["-0","-1","-2","-3","-4","-5","-6","-7","-8","-9","-10","-11","-12","-13","-14"] :
                info_list = info_list + ";CorrectType=" + res_i[9].split(",")[8][1:]
            if res_i[9].split(",")[8] in ["1","2","3"] :
                info_list = info_list + ";Denovo=" + res_i[9].split(",")[8]
            info_list = info_list + ";QUALLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                info_list = info_list + str(res_j[11]) + ","
            info_list = info_list[:-1] + ";FILTERLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                if res_j[11] == "." or res_j[11] == None:
                    filter_lable = "PASS"
                else:
                    filter_lable = "PASS" if float(res_j[11]) >= 5.0 else "q5"
                info_list = info_list + filter_lable + ","
            info_list = info_list[:-1]
            if res_i[11] == "." or res_i[11] == None:
                filter_lable = "PASS"
            else:
                filter_lable = "PASS" if float(res_i[11]) >= 5.0 else "q5"
            if args.performing_phasing :
                GT_FORMAT = "GT:HP_GT:DR:DV:PL:AP:GQ"
            else :
                GT_FORMAT = "GT:DR:DV:PL:AP:GQ"
            vcf_line = "{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}".format(
                CHR = res_i[0], 
                POS = str(int(res_i[2]) + 1), 
                ID = "cuteSVTrio.%s.<SVID>"%(res_i[1]),
                REF = ref_chrom[int(res_i[2])],
                ALT = "<%s>"%(res_i[1]), 
                INFO = info_list, 
                FORMAT = GT_FORMAT, 
                QUAL = res_i[11],
                PASS = filter_lable)
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                PL = "%s,%s,%s"%(res_j[9].split(",")[0],res_j[9].split(",")[1],res_j[9].split(",")[2])
                AP = "%s,%s"%(res_j[9].split(",")[3],res_j[9].split(",")[4])
                vcf_line = vcf_line + "\t{GT}:{DR}:{RE}:{PL}:{AP}:{GQ}".format(
                                        GT = res_j[8],
                                        DR = res_j[7],
                                        RE = res_j[4],
                                        PL = PL,
                                        AP = AP,
                                        GQ = res_j[10])
            lines.append((res_i[1],vcf_line+"\n"))
        else:
            # BND
            info_list = "{PRECISION};SVTYPE={SVTYPE};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if res_i[8] == "0/0" else "PRECISE", 
                SVTYPE = "BND", 
                RE = res_i[5],
                RNAMES = res_i[12] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(int(i[4]) / (int(i[4]) + int(i[7])), 4))
            except:
                info_list += ";AF=."
            if res_i[9].split(",")[8] in ["-1","-2","-3"] :
                info_list = info_list + ";CorrectType=" + res_i[9].split(",")[8][1]
            if res_i[9].split(",")[8] in ["1","2","3"] :
                info_list = info_list + ";Denovo=" + res_i[9].split(",")[8]
            info_list = info_list + ";QUALLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                info_list = info_list + str(res_j[11]) + ","
            info_list = info_list[:-1] + ";FILTERLIST="
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                if res_j[11] == "." or res_j[11] == None:
                    filter_lable = "PASS"
                else:
                    filter_lable = "PASS" if float(res_j[11]) >= 5.0 else "q5"
                info_list = info_list + filter_lable + ","
            info_list = info_list[:-1]
            if res_i[11] == "." or res_i[11] == None:
                filter_lable = "PASS"
            else:
                filter_lable = "PASS" if float(res_i[11]) >= 5.0 else "q5"
            try:
                ref_bnd = ref_chrom[int(res_i[3])]
            except:
                ref_bnd = 'N'
            if args.performing_phasing :
                GT_FORMAT = "GT:HP_GT:DR:DV:PL:AP:GQ"
            else :
                GT_FORMAT = "GT:DR:DV:PL:AP:GQ"
            
            vcf_line = "{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}".format(
                CHR = res_i[0], 
                POS = str(int(res_i[2]) + 1), 
                ID = "cuteSVTrio.%s.<SVID>"%("BND"),
                REF = ref_bnd,
                ALT = res_i[1], 
                INFO = info_list, 
                FORMAT = GT_FORMAT, 
                QUAL = res_i[11],
                PASS = filter_lable)
            for j in range(len(semi_result_ls)) :
                res_j = semi_result_ls[j][i]
                PL = "%s,%s,%s"%(res_j[9].split(",")[0],res_j[9].split(",")[1],res_j[9].split(",")[2])
                AP = "%s,%s"%(res_j[9].split(",")[3],res_j[9].split(",")[4])
                vcf_line = vcf_line + "\t{GT}:{DR}:{RE}:{PL}:{AP}:{GQ}".format(
                                        GT = res_j[8],
                                        DR = res_j[7],
                                        RE = res_j[4],
                                        PL = PL,
                                        AP = AP,
                                        GQ = res_j[10])
            lines.append(("BND",vcf_line+"\n"))
        if len(lines)>BATCH_SIZE:
            pickle.dump(lines,f)
            lines=[]
    if len(lines)!=0:
        pickle.dump(lines,f)


# def generate_pvcf(args, result, contigINFO, argv, ref_g):
def generate_pvcf(args, result, reference, chrom):
    fa_file = pysam.FastaFile(reference)
    try:
        ref_chrom=fa_file.fetch(chrom)
    except:
        raise Exception("No corresponding contig in reference with %s."%(chrom))
    fa_file.close()
    lines=[]
    for i in result:
        if i == []:
            continue
        if i[7][5] == "." or i[7][5] == None:
            filter_lable = "PASS"
        else:
            filter_lable = "PASS" if float(i[7][5]) >= 2.5 else "q5"
        if i[3] == 'INS':
            if abs(i[14]) > args.max_size and args.max_size != -1:
                continue
            ref = str(ref_chrom[max(i[1]-1, 0)])
            alt = i[11]
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if i[2] == "0/0" else "PRECISE", 
                SVTYPE = i[3], 
                SVLEN = i[14], 
                END = i[1], 
                CIPOS = i[5], 
                CILEN = i[6], 
                RE = i[7][0],
                RNAMES = i[8] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(i[7][0] / (i[7][0] + i[7][1]), 4))
            except:
                info_list += ";AF=."
            lines.append("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = i[1], 
                ID = i[9],
                REF = ref,
                ALT = alt, 
                QUAL = i[7][5],
                PASS = filter_lable,
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[2],
                DR = i[7][1],
                RE = i[7][0],
                PL = i[7][3],
                GQ = i[7][4]
                ))
        elif i[3] == 'DEL':
            if abs(i[14]) > args.max_size and args.max_size != -1:
                continue
            if i[12] == '<DEL>':
                ref = str(ref_chrom[max(int(i[1])-1, 0):int(i[1])-int(i[4])])
                alt = str(ref_chrom[max(int(i[1])-1, 0)])
            else:
                ref = i[10]
                alt = i[11]
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};RNAMES={RNAMES};STRAND=+-".format(
                PRECISION = "IMPRECISE" if i[2] == "0/0" else "PRECISE", 
                SVTYPE = i[3], 
                SVLEN = -abs(i[14]), 
                END = i[1] + abs(i[14]), 
                CIPOS = i[5], 
                CILEN = i[6], 
                RE = i[7][0],
                RNAMES = i[8] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(i[7][0] / (i[7][0] + i[7][1]), 4))
            except:
                info_list += ";AF=."
            lines.append("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = i[1], 
                ID = i[9],
                REF = ref,
                ALT = alt, 
                QUAL = i[7][5],
                PASS = filter_lable,
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[2],
                DR = i[7][1],
                RE = i[7][0],
                PL = i[7][3],
                GQ = i[7][4]
                ))
        elif i[3] == 'DUP':
            if abs(i[4] - i[1]) > args.max_size and args.max_size != -1:
                continue
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};RNAMES={RNAMES};STRAND=-+".format(
                PRECISION = "IMPRECISE" if i[2] == "0/0" else "PRECISE", 
                SVTYPE = i[3], 
                SVLEN = abs(i[4] - i[1]), 
                END = i[4], 
                RE = i[7][0],
                RNAMES = i[8] if args.report_readid else "NULL")
            try:
                info_list += ";AF=" + str(round(i[7][0] / (i[7][0] + i[7][1]), 4))
            except:
                info_list += ";AF=."
            lines.append("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = i[1], 
                ID = i[9],
                REF = i[10],
                ALT = i[11], 
                QUAL = i[7][5],
                PASS = filter_lable,
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[2],
                DR = i[7][1],
                RE = i[7][0],
                PL = i[7][3],
                GQ = i[7][4]
                ))
        elif i[3] == 'INV':
            if abs(i[4] - i[1]) > args.max_size and args.max_size != -1:
                continue
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if i[2] == "0/0" else "PRECISE", 
                SVTYPE = i[3], 
                SVLEN = i[4] - i[1], 
                END = i[4], 
                RE = i[7][0],
                RNAMES = i[8] if args.report_readid else "NULL")
            if i[12] != '.':
                info_list += ';STRAND=' + i[12]
            try:
                info_list += ";AF=" + str(round(i[7][0] / (i[7][0] + i[7][1]), 4))
            except:
                info_list += ";AF=."
            lines.append("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = i[1], 
                ID = i[9],
                REF = i[10],
                ALT = i[11], 
                QUAL = i[7][5],
                PASS = filter_lable,
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[2],
                DR = i[7][1],
                RE = i[7][0],
                PL = i[7][3],
                GQ = i[7][4]
                ))
        else:
            # BND
            info_list = "{PRECISION};SVTYPE={SVTYPE};RE={RE};RNAMES={RNAMES}".format(
                    PRECISION = "IMPRECISE" if i[2] == "0/0" else "PRECISE", 
                    SVTYPE = i[3], 
                    RE = i[7][0],
                    RNAMES = i[8] if args.report_readid else "NULL")
            if i[14] != 0:
                info_list += ';SVLEN=%d'%(i[14])
            try:
                info_list += ";AF=" + str(round(i[7][0] / (i[7][0] + i[7][1]), 4))
            except:
                info_list += ";AF=."
            '''
            if ':' in i[15]:
                info_list += ";CHR2={CHR2};END={END}".format(
                    CHR2 = i[15].split(':')[0],
                    END = i[15].split(':')[1])
            '''
            lines.append("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = str(i[1]), 
                ID = i[9],
                REF = i[10],
                ALT = i[11], 
                QUAL = i[7][5],
                PASS = filter_lable,
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[2],
                DR = i[7][1],
                RE = i[7][0],
                PL = i[7][3],
                GQ = i[7][4]
                ))
    return lines

            
def load_valuable_chr(path):
    valuable_chr = dict()
    valuable_chr["DEL"] = list()
    valuable_chr["DUP"] = list()
    valuable_chr["INS"] = list()
    valuable_chr["INV"] = list()
    valuable_chr["TRA"] = dict()

    for svtype in ["DEL", "DUP", "INS", "INV"]:
        file = open("%s%s.sigs"%(path, svtype), 'r')
        for line in file:
            chr = line.strip('\n').split('\t')[1]
            if chr not in valuable_chr[svtype]:
                valuable_chr[svtype].append(chr)
        file.close()
        valuable_chr[svtype].sort()

    file = open("%s%s.sigs"%(path, "TRA"), 'r')
    for line in file:
        chr1 = line.strip('\n').split('\t')[1]
        chr2 = line.strip('\n').split('\t')[4]
        
        if chr1 not in valuable_chr["TRA"]:
            valuable_chr["TRA"][chr1] = list()
        if chr2 not in valuable_chr["TRA"][chr1]:
            valuable_chr["TRA"][chr1].append(chr2)

    file.close()
    for chr1 in valuable_chr["TRA"]:
        valuable_chr["TRA"][chr1].sort()

    return valuable_chr

def load_bed(bed_file, Task_list):
    bed_regions = dict()
    if bed_file != None:
        # only consider regions in BED file
        with open(bed_file, 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[0] not in bed_regions:
                    bed_regions[seq[0]] = list()
                bed_regions[seq[0]].append((int(seq[1]) - 1000, int(seq[2]) + 1000))
        region_list = [[] for i in range(len(Task_list))]
        for chrom in bed_regions:
            bed_regions[chrom].sort(key = lambda x:(x[0], x[1]))
            for item in bed_regions[chrom]:
                for i in range(len(Task_list)):
                    if chrom == Task_list[i][0]:
                        if (Task_list[i][1] <= item[0] and Task_list[i][2] > item[0]) or item[0] <= Task_list[i][1] < item[1]:
                            region_list[i].append(item)
        assert len(region_list) == len(Task_list), "parse bed file error"
        return region_list
    else:
        return None

def unsolvable_correction(candidate_single_SV_gt_fam_ls, svtype, family_mode) :
    gt_tag = ["0/0","0/1","1/1"]
    gl_index = 9
    gt_index = 8
    if family_mode == "M1" :
        unsolvable_num = 0
        for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
            candidate_SV = candidate_single_SV_gt_fam_ls[0][i]
            gl0,gl1,gl2 = [pow(10,int(x)/-10) for x in candidate_SV[gl_index].split(",")[0:3]]
            x1,x2,p_data,c0,c1 = [float(x) for x in candidate_SV[gl_index].split(",")[3:8]]
            if int(c0) == int(c1) == 0 :
                continue
            p_data = pow(10,-1*p_data)
            a = 1
            b = gl0 - gl2 - 1
            c = gl2
            if b**2 - 4*a*c >= deta_limit :
                continue
            p_data_f_m = np.float64(pow(err, c0)*pow((1-err), c1))
            p_data_notf_notm = np.float64(pow(0.5, c0+c1))
            gl_f = [pow(10,int(x)/-10) for x in candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[0:3]]
            gl_m = [pow(10,int(x)/-10) for x in candidate_single_SV_gt_fam_ls[2][i][gl_index].split(",")[0:3]]
            p_f = 0.5*gl_f[1] + gl_f[2]
            p_m = 0.5*gl_m[1] + gl_m[2]
            p_f_data = (p_data_f_m*p_f*p_m + p_data_notf_notm*p_f*(1-p_m)) / p_data
            p_m_data = (p_data_f_m*p_f*p_m + p_data_notf_notm*(1-p_f)*p_m) / p_data
            if p_f_data > 1 or p_m_data > 1 :
                p_f_data,p_m_data = p_f_data / max(p_f_data,p_m_data),p_m_data / max(p_f_data,p_m_data)
            t0 = (1-p_f_data)*(1-p_m_data)
            t1 = (1-p_f_data)*p_m_data + p_f_data*(1-p_m_data)
            t2 = p_f_data*p_m_data
            gl0 = (1-p_f_data)*(1-p_m_data)
            gl1 = (1-p_f_data)*p_m_data + p_f_data*(1-p_m_data)
            gl2 = p_f_data*p_m_data
            p_data = int(np.around(-10*log10(p_data)))
            candidate_single_SV_gt_fam_ls[0][i][gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(-1, -1, -1, min(1,p_f_data), min(1,p_m_data), p_data, c0, c1)
            modify_gl_string(candidate_single_SV_gt_fam_ls[0][i],svtype)
            unsolvable_num += 1
        for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
            candidate_SV_s = candidate_single_SV_gt_fam_ls[0][i]
            candidate_SV_f = candidate_single_SV_gt_fam_ls[1][i]
            candidate_SV_m = candidate_single_SV_gt_fam_ls[2][i]
            gl_s_0,gl_s_1,gl_s_2 = [pow(10,int(x)/-10) for x in candidate_SV_s[gl_index].split(",")[0:3]]
            x_s_1,x_s_2,p_data_s,c_s_0,c_s_1 = [float(x) for x in candidate_SV_s[gl_index].split(",")[3:8]]
            gl_f_0,gl_f_1,gl_f_2 = [pow(10,int(x)/-10) for x in candidate_SV_f[gl_index].split(",")[0:3]]
            x_f_1,x_f_2,p_data_f,c_f_0,c_f_1 = [float(x) for x in candidate_SV_f[gl_index].split(",")[3:8]]
            gl_m_0,gl_m_1,gl_m_2 = [pow(10,int(x)/-10) for x in candidate_SV_m[gl_index].split(",")[0:3]]
            x_m_1,x_m_2,p_data_m,c_m_0,c_m_1 = [float(x) for x in candidate_SV_m[gl_index].split(",")[3:8]]
            if int(c_s_0) == int(c_s_1) == 0 or int(c_f_0) == int(c_f_1) == 0 or int(c_m_0) == int(c_m_1) == 0:
                continue
            a_f = 1
            b_f = gl_f_0 - gl_f_2 - 1
            c_f = gl_f_2
            a_m = 1
            b_m = gl_m_0 - gl_m_2 - 1
            c_m = gl_m_2
            if b_f**2 - 4*a_f*c_f >= deta_limit and b_m**2 - 4*a_m*c_m >= deta_limit:
                continue
            if b_f**2 - 4*a_f*c_f >= deta_limit and b_m**2 - 4*a_m*c_m < deta_limit:
                dis_ls = [abs(x_f_1-x_s_1),abs(x_f_1-x_s_2),abs(x_f_2-x_s_1),abs(x_f_2-x_s_2)]
                min_index = dis_ls.index(min(dis_ls))
                f_index = min_index // 2
                s_index = min_index % 2
                if s_index == 0 :
                    x_m_1 = x_s_2
                else :
                    x_m_1 = x_s_1
                if x_m_1 == 0 :
                    x_m_2 = 1-gl_m_0
                else :
                    x_m_2 = gl_m_2 / x_m_1
                gl_m_0 = (1-x_m_1) * (1-x_m_2)
                gl_m_1 = (1-x_m_1) * x_m_2 + x_m_1 * (1-x_m_2)
                gl_m_2 = x_m_1 * x_m_2
                candidate_single_SV_gt_fam_ls[2][i][gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(-1, -1, -1, min(1,x_m_1), min(1,x_m_2), p_data_m, c_m_0, c_m_1)
                modify_gl_string(candidate_single_SV_gt_fam_ls[2][i],svtype)
                continue
            if b_f**2 - 4*a_f*c_f < deta_limit and b_m**2 - 4*a_m*c_m >= deta_limit:
                dis_ls = [abs(x_m_1-x_s_1),abs(x_m_1-x_s_2),abs(x_m_2-x_s_1),abs(x_m_2-x_s_2)]
                min_index = dis_ls.index(min(dis_ls))
                m_index = min_index // 2
                s_index = min_index % 2
                if s_index == 0 :
                    x_f_1 = x_s_2
                else :
                    x_f_1 = x_s_1
                if x_f_1 == 0 :
                    x_f_2 = 1-gl_f_0
                else :
                    x_f_2 = gl_f_2 / x_f_1
                gl_f_0 = (1-x_f_1) * (1-x_f_2)
                gl_f_1 = (1-x_f_1) * x_f_2 + x_f_1 * (1-x_f_2)
                gl_f_2 = x_f_1 * x_f_2
                candidate_single_SV_gt_fam_ls[1][i][gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(-1,-1,-1, min(1,x_f_1), min(1,x_f_2), p_data_f, c_f_0, c_f_1)
                modify_gl_string(candidate_single_SV_gt_fam_ls[1][i],svtype)
                continue
            if b_f**2 - 4*a_f*c_f < deta_limit and b_m**2 - 4*a_m*c_m < deta_limit:
                if random.randint(1,2) == 1 :
                    x_f_1 = x_s_1
                    x_m_1 = x_s_2
                else :
                    x_f_1 = x_s_2
                    x_m_1 = x_s_1
                if x_f_1 == 0 :
                    x_f_2 = 1-gl_f_0
                else :
                    x_f_2 = gl_f_2 / x_f_1
                gl_f_0 = (1-x_f_1) * (1-x_f_2)
                gl_f_1 = (1-x_f_1) * x_f_2 + x_f_1 * (1-x_f_2)
                gl_f_2 = x_f_1 * x_f_2
                candidate_single_SV_gt_fam_ls[1][i][gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(-1,-1,-1, min(1,x_f_1), min(1,x_f_2), p_data_f, c_f_0, c_f_1)
                modify_gl_string(candidate_single_SV_gt_fam_ls[1][i],svtype)
                if x_m_1 == 0 :
                    x_m_2 = 1-gl_m_0
                else :
                    x_m_2 = gl_m_2 / x_m_1
                gl_m_0 = (1-x_m_1) * (1-x_m_2)
                gl_m_1 = (1-x_m_1) * x_m_2 + x_m_1 * (1-x_m_2)
                gl_m_2 = x_m_1 * x_m_2
                candidate_single_SV_gt_fam_ls[2][i][gl_index] = "%d,%d,%d,%f,%f,%d,%f,%f,0"%(-1,-1,-1, min(1,x_m_1), min(1,x_m_2), p_data_m, c_m_0, c_m_1)
                modify_gl_string(candidate_single_SV_gt_fam_ls[2][i],svtype)
                continue

def modify_genotype(candidate_single_SV_gt_fam_ls, svtype) :
    gl_index = 9
    gt_index = 8
    gq_index = 10
    qual_index = 11
    for i in range(len(candidate_single_SV_gt_fam_ls)) :
        for j in range(len(candidate_single_SV_gt_fam_ls[i])) :
            gl_split_ls = candidate_single_SV_gt_fam_ls[i][j][gl_index].split(",") 
            ori_GL00 = max((1-float(gl_split_ls[3]))*(1-float(gl_split_ls[4])),1e-100)
            ori_GL01 = max((1-float(gl_split_ls[3]))*float(gl_split_ls[4])+float(gl_split_ls[3])*(1-float(gl_split_ls[4])),1e-100)
            ori_GL11 = max(float(gl_split_ls[3])*float(gl_split_ls[4]),1e-100)
            prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
            GL_P = [pow(10, i) for i in prob]
            PL = [int(np.around(-10*log10(i))) for i in GL_P]
            gl_split_ls[0] = PL[0]
            gl_split_ls[1] = PL[1]
            gl_split_ls[2] = PL[2]
            candidate_single_SV_gt_fam_ls[i][j][gt_index] = Genotype[prob.index(max(prob))]
            candidate_single_SV_gt_fam_ls[i][j][gl_index] = ",".join(gl_split_ls)
            candidate_single_SV_gt_fam_ls[i][j][gq_index] = max([int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))])
            candidate_single_SV_gt_fam_ls[i][j][qual_index] = abs(np.around(-10*log10(GL_P[0]), 1))

def modify_gl_string(candidate_single_SV,svtype) :
    gl_index = 9
    gt_index = 8
    gq_index = 10
    qual_index = 11
    
    gl_split_ls = candidate_single_SV[gl_index].split(",")
    ori_GL00 = max((1-float(gl_split_ls[3]))*(1-float(gl_split_ls[4])),1e-100)
    ori_GL01 = max((1-float(gl_split_ls[3]))*float(gl_split_ls[4])+float(gl_split_ls[3])*(1-float(gl_split_ls[4])),1e-100)
    ori_GL11 = max(float(gl_split_ls[3])*float(gl_split_ls[4]),1e-100)
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    gl_split_ls[0] = str(PL[0])
    gl_split_ls[1] = str(PL[1])
    gl_split_ls[2] = str(PL[2])
    candidate_single_SV[gt_index] = Genotype[prob.index(max(prob))]
    candidate_single_SV[gl_index] = ",".join(gl_split_ls)
    candidate_single_SV[gq_index] = str(max([int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]))
    candidate_single_SV[qual_index] = str(abs(np.around(-10*log10(GL_P[0]), 1)))

def increase_sigs_through_pedigree(candidate_single_SV_gt_fam_ls, svtype, minimum_support_reads_list, family_mode) :
    sv_len_distance_threshold = 0.2
    gl_index = 9
    gt_index = 8
    gq_index = 10
    qual_index = 11
    sv_len_index = 14
    non_read_nameds_index = 16
    non_read_poss_index = 17
    process_thres = 0
    modify_rate = 1
    if svtype not in ["DEL","INS","INV","DUP"] :
        return
    length_limit = 1500
    if family_mode == 'M1' :
        # effective signal supplement
        for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
            sv_gl_fam = []
            sv_gl_fam.append([int(float(v)) for v in candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[6:8]])
            sv_gl_fam.append([int(float(v)) for v in candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[6:8]])
            sv_gl_fam.append([int(float(v)) for v in candidate_single_SV_gt_fam_ls[2][i][gl_index].split(",")[6:8]])
            sv_gp_fam = []
            member_gp = [pow(10,int(float(v)/-10)) for v in candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[0:3]]
            sv_gp_fam.append([v/sum(member_gp) for v in member_gp])
            member_gp = [pow(10,int(float(v)/-10)) for v in candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[0:3]]
            sv_gp_fam.append([v/sum(member_gp) for v in member_gp])
            member_gp = [pow(10,int(float(v)/-10)) for v in candidate_single_SV_gt_fam_ls[2][i][gl_index].split(",")[0:3]]
            sv_gp_fam.append([v/sum(member_gp) for v in member_gp])
            fam_sv_len = [int(x) for x in candidate_single_SV_gt_fam_ls[0][i][sv_len_index].split(",")]
            
            # if the total number of sv sigs exceeds minimum_support_reads, or if the total number of sigs for a member is less than minimum_support_reads but the sv length exceeds 1500, it will no longer be forced to be 0/0
            for j in range(3) :
                if minimum_support_reads_list[j] < process_thres :
                    continue
                if (sv_gl_fam[j][0]+sv_gl_fam[j][1]) != 0 and sv_gl_fam[j][1] < minimum_support_reads_list[j] and fam_sv_len[j] > length_limit :
                    _,gl_str,GQ,QUAL = cal_GL_3(sv_gl_fam[j][0], sv_gl_fam[j][1], svtype, minimum_support_reads_list[j], fam_sv_len[j])
                    gl_str_split = gl_str.split(",")
                    gl_str_split[8] = "-2"
                    candidate_single_SV_gt_fam_ls[j][i][gl_index] = ",".join(gl_str_split)
                    candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(GQ)
                    candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(QUAL)
                    gt1 = float(gl_str.split(",")[3])
                    gt2 = float(gl_str.split(",")[4])
                    if gt1 >= 0.5 and gt2 >= 0.5 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/1'
                    elif gt1 < 0.5 and gt2 < 0.5 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '0/0'
                    else :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/0'
                    
            
            # c0=c1=0
            # since no read information, the theoretical Mendel correction is used.
            for j in range(3) :
                if int(candidate_single_SV_gt_fam_ls[0][i][3]) < 1500 :
                    break
                if sv_gl_fam[j][0]+sv_gl_fam[j][1] == 0 :
                    # if it's the child, genotype is determined using the parents' genotypes and Mendel's laws; if it's the parents' genotypes, it's consistent with the offspring's.
                    if j == 0 :
                        gpf = sv_gp_fam[1][1]/2 + sv_gp_fam[1][2]
                        gpm = sv_gp_fam[2][1]/2 + sv_gp_fam[2][2]
                        GL_P = [(1-gpf)*(1-gpm),gpf*(1-gpm)+(1-gpf)*gpm,gpf*gpm]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                    else :
                        GL_P = [v/sum(sv_gp_fam[0]) for v in sv_gp_fam[0]]
                    if GL_P.index(max(GL_P)) == 2 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/1'
                        candidate_single_SV_gt_fam_ls[j][i][gl_index] = "100,100,0,1,1,0,0,0,-1"
                        candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(996)
                        candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(255)
                    else :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '0/1'
                        candidate_single_SV_gt_fam_ls[j][i][gl_index] = "100,0,100,0,1,0,0,0,-1"
                        candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(996)
                        candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(255)

            # the number of sigs is sufficient, but the confidence is insufficient in determining the gt.
            # child
            if minimum_support_reads_list[0] < process_thres :
                pass
            else :
                supple_support_reads_threshold = max(ceil(minimum_support_reads_list[0]/2),2)
                if (candidate_single_SV_gt_fam_ls[0][i][gt_index] == '0/0' or (svtype in ["INS","DEL","INV","DUP"] and float(candidate_single_SV_gt_fam_ls[0][i][qual_index]) < 5)) and sv_gl_fam[0][1] >= supple_support_reads_threshold:
                    new_c0 = sv_gl_fam[0][0]
                    new_c1 = sv_gl_fam[0][1]
                    if candidate_single_SV_gt_fam_ls[1][i][gt_index] != '0/0' and (abs(fam_sv_len[0]-fam_sv_len[1]) / max(fam_sv_len[0],fam_sv_len[1])) < sv_len_distance_threshold:
                        new_c0 += sv_gl_fam[1][0]
                        new_c1 += sv_gl_fam[1][1]
                    if candidate_single_SV_gt_fam_ls[2][i][gt_index] != '0/0' and (abs(fam_sv_len[0]-fam_sv_len[2]) / max(fam_sv_len[0],fam_sv_len[2])) < sv_len_distance_threshold:
                        new_c0 += sv_gl_fam[2][0]
                        new_c1 += sv_gl_fam[2][1]
                    if new_c0 == sv_gl_fam[0][0] and new_c1 == sv_gl_fam[0][1] :
                        pass
                    else :
                        gl_str = cal_Gl_3_sim(new_c0,new_c1)
                        modify_fam_gp = [pow(10,float(v)) for v in gl_str.split(",")[0:3]]
                        modify_fam_gp = [v/sum(modify_fam_gp) for v in modify_fam_gp]
                        old_gp = sv_gp_fam[0]
                        GL_P = [old_gp[v]*(1-modify_rate)+modify_fam_gp[v]*modify_rate for v in range(len(old_gp))]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                        PL = [int(np.around(-10*log10(i))) for i in GL_P]
                        GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
                        QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
                        a = 1
                        b = GL_P[0] - GL_P[2] - 1
                        c = GL_P[2]
                        roots = np.roots([a, b, c])
                        if GL_P.index(max(GL_P)) > 0 :
                            gl_str = "%d,%d,%d,%f,%f,%d,%f,%f,-3"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, int(candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[5]), sv_gl_fam[0][0], sv_gl_fam[0][1])
                            candidate_single_SV_gt_fam_ls[0][i][gl_index] = gl_str
                            candidate_single_SV_gt_fam_ls[0][i][gt_index] = Genotype[GL_P.index(max(GL_P))]
                            candidate_single_SV_gt_fam_ls[0][i][gq_index] = str(max(GQ))
                            candidate_single_SV_gt_fam_ls[0][i][qual_index] = str(QUAL)

            # parents
            if minimum_support_reads_list[1] < process_thres :
                pass
            else :
                supple_support_reads_threshold = max(ceil(minimum_support_reads_list[1]/2),2)
                if (candidate_single_SV_gt_fam_ls[1][i][gt_index] == '0/0' or (svtype in ["INS","DEL","INV","DUP"] and float(candidate_single_SV_gt_fam_ls[1][i][qual_index]) < 5)) and sv_gl_fam[1][1] >= supple_support_reads_threshold:
                    if candidate_single_SV_gt_fam_ls[0][i][gt_index] != '0/0' :
                        gl_str = cal_Gl_3_sim(sv_gl_fam[0][0]+sv_gl_fam[1][0],sv_gl_fam[0][1]+sv_gl_fam[1][1])
                        modify_fam_gp = [pow(10,float(v)) for v in gl_str.split(",")[0:3]]
                        modify_fam_gp = [v/sum(modify_fam_gp) for v in modify_fam_gp]
                        old_gp = sv_gp_fam[1]
                        GL_P = [old_gp[v]*(1-modify_rate)+modify_fam_gp[v]*modify_rate for v in range(len(old_gp))]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                        PL = [int(np.around(-10*log10(i))) for i in GL_P]
                        GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
                        QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
                        a = 1
                        b = GL_P[0] - GL_P[2] - 1
                        c = GL_P[2]
                        roots = np.roots([a, b, c])
                        if GL_P.index(max(GL_P)) > 0 and (abs(fam_sv_len[0]-fam_sv_len[1]) / max(fam_sv_len[0],fam_sv_len[1])) < sv_len_distance_threshold:
                            gl_str = "%d,%d,%d,%f,%f,%d,%f,%f,-3"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, int(candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[5]), sv_gl_fam[1][0], sv_gl_fam[1][1])
                            candidate_single_SV_gt_fam_ls[1][i][gl_index] = gl_str
                            candidate_single_SV_gt_fam_ls[1][i][gt_index] = Genotype[GL_P.index(max(GL_P))]
                            candidate_single_SV_gt_fam_ls[1][i][gq_index] = str(max(GQ))
                            candidate_single_SV_gt_fam_ls[1][i][qual_index] = str(QUAL)
            if minimum_support_reads_list[2] < process_thres :
                pass
            else :
                supple_support_reads_threshold = max(ceil(minimum_support_reads_list[2]/2),2)
                if (candidate_single_SV_gt_fam_ls[2][i][gt_index] == '0/0' or (svtype in ["INS","DEL","INV","DUP"] and float(candidate_single_SV_gt_fam_ls[2][i][qual_index]) < 5)) and sv_gl_fam[2][1] >= supple_support_reads_threshold:
                    if candidate_single_SV_gt_fam_ls[0][i][gt_index] != '0/0' :
                        gl_str = cal_Gl_3_sim(sv_gl_fam[0][0]+sv_gl_fam[2][0],sv_gl_fam[0][1]+sv_gl_fam[2][1])
                        modify_fam_gp = [pow(10,float(v)) for v in gl_str.split(",")[0:3]]
                        modify_fam_gp = [v/sum(modify_fam_gp) for v in modify_fam_gp]
                        old_gp = sv_gp_fam[2]
                        GL_P = [old_gp[v]*(1-modify_rate)+modify_fam_gp[v]*modify_rate for v in range(len(old_gp))]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                        PL = [int(np.around(-10*log10(i))) for i in GL_P]
                        GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
                        QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
                        a = 1
                        b = GL_P[0] - GL_P[2] - 1
                        c = GL_P[2]
                        roots = np.roots([a, b, c])
                        if GL_P.index(max(GL_P)) > 0 and (abs(fam_sv_len[0]-fam_sv_len[2]) / max(fam_sv_len[0],fam_sv_len[2])) < sv_len_distance_threshold:
                            gl_str = "%d,%d,%d,%f,%f,%d,%f,%f,-3"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, int(candidate_single_SV_gt_fam_ls[2][i][gl_index].split(",")[5]), sv_gl_fam[2][0], sv_gl_fam[2][1])
                            candidate_single_SV_gt_fam_ls[2][i][gl_index] = gl_str
                            candidate_single_SV_gt_fam_ls[2][i][gt_index] = Genotype[GL_P.index(max(GL_P))]
                            candidate_single_SV_gt_fam_ls[2][i][gq_index] = str(max(GQ))
                            candidate_single_SV_gt_fam_ls[2][i][qual_index] = str(QUAL)
            for j in range(len(candidate_single_SV_gt_fam_ls)) :
                gl_ls = [round(float(x)) for x in candidate_single_SV_gt_fam_ls[j][i][gl_index].split(",")[3:5]]
                candidate_single_SV_gt_fam_ls[j][i][gt_index] = str(max(gl_ls))+"/"+str(min(gl_ls))
    elif family_mode == 'M2' :
        for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
            sv_gl_fam = []
            sv_gl_fam.append([int(float(v)) for v in candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[6:8]])
            sv_gl_fam.append([int(float(v)) for v in candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[6:8]])
            sv_gp_fam = []
            member_gp = [pow(10,int(float(v)/-10)) for v in candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[0:3]]
            sv_gp_fam.append([v/sum(member_gp) for v in member_gp])
            member_gp = [pow(10,int(float(v)/-10)) for v in candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[0:3]]
            sv_gp_fam.append([v/sum(member_gp) for v in member_gp])
            fam_sv_len = [int(x) for x in candidate_single_SV_gt_fam_ls[0][i][sv_len_index].split(",")]
            
            for j in range(2) :
                if minimum_support_reads_list[j] < process_thres :
                    continue
                if (sv_gl_fam[j][0]+sv_gl_fam[j][1]) != 0 and sv_gl_fam[j][1] < minimum_support_reads_list[j] and fam_sv_len[j] > length_limit :
                    _,gl_str,GQ,QUAL = cal_GL_3(sv_gl_fam[j][0], sv_gl_fam[j][1], svtype, minimum_support_reads_list[j], fam_sv_len[j])
                    gl_str_split = gl_str.split(",")
                    gl_str_split[8] = "-2"
                    candidate_single_SV_gt_fam_ls[j][i][gl_index] = ",".join(gl_str_split)
                    candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(GQ)
                    candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(QUAL)
                    gt1 = float(gl_str.split(",")[3])
                    gt2 = float(gl_str.split(",")[4])
                    if gt1 >= 0.5 and gt2 >= 0.5 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/1'
                    elif gt1 < 0.5 and gt2 < 0.5 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '0/0'
                    else :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/0'
            
            for j in range(2) :
                if int(candidate_single_SV_gt_fam_ls[0][i][3]) < 1500 :
                    break
                if sv_gl_fam[j][0]+sv_gl_fam[j][1] == 0 :
                    if j == 0 :
                        GL_P = [v/sum(sv_gp_fam[1]) for v in sv_gp_fam[1]]
                    else :
                        GL_P = [v/sum(sv_gp_fam[0]) for v in sv_gp_fam[0]]
                    if GL_P.index(max(GL_P)) == 2 :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '1/1'
                        candidate_single_SV_gt_fam_ls[j][i][gl_index] = "100,100,0,1,1,0,0,0,-1"
                        candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(996)
                        candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(255)
                    else :
                        candidate_single_SV_gt_fam_ls[j][i][gt_index] = '0/1'
                        candidate_single_SV_gt_fam_ls[j][i][gl_index] = "100,0,100,0,1,0,0,0,-1"
                        candidate_single_SV_gt_fam_ls[j][i][gq_index] = str(996)
                        candidate_single_SV_gt_fam_ls[j][i][qual_index] = str(255)

            if minimum_support_reads_list[0] < process_thres :
                pass
            else :
                supple_support_reads_threshold = max(ceil(minimum_support_reads_list[0]/2),2)
                if (candidate_single_SV_gt_fam_ls[0][i][gt_index] == '0/0' or (svtype in ["INS","DEL","INV","DUP"] and float(candidate_single_SV_gt_fam_ls[0][i][qual_index]) < 5)) and sv_gl_fam[0][1] >= supple_support_reads_threshold:
                    new_c0 = sv_gl_fam[0][0]
                    new_c1 = sv_gl_fam[0][1]
                    if candidate_single_SV_gt_fam_ls[1][i][gt_index] != '0/0' and (abs(fam_sv_len[0]-fam_sv_len[1]) / max(fam_sv_len[0],fam_sv_len[1])) < sv_len_distance_threshold:
                        new_c0 += sv_gl_fam[1][0]
                        new_c1 += sv_gl_fam[1][1]
                    if new_c0 == sv_gl_fam[0][0] and new_c1 == sv_gl_fam[0][1] :
                        pass
                    else :
                        gl_str = cal_Gl_3_sim(new_c0,new_c1)
                        modify_fam_gp = [pow(10,float(v)) for v in gl_str.split(",")[0:3]]
                        modify_fam_gp = [v/sum(modify_fam_gp) for v in modify_fam_gp]
                        old_gp = sv_gp_fam[0]
                        GL_P = [old_gp[v]*(1-modify_rate)+modify_fam_gp[v]*modify_rate for v in range(len(old_gp))]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                        PL = [int(np.around(-10*log10(i))) for i in GL_P]
                        GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
                        QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
                        a = 1
                        b = GL_P[0] - GL_P[2] - 1
                        c = GL_P[2]
                        roots = np.roots([a, b, c])
                        if GL_P.index(max(GL_P)) > 0 :
                            gl_str = "%d,%d,%d,%f,%f,%d,%f,%f,-3"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, int(candidate_single_SV_gt_fam_ls[0][i][gl_index].split(",")[5]), sv_gl_fam[0][0], sv_gl_fam[0][1])
                            candidate_single_SV_gt_fam_ls[0][i][gl_index] = gl_str
                            candidate_single_SV_gt_fam_ls[0][i][gt_index] = Genotype[GL_P.index(max(GL_P))]
                            candidate_single_SV_gt_fam_ls[0][i][gq_index] = str(max(GQ))
                            candidate_single_SV_gt_fam_ls[0][i][qual_index] = str(QUAL)

            if minimum_support_reads_list[1] < process_thres :
                pass
            else :
                supple_support_reads_threshold = max(ceil(minimum_support_reads_list[1]/2),2)
                if (candidate_single_SV_gt_fam_ls[1][i][gt_index] == '0/0' or (svtype in ["INS","DEL","INV","DUP"] and float(candidate_single_SV_gt_fam_ls[1][i][qual_index]) < 5)) and sv_gl_fam[1][1] >= supple_support_reads_threshold:
                    if candidate_single_SV_gt_fam_ls[0][i][gt_index] != '0/0' :
                        gl_str = cal_Gl_3_sim(sv_gl_fam[0][0]+sv_gl_fam[1][0],sv_gl_fam[0][1]+sv_gl_fam[1][1])
                        modify_fam_gp = [pow(10,float(v)) for v in gl_str.split(",")[0:3]]
                        modify_fam_gp = [v/sum(modify_fam_gp) for v in modify_fam_gp]
                        old_gp = sv_gp_fam[1]
                        GL_P = [old_gp[v]*(1-modify_rate)+modify_fam_gp[v]*modify_rate for v in range(len(old_gp))]
                        GL_P = [v/sum(GL_P) for v in GL_P]
                        PL = [int(np.around(-10*log10(i))) for i in GL_P]
                        GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
                        QUAL = abs(np.around(-10*log10(GL_P[0]), 1))
                        a = 1
                        b = GL_P[0] - GL_P[2] - 1
                        c = GL_P[2]
                        roots = np.roots([a, b, c])
                        if GL_P.index(max(GL_P)) > 0 and (abs(fam_sv_len[0]-fam_sv_len[1]) / max(fam_sv_len[0],fam_sv_len[1])) < sv_len_distance_threshold:
                            gl_str = "%d,%d,%d,%f,%f,%d,%f,%f,-3"%(PL[0], PL[1], PL[2], roots[0].real, roots[1].real, int(candidate_single_SV_gt_fam_ls[1][i][gl_index].split(",")[5]), sv_gl_fam[1][0], sv_gl_fam[1][1])
                            candidate_single_SV_gt_fam_ls[1][i][gl_index] = gl_str
                            candidate_single_SV_gt_fam_ls[1][i][gt_index] = Genotype[GL_P.index(max(GL_P))]
                            candidate_single_SV_gt_fam_ls[1][i][gq_index] = str(max(GQ))
                            candidate_single_SV_gt_fam_ls[1][i][qual_index] = str(QUAL)

            for j in range(len(candidate_single_SV_gt_fam_ls)) :
                gl_ls = [round(float(x)) for x in candidate_single_SV_gt_fam_ls[j][i][gl_index].split(",")[3:5]]
                candidate_single_SV_gt_fam_ls[j][i][gt_index] = str(max(gl_ls))+"/"+str(min(gl_ls))

def check_gt_consistent(candidate_single_SV_gt_fam_ls) :
    gl_index = 9
    gt_index = 8
    for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
        for j in range(len(candidate_single_SV_gt_fam_ls)) :
            gl_ls = [round(float(x)) for x in candidate_single_SV_gt_fam_ls[j][i][gl_index].split(",")[3:5]]
            gt_ls = [int(x) for x in candidate_single_SV_gt_fam_ls[j][i][gt_index].split("/")]
            if max(gl_ls) != max(gt_ls) or min(gl_ls) != min(gt_ls) :
                logging.info(candidate_single_SV_gt_fam_ls[j][i])

def inconformity_mendel_modify(candidate_single_SV_gt_fam_ls, svtype, minimum_support_reads, family_mode) :
    family_mode_index_ls = ["M1","M2"]
    family_member_set = [["1","2","3"],["1","2"]]
    family_member_ls = family_member_set[family_mode_index_ls.index(family_mode)]
    gl_index = 9
    gt_index = 8
    fp_support_threshold = round(minimum_support_reads/2)
    for i in range(len(candidate_single_SV_gt_fam_ls[0])) :
        if abs(int(candidate_single_SV_gt_fam_ls[0][i][3])) > 100 :
            continue
        sv_gt_fam = []
        for j in range(len(family_member_ls)) :
            sv_gt_fam.append([int(v) for v in candidate_single_SV_gt_fam_ls[j][i][gt_index].split("/")])
        for j in range(len(family_member_ls)) :
            if sum(sv_gt_fam[j]) == 0 :
                continue
            fam_others = [int(x)-1 for x in family_member_ls]
            fam_others.remove(j)
            is_modify_ls = True
            for f_o in fam_others :
                if sum(sv_gt_fam[f_o]) == 0 and float(candidate_single_SV_gt_fam_ls[f_o][i][gl_index].split(",")[7]) <= fp_support_threshold :
                    pass
                else :
                    is_modify_ls = False
            if is_modify_ls :
                candidate_single_SV_gt_fam_ls[j][i][gt_index] = "0/0"
                break
        
def allele_correction(chr,member,candidate_single_SV_gt_ls, minimum_support_reads) :
    if chr not in chr_ls :
        return candidate_single_SV_gt_ls
    pre_index = -1
    gl_index = 9
    gt_index = 8
    gq_index = 10
    qual_index = 11
    for i in range(len(candidate_single_SV_gt_ls)) :
        if candidate_single_SV_gt_ls[i][1] in ["INS","DEL"] :
            if candidate_single_SV_gt_ls[i][gt_index] in ["0/1","1/0","1/1"]:
                pre_index = i
                continue
            if pre_index == -1 :
                continue
            if candidate_single_SV_gt_ls[i][1] != candidate_single_SV_gt_ls[pre_index][1] :
                if int(candidate_single_SV_gt_ls[i][2]) > int(candidate_single_SV_gt_ls[pre_index][2]) and int(candidate_single_SV_gt_ls[i][2]) + abs(int(candidate_single_SV_gt_ls[i][3])) < int(candidate_single_SV_gt_ls[pre_index][2]) + abs(int(candidate_single_SV_gt_ls[pre_index][3])) :
                    if candidate_single_SV_gt_ls[i][gt_index] == "0/0" and candidate_single_SV_gt_ls[pre_index][gt_index] in ["0/1","1/0","1/1"] :
                        c1 = int(float(candidate_single_SV_gt_ls[i][gl_index].split(",")[7]))
                        c0 = int(float(candidate_single_SV_gt_ls[i][gl_index].split(",")[6])) - int(float(candidate_single_SV_gt_ls[pre_index][gl_index].split(",")[7]))
                        if c1 < minimum_support_reads or c0 < 0 :
                            continue
                        GT, GL, GQ, QUAL = cal_GL_3(c0, c1, candidate_single_SV_gt_ls[i][1], minimum_support_reads, abs(int(candidate_single_SV_gt_ls[i][3])))
                        if GT in ["0/1","1/1"] :
                            if member == 0 and candidate_single_SV_gt_ls[i][gt_index] == "0/0" and GT == "1/1" :
                                sv_len_str_ls = candidate_single_SV_gt_ls[i][14].split(",")
                                if sv_len_str_ls[1] == sv_len_str_ls[2] == "0" :
                                    continue
                            candidate_single_SV_gt_ls[i][gt_index] = GT
                            candidate_single_SV_gt_ls[i][gl_index] = GL
                            candidate_single_SV_gt_ls[i][gq_index] = str(GQ)
                            candidate_single_SV_gt_ls[i][qual_index] = str(QUAL)
                            pre_index = i
    
    return candidate_single_SV_gt_ls
