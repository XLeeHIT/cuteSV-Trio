# -*- coding: utf-8 -*-

import sys
import random
import copy
import gzip
import numpy as np

chr_names_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
chrs_len_ls = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
    
def read_indel_answer_vcf(vcf_file,threshold_len) :
    header_vcf_lines = []
    indel_vcf_lines = []
    ins_len_ls = []
    del_len_ls = []
    with gzip.open(vcf_file, 'rt') as f:
        while(True) :
            f_line = f.readline().strip()
            if f_line == "" :
                break
            if f_line[0] == '#':
                header_vcf_lines.append(f_line)
            else :
                line_split_ls = f_line.split("\t")
                if line_split_ls[0] in chr_names_ls :
                    if "SVTYPE=INS" in line_split_ls[7] :
                        if len(line_split_ls[4]) >= 50 and len(line_split_ls[4]) >= threshold_len:
                            indel_vcf_lines.append(f_line)
                            ins_len_ls.append(len(line_split_ls[4]))
                    if "SVTYPE=DEL" in line_split_ls[7] :
                        if len(line_split_ls[3]) >= 50 and len(line_split_ls[3]) >= threshold_len:
                            indel_vcf_lines.append(f_line)
                            del_len_ls.append(len(line_split_ls[3]))
        f.close()
    return header_vcf_lines,indel_vcf_lines,ins_len_ls,del_len_ls

def read_snp_answer_vcf(vcf_file) :
    header_vcf_lines = []
    snp_vcf_lines = []
    with gzip.open(vcf_file, 'rt') as f:
        while(True) :
            f_line = f.readline().strip()
            if f_line == "" :
                break
            if f_line[0] == '#':
                header_vcf_lines.append(f_line)
            else :
                line_split_ls = f_line.split("\t")
                if line_split_ls[0] in chr_names_ls :
                    if len(line_split_ls[3]) == 1 and len(line_split_ls[4]) == 1:
                        snp_vcf_lines.append(f_line)
        f.close()
    return header_vcf_lines,snp_vcf_lines

def make_header() :
    header_vcf_lines,_,_,_ = read_indel_answer_vcf("/home/user/lixin/cutesv_trio/answer_vcf/HG002_T2T/GRCh37_HG2-T2TQ100-V1.1_stvar.sv.INDEL.vcf.gz",0)
    new_header_lines = []
    for header_lines in header_vcf_lines[:-1] :
        if header_lines[0:8] == "##contig" :
            pass
        else :
            new_header_lines.append(header_lines)
    new_header_lines.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV described in this region\">")
    new_header_lines.append("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Mate chromsome for BND SVs\">")
    new_header_lines.append("##INFO=<ID=TE,Number=1,Type=String,Description=\"TE class of insertion\">")
    for i in range(len(chr_names_ls)) :
        new_header_lines.append("##contig=<ID="+chr_names_ls[i]+",length="+str(chrs_len_ls[i])+">")
    new_header_lines.append(header_vcf_lines[-1])
    return new_header_lines

def make_simu_data_with_snp(file_base, with_snp) :
    answer_vcf_file_1 = file_base + "/answer.1.vcf"
    answer_vcf_file_2 = file_base + "/answer.2.vcf"
    answer_vcf_file_3 = file_base + "/answer.3.vcf"
    answer_vcf_file_fam = file_base + "/answer.fam.vcf"
    denovo_vcf_file = file_base + "/denovo.vcf"
    no_mendel_vcf_file = file_base + "/no_mendel.vcf"
    hap_choice_txt = file_base + "/hap_choice.txt"
    denovo_change_txt = file_base + "/denovo_change.txt"
    fam_bed_1 = file_base + "/HACK.h1.bed"
    fam_bed_2 = file_base + "/HACK.h2.bed"
    fam_bed_3 = file_base + "/HACK.h3.bed"
    fam_bed_4 = file_base + "/HACK.h4.bed"
    fam_bed_5 = file_base + "/HACK.h5.bed"
    fam_bed_6 = file_base + "/HACK.h6.bed"
    
    ans_vcf_file = "/home/user/lixin/cutesv_trio/answer_vcf/HG002_T2T/GRCh37_HG2-T2TQ100-V1.1_stvar.sv.INDEL.vcf.gz"
    _,indel_vcf_lines_total,ins_len_ls,del_len_ls = read_indel_answer_vcf(ans_vcf_file)

    if with_snp :
        snp_vcf_file_1 = file_base + "/snp.1.vcf"
        snp_vcf_file_2 = file_base + "/snp.2.vcf"
        snp_vcf_file_3 = file_base + "/snp.3.vcf"
        snp_vcf_file_fam = file_base + "/snp.fam.vcf"
        snp_vcf_file = "/home/user/lixin/cutesv_trio/answer_vcf/HG002_T2T/GRCh37_HG2-T2TQ100-V1.1_stvar.snp.vcf.gz"
        _,snp_vcf_lines_total = read_snp_answer_vcf(snp_vcf_file)

    header_vcf_lines = make_header()
    dna_ls = ["A","C","G","T"]
    indel_vcf_lines = []
    for i in indel_vcf_lines_total :
        if random.random() <= 0.75 :
            indel_vcf_lines.append(i)
    sv_num = 0
    chrs_sv_ls = [[] for i in chr_names_ls]
    # 同分布长度+原始位点(20000范围内扰动)+随机序列
    for sv_line in indel_vcf_lines :
        sv_split_ls = sv_line.split("\t")
        if "SVTYPE=INS" in sv_split_ls[7] :
            sv_len = random.choice(ins_len_ls)
            if int(sv_split_ls[1]) + sv_len > chrs_len_ls[chr_names_ls.index(sv_split_ls[0])] :
                continue
            sv_str = ""
            for i in range(sv_len) :
                sv_str = sv_str + dna_ls[random.randint(0,3)]
            pos_random_move = random.randint(-5000,5000)
            chrs_sv_ls[chr_names_ls.index(sv_split_ls[0])].append([sv_split_ls[0],"INS",int(sv_split_ls[1])+pos_random_move,int(sv_split_ls[1])+pos_random_move+1,sv_str])
        if "SVTYPE=DEL" in sv_split_ls[7] :
            sv_len = random.choice(del_len_ls)
            if int(sv_split_ls[1]) + sv_len > chrs_len_ls[chr_names_ls.index(sv_split_ls[0])] :
                continue
            pos_random_move = random.randint(-5000,5000)
            chrs_sv_ls[chr_names_ls.index(sv_split_ls[0])].append([sv_split_ls[0],"DEL",int(sv_split_ls[1])+pos_random_move,int(sv_split_ls[1])+pos_random_move+len(sv_split_ls[3]),"None"])
    
    if with_snp :
        snp_vcf_lines = []
        for i in snp_vcf_lines_total :
            if random.random() <= 0.75 :
                snp_vcf_lines.append(i)

        # 同分布长度+原始位点(200范围内扰动)+随机序列
        for snp_line in snp_vcf_lines :
            snp_split_ls = snp_line.split("\t")
            ref_len = random.choice(ref_len_ls)
            if ref_len > 1 :
                alt_len = 1
            else :
                alt_len = random.choice(alt_len_ls)
            snp_str = ""
            for i in range(alt_len) :
                snp_str = snp_str + dna_ls[random.randint(0,3)]
            pos_random_move = random.randint(-100,100)
            if ref_len == 1 and alt_len == 1 :
                chrs_sv_ls[chr_names_ls.index(snp_split_ls[0])].append([snp_split_ls[0],"SNP",int(snp_split_ls[1])+pos_random_move,int(snp_split_ls[1])+pos_random_move+1,snp_str])
            elif ref_len > 1 :
                chrs_sv_ls[chr_names_ls.index(snp_split_ls[0])].append([snp_split_ls[0],"SNP-DEL",int(snp_split_ls[1])+pos_random_move,int(snp_split_ls[1])+pos_random_move+ref_len,"None"])
            else :
                chrs_sv_ls[chr_names_ls.index(snp_split_ls[0])].append([snp_split_ls[0],"SNP-INS",int(snp_split_ls[1])+pos_random_move,int(snp_split_ls[1])+pos_random_move+1,snp_str])
    for i in range(len(chrs_sv_ls)) :
        chrs_sv_ls[i] = sorted(chrs_sv_ls[i], key = lambda x:x[2])
    chrs_sv_ls = remove_overlap(chrs_sv_ls)
    for i in chrs_sv_ls :
        sv_num += len(i)

    father_sv_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    mother_sv_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if random.random() <= 0.7 :
                father_sv_ls[i][0][j] = 1
                father_sv_ls[i][1][j] = 1
                mother_sv_ls[i][0][j] = 1
                mother_sv_ls[i][1][j] = 1
            else :
                if random.random() <= 0.5 :
                    father_sv_ls[i][0][j] = 1
                    father_sv_ls[i][1][j] = 1
                else :
                    mother_sv_ls[i][0][j] = 1
                    mother_sv_ls[i][1][j] = 1

    threshold_1 = 0.1
    threshold_2 = 0.25
    threshold_3 = 0.4
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if father_sv_ls[i][0][j] == 0 :
                continue
            random_value = random.random()
            if random_value <= threshold_1 :
                father_sv_ls[i][0][j] = 0
                father_sv_ls[i][1][j] = 0
            elif random_value > threshold_1 and random_value <= threshold_2 :
                father_sv_ls[i][0][j] = 0
            elif random_value > threshold_2 and random_value <= threshold_3 :
                father_sv_ls[i][1][j] = 0
            else :
                pass
    for i in range(len(mother_sv_ls)) :
        for j in range(len(mother_sv_ls[i][0])) :
            if mother_sv_ls[i][0][j] == 0 :
                continue
            random_value = random.random()
            if random_value <= threshold_1 :
                mother_sv_ls[i][0][j] = 0
                mother_sv_ls[i][1][j] = 0
            elif random_value > threshold_1 and random_value <= threshold_2 :
                mother_sv_ls[i][0][j] = 0
            elif random_value > threshold_2 and random_value <= threshold_3 :
                mother_sv_ls[i][1][j] = 0
            else :
                pass
    mendel_hap_choices = []
    child_sv_ls = [[[],[]] for chrs in chrs_sv_ls]
    for i in range(len(child_sv_ls)) :
        mendel_choice = []
        mendel_choice.append(random.randint(0,1))
        mendel_choice.append(random.randint(0,1))
        mendel_hap_choices.append(mendel_choice)
        child_sv_ls[i][0] = copy.deepcopy(father_sv_ls[i][mendel_choice[0]])
        child_sv_ls[i][1] = copy.deepcopy(mother_sv_ls[i][mendel_choice[1]])

    child_sv_ls_denovo_pre = copy.deepcopy(child_sv_ls)
    denovo_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            if random.random() <= 0.05 :
                denovo_record_ls[i][0][j] = 1
                if random.random() <= 0.5 :
                    child_sv_ls[i][0][j] = 1 - child_sv_ls[i][0][j]
                else :
                    child_sv_ls[i][1][j] = 1 - child_sv_ls[i][1][j]
    no_mendel_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    for i in range(len(denovo_record_ls)) :
        for j in range(len(denovo_record_ls[i][0])) :
            if denovo_record_ls[i][0][j] == 1 :
                mendel_distance = []
                mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][0][j]))
                mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][1][j]))
                mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][0][j]))
                mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][1][j]))
                mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][0][j]))
                mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][1][j]))
                mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][0][j]))
                mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][1][j]))
                if min(mendel_distance) > 0 :
                    no_mendel_record_ls[i][0][j] = 1
    
    write_vcf_lines = []
    denovo_write_vcf_lines = []
    no_mendel_write_vcf_lines = []
    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            if child_sv_ls[i][0][j] == 0 and child_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            else :
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            vcf_line.append("GT")
            vcf_line.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
            if child_sv_ls[i][0][j] + child_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
            if denovo_record_ls[i][0][j] == 1 :
                denovo_write_vcf_lines.append("\t".join(vcf_line))
            if no_mendel_record_ls[i][0][j] == 1 :
                no_mendel_write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    denovo_write_vcf_lines = header_vcf_lines + denovo_write_vcf_lines
    no_mendel_write_vcf_lines = header_vcf_lines + no_mendel_write_vcf_lines
    with open(answer_vcf_file_1,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_vcf_file,'w') as f:
        for line in denovo_write_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(no_mendel_vcf_file,'w') as f:
        for line in no_mendel_write_vcf_lines :
            f.write(line+"\n")
        f.close()

    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    write_vcf_lines = []
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if father_sv_ls[i][0][j] == 0 and father_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            else :
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            vcf_line.append("GT")
            vcf_line.append(str(father_sv_ls[i][0][j])+"|"+str(father_sv_ls[i][1][j]))
            if father_sv_ls[i][0][j] + father_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_2,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    write_vcf_lines = []
    for i in range(len(mother_sv_ls)) :
        for j in range(len(mother_sv_ls[i][0])) :
            if mother_sv_ls[i][0][j] == 0 and mother_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            else :
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            vcf_line.append("GT")
            vcf_line.append(str(mother_sv_ls[i][0][j])+"|"+str(mother_sv_ls[i][1][j]))
            if mother_sv_ls[i][0][j] + mother_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_3,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    write_vcf_lines = []
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                continue
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            else :
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            vcf_line.append("GT")
            vcf_line.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
            vcf_line.append(str(father_sv_ls[i][0][j])+"|"+str(father_sv_ls[i][1][j]))
            vcf_line.append(str(mother_sv_ls[i][0][j])+"|"+str(mother_sv_ls[i][1][j]))
            write_vcf_lines.append("\t".join(vcf_line))
    header_lines_old = copy.deepcopy(header_vcf_lines)
    header_vcf_lines[-1] = header_vcf_lines[-1] + "\tFather" + "\tMother"
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_fam,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()
    header_vcf_lines = header_lines_old

    with open(hap_choice_txt,'w') as f:
        for i in range(len(mendel_hap_choices)) :
            f.write(str(chr_names_ls[i])+"\t"+str(mendel_hap_choices[i][0])+"\t"+str(mendel_hap_choices[i][1])+"\n")
        f.close()
    
    with open(denovo_change_txt,'w') as f:
        for i in range(len(denovo_record_ls)) :
            for j in range(len(denovo_record_ls[i][0])) :
                if denovo_record_ls[i][0][j] == 0 :
                    continue
                if chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                else :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                bed_lins_ls.append("0")
                bed_lins_ls.append(str(child_sv_ls_denovo_pre[i][0][j])+"|"+str(child_sv_ls_denovo_pre[i][1][j]))
                bed_lins_ls.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_1,'w') as f:
        for i in range(len(child_sv_ls)) :
            for j in range(len(child_sv_ls[i][0])) :
                if child_sv_ls[i][0][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_2,'w') as f:
        for i in range(len(child_sv_ls)) :
            for j in range(len(child_sv_ls[i][0])) :
                if child_sv_ls[i][1][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_3,'w') as f:
        for i in range(len(father_sv_ls)) :
            for j in range(len(father_sv_ls[i][0])) :
                if father_sv_ls[i][0][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_4,'w') as f:
        for i in range(len(father_sv_ls)) :
            for j in range(len(father_sv_ls[i][0])) :
                if father_sv_ls[i][1][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_5,'w') as f:
        for i in range(len(mother_sv_ls)) :
            for j in range(len(mother_sv_ls[i][0])) :
                if mother_sv_ls[i][0][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_6,'w') as f:
        for i in range(len(mother_sv_ls)) :
            for j in range(len(mother_sv_ls[i][0])) :
                if mother_sv_ls[i][1][j] == 0 :
                    continue
                if not with_snp and chrs_sv_ls[i][j][1] not in ["INS","DEL"] :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if "INS" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                    bed_lins_ls.append("None")
                elif "DEL" in chrs_sv_ls[i][j][1] :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                    bed_lins_ls.append("None")
                else :
                    bed_lins_ls.append("SNP")
                    bed_lins_ls.append(".")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    if not with_snp :
        return

    write_vcf_lines = []
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            if child_sv_ls[i][0][j] == 0 and child_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "SNP-INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "SNP-DEL" :
                vcf_line.append("None")
                vcf_line.append("None")
            else :
                vcf_line.append("None")
                vcf_line.append(chrs_sv_ls[i][j][4])
            vcf_line.append(".")
            vcf_line.append("PASS")
            vcf_line.append("")
            vcf_line.append("GT")
            vcf_line.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
            write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(snp_vcf_file_1,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    write_vcf_lines = []
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if father_sv_ls[i][0][j] == 0 and father_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "SNP-INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "SNP-DEL" :
                vcf_line.append("None")
                vcf_line.append("None")
            else :
                vcf_line.append("None")
                vcf_line.append(chrs_sv_ls[i][j][4])
            vcf_line.append(".")
            vcf_line.append("PASS")
            vcf_line.append("")
            vcf_line.append("GT")
            vcf_line.append(str(father_sv_ls[i][0][j])+"|"+str(father_sv_ls[i][1][j]))
            write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(snp_vcf_file_2,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    write_vcf_lines = []
    for i in range(len(mother_sv_ls)) :
        for j in range(len(mother_sv_ls[i][0])) :
            if mother_sv_ls[i][0][j] == 0 and mother_sv_ls[i][1][j] == 0 :
                continue
            if chrs_sv_ls[i][j][1] in ["INS","DEL"] :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "SNP-INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "SNP-DEL" :
                vcf_line.append("None")
                vcf_line.append("None")
            else :
                vcf_line.append("None")
                vcf_line.append(chrs_sv_ls[i][j][4])
            vcf_line.append(".")
            vcf_line.append("PASS")
            vcf_line.append("")
            vcf_line.append("GT")
            vcf_line.append(str(mother_sv_ls[i][0][j])+"|"+str(mother_sv_ls[i][1][j]))
            write_vcf_lines.append("\t".join(vcf_line))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(snp_vcf_file_3,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

def make_simu_data_sv(file_base,only_indel,only_invduptra) :
    answer_vcf_file_1 = file_base + "/answer.1.vcf"
    answer_vcf_file_2 = file_base + "/answer.2.vcf"
    answer_vcf_file_3 = file_base + "/answer.3.vcf"
    answer_vcf_file_fam = file_base + "/answer.fam.vcf"
    denovo_vcf_file = file_base + "/denovo.vcf"
    denovo_1_file = file_base + "/denovo.1.vcf"
    denovo_2_file = file_base + "/denovo.2.vcf"
    denovo_3_file = file_base + "/denovo.3.vcf"
    denovo_4_file = file_base + "/denovo.4.vcf"
    no_mendel_vcf_file = file_base + "/no_mendel.vcf"
    hap_choice_txt = file_base + "/hap_choice.txt"
    denovo_change_txt = file_base + "/denovo_change.txt"
    fam_bed_1 = file_base + "/HACK.h1.bed"
    fam_bed_2 = file_base + "/HACK.h2.bed"
    fam_bed_3 = file_base + "/HACK.h3.bed"
    fam_bed_4 = file_base + "/HACK.h4.bed"
    fam_bed_5 = file_base + "/HACK.h5.bed"
    fam_bed_6 = file_base + "/HACK.h6.bed"
    
    ans_vcf_file = "/home/user/lixin/cutesv_trio/answer_vcf/HG002_T2T/GRCh37_HG2-T2TQ100-V1.1_stvar.sv.INDEL.vcf.gz"
    _,indel_vcf_lines_total,ins_len_ls,del_len_ls = read_indel_answer_vcf(ans_vcf_file,50)
    header_vcf_lines = make_header()
    dna_ls = ["A","C","G","T"]
    indel_vcf_lines = []
    for i in indel_vcf_lines_total :
        if random.random() <= 0.9 :
            indel_vcf_lines.append(i)

    sv_num = 0
    chrs_sv_ls = [[] for i in chr_names_ls]
    if not only_invduptra :
        # 同分布长度+原始位点(20000范围内扰动)+随机序列
        for sv_line in indel_vcf_lines :
            sv_split_ls = sv_line.split("\t")
            if "SVTYPE=INS" in sv_split_ls[7] :
                sv_len = random.choice(ins_len_ls)
                if int(sv_split_ls[1]) + sv_len > chrs_len_ls[chr_names_ls.index(sv_split_ls[0])] :
                    continue
                sv_str = ""
                for i in range(sv_len) :
                    sv_str = sv_str + dna_ls[random.randint(0,3)]
                pos_random_move = random.randint(-5000,5000)
                chrs_sv_ls[chr_names_ls.index(sv_split_ls[0])].append([sv_split_ls[0],"INS",int(sv_split_ls[1])+pos_random_move,int(sv_split_ls[1])+pos_random_move+1,sv_str])
            if "SVTYPE=DEL" in sv_split_ls[7] :
                sv_len = random.choice(del_len_ls)
                if int(sv_split_ls[1]) + sv_len > chrs_len_ls[chr_names_ls.index(sv_split_ls[0])] :
                    print(int(sv_split_ls[1]) + sv_len,"/",chrs_len_ls[chr_names_ls.index(sv_split_ls[0])])
                    continue
                pos_random_move = random.randint(-5000,5000)
                chrs_sv_ls[chr_names_ls.index(sv_split_ls[0])].append([sv_split_ls[0],"DEL",int(sv_split_ls[1])+pos_random_move,int(sv_split_ls[1])+pos_random_move+len(sv_split_ls[3]),"None"])
    
    if not only_indel :
        
        # 插入INV
        INV_num = 80*3
        # 中位数1000，小于50的概率小于0.01的正态分布
        INV_len_ls = np.random.normal(loc=1000, scale=407.75, size=INV_num)
        INV_len_ls = [max(50,int(x)) for x in INV_len_ls]
        for i in range(INV_num) :
            j = random.randint(0,len(chr_names_ls)-1)
            break_point = random.randint(1,chrs_len_ls[j])
            if break_point+INV_len_ls[i] > chrs_len_ls[j] :
                continue
            chrs_sv_ls[j].append([chr_names_ls[j],"INV",break_point,break_point+INV_len_ls[i],"None"])
        
        
        # 插入DUP
        DUP_num = 60*5
        # 中位数1000，小于50的概率小于0.01的正态分布
        DUP_len_ls = np.random.normal(loc=15000, scale=1000, size=DUP_num)
        DUP_len_ls = [max(50,int(x)) for x in DUP_len_ls]
        for i in range(DUP_num) :
            j = random.randint(0,len(chr_names_ls)-1)
            break_point = random.randint(1,chrs_len_ls[j])
            if break_point+DUP_len_ls[i] > chrs_len_ls[j] :
                continue
            chrs_sv_ls[j].append([chr_names_ls[j],"DUP",break_point,break_point+DUP_len_ls[i],"2"])
        
        
        # 插入TRA
        TRA_num = 430
        # 中位数1000，小于50的概率小于0.01的正态分布
        TRA_len_ls = np.random.normal(loc=15000, scale=5000, size=TRA_num)
        TRA_len_ls = [max(50,int(x)) for x in TRA_len_ls]
        for i in range(TRA_num) :
            chr1 = random.randint(0,len(chr_names_ls)-1)
            while(True) :
                chr2 = random.randint(0,len(chr_names_ls)-1)
                if chr1 != chr2 :
                    break
            chr1_break_point = random.randint(1,chrs_len_ls[chr1])
            chr2_break_point = random.randint(1,chrs_len_ls[chr2])
            if chr1_break_point+TRA_len_ls[i] > chrs_len_ls[chr1] or chr2_break_point+TRA_len_ls[i] > chrs_len_ls[chr2] :
                continue
            sv_str = "h1:"+chr_names_ls[chr2]+":"+str(chr2_break_point)+":"
            direct_ls = ["forward:forward","forward:reverse","reverse:forward","reverse:reverse"]
            sv_str = sv_str + random.choice(direct_ls)
            chrs_sv_ls[chr1].append([chr_names_ls[chr1],"BND",chr1_break_point,chr1_break_point+TRA_len_ls[i],sv_str])
        
    for i in range(len(chrs_sv_ls)) :
        chrs_sv_ls[i] = sorted(chrs_sv_ls[i], key=lambda x: x[2])

    chrs_sv_ls = remove_overlap_sv(chrs_sv_ls)
    for i in chrs_sv_ls :
        sv_num += len(i)
    print("sv总数:"+str(sv_num))

    father_sv_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    mother_sv_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if random.random() <= 1 :
                father_sv_ls[i][0][j] = 1
                father_sv_ls[i][1][j] = 1
                mother_sv_ls[i][0][j] = 1
                mother_sv_ls[i][1][j] = 1
            else :
                if random.random() <= 0.5 :
                    father_sv_ls[i][0][j] = 1
                    father_sv_ls[i][1][j] = 1
                else :
                    mother_sv_ls[i][0][j] = 1
                    mother_sv_ls[i][1][j] = 1

    threshold_1 = 0.1
    threshold_2 = 0.3
    threshold_3 = 0.6
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if father_sv_ls[i][0][j] == 0 :
                continue
            random_value = random.random()
            if random_value <= threshold_1 :
                father_sv_ls[i][0][j] = 0
                father_sv_ls[i][1][j] = 0
            elif random_value > threshold_1 and random_value <= threshold_2 :
                father_sv_ls[i][0][j] = 0
            elif random_value > threshold_2 and random_value <= threshold_3 :
                father_sv_ls[i][1][j] = 0
            else :
                pass
    for i in range(len(mother_sv_ls)) :
        for j in range(len(mother_sv_ls[i][0])) :
            if mother_sv_ls[i][0][j] == 0 :
                continue
            random_value = random.random()
            if random_value <= threshold_1 :
                mother_sv_ls[i][0][j] = 0
                mother_sv_ls[i][1][j] = 0
            elif random_value > threshold_1 and random_value <= threshold_2 :
                mother_sv_ls[i][0][j] = 0
            elif random_value > threshold_2 and random_value <= threshold_3 :
                mother_sv_ls[i][1][j] = 0
            else :
                pass
    mendel_hap_choices = []
    child_sv_ls = [[[],[]] for chrs in chrs_sv_ls]
    for i in range(len(child_sv_ls)) :
        mendel_choice = []
        mendel_choice.append(random.randint(0,1))
        mendel_choice.append(random.randint(0,1))
        mendel_hap_choices.append(mendel_choice)
        child_sv_ls[i][0] = copy.deepcopy(father_sv_ls[i][mendel_choice[0]])
        child_sv_ls[i][1] = copy.deepcopy(mother_sv_ls[i][mendel_choice[1]])
    
    if not only_invduptra :
        child_sv_ls_denovo_pre = copy.deepcopy(child_sv_ls)
        denovo_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
        
        # 以下加入固定数量突变，随机选择合适位置
        denovo_1_num = 100
        denovo_2_num = 210
        denovo_3_num = 170
        denovo_4_num = 90
        for i in range(denovo_1_num) :
            while(True) :
                chr_index = random.randint(0,len(child_sv_ls)-1)
                sv_index = random.randint(0,len(child_sv_ls[chr_index][0])-1)
                if denovo_record_ls[chr_index][0][sv_index] != 0 or chrs_sv_ls[chr_index][sv_index][1] not in ["INS","DEL"] :
                    continue
                if child_sv_ls[chr_index][0][sv_index] + child_sv_ls[chr_index][1][sv_index] == 0 and father_sv_ls[chr_index][0][sv_index] + father_sv_ls[chr_index][1][sv_index] == 0 and mother_sv_ls[chr_index][0][sv_index] + mother_sv_ls[chr_index][1][sv_index] == 0 :
                    denovo_record_ls[chr_index][0][sv_index] = 1
                    if random.random() <= 0.5 :
                        child_sv_ls[chr_index][0][sv_index] = 1
                    else :
                        child_sv_ls[chr_index][1][sv_index] = 1
                    break
        for i in range(denovo_2_num) :
            while(True) :
                chr_index = random.randint(0,len(child_sv_ls)-1)
                sv_index = random.randint(0,len(child_sv_ls[chr_index][0])-1)
                if denovo_record_ls[chr_index][0][sv_index] != 0 or chrs_sv_ls[chr_index][sv_index][1] not in ["INS","DEL"] :
                    continue
                if child_sv_ls[chr_index][0][sv_index] + child_sv_ls[chr_index][1][sv_index] == 1 and (father_sv_ls[chr_index][0][sv_index] + father_sv_ls[chr_index][1][sv_index] == 0 or mother_sv_ls[chr_index][0][sv_index] + mother_sv_ls[chr_index][1][sv_index] == 0) :
                    denovo_record_ls[chr_index][0][sv_index] = 2
                    if child_sv_ls[chr_index][0][sv_index] == 0 :
                        child_sv_ls[chr_index][0][sv_index] = 1
                    else :
                        child_sv_ls[chr_index][1][sv_index] = 1
                    break
        for i in range(denovo_3_num) :
            while(True) :
                chr_index = random.randint(0,len(child_sv_ls)-1)
                sv_index = random.randint(0,len(child_sv_ls[chr_index][0])-1)
                if denovo_record_ls[chr_index][0][sv_index] != 0 or chrs_sv_ls[chr_index][sv_index][1] not in ["INS","DEL"] :
                    continue
                if child_sv_ls[chr_index][0][sv_index] + child_sv_ls[chr_index][1][sv_index] == 2 and father_sv_ls[chr_index][0][sv_index] + father_sv_ls[chr_index][1][sv_index] == 2 and mother_sv_ls[chr_index][0][sv_index] + mother_sv_ls[chr_index][1][sv_index] == 2 :
                    denovo_record_ls[chr_index][0][sv_index] = 3
                    if random.random() <= 0.5 :
                        child_sv_ls[chr_index][0][sv_index] = 0
                    else :
                        child_sv_ls[chr_index][1][sv_index] = 0
                    break
        for i in range(denovo_4_num) :
            while(True) :
                chr_index = random.randint(0,len(child_sv_ls)-1)
                sv_index = random.randint(0,len(child_sv_ls[chr_index][0])-1)
                if denovo_record_ls[chr_index][0][sv_index] != 0 or chrs_sv_ls[chr_index][sv_index][1] not in ["INS","DEL"] :
                    continue
                if child_sv_ls[chr_index][0][sv_index] + child_sv_ls[chr_index][1][sv_index] == 1 and (father_sv_ls[chr_index][0][sv_index] + father_sv_ls[chr_index][1][sv_index] == 2 or mother_sv_ls[chr_index][0][sv_index] + mother_sv_ls[chr_index][1][sv_index] == 2) :
                    denovo_record_ls[chr_index][0][sv_index] = 4
                    if child_sv_ls[chr_index][0][sv_index] == 1 :
                        child_sv_ls[chr_index][0][sv_index] = 0
                    else :
                        child_sv_ls[chr_index][1][sv_index] = 0
                    break
        
        no_mendel_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
        for i in range(len(denovo_record_ls)) :
            for j in range(len(denovo_record_ls[i][0])) :
                if denovo_record_ls[i][0][j] != 0 :
                    mendel_distance = []
                    mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][0][j]))
                    mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][1][j]))
                    mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][0][j]))
                    mendel_distance.append(abs(child_sv_ls[i][0][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][1][j]-mother_sv_ls[i][1][j]))
                    mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][0][j]))
                    mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][0][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][1][j]))
                    mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][0][j]))
                    mendel_distance.append(abs(child_sv_ls[i][1][j]-father_sv_ls[i][1][j])+abs(child_sv_ls[i][0][j]-mother_sv_ls[i][1][j]))
                    if min(mendel_distance) > 0 :
                        no_mendel_record_ls[i][0][j] = 1
    else :
        denovo_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
        no_mendel_record_ls = [[[0 for i in chrs],[0 for i in chrs]] for chrs in chrs_sv_ls]
    
    write_vcf_lines = []
    denovo_write_vcf_lines = []
    denovo_1_vcf_lines = []
    denovo_2_vcf_lines = []
    denovo_3_vcf_lines = []
    denovo_4_vcf_lines = []
    no_mendel_write_vcf_lines = []
    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    INS_num = 0
    DEL_num = 0
    INV_num = 0
    DUP_num = 0
    TRA_num = 0
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            if child_sv_ls[i][0][j] == 0 and child_sv_ls[i][1][j] == 0 :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "BND" :
                vcf_line.append("N")
                vcf_line.append("["+str(chrs_sv_ls[i][j][4].split(":")[1])+":"+str(chrs_sv_ls[i][j][4].split(":")[2])+"[N")
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                INS_num += 1
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            elif chrs_sv_ls[i][j][1] == "DEL" :
                DEL_num += 1
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "INV" :
                INV_num += 1
                vcf_line.append("SVTYPE=INV;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "DUP" :
                DUP_num += 1
                vcf_line.append("SVTYPE=DUP;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            else :
                TRA_num += 1
                vcf_line.append("SVTYPE=BND;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3])+";CHR2="+chrs_sv_ls[i][j][4].split(":")[1])
            vcf_line.append("GT")
            vcf_line.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
            if child_sv_ls[i][0][j] + child_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
            if denovo_record_ls[i][0][j] == 1 :
                denovo_write_vcf_lines.append("\t".join(vcf_line))
                denovo_1_vcf_lines.append("\t".join(vcf_line))
            if denovo_record_ls[i][0][j] == 2 :
                denovo_write_vcf_lines.append("\t".join(vcf_line))
                denovo_2_vcf_lines.append("\t".join(vcf_line))
            if denovo_record_ls[i][0][j] == 3 :
                denovo_write_vcf_lines.append("\t".join(vcf_line))
                denovo_3_vcf_lines.append("\t".join(vcf_line))
            if denovo_record_ls[i][0][j] == 4 :
                denovo_write_vcf_lines.append("\t".join(vcf_line))
                denovo_4_vcf_lines.append("\t".join(vcf_line))
            if no_mendel_record_ls[i][0][j] == 1 :
                no_mendel_write_vcf_lines.append("\t".join(vcf_line))
    print("子代INS数量:"+str(INS_num))
    print("子代DEL数量:"+str(DEL_num))
    print("子代INV数量:"+str(INV_num))
    print("子代DUP数量:"+str(DUP_num))
    print("子代TRA数量:"+str(TRA_num))
    print("子代杂合比例:"+str(heterozygous_num/sv_num))
    print("子代纯和比例:"+str(homozygous_num/sv_num))
    print("class 1 denovo数量:"+str(len(denovo_1_vcf_lines)))
    print("class 2 denovo数量:"+str(len(denovo_2_vcf_lines)))
    print("class 3 denovo数量:"+str(len(denovo_3_vcf_lines)))
    print("class 4 denovo数量:"+str(len(denovo_4_vcf_lines)))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    denovo_write_vcf_lines = header_vcf_lines + denovo_write_vcf_lines
    denovo_1_vcf_lines = header_vcf_lines + denovo_1_vcf_lines
    denovo_2_vcf_lines = header_vcf_lines + denovo_2_vcf_lines
    denovo_3_vcf_lines = header_vcf_lines + denovo_3_vcf_lines
    denovo_4_vcf_lines = header_vcf_lines + denovo_4_vcf_lines
    no_mendel_write_vcf_lines = header_vcf_lines + no_mendel_write_vcf_lines
    with open(answer_vcf_file_1,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_vcf_file,'w') as f:
        for line in denovo_write_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_1_file,'w') as f:
        for line in denovo_1_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_2_file,'w') as f:
        for line in denovo_2_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_3_file,'w') as f:
        for line in denovo_3_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(denovo_4_file,'w') as f:
        for line in denovo_4_vcf_lines :
            f.write(line+"\n")
        f.close()
    with open(no_mendel_vcf_file,'w') as f:
        for line in no_mendel_write_vcf_lines :
            f.write(line+"\n")
        f.close()

    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    INS_num = 0
    DEL_num = 0
    INV_num = 0
    DUP_num = 0
    TRA_num = 0
    write_vcf_lines = []
    for i in range(len(father_sv_ls)) :
        for j in range(len(father_sv_ls[i][0])) :
            if father_sv_ls[i][0][j] == 0 and father_sv_ls[i][1][j] == 0 :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "BND" :
                vcf_line.append("N")
                vcf_line.append("["+str(chrs_sv_ls[i][j][4].split(":")[1])+":"+str(chrs_sv_ls[i][j][4].split(":")[2])+"[N")
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                INS_num += 1
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            elif chrs_sv_ls[i][j][1] == "DEL" :
                DEL_num += 1
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "INV" :
                INV_num += 1
                vcf_line.append("SVTYPE=INV;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "DUP" :
                DUP_num += 1
                vcf_line.append("SVTYPE=DUP;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            else :
                TRA_num += 1
                vcf_line.append("SVTYPE=BND;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3])+";CHR2="+chrs_sv_ls[i][j][4].split(":")[1])
            vcf_line.append("GT")
            vcf_line.append(str(father_sv_ls[i][0][j])+"|"+str(father_sv_ls[i][1][j]))
            if father_sv_ls[i][0][j] + father_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
    print("子代INS数量:"+str(INS_num))
    print("子代DEL数量:"+str(DEL_num))
    print("子代INV数量:"+str(INV_num))
    print("子代DUP数量:"+str(DUP_num))
    print("子代TRA数量:"+str(TRA_num))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_2,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    homozygous_num = 0
    heterozygous_num = 0
    sv_num = 0
    INS_num = 0
    DEL_num = 0
    INV_num = 0
    DUP_num = 0
    TRA_num = 0
    write_vcf_lines = []
    for i in range(len(mother_sv_ls)) :
        for j in range(len(mother_sv_ls[i][0])) :
            if mother_sv_ls[i][0][j] == 0 and mother_sv_ls[i][1][j] == 0 :
                continue
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "BND" :
                vcf_line.append("N")
                vcf_line.append("["+str(chrs_sv_ls[i][j][4].split(":")[1])+":"+str(chrs_sv_ls[i][j][4].split(":")[2])+"[N")
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                INS_num += 1
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            elif chrs_sv_ls[i][j][1] == "DEL" :
                DEL_num += 1
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "INV" :
                INV_num += 1
                vcf_line.append("SVTYPE=INV;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "DUP" :
                DUP_num += 1
                vcf_line.append("SVTYPE=DUP;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            else :
                TRA_num += 1
                vcf_line.append("SVTYPE=BND;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3])+";CHR2="+chrs_sv_ls[i][j][4].split(":")[1])
            vcf_line.append("GT")
            vcf_line.append(str(mother_sv_ls[i][0][j])+"|"+str(mother_sv_ls[i][1][j]))
            if mother_sv_ls[i][0][j] + mother_sv_ls[i][1][j] == 1 :
                heterozygous_num += 1
                sv_num += 1
            else :
                homozygous_num += 1
                sv_num += 1
            write_vcf_lines.append("\t".join(vcf_line))
    print("子代INS数量:"+str(INS_num))
    print("子代DEL数量:"+str(DEL_num))
    print("子代INV数量:"+str(INV_num))
    print("子代DUP数量:"+str(DUP_num))
    print("子代TRA数量:"+str(TRA_num))
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_3,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

    write_vcf_lines = []
    for i in range(len(child_sv_ls)) :
        for j in range(len(child_sv_ls[i][0])) :
            vcf_line = []
            vcf_line.append(chrs_sv_ls[i][j][0])
            vcf_line.append(str(chrs_sv_ls[i][j][2]))
            vcf_line.append(chrs_sv_ls[i][j][0]+"."+chrs_sv_ls[i][j][1]+"_"+str(chrs_sv_ls[i][j][2]))
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append(chrs_sv_ls[i][j][4][0])
                vcf_line.append(chrs_sv_ls[i][j][4])
            elif chrs_sv_ls[i][j][1] == "BND" :
                vcf_line.append("N")
                vcf_line.append("["+str(chrs_sv_ls[i][j][4].split(":")[1])+":"+str(chrs_sv_ls[i][j][4].split(":")[2])+"[N")
            else :
                vcf_line.append("None")
                vcf_line.append("None")
            vcf_line.append(".")
            vcf_line.append("PASS")
            if chrs_sv_ls[i][j][1] == "INS" :
                vcf_line.append("SVTYPE=INS;SVLEN="+str(len(chrs_sv_ls[i][j][4]))+";END="+str(chrs_sv_ls[i][j][2]))
            elif chrs_sv_ls[i][j][1] == "DEL" :
                vcf_line.append("SVTYPE=DEL;SVLEN="+str(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3])+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "INV" :
                vcf_line.append("SVTYPE=INV;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            elif chrs_sv_ls[i][j][1] == "DUP" :
                vcf_line.append("SVTYPE=DUP;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3]))
            else :
                vcf_line.append("SVTYPE=BND;SVLEN="+str(abs(chrs_sv_ls[i][j][2]-chrs_sv_ls[i][j][3]))+";END="+str(chrs_sv_ls[i][j][3])+";CHR2="+chrs_sv_ls[i][j][4].split(":")[1])
            vcf_line.append("GT")
            vcf_line.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
            vcf_line.append(str(father_sv_ls[i][0][j])+"|"+str(father_sv_ls[i][1][j]))
            vcf_line.append(str(mother_sv_ls[i][0][j])+"|"+str(mother_sv_ls[i][1][j]))
            write_vcf_lines.append("\t".join(vcf_line))
    header_lines_old = copy.deepcopy(header_vcf_lines)
    header_vcf_lines[-1] = header_vcf_lines[-1] + "\tFather" + "\tMother"
    write_vcf_lines = header_vcf_lines + write_vcf_lines
    with open(answer_vcf_file_fam,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()
    header_vcf_lines = header_lines_old

    with open(hap_choice_txt,'w') as f:
        for i in range(len(mendel_hap_choices)) :
            f.write(str(chr_names_ls[i])+"\t"+str(mendel_hap_choices[i][0])+"\t"+str(mendel_hap_choices[i][1])+"\n")
        f.close()
    
    with open(denovo_change_txt,'w') as f:
        for i in range(len(denovo_record_ls)) :
            for j in range(len(denovo_record_ls[i][0])) :
                if denovo_record_ls[i][0][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append(str(child_sv_ls_denovo_pre[i][0][j])+"|"+str(child_sv_ls_denovo_pre[i][1][j]))
                bed_lins_ls.append(str(child_sv_ls[i][0][j])+"|"+str(child_sv_ls[i][1][j]))
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_1,'w') as f:
        for i in range(len(child_sv_ls)) :
            for j in range(len(child_sv_ls[i][0])) :
                if child_sv_ls[i][0][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_2,'w') as f:
        for i in range(len(child_sv_ls)) :
            for j in range(len(child_sv_ls[i][0])) :
                if child_sv_ls[i][1][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_3,'w') as f:
        for i in range(len(father_sv_ls)) :
            for j in range(len(father_sv_ls[i][0])) :
                if father_sv_ls[i][0][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_4,'w') as f:
        for i in range(len(father_sv_ls)) :
            for j in range(len(father_sv_ls[i][0])) :
                if father_sv_ls[i][1][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_5,'w') as f:
        for i in range(len(mother_sv_ls)) :
            for j in range(len(mother_sv_ls[i][0])) :
                if mother_sv_ls[i][0][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

    with open(fam_bed_6,'w') as f:
        for i in range(len(mother_sv_ls)) :
            for j in range(len(mother_sv_ls[i][0])) :
                if mother_sv_ls[i][1][j] == 0 :
                    continue
                bed_lins_ls = []
                bed_lins_ls.append(chrs_sv_ls[i][j][0])
                bed_lins_ls.append(str(chrs_sv_ls[i][j][2]))
                bed_lins_ls.append(str(chrs_sv_ls[i][j][3]))
                if chrs_sv_ls[i][j][1] == "INS" :
                    bed_lins_ls.append("insertion")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                elif chrs_sv_ls[i][j][1] == "DEL" :
                    bed_lins_ls.append("deletion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "INV" :
                    bed_lins_ls.append("inversion")
                    bed_lins_ls.append("None")
                elif chrs_sv_ls[i][j][1] == "DUP" :
                    bed_lins_ls.append("tandem duplication")
                    bed_lins_ls.append("2")
                else :
                    bed_lins_ls.append("reciprocal translocation")
                    bed_lins_ls.append(chrs_sv_ls[i][j][4])
                bed_lins_ls.append("0")
                bed_lins_ls.append("1|1")
                f.write("\t".join(bed_lins_ls)+"\n")
        f.close()

def remove_overlap_sv(chrs_sv_ls):
    chrs_sv_new = []
    for chr_svs in chrs_sv_ls :
        if len(chr_svs) == 0 :
            chrs_sv_new.append(chr_svs)
            continue
        chr_new_svs = []
        last_start = int(chr_svs[0][2])
        last_end = int(chr_svs[0][3])
        chr_new_svs.append(chr_svs[0])
        for sv in chr_svs[1:] :
            if int(sv[2]) > last_end :
                chr_new_svs.append(sv)
                last_start = int(sv[2])
                last_end = int(sv[3])
        chrs_sv_new.append(chr_new_svs)
    return chrs_sv_new

def main(argv):
    # 先使用VISOR插入SNP形成六个fa，再使用六个fa插入SV
    make_simu_data_sv("xxx",True,False)

    # 额外插入SNP变异
    make_simu_data_with_snp("xxx",True)
    
if __name__ == "__main__":
    main(sys.argv[1:])