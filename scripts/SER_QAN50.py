# -*- coding: utf-8 -*-

import sys
import gzip
import csv
import math
import numpy as np
     
def cal_ser(file_base,phase_type) :
    base_gt_ls = []
    base_line_ls = []
    file_name = file_base + "/tp-base.vcf.gz"
    fh = gzip.open(file_name, 'rt') if file_name.endswith('.gz') else open(file_name, 'rt')
    for line in fh:
        if line.startswith('#'):
            continue
        line_split_ls = line.split('\t')
        gt_index = line_split_ls[-2].split(":").index("GT")
        ph_index = line_split_ls[-2].split(":").index("GT")
        base_gt_ls.append([line_split_ls[-1].split(":")[gt_index],line_split_ls[-1].split(":")[ph_index],line_split_ls[0]])
        base_line_ls.append(line)
    fh.close()
    comp_gt_ls = []
    comp_line_ls = []
    file_name = file_base + "/tp-comp.vcf.gz"
    fh = gzip.open(file_name, 'rt') if file_name.endswith('.gz') else open(file_name, 'rt')
    for line in fh:
        if line.startswith('#'):
            continue
        line_split_ls = line.split('\t')
        gt_index = line_split_ls[-2].split(":").index("GT")
        if phase_type == "trio" :
            ph_index = line_split_ls[-2].split(":").index("HP_GT")
        else :
            ph_index = line_split_ls[-2].split(":").index("GT")
        if "SVLEN=" in line_split_ls[7] :
            comp_gt_ls.append([line_split_ls[-1].split(":")[gt_index],line_split_ls[-1].split(":")[ph_index],line_split_ls[0],int(line_split_ls[1]),int(line_split_ls[1])+max(len(line_split_ls[4]),len(line_split_ls[3]),abs(int(line_split_ls[7].split("SVLEN=")[1].split(";")[0])))])
        else :
            comp_gt_ls.append([line_split_ls[-1].split(":")[gt_index],line_split_ls[-1].split(":")[ph_index],line_split_ls[0],int(line_split_ls[1]),int(line_split_ls[1])+max(len(line_split_ls[4]),len(line_split_ls[3]))])
        comp_line_ls.append(line)
    fh.close()
    total_heter_num = 0
    phased_heter_num = 0
    wrong_num = 0
    right_num = 0
    switch_num = 0
    switch_flag = -1
    N50_ls = []
    N50_records = []
    N50_flag = -1
    chr = "1"
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    for i in range(len(base_gt_ls)) :
        if base_gt_ls[i][2] not in chr_ls :
            continue
        if base_gt_ls[i][2] != chr :
            chr = base_gt_ls[i][2]
            switch_flag = -1
            if N50_flag == 1 :
                N50_ls.append(N50_records[-1][2]*(N50_records[-1][1]-N50_records[-1][0]))
                N50_flag = -1
        if len(base_gt_ls[i][0]) < 3 or len(base_gt_ls[i][1]) < 3 :
            continue
        if "." in base_gt_ls[i][0] or int(base_gt_ls[i][0][0])+int(base_gt_ls[i][0][2]) != 1 or int(comp_gt_ls[i][0][0])+int(comp_gt_ls[i][0][2]) != 1 :
            continue
        if "." in base_gt_ls[i][1] :
            continue
        if comp_gt_ls[i][1] != ".|." :
            total_heter_num += 1
        else :
            continue
        if comp_gt_ls[i][1][1] == "/" :
            switch_num += 1
            if N50_flag == 1 :
                N50_records[-1][1] = comp_gt_ls[i][4]
                N50_ls.append(N50_records[-1][2]*(N50_records[-1][1]-N50_records[-1][0]))
                N50_flag = -1
        elif base_gt_ls[i][1][0] != comp_gt_ls[i][1][0] or base_gt_ls[i][1][2] != comp_gt_ls[i][1][2] :
            phased_heter_num += 1
            wrong_num += 1
            if switch_flag == -1 :
                switch_flag = 0
            elif switch_flag == 1 :
                switch_num += 1
                switch_flag = 0
            if N50_flag == 1 :
                N50_records[-1][1] = comp_gt_ls[i][4]
                N50_ls.append(N50_records[-1][2]*(N50_records[-1][1]-N50_records[-1][0]))
                N50_flag = -1
        else :
            phased_heter_num += 1
            right_num += 1
            if switch_flag == -1 :
                switch_flag = 1
            elif switch_flag == 0 :
                switch_num += 1
                switch_flag = 1
            if N50_flag == -1 :
                N50_records.append([comp_gt_ls[i][3],comp_gt_ls[i][4],1])
                N50_flag = 1
            else :
                N50_records[-1][1] = comp_gt_ls[i][4]
                N50_records[-1][2] += 1
    heterozygous_num = 0
    no_phased_num = 0
    for i in range(len(comp_gt_ls)) :
        if int(comp_gt_ls[i][0][0])+int(comp_gt_ls[i][0][2]) == 1 :
            heterozygous_num += 1
            if comp_gt_ls[i][1][1] == "/" :
                no_phased_num += 1
    N50_ls = sorted(N50_ls, reverse=True)
    N50_sum = 0
    N50_index = -1
    for i in range(len(N50_ls)) :
        N50_sum += N50_ls[i]
        if N50_sum >= sum(N50_ls)/2 :
            N50_index = i
            break
    if N50_index == -1 :
        N50_index = len(N50_ls) - 1
    if total_heter_num == 0 :
        consistency = 0
        ser = 0
        phasing_rate = 0
    else :
        consistency = right_num/total_heter_num
        ser = switch_num/total_heter_num
        phasing_rate = 1-no_phased_num/heterozygous_num
    if len(N50_ls) == 0 :
        QAN50 = 0
    else :
        QAN50 = math.log10(N50_ls[N50_index])
    return [total_heter_num,consistency,ser,phasing_rate,QAN50]

def main(argv):
    cal_ser(argv[0],argv[1])

if __name__ == "__main__":
    main(sys.argv[1:])