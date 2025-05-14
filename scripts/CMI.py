# -*- coding: utf-8 -*-

import sys
import gzip
import csv
import math
import numpy as np

def cal_block_num_seg(file_base,phase_type,seg_interval) :
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    chrs_len_ls = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
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
        #print(line_split_ls)
        if "SVLEN=" in line_split_ls[7] :
            comp_gt_ls.append([line_split_ls[-1].split(":")[gt_index],line_split_ls[-1].split(":")[ph_index],line_split_ls[0],int(line_split_ls[1]),int(line_split_ls[1])+max(len(line_split_ls[4]),len(line_split_ls[3]),abs(int(line_split_ls[7].split("SVLEN=")[1].split(";")[0])))])
        else :
            comp_gt_ls.append([line_split_ls[-1].split(":")[gt_index],line_split_ls[-1].split(":")[ph_index],line_split_ls[0],int(line_split_ls[1]),int(line_split_ls[1])+max(len(line_split_ls[4]),len(line_split_ls[3]))])
        comp_line_ls.append(line)
    fh.close()
    N50_flag = -1
    cal_chr_num = len(chr_ls)
    block_num_ls = [[0 for x in range(y//seg_interval+1)] for y in chrs_len_ls[:cal_chr_num]]
    pre_seg_no = -1
    chr = "1"
    for i in range(len(base_gt_ls)) :
        if base_gt_ls[i][2] not in chr_ls[:cal_chr_num] :
            continue
        if base_gt_ls[i][2] != chr or comp_gt_ls[i][3] // seg_interval != pre_seg_no :
            if N50_flag == 1 :
                #print(pre_seg_no)
                block_num_ls[chr_ls.index(chr)][pre_seg_no] += 1
            chr = base_gt_ls[i][2]
            pre_seg_no = comp_gt_ls[i][3] // seg_interval
            N50_flag = -1
        if len(base_gt_ls[i][0]) < 3 or len(base_gt_ls[i][1]) < 3 :
            continue
        if "." in base_gt_ls[i][0] or int(base_gt_ls[i][0][0])+int(base_gt_ls[i][0][2]) != 1 or int(comp_gt_ls[i][0][0])+int(comp_gt_ls[i][0][2]) != 1 :
            continue
        if "." in base_gt_ls[i][1] :
            continue
        if comp_gt_ls[i][1] != ".|." :
            pass
        else :
            continue
        if comp_gt_ls[i][1][1] == "/" :
            if N50_flag == 1 :
                block_num_ls[chr_ls.index(chr)][pre_seg_no] += 1
                N50_flag = -1
        elif base_gt_ls[i][1][0] != comp_gt_ls[i][1][0] or base_gt_ls[i][1][2] != comp_gt_ls[i][1][2] :
            if N50_flag == 1 :
                block_num_ls[chr_ls.index(chr)][pre_seg_no] += 1
                N50_flag = -1
        else :
            if N50_flag == -1 :
                N50_flag = 1
            else :
                pass
    bar_data = []
    for i in range(cal_chr_num) :
        for j in range(len(block_num_ls[i])) :
            bar_data.append([chr_ls[i],j+0.5,block_num_ls[i][j]])
    #print(bar_data)
    return bar_data,block_num_ls

# 展示TrioSV和longphase在不同染色体上不同间隔下的bloack的数量
def main(argv):
    block_num_tools = []
    data = [["chr","seg_no","block_num"]]
    bar_data,block_num_ls = cal_block_num_seg("30.30.30.T2T","trio",1000000)
    data = data + bar_data
    block_num_tools.append(block_num_ls)

    data = [["chr","seg_no","block_num"]]
    bar_data,block_num_ls = cal_block_num_seg("longphase","longphase",1000000)
    data = data + bar_data
    block_num_tools.append(block_num_ls)
    
    data = [["chr","seg_no","block_num"]]
    bar_data,block_num_ls = cal_block_num_seg("whatshap","whatshap",1000000)
    data = data + bar_data
    block_num_tools.append(block_num_ls)

    data = [["chr","trio-variance","trio-abnormal","trio-continuity","longphase-variance","longphase-abnormal","longphase-continuity","whatshap-variance","whatshap-abnormal","whatshap-continuity"]]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
    tools_ls = ["cuteSV-Trio","LongPhase","WhatsHap"]
    cmi_box_ls = [["Tool","Chr","CMI"]]
    for i in range(len(chr_ls)) :
        tool_line = [chr_ls[i]]
        for j in range(len(block_num_tools)) :
            variance = np.var(block_num_tools[j][i])
            abnormal = sum([1 if x>4 else 0 for x in block_num_tools[j][i] ])
            continuity = 0
            continue_flag = -1
            for k in block_num_tools[j][i] :
                if k == 0 and continue_flag == 1:
                    continuity -= 1
                    continue_flag = -1
                elif k == 1 :
                    continuity += 1
                    continue_flag = 1
                elif k > 1 :
                    continuity = continuity + 1
                    continue_flag = 1
            tool_line += [variance,abnormal,continuity]
            cmi_box_ls.append([tools_ls[j],chr_ls[i],continuity])
        data.append(tool_line)
    
    data = [["Tool","Variance"]]
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
    tools_ls = ["cuteSV-Trio","LongPhase","WhatsHap"]
    for i in range(len(chr_ls)) :
        for j in range(len(block_num_tools)) :
            variance = np.var(block_num_tools[j][i])
            data.append([tools_ls[j],variance])

if __name__ == "__main__":
    main(sys.argv[1:])