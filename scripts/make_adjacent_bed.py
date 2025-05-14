# -*- coding: utf-8 -*-

import sys
import gzip
import csv
import math
import numpy as np

def make_bed(sv_vcf_file,bed_file,Expansion_ratio,Expansion_len) :
    line_list = []
    fh = gzip.open(sv_vcf_file, 'rt') if sv_vcf_file.endswith('.gz') else open(sv_vcf_file, 'rt')
    for line in fh:
        line = line.strip()
        if line.startswith('#'):
            continue
        line_split_ls = line.split("\t")
        if line_split_ls[6] != "PASS" :
            continue
        if "SVLEN=" not in line_split_ls[7] :
            continue
        sv_len = max(len(line_split_ls[4]),len(line_split_ls[3]),abs(int(line_split_ls[7].split("SVLEN=")[1].split(";")[0])))
        start_pos = int(max(0,int(line_split_ls[1]) - Expansion_ratio * sv_len - Expansion_len))
        end_pos = int(int(line_split_ls[1]) + sv_len + Expansion_ratio * sv_len + Expansion_len)
        line_list.append([line_split_ls[0],start_pos,end_pos])
    fh.close()
    bed_list = [line_list[0]]
    
    for i in range(1,len(line_list)) :
        if bed_list[-1][0] == line_list[i][0] and bed_list[-1][2] >= line_list[i][1] :
            bed_list[-1] = [bed_list[-1][0],bed_list[-1][1],line_list[i][2]]
        else :
            bed_list.append(line_list[i])
    for i in range(len(bed_list)) :
        bed_list[i] = [str(s) for s in bed_list[i]]
    
    with open(bed_file,'w') as f:
        for ls in bed_list :
            f.write("\t".join(ls)+"\n")
        f.close()

def main(argv):
    make_bed(argv[0],argv[1],argv[2],argv[3])

if __name__ == "__main__":
    main(sys.argv[1:])