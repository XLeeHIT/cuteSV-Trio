# -*- coding: utf-8 -*-

import sys
import os
import csv
import gzip
import copy
import math
import numpy as np
from scipy.interpolate import CubicSpline
import random

child_ls = ["HG002","HG005","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877","HG01891","HG02668","HG03371","HG03453","NA19705"]
father_ls = ["HG003","HG006","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889","HG01890","HG02666","","","NA19703"]
mother_ls = ["HG004","HG007","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890","","","HG03369","HG03452",""]

child_trio_ls = ["HG002","HG005","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877"]
father_trio_ls = ["HG003","HG006","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889"]
mother_trio_ls = ["HG004","HG007","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890"]

child_duo_ls = ["HG01891","HG02668","HG03371","HG03453","NA19705"]
father_duo_ls = ["HG01890","HG02666","","","NA19703"]
mother_duo_ls = ["","","HG03369","HG03452",""]

# 只有这7个完整家庭同时出现在1kGP的样本集中
child_1kgp_trio_ls = ["HG00405","HG00408","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877"]
father_1kgp_trio_ls = ["HG00403","HG00406","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889"]
mother_1kgp_trio_ls = ["HG00404","HG00407","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890"]

cuteSVTrio_1kgp_member_ls = ["HG00512","HG00513","HG00514","HG00731","HG00732","HG00733",
                             "HG01890","HG01891","HG02666","HG02668","HG03369","HG03371",
                             "HG03452","HG03453","HG03520","HG03521","HG03522","NA12877",
                             "NA12878","NA12889","NA12890","NA12891","NA12892","NA19238",
                             "NA19239","NA19240","NA19703","NA19705","NA19834","NA19835","NA19836"]

HGSVC_1kgp_member_ls = ["HG00171","HG00268","HG00358","HG00512","HG00513","HG00514",
                        "HG00731","HG00732","HG00733","HG01457","HG01890","HG02282","HG02587",
                        "HG02666","HG02769","HG03009","HG03371","HG03452","HG03520",
                        "HG03732","NA18534","NA19036","NA19129","NA19238","NA19239",
                        "NA19240","NA19331","NA19384","NA19705","NA19836","NA19983"]

cuteSVTrio_member_ls = ["HG002","HG003","HG004","HG005","HG006","HG007",
                        "HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877",
                        "HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889",
                        "HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890",
                        "HG01891","HG02668","HG03371","HG03453","NA19705",
                        "HG01890","HG02666","NA19703","HG03369","HG03452",
                        "HG01889","HG02667","HG03370","HG03451","NA19704"]

# 这个样本列表是从1kGP_chr1_SV.15.target.vcf.gz中提取出来的
target_member_ls = ["HG00150","HG00106","HG00105","HG00109","HG00243","HG00253","HG00155","HG00237","HG00261","HG00244","HG00185","HG00269","HG00320","HG00334","HG00382","HG00188","HG00343","HG00179","HG00266","HG00383","NA12272","NA11831","NA12829","NA11917","NA12827","NA07357","NA12813","NA12006","NA12749","NA12043","NA20785","NA20807","NA20759","NA20758","NA20539","NA20752","NA20809","NA20792","NA20763","NA20543","HG01507","HG01626","HG02235","HG01707","HG01506","HG01537","HG01623","HG02234","HG01509","HG01747","HG02016","HG02128","HG02058","HG02056","HG01596","HG02050","HG01865","HG02073","HG02522","HG02026","NA19075","NA19072","NA18959","NA18970","NA18951","NA19055","NA18948","NA18943","NA19091","NA18986","HG00728","HG00671","HG00533","HG00608","HG00473","HG00427","HG00537","HG00595","HG00623","HG00593","NA18638","NA18532","NA18579","NA18565","NA18557","NA18582","NA18605","NA18555","NA18645","NA18539","HG01046","HG02353","HG01794","HG02184","HG02395","HG01815","HG02373","HG02164","HG02351","HG02152","HG03772","HG04211","HG03781","HG03862","HG03779","HG03973","HG03871","HG04219","HG03789","HG03778","NA21127","NA20897","NA20866","NA21133","NA20881","NA20874","NA21088","NA20885","NA20890","NA21091","HG03998","HG03673","HG04033","HG04228","HG03681","HG04099","HG04106","HG03711","HG03850","HG03689","HG04177","HG03600","HG04134","HG04182","HG03904","HG03927","HG03907","HG04150","HG03802","HG03616","HG03019","HG03708","HG02790","HG02687","HG03640","HG02652","HG02737","HG03489","HG02731","HG01589","HG02150","HG01578","HG01955","HG01565","HG02253","HG02002","HG02104","HG02291","HG01976","HG01967","NA19788","NA19654","NA19651","NA19670","NA19661","NA19795","NA19779","NA19755","NA19740","NA19777","HG01272","HG01550","HG01391","HG01552","HG01250","HG01140","HG01491","HG01374","HG01479","HG01126","HG01102","HG01097","HG01191","HG01167","HG00740","HG01398","HG01182","HG01169","HG01302","HG01067","NA19923","NA19916","NA20356","NA20321","NA20289","NA20362","NA20351","NA20291","NA20126","NA20287","HG02332","HG02316","HG01886","HG01894","HG02546","HG02014","HG02419","HG02489","HG02455","HG01896","HG03343","HG03307","HG02943","HG03117","HG03135","HG03169","HG03312","HG03342","HG02968","HG03352","NA19456","NA19472","NA19466","NA19307","NA19321","NA19448","NA19310","NA19025","NA19437","NA19434","HG02595","HG02573","HG03246","HG02588","HG02869","HG02624","HG03033","HG02839","HG02568","HG03028","HG03057","HG03397","HG03085","HG03571","HG03096","HG03556","HG03064","HG03563","HG03060","HG03548","NA18489","NA18917","NA18487","NA19213","NA19147","NA19101","NA19099","NA18861","NA19191","NA19102"]

cuteSVTrio_relative_member_ls = ["HG002","HG003","HG004","HG005","HG006","HG007",
                        "HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877",
                        "HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889",
                        "HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890",
                        "HG01891","HG02668","HG03371","HG03453","NA19705",
                        "HG01890","HG02666","NA19703","HG03369","HG03452",
                        "HG01889","HG02667","HG03370","HG03451","NA19704"]

HGSVC_relative_member_ls = ["HG00171","HG00268","HG00358","HG00512","HG00513","HG00514",
                   "HG00731","HG00732","HG00733","HG01890","HG02282","HG02587",
                   "HG02666","HG02769","HG03009","HG03371","HG03452","HG03520",
                   "HG03732","NA18534","NA19036","NA19129","NA19238","NA19239",
                   "NA19240","NA19331","NA19384","NA19705","NA19836","NA19983",
                   "HG01455","HG01456","HG02281","HG02585","HG02586","HG02770","HG02768","HG02667",
                   "HG02668","HG03727","HG03721","NA19128","NA19127","NA19982","NA19713"]


area_ls = ["AFR","EAS","EAS","AMR","EUR","AFR","AFR","AFR","EUR","AFR","AFR","AFR","AFR","AFR"]

chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

# 将输入vcf文件按照比例分为两个文件，对应的是imputation过程中的target文件和answer文件
def split_vcf_file(read_vcf,target_vcf,answer_vcf,target_rate) :
    af_range_ls = [0,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1]
    af_random_ls = [0,0.75,0.75,0.6,0.6,0.6,0.5,0.5,0.45]
    fh = gzip.open(read_vcf, 'rt') if read_vcf.endswith('.gz') else open(read_vcf, 'rt')
    target_lines = []
    answer_lines = []
    sample_index = []
    for line in fh:
        if line == "" :
            continue
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                line_split_ls = line.split("\t")
                new_line_split_ls = line_split_ls[0:9]
                for sample_name in target_member_ls :
                    sample_index.append(line_split_ls.index(sample_name))
                    new_line_split_ls.append(sample_name)            
                target_lines.append("\t".join(new_line_split_ls))
                answer_lines.append("\t".join(new_line_split_ls))
            else :
                target_lines.append(line)
                answer_lines.append(line)
        else :
            line_split_ls = line.split("\t")
            new_line_split_ls = line_split_ls[0:9]
            is_sv_flag = False
            for i in sample_index :
                if line_split_ls[i] == "0" :
                    line_split_ls[i] = "0|0"
                if line_split_ls[i] == "1" :
                    line_split_ls[i] = "1|1"
                if int(line_split_ls[i][0])+int(line_split_ls[i][2]) > 1 :
                    is_sv_flag = True
                new_line_split_ls.append(line_split_ls[i])
            if is_sv_flag :
                af = float(line_split_ls[7].split("AF=")[1].split(";")[0])
                for j in range(1,len(af_range_ls)) :
                    if af > af_range_ls[j-1] and af <= af_range_ls[j] :
                        #print(j)
                        if random.random() > af_random_ls[j] :
                            target_lines.append("\t".join(new_line_split_ls))
                        else :
                            answer_lines.append("\t".join(new_line_split_ls))
    fh.close()
    #print(len(target_lines))
    #print(len(answer_lines))
    with open(target_vcf,'w') as f:
        for line in target_lines :
            f.write(line+"\n")
    f.close()
    with open(answer_vcf,'w') as f:
        for line in answer_lines :
            f.write(line+"\n")
    f.close()


def main(argv):
    split_vcf_file("1kGP_chr1_SV.vcf.gz".replace("chr1",argv[0]),"1kGP_chr1_SV.260.target.vcf".replace("chr1",argv[0]),"1kGP_chr1_SV.260.answer.vcf".replace("chr1",argv[0]),0.5)
    

if __name__ == "__main__":
    main(sys.argv[1:])