import sys
import gzip
import subprocess
import csv
import os
import json

def parse_bnd(alt):
    if alt[0] == ']':
        form = ']]N'
        chr2 = alt.split(':')[0][1:]
        pos2 = int(alt.split(':')[1][:-2])
    elif alt[0] == '[':
        form = '[[N'
        chr2 = alt.split(':')[0][1:]
        pos2 = int(alt.split(':')[1][:-2])
    else:
        if alt[1] == ']':
            form = 'N]]'
            chr2 = alt.split(':')[0][2:]
            pos2 = int(alt.split(':')[1][:-1])
        else:
            form = 'N[['
            chr2 = alt.split(':')[0][2:]
            pos2 = int(alt.split(':')[1][:-1])
    return pos2, chr2

def trans_gt(gt):
    gt1 = 0 if gt[0] == '.' else int(gt[0])
    gt2 = 0 if gt[2] == '.' else int(gt[2])
    return gt1+gt2

base = sys.argv[1]
re_dist = 1000
tra_dict = dict()
tra_id = dict()
ii = 0
f = gzip.open(base, 'rt') if base.endswith('.gz') else open(base, 'rt')
for line in f:
    if line[0] == '#':
        continue
    seq = line.strip().split('\t')
    svtype = seq[7].split('SVTYPE=')[1].split(';')[0]
    if svtype != 'TRA' and svtype != 'BND':
        continue
    ii += 1
    chr1 = seq[0]
    pos = int(seq[1])
    pos2, chr2 = parse_bnd(seq[4])
    length = int(seq[7].split(';SVLEN=')[1].split(';')[0])
    gt = trans_gt(seq[9])
    if chr1 not in tra_dict:
        tra_dict[chr1] = list()
    if chr2 not in tra_dict:
        tra_dict[chr2] = list()
    tra_dict[chr1].append([pos, chr2, pos2, pos2+length, ii, gt])
    tra_dict[chr1].append([pos+length, chr2, pos2, pos2+length, ii, gt])
    tra_dict[chr2].append([pos2, chr1, pos, pos+length, ii, gt])
    tra_dict[chr2].append([pos2+length, chr1, pos, pos+length, ii, gt])

call = sys.argv[2]
found = set()
found_gt = set()
tp_gt = 0
call_sv = 0
call_sv_gt = 0
total_sv = 0
f = gzip.open(call, 'rt') if call.endswith('.gz') else open(call, 'rt')
for line in f:
    if line[0] == '#':
        continue
    seq = line.strip().split('\t')
    svtype = seq[7].split('SVTYPE=')[1].split(';')[0]
    if svtype != 'TRA' and svtype != 'BND':
        continue
    if seq[6] != 'PASS':
        continue
    total_sv += 1
    flag = False
    chr1 = seq[0]
    pos = int(seq[1])
    if seq[4] == "<TRA>" :
        pos2 = int(seq[7].split('END=')[1].split(';')[0])
        chr2 = seq[7].split('CHR2=')[1].split(';')[0]
    else : 
        pos2, chr2 = parse_bnd(seq[4])
    gt = trans_gt(seq[9].split(':')[0])
    if chr1 not in tra_dict:
        continue
    for base_sv in tra_dict[chr1]:
        if abs(base_sv[0] - pos) < re_dist and (abs(base_sv[2] - pos2) < re_dist or abs(base_sv[3] - pos2) < re_dist) and base_sv[1] == chr2:
            found.add(base_sv[4])
            if base_sv[5] == gt:
                found_gt.add(base_sv[4])
                tp_gt += 1
                call_sv_gt += 1
            call_sv += 1
            flag = True
            break
recall = 0 if ii==0 else len(found) / ii
precision = 0 if total_sv==0 else call_sv / total_sv
f1 = 0 if precision+recall==0 else 2*precision*recall/(precision+recall)
gt_recall = 0 if ii==0 else len(found_gt) / ii
gt_precision = 0 if total_sv==0 else call_sv_gt / total_sv
gt_f1 = 0 if gt_precision+gt_recall==0 else 2*gt_precision*gt_recall/(gt_precision+gt_recall)
gt_concordance = 0 if call_sv==0 else call_sv_gt / call_sv

os.mkdir(sys.argv[3])
write_file = sys.argv[3]+"/summary.json"
new_dict = {
    "TP-base": len(found),
    "TP-comp": call_sv,
    "FP": total_sv-call_sv,
    "FN": ii-len(found),
    "precision": precision,
    "recall": recall,
    "f1": f1,
    "TP-comp_TP-gt": call_sv_gt,
    "TP-base_TP-gt": len(found_gt),
    "gt_concordance": gt_concordance
}
# 将字典写入到 JSON 文件
with open(write_file, 'w') as json_file:
    json.dump(new_dict, json_file, indent=4)  # indent 参数用于格式化输出

#rm -r /home/user/lixin/cutesv_trio/output_real/cmp_svision/simu/HG002_simu/20.5.5.1.T2T.TRA; python /home/user/lixin/cutesv_trio/cuteSV/code/bench_tra.py /home/user/lixin/cutesv_trio/simu/answer.1.TRA.vcf.gz /home/user/lixin/cutesv_trio/output_real/output_svision/simu/HG002_simu/svision_pro_20.5.5/sample1.svision_pro_v1.8.s5.1.TRA.vcf.gz /home/user/lixin/cutesv_trio/output_real/cmp_svision/simu/HG002_simu/20.5.5.1.T2T.TRA
