import sys
import gzip

def write_vcf_file(write_vcf_lines,write_vcf_file) :
    with open(write_vcf_file,'w') as f:
        for line in write_vcf_lines :
            f.write(line+"\n")
        f.close()

def divide_denovo_sv(file_base,read_file,denovo_file_ls) :
    header_vcf_lines = []
    denovo_1_vcf_lines = []
    denovo_2_vcf_lines = []
    denovo_3_vcf_lines = []
    denovo_4_vcf_lines = []
    read_vcf = file_base + read_file
    denovo_1_file = file_base + denovo_file_ls[0]
    denovo_2_file = file_base + denovo_file_ls[1]
    denovo_3_file = file_base + denovo_file_ls[2]
    denovo_4_file = file_base + denovo_file_ls[3]
    f = gzip.open(read_vcf, 'rt') if read_vcf.endswith('.gz') else open(read_vcf, 'rt')
    while(True) :
        f_line = f.readline().strip().replace("/","|")
        if f_line == "" :
            break
        if f_line[0] == '#':
            if f_line.startswith('#CHROM'):
                f_split_ls = f_line.split("\t")
                header_vcf_lines.append("\t".join(f_split_ls[0:10]))
            else :
                header_vcf_lines.append(f_line)
        else :
            f_split_ls = f_line.split("\t")
            if "SVTYPE=INS" not in f_split_ls[7] and "SVTYPE=DEL" not in f_split_ls[7] :
                continue
            if abs(int(f_split_ls[7].split("SVLEN=")[1].split(";")[0])) < 50 :
                continue
            family_genotype_ls = [x[0:3] for x in f_split_ls[9:12]]
            if ("0|1" in family_genotype_ls[0] or "1|0" in family_genotype_ls[0]) and "0|0" == family_genotype_ls[1] and "0|0" == family_genotype_ls[2] :
                denovo_1_vcf_lines.append("\t".join(f_split_ls[0:10]))
            elif "1|1" == family_genotype_ls[0] and ("0|0" == family_genotype_ls[1] or "0|0" == family_genotype_ls[2]) :
                denovo_2_vcf_lines.append("\t".join(f_split_ls[0:10]))
            elif ("0|1" in family_genotype_ls[0] or "1|0" in family_genotype_ls[0]) and "1|1" == family_genotype_ls[1] and "1|1" == family_genotype_ls[2] :
                denovo_3_vcf_lines.append("\t".join(f_split_ls[0:10]))
            elif "0|0" == family_genotype_ls[0] and ("1|1" == family_genotype_ls[1] or "1|1" == family_genotype_ls[2]) :
                denovo_4_vcf_lines.append("\t".join(f_split_ls[0:10]))
    f.close()
    
    print("denovo 1 :",len(denovo_1_vcf_lines))
    denovo_1_vcf_lines = header_vcf_lines + denovo_1_vcf_lines
    write_vcf_file(denovo_1_vcf_lines,denovo_1_file)
    print("denovo 2 :",len(denovo_2_vcf_lines))
    denovo_2_vcf_lines = header_vcf_lines + denovo_2_vcf_lines
    write_vcf_file(denovo_2_vcf_lines,denovo_2_file)
    print("denovo 3 :",len(denovo_3_vcf_lines))
    denovo_3_vcf_lines = header_vcf_lines + denovo_3_vcf_lines
    write_vcf_file(denovo_3_vcf_lines,denovo_3_file)
    print("denovo 4 :",len(denovo_4_vcf_lines))
    denovo_4_vcf_lines = header_vcf_lines + denovo_4_vcf_lines
    write_vcf_file(denovo_4_vcf_lines,denovo_4_file)

def main(argv):
    # 该函数的作用是按照mendel遗传定律和三种denovo的定义，对家系sv检测结果进行筛选
    # 因为这个函数没有phasing信息，所以只能通过孟德尔失序性检测突变
    divide_denovo_sv(argv[0],argv[1],[argv[2]+".1.vcf",argv[2]+".2.vcf",argv[2]+".3.vcf",argv[2]+".4.vcf"])


if __name__ == "__main__":
    main(sys.argv[1:])