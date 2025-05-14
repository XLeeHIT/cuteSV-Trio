# run SVision-pro
rm -r svision_pro_30.30.30 ; SVision-pro --target_path HG002.30.bam  --base_path HG003.30.bam HG004.30.bam --genome_path human_hs37d5.fasta --model_path model_liteunet_256_8_16_32_32_32.pth --out_path svision_pro_30.30.30 --sample_name sample1 --detect_mode denovo --process 32
rm -r svision_pro_HG002_30 ; SVision-pro --target_path HG002.30.bam  --genome_path human_hs37d5.fasta --model_path model_liteunet_256_8_16_32_32_32.pth --out_path svision_pro_HG002_30 --sample_name sample1 --detect_mode germline --process 32

# run Sniffles2
rm fam.1.vcf ; sniffles -i HG002.30.bam -v fam.1.vcf --reference human_hs37d5.fasta --threads 32 --output-rnames
rm fam.1.snf ; sniffles -i HG002.30.bam --snf fam.1.snf --reference human_hs37d5.fasta --threads 32
rm fam.2.snf ; sniffles -i HG003.30.bam --snf fam.2.snf --reference human_hs37d5.fasta --threads 32
rm fam.3.snf ; sniffles -i HG004.30.bam --snf fam.3.snf --reference human_hs37d5.fasta --threads 32
rm fam.trio.vcf ; sniffles -i fam.1.snf fam.2.snf fam.3.snf --vcf fam.trio.vcf --threads 32

# run cuteSV
rm -r work ; mkdir work ; cuteSV --retain_work_dir --write_old_sigs --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 HG002.30.bam human_hs37d5.fasta sv.vcf work --threads 32 --min_support 5 ;

# run DeBreak
rm -r 30x/ ; mkdir 30x/ ; debreak --bam HG002.30.bam --outpath 30x/ --rescue_large_ins --poa --ref human_hs37d5.fasta -t 32 --min_support 5 ; 

# run cuteSV-Trio
rm -r fam.trio.30.30.30/ ; mkdir fam.trio.30.30.30/ ; cuteSV --retain_work_dir --write_old_sigs --performing_phasing -p HiFI -g T2T -r human_hs37d5.fasta -o output.fam.trio.M1.30.30.30.vcf -w fam.trio.30.30.30 --family_mode M1 --input_offspring HG002.30.bam --input_parent_1 HG003.30.bam --input_parent_2 HG004.30.bam --threads 32 --execute_stage 0 --min_support_list 5,5,5; 

#truvari bench
rm -r 30.30.30.T2T; truvari bench -b GRCh37_HG2-T2TQ100-V1.1_stvar.sv.INDEL.vcf.gz -c output.fam.trio.M1.30.30.30.1.vcf.gz --includebed GRCh37_HG2-T2TQ100-V1.1_stvar.benchmark.bed -o 30.30.30.T2T -p 0 -r 1000 --passonly

#truvari refine
rm -r 30.30.30.T2T; truvari bench --reference human_hs37d5.fasta --includebed GRCh37_HG2-T2TQ100-V1.1_stvar.benchmark.bed --base GRCh37_HG2-T2TQ100-V1.1_stvar.sv.INDEL.vcf.gz --comp output.fam.trio.M1.30.30.30.1.vcf.gz --output 30.30.30.T2T --refdist 2000 -C 5000 --passonly --pick ac ; truvari refine --reference human_hs37d5.fasta --regions 30.30.30.T2T/candidate.refine.bed --use-original-vcfs --align mafft --mafft-params '--auto --thread 32' -t 32 30.30.30.T2T

# use truvari benchmark INS,DEL,INV,DUP, but use the bench_tra.py to benchmark the TRA
rm -r 30.30.30.1.T2T.TRA; python 30.30.30.bench_tra.py 30.30.30.answer.1.TRA.vcf.gz output.fam.trio.M1.30.30.30.1.TRA.vcf.gz 30.30.30.1.T2T.TRA

# run clair3
run_clair3.sh --bam_fn=HG002.30.bam --ref_fn=human_hs37d5.fasta --threads=32 --platform=hifi --model_path=r941_prom_sup_g5014 --output=without_sv_bed

# use make_adjacent_bed.py to make bed file to limit the range of SNPs around SV

# run WhatsHap
whatshap phase --ignore-read-groups --indels -r human_hs37d5.fasta merge_sv.vcf.gz HG002.30.bam -o whatshap.multi.vcf

# run LongPhase
longphase_linux-x64 phase -s merge_output.vcf.gz --sv-file sv.vcf.gz -b HG002.30.bam -r human_hs37d5.fasta -t 32 -o phased.longphase --pb

# use SER_QAN50.py to calculate SER and QAN50 in phasing benchmark
# use CMI.py to calculate CMI in phasing benchmark

# 挑选denovo变异的脚本


