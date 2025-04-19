cuteSV-Trio

安装：
cd cuteSV-Trio/src/
python setup.py install

使用样例：
mkdir tmp/ ; 
cuteSVTrio --retain_work_dir --write_old_sigs --performing_phasing -p HiFI -r ref/human_hs37d5.fasta -o tmp.vcf -w tmp/ --family_mode M1 --input_offspring HG002.bam --input_parent_1 HG003.bam --input_parent_2 HG004.bam --threads 32 --min_support_list 5,5,5; 

相比于cuteSV增加以下参数：
--performing_phasing：执行SV的phasing，不设置时只执行calling
--min_support_list：每个个体SV过滤的最小支持read数量，按顺序分别为offspring：input_parent_1：input_parent_2