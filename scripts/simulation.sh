# use make_simulation_bed.py to generate simulated SV bed file.

# simulate family data
# child
VISOR HACK -g human_hs37d5.fasta -b HACK.h1.bed -o HackOutput.h1/
VISOR HACK -g human_hs37d5.fasta -b HACK.h2.bed -o HackOutput.h2/
lrsim -e 0.01 -t 64 -m HG002_ONT_UL.lrsm -b xxxG HackOutput.h1/h1.fa HackOutput.h2/h1.fa > tmp/fam.1.fq
samtools faidx tmp/fam.1.fq
minimap2 -t 64 -ax map-ont human_hs37d5.fasta tmp/fam.1.fq | samtools sort -o tmp/fam.1.bam
samtools index tmp/fam.1.bam
pandepth -i tmp/fam.1.bam -o tmp/fam.1.coverage -t 32
samtools view -@ 32 -s xxxx -o 30x/fam.1.bam tmp/fam.1.bam
samtools index -@ 32 30x/fam.1.bam
pandepth -i 30x/fam.1.bam -o 30x/fam.1.coverage -t 32

# father
VISOR HACK -g human_hs37d5.fasta -b HACK.h3.bed -o HackOutput.h3/
VISOR HACK -g human_hs37d5.fasta -b HACK.h4.bed -o HackOutput.h4/
lrsim -e 0.01 -t 64 -m HG002_ONT_UL.lrsm -b xxxG HackOutput.h3/h1.fa HackOutput.h4/h1.fa > tmp/fam.2.fq
samtools faidx tmp/fam.2.fq
minimap2 -t 64 -ax map-ont human_hs37d5.fasta tmp/fam.2.fq | samtools sort -o tmp/fam.2.bam
samtools index tmp/fam.2.bam
pandepth -i tmp/fam.2.bam -o tmp/fam.2.coverage -t 32
samtools view -@ 32 -s xxxx -o 30x/fam.2.bam tmp/fam.2.bam
samtools index -@ 32 30x/fam.2.bam
pandepth -i 30x/fam.2.bam -o 30x/fam.2.coverage -t 32

# mother
VISOR HACK -g human_hs37d5.fasta -b HACK.h5.bed -o HackOutput.h5/
VISOR HACK -g human_hs37d5.fasta -b HACK.h6.bed -o HackOutput.h6/
lrsim -e 0.01 -t 64 -m HG002_ONT_UL.lrsm -b xxxG HackOutput.h5/h1.fa HackOutput.h6/h1.fa > tmp/fam.3.fq
samtools faidx tmp/fam.3.fq
minimap2 -t 64 -ax map-ont human_hs37d5.fasta tmp/fam.3.fq | samtools sort -o tmp/fam.3.bam
samtools index tmp/fam.3.bam
pandepth -i tmp/fam.3.bam -o tmp/fam.3.coverage -t 32
samtools view -@ 32 -s xxxx -o 30x/fam.3.bam tmp/fam.3.bam
samtools index -@ 32 30x/fam.3.bam
pandepth -i 30x/fam.3.bam -o 30x/fam.3.coverage -t 32

