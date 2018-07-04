for i in `seq 1 22`;
do
     tmux new-session -d -s "chr${i}" "python3 ../treeseq-inference/src/convert_1kg.py ALL.chr${i}_GRCh38.genotypes.20170504.vcf.gz 20130606_g1k.ped H_sap_chr${i}.samples -a homo_sapiens_chr${i}.vcf.gz -p"
done    


for i in `seq 9`;
do
     ./miniconda3/bin/tmux new-session -d -s "chr${i}" "./miniconda3/bin/python3 bgen_to_samples.py /well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr${i}_v2.bgen 1000G_GRCh38/H_sap_chr${i}.samples UK_BB/UK-BB_chr${i}.samples -p; exec bash -i"
done    
