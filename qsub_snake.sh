# 20200103
# configfile: config_test.yaml
snakemake --jobs 50 --latency-wait 120 --use-conda --rerun-incomplete --cluster "qsub -l h_rt=24:00:00 -l h_vmem=40G -l tmpspace=100G -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"