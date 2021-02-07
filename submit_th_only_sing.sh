NAME=$1
snakemake --profile slurm --configfiles configfiles/config-$NAME.yaml --snakefile tidehunter_only_sing.smk --cores 4 --use-conda --use-singularity --keep-going  --restart-times 3 --latency-wait 90 --rerun-incomplete

