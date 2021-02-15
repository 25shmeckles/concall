NAME=$1
snakemake --profile slurm --configfiles configfiles/config-$NAME.yaml --snakefile medaka.smk --cores 4 --use-conda --keep-going --restart-times 2 --latency-wait 90 --rerun-incomplete
snakemake --configfiles configfiles/config-$NAME.yaml --snakefile medaka.smk --cores 4 --report reports/$NAME.html
