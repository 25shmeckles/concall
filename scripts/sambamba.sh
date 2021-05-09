# check calculate_depth_*
name=$1
sambamba depth base -L chr17:7668402-7687550  ${name}.sorted.bam > ${name}_sambamba_output_TP53.txt
