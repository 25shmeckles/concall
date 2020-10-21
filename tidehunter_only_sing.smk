# this is a snakefile for tidehunter_only pipeline. Provide configfile for a successful run.

# User defined file prefix name
SUP_SAMPLES = config['SUP_SAMPLES']
# check input data type. Should be specified in configfiles.
if config['gz'] == True:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")
else:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq")
print(f"There are {len(SAMPLES)} samples, {SAMPLES}.")

# end products of this pipeline
rule all:
    input:
        expand("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle", SUP_SAMPLE=SUP_SAMPLES),
# bam files index for tide and tide_fl 
        expand("output/{SUP_SAMPLE}/07_stats_done/bam_index_tide.done", SUP_SAMPLE=SUP_SAMPLES),
# samtool stats
        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done", SUP_SAMPLE=SUP_SAMPLES),
# map fasta file to genome without bb
        expand("output/{SUP_SAMPLE}/07_stats_done/bwa_mem_tide_no_bb.done", SUP_SAMPLE=SUP_SAMPLES),
# stats mapped full length bam file
        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_no_bb.done", SUP_SAMPLE=SUP_SAMPLES),
# stats mapped all consensus bam file
        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_no_bb_not_fl.done", SUP_SAMPLE=SUP_SAMPLES),
# map fasta file to bb only. (check which reads contain backbones)
        expand("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_bb_only.done", SUP_SAMPLE=SUP_SAMPLES),
localrules: all, get_timestamp, gz_fastq_get_fasta, fastq_get_fasta, aggregate_tide 

# check if use singularity image for Tidehunter or not. Please specify in configfiles.
if config['sing'] == True:
    ruleorder: tidehunter_sing > tidehunter_conda
    ruleorder:  tidehunter_sing_fl > tidehunter_conda_full_length
else:
    ruleorder: tidehunter_conda > tidehunter_sing
    ruleorder:  tidehunter_conda_full_length > tidehunter_sing_fl
 
ruleorder: bwa_mem > bwa_wrapper_tide_full_length

# retrieve time_stamp for each read, store as pickle file.
rule get_timestamp:
    input:
        fasta = expand("output/{SUP_SAMPLE}/00_fasta/{SAMPLE}.fasta", SAMPLE=SAMPLES, SUP_SAMPLE=SUP_SAMPLES)
    output:
        timestamp = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle"
    params:
        fa = "output/{SUP_SAMPLE}/00_fasta/",
        name= "{SUP_SAMPLE}"
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/get_timestamp.py -i {params.fa} -o {output.timestamp} -n {params.name} --datype 'fa' " 

# convert fastq files to fasta files.
rule gz_fastq_get_fasta:
    # this rule takes > 5 sec < 60 sec to generate output files while submitted to cluster.
    input:
        gz = ancient(config['rawdir']+"/{sample}.fastq.gz")
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fastq = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fastq"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    conda:
        "envs/bt.yaml"
    shell:
        "pyfastx fq2fa {input.fastq} > -o {output.fasta}"
rule fastq_get_fasta:
    input:
        fastq  = ancient(config['rawdir']+"/{sample}.fastq")
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    conda:
        "envs/bt.yaml"
    shell:
        "pyfastx fq2fa {input.fastq} > -o {output.fasta}"
#        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"


rule aggregate_tide:
    input:
        tide_all = expand("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus.fasta", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
        tide_fasta_full_length = expand("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus_full_length.fasta", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
    output:
        tide_all = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus.fasta"),
        tide_full_length = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_full_length.fasta"),
    params:
        tide = "output/{SUP_SAMPLE}/09_tide",
    shell:
        "cat {params.tide}/*_tide_consensus.fasta > {output.tide_all};"
        "cat {params.tide}/*_tide_consensus_full_length.fasta > {output.tide_full_length};"


rule cutadapt:
    # -e error rate
    # -b adaptor sequence (bb) can be anywhere in sequence. If in middle, all down stream of read is trimmed.
    # check cutadapt program
    input:
        ins = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_ins.fasta"
    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cut_info.csv",
        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cutadapt_summary.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:"envs/bt.yaml"
    shell:
        "cutadapt -e 0.15 -b GGGCGGTATGTCATGCACACGAATCCCGAAGAnTGTTGTCCATTCATTGAATATGAGATCTCnATGGTATGATCAATATnCGGATGCGATATTGATAnCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC --action=lowercase --info-file {output.cut_info} -o {output.fasta} {input.ins} > {output.summary}"

rule bwa_wrapper_after_cutadapt:
    input:
        reads="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config['genome'],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"

rule bwa_wrapper_bb:
    input:
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.sorted.bam",
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config['genome'],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"

rule bam_index:
    input:
        tide_bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide.sorted.bam",
        tide_bam_full_length = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl.sorted.bam",
    output:
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bam_index_tide.done")
    conda:
        "envs/bt.yaml"
    shell:
        "samtools index {input.tide_bam};"
        "samtools index {input.tide_bam_full_length}"

rule tidehunter_conda_full_length:
    input:
       prime_3=config['5_prime'],
       prime_5=config['3_prime'],
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        fasta=temp("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus_full_length.fasta"),
    log:
        stdout= "output/{SUP_SAMPLE}/04_done/{sample}_resource_fl.txt"
    conda:
        "envs/tidehunter.yaml"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "TideHunter -t {threads} -5 {input.prime_5} -3 {input.prime_3} -p 20 -a 0.70 -F {input.fasta} > {output.fasta} 2>{log.stdout}"

rule tidehunter_conda:
    input:
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        fasta=temp("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus.fasta"),
    log:
        stdout= "output/{SUP_SAMPLE}/04_done/{sample}_resource.txt"    
    conda:
        "envs/tidehunter.yaml"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
       "TideHunter -t {threads}  {input.fasta} > {output.fasta} 2> {log.stdout}"

rule tidehunter_sing:
    input:
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        fasta=temp("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus.fasta"),
    threads: 4
    log:
        stdout= "output/{SUP_SAMPLE}/04_done/{sample}_resource.txt"    
    singularity:
        "tidehunter.sif"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "/TideHunter-v1.4.2/bin/TideHunter -t {threads} {input.fasta} > {output.fasta} 2> {log.stdout}"

rule tidehunter_sing_fl:
    input:
       prime_3=config['5_prime'],
       prime_5=config['3_prime'],
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        fasta=temp("output/{SUP_SAMPLE}/09_tide/{sample}_tide_consensus_full_length.fasta"),
    threads: 4
    singularity:
        "tidehunter.sif"
    log:
        stdout= "output/{SUP_SAMPLE}/04_done/{sample}_resource_fl.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "/TideHunter-v1.4.2/bin/TideHunter -t {threads} -5 {input.prime_5} -3 {input.prime_3} -p 20 -a 0.70 -F {input.fasta} > {output.fasta} 2>{log.stdout}"

rule trim_tide:
    # trim = cut too long read names into supplemental files
    # before BWA
    input:
        fasta_full_length="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_full_length.fasta",
        fasta_all="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus.fasta"
    output:
        fasta_all="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_trimmed.fasta",
        fasta_full_length="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_full_length_trimmed.fasta",
        metadata_all="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_all.txt"
    shell:
        "sed -n -e 's/^>//p' {input.fasta_all} > {output.metadata_all};" # only export metadata from flexible parameter since it is the same if the pattern is found.
        "sed 's/,.*//' {input.fasta_all} > {output.fasta_all};"
        "sed 's/,.*//' {input.fasta_full_length} > {output.fasta_full_length};"
rule cutadapt_tide:
    input:
        tide="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_full_length_trimmed.fasta",
    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_tide_cutadapt_info.csv",
        fasta_fl_cutadapt="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_fl_cutadapt_tide.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cutadapt_tide_summary.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:"envs/bt.yaml"
    shell:
        "cutadapt -e 0.15 -b GGGCGGTATGTCATGCACACGAATCCCGAAGAnTGTTGTCCATTCATTGAATATGAGATCTCnATGGTATGATCAATATnCGGATGCGATATTGATAnCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC --action=lowercase --info-file {output.cut_info} -o {output.fasta_fl_cutadapt} {input.tide} > {output.summary}"

rule bwa_wrapper_tide:
    # without cutadapt
    input:
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_trimmed.fasta",
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config["genome"],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.58.0/bio/bwa/mem"



rule bwa_wrapper_tide_no_bb:
    input:
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_trimmed.fasta",
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_no_bb.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide_no_bb.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa_no_bb.log"
    params:
        index=config['genome'],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.58.0/bio/bwa/mem"

# preprocessing. Run only if backbone is new.
rule generate_bb_only_ref:
    input:
        bb = config['backbone_fa']
    output:
        done = touch("output/{SUP_SAMPLE}/04_done/gen_bb_ref.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_bwa_gen_bb_ref.log" 
    conda:
        "envs/bt.yaml"
    threads:4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1000,
        runtime = lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "bwa index {input.bb}" 

rule bwa_wrapper_bb_only:
    input:
        ref_built = "output/{SUP_SAMPLE}/04_done/gen_bb_ref.done",
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_consensus_trimmed.fasta",
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb_only.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_bb_only.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa_bb_only.log"
    params:
        index=config['backbone_fa'],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.58.0/bio/bwa/mem"


rule bwa_wrapper_tide_full_length:
    #after cutadapt
    input:
        reads=rules.cutadapt_tide.output.fasta_fl_cutadapt
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide_full_length_reads.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa_full_length.log"
    params:
        index=config['genome'],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.58.0/bio/bwa/mem"

rule bwa_mem:
    input:
        reads=rules.cutadapt_tide.output.fasta_fl_cutadapt
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl.sam",
        bam = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl.bam"),
        sorted = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_mem_tide.done")
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:
       "envs/bt.yaml"
    params:
        ref_genome_fasta = config['genome'],
        name = "{SUP_SAMPLE}"
    shell:
        "bwa mem -t 8 -c 100 -M -R '@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:NANOPORE\\tLB:{params.name}' {params.ref_genome_fasta} {input.reads} > {output.sam};"
        "samtools view -h {output.sam} > {output.bam};"
        "samtools sort -l 7  {output.bam} > {output.sorted};"
        "samtools index {output.sorted};"

rule bwa_mem_ref_no_bb:
    input:
        reads=rules.cutadapt_tide.output.fasta_fl_cutadapt
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl_no_bb.sam",
        bam = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl_no_bb.bam"),
        sorted = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl_no_bb.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_mem_tide_no_bb.done")
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:
       "envs/bt.yaml"
    params:
        ref_genome_fasta = config['genome'],
        name = "{SUP_SAMPLE}"
    shell:
        "bwa mem -t 8 -c 100 -M -R '@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:NANOPORE\\tLB:{params.name}' {params.ref_genome_fasta} {input.reads} > {output.sam};"
        "samtools view -h {output.sam} > {output.bam};"
        "samtools sort -l 7  {output.bam} > {output.sorted};"
        "samtools index {output.sorted}"


rule plot_samtools_stats_no_bb_fl:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_fl_no_bb.sorted.bam",
    output:
        stats = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_fl_no_bb.stats",
        SN = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_SN_tag_read_mapped_fl_no_bb.txt",
        RL =  "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_RL_tag_read_length_fl_no_bb.txt",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_no_bb.done"),
    params:
        name = "{SUP_SAMPLE}_tide",
        plot = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_plot_fl_no_bb/"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {params.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^SN | cut -f 2- > {output.SN};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.RL};"


rule plot_samtools_stats_no_bb:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide_no_bb.sorted.bam",
    output:
        stats = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_no_bb.stats",
        SN = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_SN_tag_read_mapped_no_bb.txt",
        RL =  "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_RL_tag_read_length_no_bb.txt",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_no_bb_not_fl.done"),
    params:
        name = "{SUP_SAMPLE}_tide",
        plot = "output/{SUP_SAMPLE}/05_aggregated/tide_stats_no_bb/{SUP_SAMPLE}_plot_no_bb/"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {params.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^SN | cut -f 2- > {output.SN};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.RL};"


rule plot_samtools_stats:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide.sorted.bam",
    output:
        stats = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}.stats",
        SN = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_SN_tag_read_mapped.txt",
        RL =  "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_RL_tag_read_length.txt",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done"),
    params:
        name = "{SUP_SAMPLE}_tide",
        plot = "output/{SUP_SAMPLE}/05_aggregated/tide_stats/{SUP_SAMPLE}_plot/"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {params.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^SN | cut -f 2- > {output.SN};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.RL};"

if config['mail'] == True:
    onsuccess:
        print("Workflow finished, no error. Success!")
        shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
    onerror:
        print("An notice sent to Liting by mail.")
        shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
