#configfile:"configfiles/config-test2.yaml"
SUP_SAMPLES = config['SUP_SAMPLES']
TYPES = ["bb","ins"]

# gz or not gz
if config['gz'] == True:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")
else:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq")
print(f"There are {len(SAMPLES)} samples, {SAMPLES}.")
rule all:
    input:
        # timestamp, tide results, medaka results, cutadapt_info, bb_barcode(to be implement), tagged_bam_files (to be implement), samtools stats
# AFTER TEMP: enable this!!
#        expand("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle", SUP_SAMPLE=SUP_SAMPLES),
# AFTER TEMP: enable this!! bwa index for both tide and medaka
#        expand("output/{SUP_SAMPLE}/07_stats_done/bwa_index.done", SUP_SAMPLE=SUP_SAMPLES),
# Temp modification for plotting structure file
        expand("output/{SUP_SAMPLE}/02_split/stats/{sample}.pickle", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
# samtool stats for tide
#        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done", SUP_SAMPLE=SUP_SAMPLES),
#        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_medaka.done", SUP_SAMPLE=SUP_SAMPLES),
# align bb sequence for extracting barcode (can also use tide to extract barcode)
#        expand( "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.bam", SUP_SAMPLE=SUP_SAMPLES),
# per repeat bwa + sambamba result -- only for targeted
        #expand("output/{SUP_SAMPLE}/04_done/{type}_sambamba.done", SUP_SAMPLE=SUP_SAMPLES, type = TYPES)
localrules: all, bwasw, bwa_mem, get_timestamp, bedtool_getfasta, gz_fastq_get_fasta, fastq_get_fasta, aggregate_python, aggregate_tide, count_repeat, sambamba

rule get_timestamp:
    input:
        fasta = expand("output/{SUP_SAMPLE}/00_fasta/{SAMPLE}.fasta", SAMPLE=SAMPLES, SUP_SAMPLE=SUP_SAMPLES)
        #fastq = ancient(config['rawdir']+"/{sample}.fastq.gz"),
    output:
        timestamp = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle"
    params:
        fq = config['rawdir'],
        name= "{SUP_SAMPLE}"
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/get_timestamp.py -i {params.fq} -o {output.timestamp} -n {params.name} --datype 'fq' "



rule bedtool_getfasta:
#    group: "bowtie_split"
    input:
        seg = "data/seg/{seg}.bed",
        ref = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
	#ref = config["genome"]
    output:
	"data/seg/{seg}.fa"
#        expand("data/seg/{seg}.fa",seg = SEG)
    shell:
        "bedtools getfasta -fi /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -bed {input.seg} -fo {output}"

# convert fastq files to fasta files.
rule gz_fastq_get_fasta:
    # this rule takes > 5 sec < 60 sec to generate output files while submitted to cluster.
    input:
        gz = ancient(config['rawdir']+"/{sample}.fastq.gz")
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"),
    conda:
        "envs/bt.yaml"
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_{sample}_pyfastx_gz_fastq_to_fasta.log"
    shell:
        "pyfastx fq2fa {input.gz} -o {output.fasta} > {log}"

rule fastq_get_fasta:
    input:
        fastq  = ancient(config['rawdir']+"/{sample}.fastq")
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    conda:
        "envs/bt.yaml"
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_{sample}_pyfastx_fastq_to_fasta.log"
    shell:
        "pyfastx fq2fa {input.fastq} -o {output.fasta} > {log}"
#        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"


#rule batch_fasta:
#    input:
#        fasta = "output/{SUP_SAMPLE}/00_fasta/{meta_sample}.fasta"
#    output:
#        fasta = dynamic("output/{SUP_SAMPLE}/001_split_fasta/{sample}.fasta")
#    params:
#        prefix = "output/{SUP_SAMPLE}/001_split_fasta/{meta_sample}",
#    shell:
#        "split -l 8000 --numeric-suffixes {input.fasta} {params.prefix} --additional-suffix=.fasta"

rule bowtie_build:
    # How to write output? answer here: https://www.biostars.org/p/342988/
#    group: "bowtie_split"
    input:
        fasta = ancient("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"),
        done =  ancient("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done")
    params:
        re_path="output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    log:
        "log/{SUP_SAMPLE}/bt-build_{sample}.log"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_bowtie_build.txt"
    threads: 8
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)    
    output:
        touch('output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done')
    shell:
        "bowtie2-build --threads {threads} -f {input.fasta} {params.re_path} > {log} 2>&1"

#rule bowtie_wrapper_map:
#    group:"bowtie_split"
#    input:
#        split_by = "data/seg/bb4.fa",
#        done = "output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done"
#    output:
#        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
#    log:
#        "log/bt-split_{SUP_SAMPLE}_{sample}.log"    
#    params:
#        index = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference",
#        extra="--local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -f"
#    threads : 4
#    wrapper:
#        "0.42.0/bio/bowtie2/align"
        
    
rule bowtie_map_backbone_to_read:
    ## --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.5
    # we use -L 15: even more sensitive
    ## -i the interval between extracted seeds. (eg S, 1, 0.5)
    # interval function f(x) = <init1> + <init2> * sqrt(x). x = read length
    ## -rdg: read gap open, -rfg: ref gap open.
    # A gap of N gets a penalty of var1 + var2*N. default: 5,3
    ## -mp MX,MN (maximum and minimum mismatch penalties): modify alignment score considering Q(Phred score)
    ## -ma match bonus
    # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
#    group: "bowtie_split"
    input:
        split_by = config["backbone_fa"],
        done = ancient("output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done")
    params:
        # -x expect to find index files at current folder, then in the
        # directory specified in the BOWTIE2_INDEXES environment variable.
        # if this doesn't work, try:
        # export BOWTIE2_INDEXES=/path/to/my/bowtie2/databases/
        basename = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference",
        #runtime="6h"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 40000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1),
    conda:
        "envs/bt.yaml"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_map_backbone_to_read_time.txt"
    log:
        "log/{SUP_SAMPLE}/bt-split_{sample}.log"
    output:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
    shell: "bowtie2 --local --very-sensitive -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"
    #shell: "bowtie2 --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"

#rule clean_bowtie_ref:
#    input:
#        done = "output/02_split/{sample}_split.done",
#    output:
#        touch("output/02_split/{sample}_clean.done")
#    shell:
#        "rm -rf output/01_bowtie/{sample}/"

rule split_by_backbone:
#    group: "bowtie_split"
    input:
        sam = ancient("output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"),
        fasta = ancient("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    params:
        min_insert_length = config['min_insert_length'],
        max_insert_length = config['max_insert_length']
    benchmark:
# specify wildcard!!!
        "log/benchmark/{SUP_SAMPLE}xxx{sample}.sam2_split_time.txt"
    output:
        bb = "output/{SUP_SAMPLE}/02_split/bb/{sample}.fasta",
        ins = "output/{SUP_SAMPLE}/02_split/ins/{sample}.fasta",
        stats = "output/{SUP_SAMPLE}/02_split/stats/{sample}.csv",
        stats_read_dict = "output/{SUP_SAMPLE}/02_split/stats/{sample}.pickle"
    log:
        "log/{SUP_SAMPLE}/sam2_split_{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        runtime=1
    conda:
        "envs/pysam.yaml"
    shell:
        #"python scripts/sam2.py {input.sam} {input.fastq} {output}"
        "python3 scripts/sam2.py {input.sam} {input.fasta} {output.bb} {output.ins} {output.stats} {params.min_insert_length} {params.max_insert_length} {output.stats_read_dict}"

rule smolecule_ins:
#    group: "smolecule"
    input:
#        venv = ancient(IN_MEDAKA),
        fasta = "output/{SUP_SAMPLE}/02_split/ins/{sample}.fasta"
    params:
        path = directory("output/{SUP_SAMPLE}/03_consensus/ins/{sample}/"),
    output:
        con_folder = directory("output/{SUP_SAMPLE}/03_consensus/ins/{sample}/"),
#        consensus = "output/{SUP_SAMPLE}/03_consensus/ins/{sample}/consensus.fasta",
#        consensus_folder = "output/{SUP_SAMPLE}/03_consensus/ins/{sample}",
        done = touch("output/{SUP_SAMPLE}/04_done/{sample}_ins.done"),
#        consensus = "output/{SUP_SAMPLE}/03_consensus/ins/{sample}/consensus.fasta"
    log:
        "log/{SUP_SAMPLE}/smol_ins_{sample}.log"
    threads: 4
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}.smolecule_ins_time.txt"
    resources:
        attempt = lambda wildcards, attempt: attempt,        
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    conda:
        "envs/smolecule-env.yaml"
    shell:
         "ulimit -c 0; medaka smolecule --length 30 --depth 1 --threads {threads} {output.con_folder} {input.fasta} > {log}{resources.attempt} 2>&1"
         #set +u; source activate snake-mdk; set -u;
         #medaka smolecule --length 50 --depth 5 --threads {threads} {input.fasta} {output.path} > {log} 2>&1

rule bwa:
    input:
        fa = config['rawdir']+"{SAMPLES}.bam.bai",
    params:
        ref = config['genome']
    output:
        done = touch("output/{SAMPLES}/bwa.done"),
        sam = "output/{SAMPLES}.sam"
    resources:
        attempt = lambda wildcards, attempt: attempt,
        mem_mb=lambda wildcards, attempt: attempt * 4,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    log: '{SAMPLES}.log'

rule smolecule_bb:
#    group: "smolecule"
    input:
        fasta = ancient("output/{SUP_SAMPLE}/02_split/bb/{sample}.fasta")
    output:
        con_folder = directory("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/"),
        done = touch("output/{SUP_SAMPLE}/04_done/{sample}_bb.done"),
#        consensus = "output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta"
#    params:
#        path = directory("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/")
    log:
        "log/{SUP_SAMPLE}/smol_bb_{sample}.log"
    threads: 4
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}.smolecule_bb_time.txt"
    resources:
        attempt = lambda wildcards, attempt: attempt,
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    conda:
        "envs/smolecule-env.yaml"
    shell:
         "ulimit -c 0; medaka smolecule --length 30 --depth 1 --threads {threads} {input.fasta} {output.con_folder} > {log}{resources.attempt} 2>&1"

# how to make sure all samples are done?
#ALL_BB, = glob_wildcards(output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta")


#checkpoint clustering:
#    input:
#        expand("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta", SUP_SAMPLE=SUP_SAMPLES,  sample = SAMPLES)
#    output:
#        clusters=directory("output/{SUP_SAMPLE}/05_aggregated")
#    shell:
#        "mkdir {output}"

#rule intermediate:
#    # everything post processing goes in here.
#    input:
#        "output/aggregated/{SUP_SAMPLE}/{sample}.csv"
#    output:
#        "output/post_processing/{SUP_SAMPLE}/{sample}.csv"
#    shell:
#        "cp {input} {output}"

#def aggregate_input(wildcards):
#    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
#    return expand("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta",
#           SUP_SAMPLE=SUP_SAMPLES,
#           sample=SAMPLES)
#           sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}/consensus.fasta")).sample)

#rule aggregate_csv:
#    input:
#        aggregate_input
#    output:
#        csv = "output/{SUP_SAMPLE}/05_aggregated/stats.csv"
#    shell:
#        "cat {input} > {output.csv}"

rule aggregate_python:
    input:
        mdk_bb = expand("output/{SUP_SAMPLE}/04_done/{sample}_bb.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
        mdk_ins = expand("output/{SUP_SAMPLE}/04_done/{sample}_ins.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
    output:
        csv = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_medaka_stats.csv",
        bb = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_bb.fasta",
        ins = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_ins.fasta",
        done = touch("output/{SUP_SAMPLE}/04_done/aggregate.done")
    params:
        bb = "output/{SUP_SAMPLE}/03_consensus/bb",
        ins = "output/{SUP_SAMPLE}/03_consensus/ins",
        stats = "output/{SUP_SAMPLE}/02_split/stats"
    shell:
        "cat {params.bb}/*/consensus.fasta > {output.bb};"
        "cat {params.ins}/*/consensus.fasta > {output.ins};"
        "cat {params.stats}/*.csv > {output.csv};"

rule aggregate_tide:
    input:
        tide_fasta = expand("output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
        tsv=expand("output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
    output:
        tide = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta"),
        tsv = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide.tsv",
        done = touch("output/{SUP_SAMPLE}/04_done/aggregate_tide.done")
    params:
        tide = "output/{SUP_SAMPLE}/05_aggregated/tide",
    shell:
        "cat {params.tide}/*_tide_consensus.fasta > {output.tide};"
        "cat {params.tide}/*_tide_consensus.tsv > {output.tsv}"


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


#rule bwa_after_cutadapt:
#    input:
#        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta",        
#    output:
#        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa-whole-{SUP_SAMPLE}_clean_bb.done"),
#        sam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sam",
#        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.bam",
#        sorted = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam"
#    threads: 1
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 4000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 4)
#    conda:
#       "envs/bt.yaml"
#    params:
#        ref_genome_fasta = config["ref_genome_final"],
#        name = "{SUP_SAMPLE}"
#    shell:
#        "bash scripts/bwa_mem.sh {input.fasta} {output.sam} {params.ref_genome_fasta} {params.name} {threads} {output.bam} {output.sorted}"

rule bwa_wrapper_after_cutadapt:
    input:
        reads="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper.done")
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
        "0.50.0/bio/bwa/mem"

rule bwa_wrapper_bb:
    input:
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.sorted.bam",
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
        "0.50.0/bio/bwa/mem"

rule bwa_index:
    input:
        bam="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam",
#        tide_bam= "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam",
        bb_bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.sorted.bam",
    output:
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_index.done")
    conda:
        "envs/bt.yaml"
    shell:
        "samtools index {input.bam};"
#        "samtools index {input.tide_bam};"
        "samtools index {input.bb_bam};"
#rule bwa_whole:
#    input:
##        csv = "output/{SUP_SAMPLE}/05_aggregated/stats.csv",
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_{type}.fasta",
##        ins = "output/{SUP_SAMPLE}/05_aggregated/all_consensus_ins.fasta"
#    output:
#        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa-whole-{SUP_SAMPLE}-{type}.done"),
#        sam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-{type}.sam",
#        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-{type}.bam",
#        sorted = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-{type}.sorted.bam"
#    threads: 1
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 8000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 4)
#    conda:
#       "envs/bt.yaml"
#    params:
#        ref_genome_fasta = config["ref_genome_final"],
#        name = "{SUP_SAMPLE}"
#    benchmark:
#        "log/benchmark/{SUP_SAMPLE}-bwa-whole-time-{type}.txt"
#    shell:
#        "bash scripts/bwa_mem.sh {input.fasta} {output.sam} {params.ref_genome_fasta} {params.name} {threads} {output.bam} {output.sorted}" 

rule count_repeat:
    input:
        csv = rules.aggregate_python.output.csv
    output:
        folder = directory("output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}"),
        done = touch("output/{SUP_SAMPLE}/04_done/bin_name_{type}.done")
    conda:
        "envs/bt.yaml"
    log:
        "log/{SUP_SAMPLE}/{type}_count_repeat.log"
    shell:
        "python scripts/fastq_splitby_consensus.py {output.folder} {input.csv} > {log};"

#rule split_fasta:
#    input:
#         donefile = "output/{SUP_SAMPLE}/04_done/{type}_bin_name.done",
#         fasta = "output/{SUP_SAMPLE}/05_aggregated/all_consensus_{type}.fasta"
#    output:
#        split = directory("output/{SUP_SAMPLE}/05_aggregated/02_split_{type}"),
#        done = touch("output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done")
#    conda:
#        "envs/bt-new.yaml"
#    params:
#        txt = "output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}",
#        split = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type}"
#    shell:
#        "python scripts/split_fasta.py {output.split} {params.txt} {input.fasta}"

# this rule is broken
rule postprocessing:
    input:
        all_fasta = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_{type}.fasta",
# in params        txt = "output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}"
        done= "output/{SUP_SAMPLE}/04_done/bin_name_{type}.done"
    output:
        done = touch("output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done"),
        split_folder = directory("output/{SUP_SAMPLE}/05_aggregated/02_split_{type}"),
        bwa_folder = directory("output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}")
    threads: 2
    params:
        txt = "output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}",
        ref_genome_fasta = config["genome"],
        type= "{type}"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)    
    shell:
        "bash scripts/postprocessing.sh {input.all_fasta} {params.txt} {output.split_folder} {output.bwa_folder} {params.type} {params.ref_genome_fasta}"


rule bwasw:
    input:
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 1
    params:
        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["genome"],
        name = "{SUP_SAMPLE}"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:
       "envs/bt.yaml"
    shell:
        "bwasw \
        -b 5 \
        -q 2 \
        -r 1 \
        -z 10 \
        -T 15 \
        -t 4 \
        -f {output.sam} {params.ref_genome_fasta} {params.fasta}"

rule sambamba:
    input:
        "output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done"
    params:
        bamfiles_folder = directory("output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/")
    output:
        done = touch("output/{SUP_SAMPLE}/04_done/{type}_sambamba.done")
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "python scripts/calculate_depth_new_py37.py -i {params.bamfiles_folder}"

#rule bwasw:
#    input:
#        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
#    output:
#        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
#        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
#        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
#    threads: 1
#    params:
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
#        sorted = "output/{SUP_SAMPLE}/06_sorted",
#        ref_genome_fasta = config["ref_genome_final"],
#        name = "{SUP_SAMPLE}"
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 4000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 4)
#    conda:
#       "envs/bt.yaml"
#    shell:
#        "bwasw \
#        -b 5 \
#        -q 2 \
#        -r 1 \
#        -z 10 \
#        -T 15 \
#        -t 4 \
#        -f {output.sam} {params.ref_genome_fasta} {params.fasta}"
#
#rule sambamba:
#    input:
#        "output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done"
#    params:
#        bamfiles_folder = directory("output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/")
#    output:
#        done = touch("output/{SUP_SAMPLE}/04_done/{type}_sambamba.done")
#    conda:
#        "envs/bt.yaml"
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 2000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 1)
#    shell:
#        "python scripts/calculate_depth_new_py37.py -i {params.bamfiles_folder}"

### WOrking to replace last part of postprocessing

#rule bwa_wrapper_per_repeat:
#    input:
#        reads="output/{SUP_SAMPLE}/05_aggregated/02_split_fasta/{rep}_consensus_{type}.fasta"
#    output:
#        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{rep}.bam",
#        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_per_repeat.done")
#    log:
#        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa_per_repeat.log"
#    params:
#        index=config["ref_genome_final"],
#        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
#        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
#        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
#        sort_extra="-l 9"            # Extra args for samtools/picard.
#    threads: 4
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 10000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 4)
#    wrapper:
#        "0.50.0/bio/bwa/mem"


rule bwa_mem:
    input:
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:
       "envs/bt.yaml"
    params:
        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["genome"],
        name = "{SUP_SAMPLE}"
    shell:
        "bash scripts/bwa_mem.sh {params.fasta} {output.sam} {params.ref_genome_fasta} {params.name} {threads} {output.bam} {output.sorted}"
#        "bwa mem -t {threads} -c 100 -M -R’@RG\tID:{params.name}\tSM:{params.name}\tPL:NANOPORE\tLB:{params.name}’ {params.ref_genome_fasta} {input.fasta} > {output.sam} "
#FILE_FASTA=$1
#FILE_SAM=$2
#FILE_REF=$3
#NAME=$4
#THREADS=$5
#FILE_BAM=$6
#FILE_SORTED_BAM=$7
#

rule final_stats_aggregation_bb:
    input:
        sorted = expand("output/{SUP_SAMPLE}/05_aggregated/03_bwa_bb/{count_name}.sorted.bam", SUP_SAMPLE=SUP_SAMPLES, count_name=range(1,41))
    output:
        done=touch("output/{SUP_SAMPLE}/07_stats_done/stats_bb.done"),
        out_dir = directory("output/{SUP_SAMPLE}/06_stats_bb/")
    params:
        base_folder = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_bb/",
#        out_folder = "output/{SUP_SAMPLE}/06_stats_bb/",
        type = "bb"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {output.out_dir}"
#        "python scripts/calculate_depth_new_py37.py -d {params.base_folder}; mkdir -p {params.base_folder}/04_sam {params.base_folder}/05_bam {params.base_folder}/06_sorted_bam {params.base_folder}/../../06_stats_{params.type};mv {params.base_folder}/*.sorted.bam {params.base_folder}/06_sorted_bam; mv {params.base_folder}/*.sam {params.base_folder}/04_sam; mv {params.base_folder}/*.bam* {params.base_folder}/05_bam; mv {params.base_folder}/*sambamba_out* {params.base_folder}/../../06_stats_{params.type}"


rule final_stats_aggregation_ins:
    input:
        sorted = expand("output/{SUP_SAMPLE}/05_aggregated/03_bwa_ins/{count_name}.sorted.bam", SUP_SAMPLE=SUP_SAMPLES, count_name=range(1,41))
    output:
        done=touch("output/{SUP_SAMPLE}/07_stats_done/stats_ins.done"),
        out_dir = directory("output/{SUP_SAMPLE}/06_stats_ins/")
    params:
        base_folder = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_ins/",
#        out_folder = "output/{SUP_SAMPLE}/06_stats_ins/",
        type = "ins"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {output.out_dir}"
#        "mkdir -p {params.base_folder}/../../06_stats_{params.type}; python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {params.out_folder}"
#        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {params.out_folder}; mkdir -p {params.base_folder}/../../06_stats_{params.type};cp {params.base_folder}/*sambamba_out* {params.base_folder}../../06_stats_{params.type}"
#        "python scripts/calculate_depth_new_py37.py -d {params.base_folder}; mkdir -p {params.base_folder}/04_sam {params.base_folder}/05_bam {params.base_folder}/06_sorted_bam {params.base_folder}/../../06_stats_{params.type};mv {params.base_folder}/*.sorted.bam {params.base_folder}/06_sorted_bam; mv {params.base_folder}/*.sam {params.base_folder}/04_sam; mv {params.base_folder}/*.bam* {params.base_folder}/05_bam; mv {params.base_folder}/*sambamba_out* {params.base_folder}/../../06_stats_{params.type}"



rule tidehunter_conda:
    input:
       prime_3=config['5_prime'],
       prime_5=config['3_prime'],       
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        tsv="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv",
        done=touch("output/{SUP_SAMPLE}/04_done/{sample}_tide.done")
    threads: 4
    conda:
        "envs/tidehunter.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "TideHunter -f 2 -t {threads} -5 {input.prime_5} -3 {input.prime_3} {input.fasta} > {output.tsv}"

rule trim_tide:
    # before BWA
    input:
        tsv="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv"
    output:
        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
#        pickle="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus_tsv.pickle.gz"
    conda:
       "envs/bt.yaml"
    shell:
        "python3 scripts/trim_tide_fasta_long_readname.py {input.tsv} {output.fasta}"

rule cutadapt_tide:
    input:
        tide="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta",
        done="output/{SUP_SAMPLE}/04_done/aggregate_tide.done"

    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_tide_cut_info.csv",
        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cutadapt_tide_summary.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:"envs/bt.yaml"
    shell:
        "cutadapt -e 0.15 -b GGGCGGTATGTCATGCACACGAATCCCGAAGAnTGTTGTCCATTCATTGAATATGAGATCTCnATGGTATGATCAATATnCGGATGCGATATTGATAnCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC --action=lowercase --info-file {output.cut_info} -o {output.fasta} {input.tide} > {output.summary}"

rule bwa_wrapper_tide:
    input:
        reads="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam",
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
        "0.50.0/bio/bwa/mem"

rule plot_samtools_stats:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam", 
    output:
        stats = "output/{SUP_SAMPLE}/08_samtools_stats/tide/{SUP_SAMPLE}.stats",
        plot = directory("output/{SUP_SAMPLE}/08_samtools_stats/tide/{SUP_SAMPLE}_plot/"),
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done"),
    params:
        name = "{SUP_SAMPLE}_tide"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {output.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.plot}{params.name}_RL.txt;"

rule plot_samtools_stats_medaka:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam",
    output:
        stats = "output/{SUP_SAMPLE}/08_samtools_stats/medaka/{SUP_SAMPLE}.stats",
        plot = directory("output/{SUP_SAMPLE}/08_samtools_stats/medaka/{SUP_SAMPLE}_plot/"),
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_medaka.done"),
    params:
        name = "{SUP_SAMPLE}_medaka"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {output.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.plot}{params.name}_RL.txt;"
#rule medaka:
#    input:
#        read = dynamic("data/output/03_read_repeats/{type}/{sample}/{readname}.fastq"),
#        draft = "data/seg/{seg}.fa",
#        ref = config["genome"]
#    params:
#        model = "r941_min_high"
#    output:
#        "data/output/04_consensus/{seg}/{sample}/{readname}"
#    shell:
#        "medaka_consensus -i {input.ref} -d {input.draft} -o {output} -m {params.model}"

#rule clean_fastq:
#    input:
#        "output/04_done/{sample}_bb.done"
#    output:
#        done = touch("output/04_done/{sample}_cleanup.done")
#    shell:
#        "rm -rf output/00_fasta/{sample}/* output/01_bowtie/{sample}*"
if config["mail"]:
    onsuccess:
        print("Workflow finished, no error. Success!")
        shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
    onerror:
        print("An notice sent to Liting by mail.")
        shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
