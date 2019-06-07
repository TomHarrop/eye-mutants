#!/usr/bin/env python3

import csv
import multiprocessing
import re
import snakemake
import os
import pathlib2


#############
# FUNCTIONS #
#############

def map_indiv_to_raw_reads(wildcards):
    cs_id = indiv_to_csid[wildcards.indiv]
    return({
        'r1': os.path.join(read_dir, f'{cs_id}/{cs_id}_R1.fq.gz'),
        'r2': os.path.join(read_dir, f'{cs_id}/{cs_id}_R2.fq.gz')
        })


def resolve_path(x):
    return(str(pathlib2.Path(x).resolve()))


###########
# GLOBALS #
###########

sample_key = 'data/sample_key.csv'
read_dir = 'data/raw'
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'
honeybee_ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
vcflib_container = 'shub://TomHarrop/singularity-containers:vcflib_1.0.0-rc2'
biopython_container = 'shub://TomHarrop/singularity-containers:biopython_1.73'
bioconductor_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.90beta5'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'


########
# MAIN #
########

# map csid to our id
csid_to_indiv = {}
with open(sample_key, 'rt') as f:
    csv_reader = csv.reader(f)
    next(csv_reader)
    for row in csv_reader:
        csid_to_indiv[row[0]] = re.sub('\s', '_', row[1])
indiv_to_csid = {v: k for k, v in csid_to_indiv.items()}

# find all available csids
all_csid = list(set(snakemake.io.glob_wildcards(os.path.join(
    read_dir,
    '{cs_id1}/{cs_id2}_R{r}.fq.gz')).cs_id1))
all_indivs = [csid_to_indiv[x] for x in all_csid]

#########
# RULES #
#########

rule target:
    input:
        expand('output/050_variant-annotation/{vcf}_reheadered.vcf',
               vcf=['csd', 'goi']),
        expand('output/060_plink/goi.{suffix}',
               suffix=['map', 'ped', 'assoc']),
        expand('output/025_pileup/basecov/{indiv}.txt.gz',
               indiv=all_indivs),
        expand('output/070_regions/{indiv}_consensus.fa',
               indiv=all_indivs)

# extract aa sequences
rule extract_derived_cds:
    input:
        fa = honeybee_ref,
        regions = 'output/070_regions/regions.txt',
        vcf = 'output/050_variant-annotation/goi_reheadered.vcf.gz'
    output:
        'output/070_regions/{indiv}_consensus.fa'
    log:
        'output/logs/070_regions/{indiv}.log'
    singularity:
        samtools_container
    shell:
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '2> {log} '
        '| '
        'bcftools consensus '
        '-s Red1 '
        '-H 1 '
        '{input.vcf} '
        '> {output} '
        '2>> {log}'

rule extract_goi_cds_regions:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
        coding = 'output/050_variant-annotation/coding.Rds',
        assoc = 'output/060_plink/goi.assoc'
    output:
        regions = 'output/070_regions/regions.txt'
    log:
        'output/logs/070_regions/extract_goi_cds_regions.log'
    singularity:
        bioconductor_container
    script:
        'src/extract_goi_cds_regions.R'

# run association tests
rule plink_association:
    input:
        vcf = 'output/050_variant-annotation/goi.vcf',
        pheno = 'output/060_plink/pheno.txt'
    output:
        expand('output/060_plink/goi.{suffix}',
               suffix=['map', 'ped', 'assoc'])
    params:
        prefix = 'goi',
        wd = 'output/060_plink',
        vcf = lambda wildcards,input: resolve_path(input.vcf),
        pheno = lambda wildcards,input: resolve_path(input.pheno)
    log:
        resolve_path('output/logs/060_plink/plink_goi.log')
    singularity:
        plink_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'plink --recode '
        '--vcf {params.vcf} '
        '--allow-no-sex --allow-extra-chr --1 '
        '--out {params.prefix} '
        '--pheno {params.pheno} '
        '&>> {log} ; '
        'plink --assoc '
        '--allow-no-sex --allow-extra-chr --1 '
        '--out {params.prefix} '
        '--file {params.prefix} '
        '--pheno {params.pheno} '
        '&>> {log}'

rule recode_sample_key:
    input:
        key = 'data/sample_key.csv'
    output:
        pheno = 'output/060_plink/pheno.txt'
    log:
        'output/logs/060_plink/recode_sample_key.log'
    singularity:
        bioconductor_container
    script:
        'src/recode_sample_key.R'

# run variant annotation
rule fix_vcf:
    input:
        'output/050_variant-annotation/{vcf}.vcf'
    output:
        'output/050_variant-annotation/{vcf}_reheadered.vcf'
    singularity:
        bbduk_container
    shell:
        'grep -v "^##FILTER=All filters passed" {input} > {output}'


rule annotate_variants:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
        vcf = 'output/040_freebayes/variants_filtered.vcf.gz',
        tbi = 'output/040_freebayes/variants_filtered.vcf.gz.tbi',
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        goi = 'output/050_variant-annotation/dros_genes.csv'
    output:
        coding = 'output/050_variant-annotation/coding.Rds',
        csd = 'output/050_variant-annotation/csd.vcf',
        goi = 'output/050_variant-annotation/goi.vcf'
    log:
        'output/logs/050_variant-annotation/annotate_variants.log'
    singularity:
        bioconductor_container
    script:
        'src/annotate_variants.R'


# get genes of interest
rule get_apis_genes:
    input:
        'data/dros_eye_genes.txt'
    output:
        'output/050_variant-annotation/dros_genes.csv'
    log:
        'output/logs/050_variant-annotation/get_apis_genes.log'
    singularity:
        biopython_container
    script:
        'src/get_apis_genes.py'


# run freebayes
rule filter_vcf:
    input:
        'output/040_freebayes/variants.vcf'
    output:
        'output/040_freebayes/variants_filtered.vcf'
    params:
        filter = 'QUAL > 20'
    log:
        'output/logs/040_freebayes/freebayes_filter.log'
    singularity:
        vcflib_container
    shell:
        'vcffilter -f \'{params.filter}\' {input} > {output} 2> {log}'


rule freebayes:
    input:
        bam = expand('output/030_process-aln/{indiv}_marked.bam',
                     indiv=all_indivs),
        bai = expand('output/030_process-aln/{indiv}_marked.bam.bai',
                     indiv=all_indivs),
        fa = honeybee_ref
    output:
        vcf = 'output/040_freebayes/variants.vcf'
    params:
        ploidy = '2'
    log:
        'output/logs/040_freebayes/freebayes.log'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        # '{params.region} '
        '--ploidy {params.ploidy} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'


# process samfiles
rule index_bamfile:
    input:
        'output/030_process-aln/{indiv}_marked.bam'
    output:
        'output/030_process-aln/{indiv}_marked.bam.bai'
    log:
        'output/logs/030_process-aln/{indiv}_index.log',
    threads:
        2
    singularity:
        samtools_container
    shell:
        'samtools index -@ {threads} {input} 2> {log}'


rule markdup:
    input:
        'output/020_bwa/{indiv}.sam'
    output:
        sorted = temp('output/030_process-aln/{indiv}_sorted.bam'),
        marked = 'output/030_process-aln/{indiv}_marked.bam'
    threads:
        multiprocessing.cpu_count()
    log:
        s = 'output/logs/030_process-aln/{indiv}_sort.log',
        m = 'output/logs/030_process-aln/{indiv}_markdup.log'
    singularity:
        samtools_container
    shell:
        'samtools fixmate '
        '-m '
        '-O BAM '
        '-@ {threads} '
        '{input} '
        '- '
        '2> {log.s} '
        '| '
        'samtools sort '
        '-o {output.sorted} '
        '-O BAM '
        '-l 0 '
        '-@ {threads} '
        '- '
        '2>> {log.s} '
        '; '
        'samtools markdup '
        '-@ {threads} '
        '-s '
        '{output.sorted} '
        '{output.marked} '
        '2> {log.m}'

# calculate coverage from bamfiles
rule pileup:
    input:
        sam = 'output/020_bwa/{indiv}.sam',
        fa = 'output/000_data/Amel_HAv3.1_reheadered.fa'
    output:
        covstats = 'output/025_pileup/covstats/{indiv}.txt',
        hist = 'output/025_pileup/hist/{indiv}.txt',
        basecov = temp('output/025_pileup/basecov/{indiv}.txt'),
        bincov = 'output/025_pileup/bincov/{indiv}.txt',
        normcov = 'output/025_pileup/normcov/{indiv}.txt',
        normcovo = 'output/025_pileup/normcovo/{indiv}.txt',
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/025_pileup/{indiv}.log'
    singularity:
        bbduk_container
    shell:
        'pileup.sh '
        'in={input.sam} '
        'ref={input.fa} '
        'out={output.covstats} '
        'hist={output.hist} '
        'basecov={output.basecov} '
        'bincov={output.bincov} '
        'normcov={output.normcov} '
        'normcovo={output.normcovo} '
        'secondary=f '
        '2> {log}'

rule reheader_fasta:
    input:
        honeybee_ref
    output:
        'output/000_data/Amel_HAv3.1_reheadered.fa'
    threads:
        1
    log:
        'output/logs/000_data/reheader_fasta.log'
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input} '
        'out={output} '
        'trimreaddescription=t '
        '2> {log}'

# map individuals
rule bwa:
    input:
        fq = 'output/010_trim-decon/{indiv}/{indiv}.fq.gz',
        index = expand('output/020_bwa/honeybee_ref.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/020_bwa/{indiv}.sam'
    params:
        prefix = 'output/020_bwa/honeybee_ref.fasta',
        rg = '\'@RG\\tID:{indiv}\\tSM:{indiv}\''
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/020_bwa/{indiv}.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '-p '
        '-R {params.rg} '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


# prepare ref
rule index:
    input:
        honeybee_ref
    output:
        expand('output/020_bwa/honeybee_ref.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/020_bwa/honeybee_ref.fasta'
    threads:
        1
    log:
        'output/logs/020_bwa/index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'


# process reads
rule trim_decon:
    input:
        unpack(map_indiv_to_raw_reads)
    output:
        fq = 'output/010_trim-decon/{indiv}/{indiv}.fq.gz',
        f_stats = 'output/010_trim-decon/{indiv}/filter-stats.txt',
        t_stats = 'output/010_trim-decon/{indiv}/trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/{indiv}_filter.log',
        trim = 'output/logs/010_trim-decon/{indiv}_trim.log',
        repair1 = 'output/logs/010_trim-decon/{indiv}_repair1.log',
        repair2 = 'output/logs/010_trim-decon/{indiv}_repair2.log'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '2> {log.repair1} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq '
        'out={output.fq} '
        '2> {log.repair2} '

# generic compression rule
rule compress:
    input:
        '{folder}/{file}.{ext}'
    output:
        '{folder}/{file}.{ext}.gz'
    wildcard_constraints:
        ext = '(?!vcf)'     # not vcf
    singularity:
        pigz_container
    shell:
        'pigz '
        '--best '
        '--keep '
        '--stdout '
        '{input} '
        '> {output}'

# generic index rule
rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    singularity:
        samtools_container
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'
