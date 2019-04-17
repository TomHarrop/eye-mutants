#!/usr/bin/env python3

import csv
import multiprocessing
import re
import snakemake
import os

#############
# FUNCTIONS #
#############

def map_indiv_to_raw_reads(wildcards):
    cs_id = indiv_to_csid[wildcards.indiv]
    return({
        'r1': os.path.join(read_dir, f'{cs_id}/{cs_id}_R1.fq.gz'),
        'r2': os.path.join(read_dir, f'{cs_id}/{cs_id}_R2.fq.gz')
        })

###########
# GLOBALS #
###########

sample_key = 'data/sample_key.csv'
read_dir = 'data/raw'
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'


########
# MAIN #
########

# map csid to our id
csid_to_indiv = {}
with open(sample_key, 'rt') as f:
    csv_reader = csv.reader(f)
    next(csv_reader)
    for row in csv_reader:
        csid_to_indiv[row[0]] = re.sub("\s", "_", row[1])
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
        expand('output/010_trim-decon/{indiv}/reads.fq.gz',
               indiv=all_indivs)
rule rename:
    input:
        unpack(map_indiv_to_raw_reads)
    output:
        'output/001_test/{indiv}.fq.gz'
    shell:
        'cp {input} {output}'


rule trim_decon:
    input:
        unpack(map_indiv_to_raw_reads)
    output:
        fq = 'output/010_trim-decon/{indiv}/reads.fq.gz',
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
