include: 'config.py'

from os.path import isfile

rule init:
    output:
        'metadata/sample_list.txt',
        'metadata/exp_info.txt'
    shell:
        source_r('R', 'metadata.R')

SAMPLES_FILE = 'metadata/sample_list.txt'
SAMPLES = []
if isfile(SAMPLES_FILE):
    with open('metadata/sample_list.txt', 'r') as input:
        for line in input:
            SAMPLES.append(line.strip())
else:
    print('Warning: the sample_list.txt is missing or init script has not yet been run. Running now.')
    rule:
       input: rules.init.output

rule all:
    input:
        'metadata/sample_list.txt',
        'metadata/exp_info.txt',
        expand('data/{id}.fastq.gz', id = SAMPLES),
        expand('results/kallisto/{id}/abundance.h5', id = SAMPLES),
        'results/sleuth_kal/L2KD_alr_actb.rds',
        'results/sleuth_kal/L2KD_alr_hprt1.rds',
        'results/sleuth_kal/L2KD_alr_both.rds',
        'results/sleuth_kal/L2KD_alr_actb_gene_aggregate_allResults.txt',
        'results/sleuth_kal/L2KD_alr_hprt1_gene_aggregate_allResults.txt',
        'results/sleuth_kal/L2KD_alr_both_gene_aggregate_allResults.txt',
        'results/go_analysis/L2KD_kal_alr_hprt1_sig_go_res.txt',
        'results/go_analysis/L2KD_kal_alr_actb_sig_go_res.txt',
        'results/go_analysis/L2KD_kal_alr_both_sig_go_res.txt'


def sample_input(wildcards):
    id = wildcards['id']
    return expand('data/{id}.fastq.gz', id = id)

rule kallisto:
    input:
        sample_input
    params:
        dir = 'results/kallisto/{id}'
    output:
        'results/kallisto/{id}/abundance.h5'
    threads: 5
    shell:
        '{UPDATED_PATH} '
        'kallisto quant'
        ' -i "{KALLISTO_INDEX}"'
        ' --rf-stranded'
        ' -b 100'
        ' --single'
        ' --bias'
        ' -l 150'
        ' -s 20'
        ' -o {params.dir}'
        ' -t {threads}'
        ' {input[0]}'

def get_sleuth_kal(wildcards):
    return expand('results/kallisto/{id}/abundance.h5', id = SAMPLES)

rule sleuth_kal:
    input:
        get_sleuth_kal
    output:
        'results/sleuth_kal/L2KD_alr_actb.rds',
        'results/sleuth_kal/L2KD_alr_hprt1.rds',
        'results/sleuth_kal/L2KD_alr_both.rds'
    shell:
        source_r('R', 'sleuth_kal.R')

rule sleuth_agg:
    input:
        'results/sleuth_{tool}/L2KD_alr_{gene}.rds'
    output:
        'results/sleuth_{tool}/L2KD_alr_{gene}_gene_aggregate_allResults.txt'
    shell:
        source_r('R', 'sleuth_aggregate.R') + ' ../{input[0]} ../results/sleuth_{wildcards.tool}/L2KD_alr_{wildcards.gene}'

rule go_analysis:
    input:
        'results/sleuth_{tool}/L2KD_alr_{gene}_gene_aggregate_allResults.txt'
    output:
        'results/go_analysis/L2KD_{tool}_alr_{gene}_sig_go_res.txt'
    shell:
        source_r('R', 'go_analysis.R') + ' ../{input[0]} ../results/go_analysis/L2KD_{wildcards.tool}_alr_{wildcards.gene}'
