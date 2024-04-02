'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule get_fasta:
    input: config['out_dir']+'/bambu/BAM/{sample}.bam'
    output: config['out_dir']+'/FASTA/{sample}.fasta'
    conda: '../envs/samtools.yml'
    shell: 'samtools fasta {input} > {output}'

rule nanopolish_index:
    input:
        fast5=config['out_dir']+'/FAST5/{sample}',
        fasta=config['out_dir']+'/FASTA/{sample}.fasta',
    output:
        index=temp(config['out_dir']+'/FASTA/{sample}.fasta.index'),
        fai=temp(config['out_dir']+'/FASTA/{sample}.fasta.index.fai'),
        gzi=temp(config['out_dir']+'/FASTA/{sample}.fasta.index.gzi'),
        readdb=temp(config['out_dir']+'/FASTA/{sample}.fasta.index.readdb')
    conda: '../envs/nanopolish.yml'
    shell: 'nanopolish index -d {input}'

rule nanopolish_eventalign_tx:
    input:
        fast5=config['out_dir']+'/FAST5/{sample}',
        fastq=config['out_dir']+'/FASTQ/{sample}.fq',
        fasta=config['out_dir']+'/FASTA/{sample}.fasta',
        bam=config['out_dir']+'/bambu/BAM/{sample}.bam',
        bai=config['out_dir']+'/bambu/BAM/{sample}.bam.bai',
        transcriptome=config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.fasta',
        index=config['out_dir']+'/FASTA/{sample}.fasta.index',
        fai=config['out_dir']+'/FASTA/{sample}.fasta.index.fai',
        gzi=config['out_dir']+'/FASTA/{sample}.fasta.index.gzi',
        readdb=config['out_dir']+'/FASTA/{sample}.fasta.index.readdb',
        summary=config['input_dir']+'/{sample}'
    output: temp(config['out_dir']+'/eventalign/{sample}.txt')
    threads: workflow.cores/4
    conda: '../envs/nanopolish.yml'
    shell:
        '''
        nanopolish eventalign \
            --reads {input.fasta} \
            --bam {input.bam} \
            --genome {input.transcriptome} \
            --scale-events \
            --signal-index \
            --summary {input.summary}/sequencing_summary*.txt \
            --threads {threads} > {output}
        '''
