'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule nanopolish_index:
    input:
        fast5='results/FAST5/{sample}',
        fasta='results/FASTA/{sample}.fasta',
    output:
        index=temp('results/FASTA/{sample}.fasta.index'),
        fai=temp('results/FASTA/{sample}.fasta.index.fai'),
        gzi=temp('results/FASTA/{sample}.fasta.index.gzi'),
        readdb=temp('results/FASTA/{sample}.fasta.index.readdb')
    conda: 'envs/nanopolish.yml'
    shell: 'nanopolish index -d {input}'

rule get_fasta:
    input: 'results/bambu/BAM/{sample}.bam'
    output: 'results/FASTA/{sample}.fasta'
    conda: 'envs/samtools.yml'
    shell: 'samtools fasta {input} > {output}'

rule nanopolish_eventalign_tx:
    input:
        fast5='results/FAST5/{sample}',
        fastq='results/FASTQ/{sample}.fq',
        fasta='results/FASTA/{sample}.fasta',
        bam='results/bambu/BAM/{sample}.bam',
        bai='results/bambu/BAM/{sample}.bam.bai',
        transcriptome='results/bambu/sqanti_filter/extended_annotations.filtered.fasta',
        index='results/FASTA/{sample}.fasta.index',
        fai='results/FASTA/{sample}.fasta.index.fai',
        gzi='results/FASTA/{sample}.fasta.index.gzi',
        readdb='results/FASTA/{sample}.fasta.index.readdb',
        summary=config['summaries_dir']+'/{sample}'
    output: temp('results/eventalign/{sample}.txt')
    threads: workflow.cores/4
    conda: 'envs/nanopolish.yml'
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
