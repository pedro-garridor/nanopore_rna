'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule nanopolish_polya:
    input:
        fast5=config['out_dir']+'/FAST5/{sample}',
        fastq=config['out_dir']+'/FASTQ/{sample}.fq',
        bam=config['out_dir']+'/bambu/BAM/{sample}.bam',
        bai=config['out_dir']+'/bambu/BAM/{sample}.bam.bai',
        transcriptome=config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.fasta',
        index=config['out_dir']+'/FASTA/{sample}.fasta.index',
        fai=config['out_dir']+'/FASTA/{sample}.fasta.index.fai',
        gzi=config['out_dir']+'/FASTA/{sample}.fasta.index.gzi',
        readdb=config['out_dir']+'/FASTA/{sample}.fasta.index.readdb'
    output: protected(config['out_dir']+'/polyA/{sample}.tsv')
    threads: workflow.cores/4
    conda: '../envs/nanopolish.yml'
    shell:
        '''
        nanopolish polya \
            -t {threads} \
            -r {input.fastq} \
            -b {input.bam} \
            -g {input.transcriptome} > {output}
        '''