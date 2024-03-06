'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule samtools_fastq:
    input: config['input_dir']+'/uBAM/{sample}.bam'
    output: temp('results/FASTQ/{sample}.fq')
    conda: 'envs/samtools.yml'
    shell: "samtools fastq -T '*' {input} > {output}"

rule minimap2_genome:
    input:
        genome=config['input_dir']+'/hg38/hg38.fa',
        fastq='results/FASTQ/{sample}.fq'
    output: temp('results/BAM/{sample}.sam')
    threads: workflow.cores/2
    conda: 'envs/minimap2.yml'
    shell:
        '''
        minimap2 \
            -ax splice \
            -uf \
            -k 14 \
            -t {threads} \
            -y \
            --secondary=no \
            {input} > {output}
        '''

rule minimap2_tx:
    input:
        tx='results/bambu/sqanti_filter/extended_annotations.filtered.fasta',
        fastq='results/FASTQ/{sample}.fq'
    output: temp('results/bambu/BAM/{sample}.sam')
    threads: workflow.cores/2
    conda: 'envs/minimap2.yml'
    shell:
        '''
        minimap2 \
            -ax map-ont \
            -uf \
            -t {threads} \
            --secondary=no \
            {input} > {output}
        '''

rule samtools_sort_genome:
    input: 'results/BAM/{sample}.sam'
    output: protected('results/BAM/{sample}.bam')
    threads: workflow.cores/4
    conda: 'envs/samtools.yml'
    shell:'samtools sort {input} -@ {threads} -o {output}'

rule samtools_sort_tx:
    input: 'results/bambu/BAM/{sample}.sam'
    output: protected('results/bambu/BAM/{sample}.bam')
    threads: workflow.cores/4
    conda: 'envs/samtools.yml'
    shell:'samtools sort {input} -@ {threads} -o {output}'

rule samtools_index_genome:
    input: 'results/BAM/{sample}.bam',
    output: protected('results/BAM/{sample}.bam.bai')
    conda: 'envs/samtools.yml'
    shell: 'samtools index {input}'

rule samtools_index_tx:
    input: 'results/bambu/BAM/{sample}.bam',
    output: protected('results/bambu/BAM/{sample}.bam.bai')
    conda: 'envs/samtools.yml'
    shell: 'samtools index {input}'
