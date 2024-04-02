'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule fastq:
    input: config['input_dir']+'/{sample}/fastq_pass/'
    output: temp(config['out_dir']+'/FASTQ/{sample}.fq')
    shell: 'cat {input} | gunzip -c > {output}'

rule minimap2_genome:
    input:
        genome=config['genome_dir']+'/hg38.fa',
        fastq=config['out_dir']+'/FASTQ/{sample}.fq'
    output: temp(config['out_dir']+'/BAM/{sample}.sam')
    threads: workflow.cores/2
    conda: '../envs/minimap2.yml'
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
        tx=config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.fasta',
        fastq=config['out_dir']+'/FASTQ/{sample}.fq'
    output: temp(config['out_dir']+'/bambu/BAM/{sample}.sam')
    threads: workflow.cores/2
    conda: '../envs/minimap2.yml'
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
    input: config['out_dir']+'/BAM/{sample}.sam'
    output: protected(config['out_dir']+'/BAM/{sample}.bam')
    threads: workflow.cores/4
    conda: '../envs/samtools.yml'
    shell:'samtools sort {input} -@ {threads} -o {output}'

rule samtools_sort_tx:
    input: config['out_dir']+'/bambu/BAM/{sample}.sam'
    output: protected(config['out_dir']+'/bambu/BAM/{sample}.bam')
    threads: workflow.cores/4
    conda: '../envs/samtools.yml'
    shell:'samtools sort {input} -@ {threads} -o {output}'

rule samtools_index_genome:
    input: config['out_dir']+'/BAM/{sample}.bam',
    output: protected(config['out_dir']+'/BAM/{sample}.bam.bai')
    conda: '../envs/samtools.yml'
    shell: 'samtools index {input}'

rule samtools_index_tx:
    input: config['out_dir']+'/bambu/BAM/{sample}.bam',
    output: protected(config['out_dir']+'/bambu/BAM/{sample}.bam.bai')
    conda: '../envs/samtools.yml'
    shell: 'samtools index {input}'

rule bam_slam:
    input:
        bam=config['out_dir']+'/bambu/BAM/{sample}.bam',
        bai=config['out_dir']+'/bambu/BAM/{sample}.bam.bai'
    output:
        coverage_fraction=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_coverage_fraction.pdf'),
        csv=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_data.csv'),
        density=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_density.pdf'),
        read_accuracy=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_read_accuracy.pdf'),
        sec_alns=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_sec_alns.pdf'),
        stats=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_stats.csv'),
        tx_length_distr_per_tx=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_transcript_length_distribution_per_distinct_transcript.pdf'),
        tx_length_distr_per_read=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_transcript_length_distribution_per_read.pdf'),
        tx_level_data=protected(config['out_dir']+'/bambu/BAM/BamSlam/{sample}_transcript_level_data.csv')
    conda: '../envs/bamslam.yml'
    params:
        bamslam_script='workflow/scripts/BamSlam.R'
    shell: 'Rscript {params.bamslam_script} rna {input.bam} results/bambu/BAM/BamSlam/{wildcards.sample}'
    
