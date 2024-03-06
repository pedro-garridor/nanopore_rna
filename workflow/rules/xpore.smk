'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule xpore_dataprep:
    input:
        eventalign='results/eventalign/{sample}.txt',
        gtf='results/bambu/sqanti_filter/extended_annotations.filtered.gtf',
        fasta='results/bambu/sqanti_filter/extended_annotations.filtered.fasta'
    output: temp(directory('results/xpore/dataprep/{sample}'))
    threads: workflow.cores/4
    conda: 'envs/xpore.yml'
    shell:
        '''
        xpore dataprep \
            --eventalign {input.eventalign} \
            --gtf_or_gff {input.gtf} \
            --transcript_fasta {input.fasta} \
            --out_dir {output} \
            --n_processes {threads}
        '''

rule xpore_diffmod_postprocessing:
    input:
        samples=expand('results/xpore/dataprep/{sample}', sample=SAMPLES),
        config=config['xpore_diffmod']
    output: protected(directory('results/xpore/diffmod'))
    threads: workflow.cores/2
    conda: '/home/pgarrido/mamba_envs_yaml/xpore.yml'
    shell:
        '''
        xpore diffmod \
            --config {input.config} \
            --n_processes {threads}
        xpore postprocessing \
            --diffmod_dir {output}
        '''