'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule xpore_dataprep:
    input:
        eventalign=config['out_dir']+'/eventalign/{sample}.txt',
        gtf=config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.gtf',
        fasta=config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.fasta'
    output: temp(directory(config['out_dir']+'/xpore/dataprep/{sample}'))
    threads: workflow.cores/4
    conda: '../envs/xpore.yml'
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
        samples=expand(config['out_dir']+'/xpore/dataprep/{sample}', sample=SAMPLES),
        config=config['xpore_diffmod']
    output: protected(directory(config['out_dir']+'/xpore/diffmod'))
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