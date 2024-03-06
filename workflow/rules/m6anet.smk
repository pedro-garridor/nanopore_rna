'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule m6anet_dataprep:
    input: 'results/eventalign/{sample}.txt'
    output: temp(directory('results/m6anet/dataprep/{sample}'))
    threads: workflow.cores/4
    conda: 'envs/m6anet.yml'
    shell:
        '''
        m6anet dataprep \
            --eventalign {input} \
            --out_dir {output} \
            --n_processes {threads}
        '''

rule m6anet_inference:
    input: 'results/m6anet/dataprep/{sample}'
    output: protected(directory('results/m6anet/inference/{sample}'))
    threads: workflow.cores/4
    conda: 'envs/m6anet.yml'
    shell:
        '''
        m6anet inference \
            --input_dir {input} \
            --out_dir {output} \
            --n_processes {threads}
        '''
