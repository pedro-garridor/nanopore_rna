'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

rule get_fast5:
    input: config['input_dir']+'/{sample}/'
    output: temp(directory(config['out_dir']+'/FAST5/{sample}'))
    threads: workflow.cores/4
    conda: '../envs/pod5.yml'
    shell:
        '''
        pod5 convert to_fast5 \
            -o {output} \
            -t {threads} \
            {input}
        '''
