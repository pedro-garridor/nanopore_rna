'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

SAMPLES = glob_wildcards(config['input_dir']+'/{sample,[^/]+}').sample

rule bambu:
    input:
        genome=config['genome_dir']+'/hg38.fa',
        ref_tx=config['ref_tx'],
        bam=expand(config['out_dir']+'/BAM/{sample}.bam', sample=SAMPLES),
        bai=expand(config['out_dir']+'/BAM/{sample}.bam.bai', sample=SAMPLES)
    output:
        counts_gene=protected(config['out_dir']+'/bambu/counts_gene.txt'),
        counts_tx=protected(config['out_dir']+'/bambu/counts_transcript.txt'),
        CPM_tx=protected(config['out_dir']+'/bambu/CPM_transcript.txt'),
        gtf=protected(config['out_dir']+'/bambu/extended_annotations.gtf'),
        full_length_tx=protected(config['out_dir']+'/bambu/fullLengthCounts_transcript.txt'),
        unique_tx=protected(config['out_dir']+'/bambu/uniqueCounts_transcript.txt'),
        rdata=protected(config['out_dir']+'/bambu/bambu.RData')
    threads: workflow.cores
    conda: '../envs/bambu.yml'
    params:
        bambu_script='workflow/scripts/bambu.R'
    shell:
        '''
        Rscript {params.bambu_script} \
                {threads} \
                {input.genome} \
                {input.ref_tx} \
                {input.bam}
        '''

rule sqanti_qc:
    input: 
        bambu_gtf=config['out_dir']+'/bambu/extended_annotations.gtf',
        ref_tx=config['ref_tx'],
        genome=config['genome_dir']+'/hg38.fa'
    output:
        GMST=protected(directory(config['out_dir']+'/bambu/sqanti_qc/GMST')),
        RTS=protected(directory(config['out_dir']+'/bambu/sqanti_qc/RTS')),
        classification=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_classification.txt'),
        faa=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.faa'),
        fasta=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.fasta'),
        genepred=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.genePred'),
        gtf=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.gtf'),
        cds=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.gtf.cds.gff'),
        junctions=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_junctions.txt'),
        report=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations_SQANTI3_report.html'),
        params=protected(config['out_dir']+'/bambu/sqanti_qc/extended_annotations.params.txt'),
        refann=protected(config['out_dir']+'/bambu/sqanti_qc/refAnnotation_extended_annotations.genePred')
    conda: '../envs/SQANTI3_env_fixed.yml'
    params:
        sqanti_dir=config['sqanti_dir'],
        cdnacupcake_dir=config['cdnacupcake_dir']
    shell:
        '''
        export PYTHONPATH={params.cdnacupcake_dir}/
        export PYTHONPATH=$PYTHONPATH:{params.cdnacupcake_dir}/sequence/
        {params.sqanti_dir}/sqanti3_qc.py \
            {input.bambu_gtf} \
            {input.ref_tx} \
            {input.genome} \
            -d $(dirname {output.report})
        '''

rule sqanti_filter:
    input: 
        classification=config['out_dir']+'/bambu/sqanti_qc/extended_annotations_classification.txt',
        fasta=config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.fasta',
        faa=config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.faa',
        gtf=config['out_dir']+'/bambu/sqanti_qc/extended_annotations_corrected.gtf',
    output:
        filter_reasons=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations_filtering_reasons.txt'),
        inclusion_list=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations_inclusion-list.txt'),
        classification=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations_RulesFilter_result_classification.txt'),
        faa=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.faa'),
        fasta=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.fasta'),
        gtf=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations.filtered.gtf'),
        params=protected(config['out_dir']+'/bambu/sqanti_filter/extended_annotations.params.txt')
    conda: '../envs/SQANTI3_env_fixed.yml'
    params:
        sqanti_dir=config['sqanti_dir'],
        cdnacupcake_dir=config['cdnacupcake_dir']
    shell:
        '''
        export PYTHONPATH={params.cdnacupcake_dir}/
        export PYTHONPATH=$PYTHONPATH:{params.cdnacupcake_dir}/sequence/
        {params.sqanti_dir}/sqanti3_filter.py rules \
            {input.classification} \
            --isoforms  {input.fasta} \
            --gtf {input.gtf} \
            --faa {input.faa} \
            -d $(dirname {output.classification}) \
            --skip_report
        ls -lah $(dirname {output.classification})
        '''
    