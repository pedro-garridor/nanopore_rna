'''
Nanopore RNA-Seq workflow
Copyright (C) 2024, Pedro Garrido Rodríguez
'''

input_table = pd.read_csv(config['manifest'], header=None)
IDS = input_table.iloc[:, 0].to_list()
SAMPLES = [re.sub(r'^.*_OM', 'OM', ID) for ID in IDS]

rule bambu:
    input:
        genome=config['genome_dir']+'/hg38.fa',
        ref_tx=config['ref_tx'],
        bam=expand('results/BAM/{sample}.bam', sample=SAMPLES),
        bai=expand('results/BAM/{sample}.bam.bai', sample=SAMPLES)
    output:
        counts_gene=protected('results/bambu/counts_gene.txt'),
        counts_tx=protected('results/bambu/counts_transcript.txt'),
        CPM_tx=protected('results/bambu/CPM_transcript.txt'),
        gtf=protected('results/bambu/extended_annotations.gtf'),
        full_length_tx=protected('results/bambu/fullLengthCounts_transcript.txt'),
        unique_tx=protected('results/bambu/uniqueCounts_transcript.txt'),
        rdata=protected('results/bambu/bambu.RData')
    threads: workflow.cores
    conda: 'envs/bambu.yml'
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
        bambu_gtf='results/bambu/extended_annotations.gtf',
        ref_tx=config['ref_tx'],
        genome=config['genome_dir']+'/hg38.fa'
    output:
        GMST=protected(directory('results/bambu/sqanti_qc/GMST')),
        RTS=protected(directory('results/bambu/sqanti_qc/RTS')),
        classification=protected('results/bambu/sqanti_qc/extended_annotations_classification.txt'),
        faa=protected('results/bambu/sqanti_qc/extended_annotations_corrected.faa'),
        fasta=protected('results/bambu/sqanti_qc/extended_annotations_corrected.fasta'),
        genepred=protected('results/bambu/sqanti_qc/extended_annotations_corrected.genePred'),
        gtf=protected('results/bambu/sqanti_qc/extended_annotations_corrected.gtf'),
        cds=protected('results/bambu/sqanti_qc/extended_annotations_corrected.gtf.cds.gff'),
        junctions=protected('results/bambu/sqanti_qc/extended_annotations_junctions.txt'),
        report=protected('results/bambu/sqanti_qc/extended_annotations_SQANTI3_report.html'),
        params=protected('results/bambu/sqanti_qc/extended_annotations.params.txt'),
        refann=protected('results/bambu/sqanti_qc/refAnnotation_extended_annotations.genePred')
    conda: 'envs/SQANTI3_env_fixed.yml'
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
        classification='results/bambu/sqanti_qc/extended_annotations_classification.txt',
        fasta='results/bambu/sqanti_qc/extended_annotations_corrected.fasta',
        faa='results/bambu/sqanti_qc/extended_annotations_corrected.faa',
        gtf='results/bambu/sqanti_qc/extended_annotations_corrected.gtf',
    output:
        filter_reasons=protected('results/bambu/sqanti_filter/extended_annotations_filtering_reasons.txt'),
        inclusion_list=protected('results/bambu/sqanti_filter/extended_annotations_inclusion-list.txt'),
        classification=protected('results/bambu/sqanti_filter/extended_annotations_RulesFilter_result_classification.txt'),
        faa=protected('results/bambu/sqanti_filter/extended_annotations.filtered.faa'),
        fasta=protected('results/bambu/sqanti_filter/extended_annotations.filtered.fasta'),
        gtf=protected('results/bambu/sqanti_filter/extended_annotations.filtered.gtf'),
        params=protected('results/bambu/sqanti_filter/extended_annotations.params.txt')
    conda: 'envs/SQANTI3_env_fixed.yml'
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
    