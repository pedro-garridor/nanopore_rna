# nanopore_rna
Identification of novel transcripts and epitranscriptome profiling for nanopore RNA-Seq experiments

`nanopore_rna` is a transcriptome and epitranscriptome pipeline for nanopore RNA-Seq experiments.

This implementation allows the analysis of nanopore RNA-Seq raw data (from FASTQ and POD5) either on local workstations, as well as in shared servers or clusters, **without requiring admin permissions** to run it.

## Get started

To use this pipeline, you will need to install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), [snakemake](https://snakemake.readthedocs.io/en/stable/) and [SingularityCE](https://github.com/sylabs/singularity/releases/). To do this:

**1.** Download in your home the [latest Miniforge distribution](https://github.com/conda-forge/miniforge/releases/latest/) selecting your OS and architecture. For example, for Linux x64 it would be:
   
   ```
   wget https://github.com/conda-forge/miniforge/releases/download/<MINIFORGE_VERSION>/Mambaforge-<MINIFORGE_VERSION>-Linux-x86_64.sh
   ```
   

Changing `<MINIFORGE_VERSION>` with the latest version available.

   If you already have mamba installed in your account, go to step **2**.

**2.** Clone this repository in your system with:

    
    git clone https://github.com/pedro-garridor/nanopore_rna.git
    

**3.** For this step, install SingularityCE, move into `nanopore_rna` and build the pipeline image. You will need to do this last thing as a privileged user in your local computer. However, once you build it you will be able to copy the image to another computer and use it without sudo permissions:

    wget https://github.com/sylabs/singularity/releases/download/v4.2.1/singularity-ce_4.2.1-noble_amd64.deb
    sudo apt install -y singularity-ce_4.2.1-noble_amd64.deb
    rm singularity-ce_4.2.1-noble_amd64.deb
    cd nanopore_rna
    sudo singularity build nanopore_rna.sif nanopore_rna.def

**NOTE**: This is the only step where you will need sudo permissions. You can run it in your local computer and then copy the `nanopore_rna` folder to another computer, where you will be able to run the pipeline without admin privileges. 


**4.** Once you've done all of the above, you can run the pipeline with

    
    bash nanopore_rna.sh \
        -i <INPUT_FOLDER> \
        -o <OUTPUT_FOLDER> \
        -r <REFERENCE_GENOME> \
        -g <REFERENCE_TRANSCRIPTOME> \
        -d <DIFFMOD_YML> \
        -t <THREADS>
    

where:

- `<INPUT_FOLDER>`: directory where all MinKNOW result files are located. One folder per sample, named with the sample ID.
- `<OUTPUT_FOLDER>`: path where you want to get the results.
- `<REFERENCE_GENOME>`: path to your genome `.fa`. `.fai` index is required to be present on the same folder.
- `<REFERENCE_TRANSCRIPTOME>`: path to your reference transcriptome `.gtf`.
- `<DIFFMOD_YML>`: YML file with the grouping of your samples for xPore.
- `<THREADS>`: threads you want the pipeline to use.

As an example, the `<DIFFMOD_YML>` file needs to have a format like the following:

    data:
        GROUP_A:
            sample_a_1: <OUTPUT_FOLDER>/xpore/dataprep/sample_a_1_
            sample_a_2: <OUTPUT_FOLDER>/xpore/dataprep/sample_a_2
            sample_a_3: <OUTPUT_FOLDER>/xpore/dataprep/sample_a_3
        GROUP_B:
            sample_b_1: <OUTPUT_FOLDER>/xpore/dataprep/sample_b_1
            sample_b_2: <OUTPUT_FOLDER>/xpore/dataprep/sample_b_2
            sample_b_3: <OUTPUT_FOLDER>/xpore/dataprep/sample_b_3

    out: <OUTPUT_FOLDER>/xpore/diffmod

If you need more info on what this file is, or how to make it, take a look [here](https://xpore.readthedocs.io/en/latest/configuration.html).

## Help

Finally, if you need more information on how to run the pipeline, you can run:

    
    bash nanopore_rna.sh -h
    

## Citations

This pipeline is part of my PhD Thesis. If you're using it, please cite it :)

> P. Garrido-Rodríguez, “Aplicación de la bioinformática en la descripción y resolución de patologías hematológicas y mecanismos biológicos relacionados,” Ph.D dissertation, Universidad de Murcia, Spain, 2025. [https://digitum.um.es/digitum/handle/10201/152160](https://digitum.um.es/digitum/handle/10201/152160)
