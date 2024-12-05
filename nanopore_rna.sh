#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "No arguments supplied."
    echo "Use 'singularity run-help nanopore_rna.sif' to get help."
    exit 2
fi

while [ $# -gt 0 ]; do
    case $1 in 
        -i|--input)
            INPUT="$2"
            shift
            shift
            ;;
        -o|--output)
            OUTDIR="$2"
            shift
            shift
            ;;
        -r|--reference)
            REF="$2"
            shift
            shift
            ;;
        -g|--gtf)
            REF_TX="$2"
            shift
            shift
            ;;
        -d|--diffmod)
            DIFFMOD="$2"
            shift
            shift
            ;;
        -n|--dryrun)
            DRYRUN=1
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            shift
            ;;
        -v|--version)
            exit
            shift
            ;;
        -h|--help)
            singularity run-help nanopore_rna.sif
            exit 0
            ;;
        -*|--*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            echo "Unknown argument $1"
            exit 1
            ;;
    esac
done

if [ -z $INPUT ] || [ -z $OUTDIR ] || [ -z $REF ] || [ -z $REF_TX ] || [ -z $DIFFMOD ]; then
    echo "Arguments -i, -o , -r, -g and -d are required."
    echo "Use 'bash nanopore_rna.sh -h' to get help."
    exit 2
fi

if [ $DRYRUN ]; then
    mkdir -p $OUTDIR
    singularity run \
        -B $INPUT,$OUTDIR,$REF,$REF_TX,$DIFFMOD \
        nanopore_rna.sif \
        -n \
        -t $THREADS \
        -i $INPUT \
        -o $OUTDIR \
        -r $REF \
        -g $REF_TX \
        -d $DIFFMOD 
else
    mkdir -p $OUTDIR
    singularity run \
        -B $INPUT,$OUTDIR,$REF,$REF_TX,$DIFFMOD \
        nanopore_rna.sif \
        -t $THREADS \
        -i $INPUT \
        -o $OUTDIR \
        -r $REF \
        -g $REF_TX \
        -d $DIFFMOD 
fi
