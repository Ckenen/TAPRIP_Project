# Snakemake workflow of seCLIP-seq datasets

To run the workflow:

    snakemake -s 1_SnakePipeline.smk -np
    snakemake -s 2_SnakePeaks -np
    snakemake -s 3_SnakeTracks.smk -np
    snakemake -s 4_SnakeNormPeaks.smk -np
    snakemake -s 5_SnakeMergePeaks.smk -np

    