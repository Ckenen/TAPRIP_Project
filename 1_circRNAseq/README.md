# Snakemake workflow of circRNA-seq datasets

To run the workflow:

    snakemake -s 1_SnakePrepare.smk -np
    snakemake -s 2_SnakeCircExplorer -np
    snakemake -s 3_SnakeCIRI.smk -np
    sh 4_PostCIRI.sh
    snakemake -s 5_SnakeTracks.smk -np
    