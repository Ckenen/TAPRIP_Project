# Snakemake workflow of circRNA-seq datasets

To run the workflow:

    snakemake -s 1_SnakePrepare.smk -np
    snakemake -s 2_SnakeMapping.smk -np
    snakemake -s 3_SnakeExpression.smk -np
    snakemake -s 4_SnakeTracks.smk -np
    snakemake -s 5_SnakeCircExplorer -np
    snakemake -s 5_SnakeTracks.smk -np
    snakemake -s 6_SnakeCIRI.smk -np
    sh 7_RunCIRI.sh