# RepExpress
<<<<<<< HEAD

A pipeline for quantifying and exploring transposable element expression

Developed by Peter A. Stockwell and Chiemi F. Lynch-Sutherland.
Department of Pathology, Otago Medical School, University of Otago,
Dunedin, New Zealand.

----------------------------------------------------------------------

## Contents of Package

README.md - this file

RepExpress_processing.pdf - documentation.

Main processing scripts:

    generate_combine_gtfs.sh - to build STAR index and build required GTFs

    map_count_reads.sh - run STAR, featureCounts & stringtie on each sample

    compare_express.sh - combine expression counts for a series of similar samples

    contrast_express.sh - take combined expression counts for two different series of samples

Example_scripts:  A series of parameter files (examples) showing the required parameters for each process.

    Readme.txt - a description of each example parameter file

    basic_params.sh - system-specific details required for all other processing

    map_params.sh - mapping details for map_count_reads.sh execution

    combine_params.sh - details for the compare_express.sh run

    genloc_file_list.txt - file with names of genloc files needed for compare_express.sh

    contrast_params.sh - details for the contrast_express.sh run
=======
A pipeline for quantifying and exploring transposable element expression
>>>>>>> d18818786de4bd04803bd882c33fc34b6800c5ab
