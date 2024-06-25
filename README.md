# Clostrioides Dificiles ncRNA target prediction

Repository for the C.difficiles target prediction with intaRNA and MAPS analysis.


## Usage 

### intaRNA


1 - clone the repo.

2 - rewrite the path/filename of `inta_smk_config.yml`, the config file for the snakefile, to your own environment.

3 - run `snakemake -c2 -s inta_smk.smk --configfile inta_smk_config.yml --use-conda`


### MAPS


1 - clone the repo

2 - download of the template_script_DESeq2_CL.r of the the SARTools package in the 07_Fletcher_and_Pruss_analyses folder (needed for the differential analysis step)

3 - rewrite the path/filename of `ftp_fQC_bwt2_ftCounts_DEseq2_annot.yml`, the config file for the snakefile, to your own environment.

4 - run `snakemake -c2 -s ftp_fQC_bwt2_ftCounts_DEseq2_annot.smk --configfile ftp_fQC_bwt2_ftCounts_DEseq2_annot.yml --use-conda`


### Circos plot

1 - clone the repo

2 - create circlizeR_ce, a conda environment with "r-circlize" (https://github.com/jokergoo/circlize?tab=readme-ov-file) package using the command `conda env create -f condaEnv_rCirclize.yml `

3 - activate environment `circlizeR_ce`

2 - run `python3 linkfiles_MapsInta.py [path_to_MAPS_output_file] [MAPS_padj_filter][MAPS_l2fc_filter] [path_to_intaRNA_output_file:passFilter_list] [intaRNA_interaction_energy_filter] [path_to_annotation_file] [ncRNA_id]`
(example: `python3 linkfiles_MapsInta.py 01820vsctrl.complete.txt 0.05 2 passFilter_list_01820.tsv -10 CD630_CDS-ncRNA.gff 01820`)

3 - run `Rscript circlize_MapsInta.r [ncRNA_id][working_directory_path][MAPS_padj_filter][MAPS_l2fc_filter][intaRNA_interaction_energy_filter]`
(example: `Rscript circlize_MapsInta.r 01820 $(pwd)0.05 2 -10`)



