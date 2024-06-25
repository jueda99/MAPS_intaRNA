# Clostrioides Dificiles ncRNA target prediction

Repository for the C.difficiles target prediction with intaRNA and MAPS analysis.

## intaRNA

### Usage

1 - clone the repo.
2 - rewrite the path/filename of `inta_smk_config.yml`, the config file for the snakefile, to your own environment.
3 - run `snakemake -c2 -s inta_smk.smk --configfile inta_smk_config.yml --use-conda`


## MAPS

1 - clone the repo
2 - download of the template_script_DESeq2_CL.r of the the SARTools package in the 07_Fletcher_and_Pruss_analyses folder (needed for the differential analysis step)
3 - rewrite the path/filename of `ftp_fQC_bwt2_ftCounts_DEseq2_annot.yml`, the config file for the snakefile, to your own environment.
4 - run `snakemake -c2 -s ftp_fQC_bwt2_ftCounts_DEseq2_annot.smk --configfile ftp_fQC_bwt2_ftCounts_DEseq2_annot.yml --use-conda`

