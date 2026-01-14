tsv_E=Sample_stdata.tsv
path_to_annotation_file=/gpfs/ycga/project/fan/yl2224/Alignment/Dropseq_Alignment_References/mm10/mm10.gtf

convertEnsemblToNames.py $tsv_E --annotation $path_to_annotation_file --output Sample_exp_matrix.tsv