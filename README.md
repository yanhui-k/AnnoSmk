Install


git clone https://github.com/yanhui-k/annotation_smk.git


conda env create -f environment.yml


conda activate annotation


mamba install -c conda-forge -c bioconda snakemake


snakemake --cluster "bsub -o log/output1 -e log/error1 -q Q104C512G_X4 -m yi02" -j 10 -p --use-conda
