# annotation_smk
This is a simple process of genome structure annotation, which can be parallelized and automated. <br>

## Install
```$ git clone https://github.com/yanhui-k/annotation_smk.git``` <br>

## Getting Started
You can use the accompanying `environment.yaml` to creat a general conda envirment <br>
```$ mamba env create -f environment.yaml``` <br>

Activate the environment <br>
```$ conda activate annotation``` <br>

Use the following code to modify the path of the input file <br>
```$ vim config/config.yaml``` <br>

To start the main pipeline, type in: <br>
```$ nohup snakemake --cluster "bsub -o log/output2 -e log/error2 -q Q104C512G_X4 -m yi02" -j 10 -p --use-conda &``` <br>

If you want to terminate a running pipline, please type in: <br>
```$ killall -TERM snakemake``` <br>
