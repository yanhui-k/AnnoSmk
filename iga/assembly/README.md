[TOC]

# I. 基因组Contig的三代组装
## 1.1组装

### 1.1.1Raw reads 组装
> 使用未纠错的下机subreads进行组装（包含纠错步骤）

#### falcon组装命令
默认情况：

python -m iga.assembly.assemble falcon  Zn.fasta.gz  464000000

修改参数：

python -m iga.assembly.assemble falcon  Zn.fasta.gz  464000000 --etc '[General]length_cutoff=35000;[General]length_cutoff_pr=34000'

**结果**： 

* p_ctg.fasta： 组装好的contig
* preads4falcon.fasta：纠错的subreads

#### Canu组装命令

python -m iga.assembly.assemble canu  Zn.fasta.gz  464000000

**结果**：

* Eu.contigs.fasta：组装好的contig
* Eu.correctedReads.fasta.gz：纠错的subreads

### 1.1.2 Corrected reads 组装

> 使用纠错好的reads进行组装，跳过纠错步骤

#### wtdbg组装命令

python -m iga.assembly.assemble wtdbg corrected.fasta 123m

**结果**: assembly.fasta

#### flye组装命令：

python -m iga.assembly.assemble flye corrected.fasta 123m

**结果**：Eu_wtdbg.cns.fa

### 1.1.3 HiFi reads 组装

#### hifiasm组装命令

python -m iga.assembly.assemble hifiasm --threads 40 output.fasta.gz

**结果**：workdir_hifiasm_output/output.p_ctg.fa

---

## 1.2 组装后处理

### 1.2.1 抛光（组装后纠错）

> 用三代和二代数据对基因组进行再次纠错

先用三代数据纠错，需要用下机未纠错的reads

python -m iga.assembly.polish gcpp contig pacbio_subreads.bam

**结果**：workdir_GCPP/output.fasta

再用二代纠错

python -m iga.assembly.polish nextpolish contig 'left.fq right.fq'

**结果**： workdir_nextpolish_canu_polished/canu_polished.1.fa

### 1.2.2 组装评估

#### 算N50

N50v1.pl contig > contig.N50

#### 算BUSCO的

python -m iga.assembly.assess busco contig

**结果文件**：`flye_polished.fasta.busco.embryophyta.v4.1.2/short_summary.specific.embryophyta_odb10.flye_polished.fasta.busco.embryophyta.v4.1.2.txt`

#### 算LAI的

python -m iga.assembly.assess lai contig.fasta

**结果文件**：`workdir_LAI_canu_polished.fasta/canu_polished.fasta.mod.out.LAI`


---

## 1.3 其它脚本

### 1.3.1 格式转换

#### bam转fasta

把subreads.bam转换成组装软件需要的fasta文件，目录需要同时存在 zn.bam.pbi

python -m iga.assembly.assemble bam2fastq zn.bam

**结果文件**： zn.fasta.gz

# II. HiC 挂载

## juicer流程

```
python -m iga.assembly.hic juicer_pipe Zn_falcon.2.ctg.fasta "hic1.fq hic2.fq"
```
