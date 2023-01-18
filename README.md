# Code for the manuscript "Machine learning reveals STAT motifs as predictors for GR-mediated gene repression"

TODO: Add citation once paper is accepted

# Table of contents
1. [Software dependencies](#software)  
    1. [Conda](#conda)
    2. [MEME suite](#meme)
    3. [Activity-by-Contact (ABC) model](#abc)
2. [Required inputs](#inputs)
3. [Snakemake workflow](#snakemake)
    1. [Workflow execution](#execution)
    2. [Workflow configuration](#config)
    3. [Workflow visualization](#visualization)
    4. [Troubleshooting](#troubleshooting)
4. [Acknowledgments](#acknowledgments)

<a id="software"></a>

# Software dependencies 

<a id="conda"></a>

## Conda

In case you do not have yet a conda installation, you can set up one following instructions in <https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation>. Once it is finished, update it if necessary.

In order to setup the main workflow environment (called "gr_env"), a wrapper `setup_env_workflow.sh` is provided. It uses the conda drop-in replacement **mamba** to resolve dependencies.
Make sure you have mamba installed in your base environment. If necessary install it by running:

```
conda install -n base conda-forge mamba
```

Then you can create the environment *gr_env* by running:

```
sh setup_env_workflow.sh
```

Several rules require packages that are not available within *gr_env*.  
The required packages are specified through .yml files and the path provided to the respective Snakemake rule definitions. 
The .yml files are located within `src/conda/` and are used by snakemake to create the envs on the fly.  
The snakemake call within `run_workflow.sh` has the flag `--use-conda` set, so that creatign environments from the provided yml files is the default behavior.

<a id="meme"></a>

## MEME suite

In order to run commands from the MEME suite (such as STREME or FIMO), the software needs to be available locally (unfortunately the conda package of meme is buggy and cannot be used).  
Depending on the user software environment, this might require the installation of perl packages. Since we lacked sudo privelege, we accomplished this through perlbrew.  
In detail, the process to install the tool was

```
\curl -L https://install.perlbrew.pl | bash
~/perl5/perlbrew/bin/perlbrew install perl-5.34.0

perlbrew use perl-5.34.0
# install perl dependencies through cpanm (the perl package manager)
cpanm HTML::Template
cpanm HTML::TreeBuilder
cpanm XML::Parser::Expat
cpanm File::Which
cpanm JSON
cpanm XML::Simple
cpanm Sys::Info

version=5.4.1
wget http://meme-suite.org/meme-software/$version/meme-$version.tar.gz
tar zxf meme-$version.tar.gz
cd meme-$version

./configure --prefix=$HOME/software/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install

# add to .bash_profile:
export PATH=$HOME/software/meme/bin:$HOME/software/meme/libexec/meme-5.4.1:$PATH

```
After that, snakemake rules containing meme tools can be run and will be executed from within a memes environment that gets created on the fly.

<a id="abc"></a>

## Activity-by-Contact (ABC) model

The ABC repository was cloned from https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction on July 6th 2021 (SHA of last commit was 74ee57e7baff7c5a584d53c36c22681620a3e491) and symlinked into `src/scripts/abcmodel/`.
Besides adding bedtools as dependencies to the conda env, no changes were made to the ABC code. 

<a id="inputs"></a>

# Required inputs

Data files are read from `/data/current/` and output files written to `/results/current/`.  
It is recommended to set these directories up as symbolic link to a data storage location that is timestamped.  

Within `data/current/` the workflow expects to find the following subdirectories and input files:

- genome/
    - GRCm38.primary_assembly.standchr.fa

- chipseq/

    - H3K27ac/fastq/
        - GAR0814_S7_L002_R1_001.fastq.gz, GAR0814_S7_L002_R2_001.fastq.gz
        - GAR0815_S8_L002_R1_001.fastq.gz, GAR0815_S8_L002_R2_001.fastq.gz  
        - GAR0816_S9_L002_R1_001.fastq.gz, GAR0816_S9_L002_R2_001.fastq.gz  
        - GAR0823_S10_L002_R1_001.fastq.gz, GAR0823_S10_L002_R2_001.fastq.gz  

    - GR_2020/fastq_merged/
        - Sample_MUC20387/MUC20387_S6_R1_001.fastq.gz, MUC20387_S6_R2_001.fastq.gz
        - Sample_MUC20388/MUC20388_S7_R1_001.fastq.gz, MUC20388_S7_R2_001.fastq.gz

    - GR_2020/input_ctrls/
        - GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R1_001.fastq.gz, GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R2_001.fastq.gz  
        - Sample_MUC9117/MUC9117_R1_merged.fastq.gz, Sample_MUC9117/MUC9117_R2_merged.fastq.gz
        - Sample_MUC9118/MUC9118_R1_merged.fastq.gz, Sample_MUC9118/MUC9118_R2_merged.fastq.gz
        - Sample_GAR1531/GAR1531_S13_L002_R1_001.fastq.gz, Sample_GAR1531/GAR1531_S13_L002_R2_001.fastq.gz

- atacseq/  

    - bam/processed/DexLPS/*.bam
    - bam/processed/LPS/*.bam
    - data/current/atacseq/footprints/DexLPSvLPS_diff_footprint/differential_statistics.txt   
    ATACseq raw data was processed by an internal pipeline that handles adapter removal, read alignment, deduplication and the footprinting analysis.

- hic/
    - GSM3181429_BMDM_HiC_B1.hic

- rnaseq_4su/  

    - Count_matrix/
        - Count.matrix.xls
        - FPKM.matrix.xls
        - Gene.length.xls
        - TPM.matrix.xls
    
    4sU-seq raw data was processed by a pipeline based on https://nf-co.re/rnaseq/3.6.

- chipms_macro/
    - Kopie von GR ChIP MS.xlsx

- motifs/custom/
    - nr3c1_simplified_fullsite.motif
    - nr3c1_simplified_halfsite.motif  

    These are custom fullsite and halfsite motifs of NR3C1 for the genomewide motifscan.

Additionally, the workflow requires a manually generated IGV snapshot named `results/current/abcmodel/igv_snapshot_cd83.png` for the generation of Figure 3. Figure 5 does not get created in its final state by the workflow but had content added manually.

To keep the workflow organized, the individual datatypes are analyzed and integrated in several **snakemake subworkflows** located within `src/snakefiles/` that are included through the Snakefile in the project root directory.  

<a id="snakemake"></a>

# Snakemake workflow

<a id="execution"></a>

## Workflow execution

To run the workflow you need to be in the workflow's root directory.

The workflow can be executed by running the wrapper script `run_workflow.sh`. If you use the wrapper script it is not necessary to activate the conda environment, since it will be activated during execution. The general command for the execution of the workflow using the wrapping script is:

```
sh run_workflow.sh [njobs] [snakemake arguments] [rule1 rule2 ... rulen]
```

There are two modes to run the workflow. To run it in test mode, just execute:

```
sh run_workflow.sh 
```

This will show all the jobs that would be executed if the workflow was actually run. To know the jobs to be executed to achieve a certain rule or output file, add the target rule or output file (e.g.: `get_remap_data` or `results/current/Figure_chipseq.png`) following snakemake's specificiations (be aware of restrictions in rules using wildcards):

```
sh run_workflow.sh get_remap_data
```

If you want to run the workflow in production mode for a certain file, you need to add the number of jobs to be simultaneously submitted as the first argument followed by the desired file:

```
sh run_workflow.sh 10 results/current/Figure_chipseq.png
```

You can run the workflow in production mode, without specifying a target:

```
sh run_workflow.sh 10
```

This will use the input files of the first rule of the workflow as target files.

```
rule all:
  input:
    "results/current/Figures/Figure_chipseq.png",
    "results/current/Figures/Figure_peakgeneannotation.png",
    "results/current/Figures/Figure_abcresults.png",
    "results/current/Figures/Figure_GLMs.png",
    "results/current/Figures/Figure_stats.png"
```

Note: Filepaths can also be expanded with arguments from the config file.  
See more on that in the [snakemake readthedocs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html).

<a id="config"></a>

## Workflow configuration

Parameters of the workflow can be defined in the configuration file `config/sample_config.yaml`.

The resources necessary to run the jobs are defined in the file `config/cluster.yaml`. The ones already defined are based on test executions, so you may need to adjust them according to your own requirements.

You may add any resource allowed by the [`sbatch` command of the SLURM Workload Manager](https://slurm.schedmd.com/sbatch.html), but it has to be using the long names (e.g.: 'partition' and not 'p'). **The only limitation is for the `--cpus-per-task` option.** Since some of the commands executed require parallelization, the specification of the number of cpus to be used is **defined in the `config/sample_config.yaml` parameter of the corresponding rule**. The definition of `--cpus-per-task` in the `config/cluster.yaml` would override the values set in the `config/sample_config.yaml`, and since the change would not propagate to the actual execution of the command it could lead to abuse/misuse of the cluster resources.

<a id="visualization"></a>

## Workflow visualization

After activating the environment with `conda activate gr_env`, an image of the directed acyclic graph (DAG) can be generated using

```
snakemake --configfile config/sample_config.yaml --forceall --dag | dot -Tpdf > <imagename>.pdf
```

Alternatively, to avoid long waittimes, visualizing the rulegraph shows rule dependencies, without including every wildcard in the output

```
snakemake --configfile config/sample_config.yaml -s Snakefile --forceall --rulegraph | dot -Tpdf > dag.pdf
```

<a id="troubleshooting"></a>

## Troubleshooting

If you manually cancel submitted jobs, the working directory may be left in a locked state. To unlock the working directory execute:

```
sh run_workflow.sh --unlock
```

<a id="acknowledgments"></a>

# Acknowledgments

This workflow was developed with the help of Xavier Pastor from the Bioinformatics Core Facility at Helmholtz Munich.
