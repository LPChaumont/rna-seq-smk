# RNA-seq Snakemake workflow

A Snakemake workflow for short-read RNA-seq. It takes FASTQ files as input, performs quality control, trimming, alignment, quantification, differential gene expression analysis and splicing analysis.

## Dependencies

This Snakemake workflow uses mamba for installation and dependency resolution, so [Miniforge](https://github.com/conda-forge/miniforge) needs to be installed. In case you don’t use Miniforge you can always install Mamba.

```bash
conda install -n base -c conda-forge mamba
```

Afterwards, create a new conda environment called "smake" with Snakemake 7.32.4. Be aware that Snakemake 8.0.0 introduced some [breaking changes](https://snakemake.readthedocs.io/en/stable/project_info/history.html#breaking-changes).

```bash
mamba create -n smake -c bioconda -c conda-forge snakemake=7.32.4
```

## Installation

Clone this github repository.
```bash
git clone https://github.com/LPChaumont/rna-seq-smk.git
```

## Usage

Before running the workflow, ensure that the configuration is complete. The [`config`](config/README.md) directory contains the required files for defining sample information, experimental conditions, and other key parameters. The [`profiles`](profiles/README.md) directory includes the necessary settings for executing the workflow locally or on a SLURM cluster. See their respective `README.md` files.

Once the configuration is set, follow these steps to execute the workflow on a SLURM cluster:

- **Activate the Snakemake environment:**

    ```bash
    mamba activate smake
    ```

- **Download required files on the head node:**

    ```bash
    snakemake download --profile profiles/local/
    ```

- **Run the workflow on the compute nodes:**

    ```bash
    snakemake --profile profiles/slurm/
    ```

## Contact

__Author__ : Louis-Philippe Chaumont

__Email__ : <louis-philippe.chaumont@usherbrooke.ca>

__Research groups__: [Choquet Lab](https://www.choquetlab.com/) and [Scott Lab](https://bioinfo-scottgroup.med.usherbrooke.ca/)