# The pangenome pipeline

The pangenome pipeline is a set of workflows and scripts to create pangenomes from the sourmash and roary software.

## Installation

The pipelines require conda and snakemake to operate:

- [conda installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- [snakemake installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Once these systems are installed, running the snakemake workflow will install the conda environment for the applicable rule.

Grab the repo with:

```
git clone git@github.com:ccbaumler/pangenome_pipeline.git
```

And `cd` into it.

```
cd pangenome_pipeline/
```

## Usage

Before running the snakefiles, configure the `config.yaml` file in the `config` directory according to your deepest desires.

The config file requires information on:

1. `output_directory`: An output directory path for the output.
2. `taxon_name_and_id`: A taxonomic name and identification number from [NCBI Datasets' Taxonomy Browser](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/).
> Note: the taxonomic name should not be split by white space, use `-` or `_`.
3. `k_values`: The word length for each kmer that sourmash generates. (Only sourmash)
4. `scale_value`: The magnitude scale to downsample the kmer space. (Only sourmash)
5. `batch_size`: The batch size to extract genome sequence files from NCBI. (Only sourmash)

The standard pangenomic pipeline for bacterial species is to use prokka and roary. You can accomplish this with the snakefiles:

1. `snakemake -s taxa.smk`
2. `snakemake -s prokka.smk`
3. `snakemake -s roary.smk`

The novel sourmash pangenome pipeline:

1. `snakemake -s directsketch.smk`
2. In - progress

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Author

Colton Baumler

[![UC Davis Email](https://img.shields.io/badge/UC_Davis-Email-blue?style=for-the-badge&colorA=blue&colorB=gold)](mailto:ccbaumler@ucdavis.edu) <a href="mailto:ccbaumler@gmail.com"><img src="https://img.shields.io/badge/gmail-%23DD0031.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
