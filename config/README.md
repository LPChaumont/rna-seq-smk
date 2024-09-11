# General settings

To configure this workflow, modify the following files to reflect your data and required settings:
* `config/config.yaml`: General workflow configuration and tool-specific settings.
* `config/samples.tsv`: Sample names, raw data paths and associated metadata.
* `config/comparisons.tsv`: Comparisons to be made between different conditions or groups.

## config.yaml
This file contains the general workflow configuration and the settings for the different tools. Configurable options should be explained in the comments above the respective entry or right here in this section.


## samples.tsv
For each sample, add a line to the sample sheet in `config/samples.tsv`. The sample sheet has to be a tab-separated file with a at least 4 columns and a header row as shown in `config/samples.tsv`. Additional columns can specify covariates (including batch effects) or any other metadata. These columns can then be used in the `diffexp: models:` section in `config/config.yaml`. The first 4 columns must match those defined in the table below.

The `sample` column is required and gives the sample name. The `condition` column is also required and specifies the condition for each sample. 

Column | Description
--- | ---
sample | Custom sample name. 
condition | The name of the condition a sample belongs to.
fastq_1 | The path to paired-end FASTQ file for read 1. 
fastq_2 | The path to paired-end FASTQ file for read 2.


## contrasts.tsv
Create a contrast sheet with information about the contrasts to analyse. For each contrast, add a line in the contrast sheet. The contrast sheet has to be a tab-separated file with at least 3 columns and a header row as shown in `config/contrats.tsv`. The reference/tested level values must be the same as in the condition column of the sample sheet. The first 3 columns must match those defined in the table below.

Column | Description
--- | ---
contrast | An arbitrary identifier, will be used to name contrast-wise output files.
reference_level | The control/base level for the comparison. 
tested_level | The treatment/target level for the comparison. 


