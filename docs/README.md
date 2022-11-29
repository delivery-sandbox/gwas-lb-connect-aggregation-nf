# Pipeline documentation

<ins>Table of contents</ins>

  - [1 - Pipeline description](#1---pipeline-description)
    - [1.1 - Pipeline overview](#11---pipeline-overview)
    - [1.2 - Input](#12---input)
    - [1.3 - Default processes](#13---default-processes)
    - [1.4 - Optional processes](#14---optional-processes)
    - [1.5 - Output](#15---output)
  - [2 - Options](#2---options)
  - [3 - Usage](#3---usage) 
    - [3.1 - Running with Docker or Singularity](#31---running-with-docker-or-singularity)
    - [3.2 - Running with cloudos-cli](#32---running-with-cloudos-cli)
    - [3.3 - Execution examples](#33---execution-examples) 
  - [4 - Additional information](#4---additional-information) 
    - [4.1 - Stress testing](#41---stress-testing)
 
## 1 - <ins>Pipeline description</ins>

### 1.1 - <ins>Pipeline overview</ins>

  - Name: gwas-sumstats-harmonisation-nf
  - Tools: are listed in the environment.yml files in each respective tool bundle in the folder containers/<tool> in the repo.
  - Dependencies: are listed in the [docs/dependencies.md](docs/dependencies.md)

This pipeline is designed to harmonise and convert GWAS summary statistics from a number of sources into [MRCIEU GWAS-VCF V1.0 format](https://github.com/MRCIEU/gwas-vcf-specification/tree/1.0.0) with a standardised schema. Compatible input sources:
- Study IDs for GWAS studies in the [EBI GWAS catalog](https://www.ebi.ac.uk/gwas/)
- Study IDs for GWAS studies in the [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/)
- GWAS summary statistics in GWAS-VCF 1.0 format
- GWAS output tables from gwas-nf pipeline

Aditionally, this pipeline can also perform:
- **BETA** coefficient conversion to **Odds Ratio**.
- **BETA** and **BETA Standard Error** standardisation.
- Convert GWAS-VCF to `Hail` format.
- Filtering of the harmonised sumstats.

### 1.2 - <ins>Input</ins>

This pipeline expects as input (indicated using `--input` parameter) a text file containing the either a list of study IDs of GWAS studies in IEU OpenGWAS, or EBI GWAS catalogues, or paths to GWAS summary statistics files.

The `--gwas_source` parameter should be used to specify which type of input is being given ("ebi", "ieu", "gwas_vcf", or "gwas_table").

> NOTE: If EBI GWAS catalog is used as a source, only the GWAS studies for which harmonised datasets are available will be processed.

Input file should only contain IDs from a single source.

Example of input file for collecting EBI studies:

```
GCST001969
GCST90077560
``` 

Example of input file for collecting IEU studies:

```
ukb-b-14043
ieu-a-297
```


Example of input file for harmonising GWAS VCF files:

```
s3://lifebit-featured-datasets/pipelines/downstream-benchmarking/gwas_vcf/GRCh38/GCST90077560.harmonised.gwas.vcf
s3://lifebit-featured-datasets/pipelines/downstream-benchmarking/gwas_vcf/GRCh38/GCST001969.harmonised.gwas.vcf
```

Example of input file for harmonising GWAS table files:

```
s3://lifebit-featured-datasets/pipelines/downstream-omics/sumstats-harmonisation/testdata/gwas_table_inputs/allancs-gwas_bin-bolt_lmm.tsv
```

> Note: GWAS tables as input that are supported are output tables from the following GWAS tools: **regenie**, **saige**, **plink2_gwas**, **bolt-lmm**, **fastgwa**, **hail-gwas**, and **metal**.

When `--gwas_source` is set to "gwas_vcf", or "gwas_table", the user can directly provide the path of a GWAS VCF file or a GWAS table file to the `--input` parameter by also setting `--input_type` to `single`. For example:
```
nextflow run main.nf --input s3://lifebit-featured-datasets/pipelines/downstream-benchmarking/gwas_vcf/GRCh38/GCST90077560.harmonised.gwas.vcf \
    --input_type single \
    --gwas_source gwas_vcf \
    -with-docker
```


### 1.3 - <ins>Default processes</ins>

- `fetch_from_ebi`, `munge_ebi` : If `--gwas_source` is set to "ebi", these processes fetch the GWAS summary statistics and the available metadata from the specified studies in the input, then standardise and harmonise this data into an expanded GWAS VCF format.
- `fetch_from_ieu`, `munge_ieu`: If `--gwas_source` is set to "ieu", these processes fetch the GWAS summary statistics and the available metadata from the specified studies in the input, then standardise and harmonise this data into an expanded GWAS VCF format.
- `munge_gwas_tables`: If `--gwas_source` is set to "gwas_tables", this process standardises and harmonises the tabular data into an expanded GWAS VCF format.
- `munge_gwas_vcf` / `expand_gwas_vcf`: If `--force_munge` is `true`, `munge_gwas_vcf` will standardise and harmonise the input GWAS VCFs into an expanded GWAS VCF format. If `--force_munge` is `false` (`false` by default) the pipeline will not attempt to further harmonise the input GWAS VCFs but will just expand them. 
- `collapse_vcf`: Collapses the final GWAS VCF (after optional trait mapping, standardisation, filtering) into a GWAS VCF with multi-allelic variants collapsed into a single row (to fully conform with MRCIEU GWAS-VCF V1.0 format).
- `make_qc_plots`: generate some plots from the harmonised VCF file.
- `create_report`: create an HTML report with the plots generated in the previous process.

### 1.4 - <ins>Optional processes</ins>

- `standardisation` [`--standardise true`]: performs BETA and SE standardisation and adds the standardised values to the output VCF.
- `beta2or` [`--coef_conversion true`]: performs BETA to Odds Ratio conversion and adds the Odds Ratio value to the output VCF.
- `map_traits`, `insert_traits` [`--map_traits true`]: These processes map the available metadata traits to the OMOP common data model, then insert the mapped trait metadata into the GWAS VCF headers.
- `filter_sumstats` [`--filter_beta_smaller_than`, `--filter_beta_greater_than`, `--filter_p_greater_than`, `--filter_freq_smaller_than`, `--filter_freq_greater_than`, `--filter_missing_info`, `--filter_info`]: performs different filterings to the harmonised VCF file. More details on each of the specific parameters can be found in the pipeline help (`nextflow run main.nf --help`, see [Options](#opt) section).
- `vcf2hail` [`--convert_to_hail true`]: converts harmonised VCF to [`Hail`](https://hail.is/docs/0.2/index.html) format.

### 1.5 - <ins>Output</ins>

The main output of this pipeline is a harmonised VCF file for each of the requested studies. They can be found in `results/harmonised` folder.

Example of harmonised VCF file:

```shell
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Variant-specific number of samples/individuals with called genotypes used to test association with specified trait">
##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">
##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">
##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">
##GenomeBuild=GRCh38
##fileDate=20220923
##phasing=unphased
##source=gwas-sumstats-harmonisation-nf
##META=<ID=TotalSamples,Number=1,Type=Integer,Description="Total number of Samples in the association study">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=TraitEfos,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubTitle,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubAuthor,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubDate,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubJournal,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitReported,Number=1,Type=String,Description="Trait name used as input for trait search and mapping">
##META=<ID=TraitSourceOrigin,Number=1,Type=String,Description="Was vocabulary specified for trait name used as input for trait search and mapping.">
##META=<ID=TraitSourceConceptId,Number=1,Type=String,Description="OMOP concept ID of trait search-hit before mapping.">
##META=<ID=TraitSourceConceptName,Number=1,Type=String,Description="OMOP concept name of trait search-hit before mapping.">
##META=<ID=TraitSourceConceptCode,Number=1,Type=String,Description="Identifier within the vocabulary of trait search-hit before mapping.">
##META=<ID=TraitSourceVocabularyId,Number=1,Type=String,Description="Vocabulary of trait search-hit before mapping.">
##META=<ID=TraitSourceDomainId,Number=1,Type=String,Description="OMOP domain ID of trait search-hit before mapping.">
##META=<ID=TraitSourceConceptClassId,Number=1,Type=String,Description="OMOP class ID of trait search-hit before mapping.">
##META=<ID=TraitSourceStandardConcept,Number=1,Type=String,Description="Is the trait search-hit a standard OMOP concept before mapping.">
##META=<ID=TraitConceptId,Number=1,Type=String,Description="OMOP concept ID after mapping to a standard OMOP concept.">
##META=<ID=TraitConceptName,Number=1,Type=String,Description="OMOP concept name after mapping to a standard OMOP concept.">
##META=<ID=TraitConceptCode,Number=1,Type=String,Description="Identifier within the vocabulary after mapping to a standard OMOP concept.">
##META=<ID=TraitVocabularyId,Number=1,Type=String,Description="Vocabulary after mapping to standard concept.">
##META=<ID=TraitDomainId,Number=1,Type=String,Description="OMOP domain ID of concept after mapping.">
##META=<ID=TraitConceptClassId,Number=1,Type=String,Description="OMOP class ID of concept after mapping.">
##META=<ID=TraitStandardConcept,Number=1,Type=String,Description="Is the trait a standard OMOP concept after mapping. (Should always be Standard)">
##META=<ID=TraitGroupConceptName,Number=1,Type=String,Description="OMOP concept name of ontology grouping after mapping.">
##META=<ID=TraitGroupConceptId,Number=1,Type=String,Description="OMOP concept ID of ontology grouping after mapping.">
##SAMPLE=<ID=GCST90077560,TraitReported="heart rate",TotalSamples=408215,TraitEfos="http://www.ebi.ac.uk/efo/EFO_0004326",PubTitle="Exome sequencing and analysis of 454,787 UK Biobank participants.",PubAuthor="Backman JD",PubDate="2021-10-18",PubJournal="Nature",TotalVariants=405602,HarmonisedVariants=273557,TraitSourceOrigin="Inferred",TraitSourceConceptId="40462158",TraitSourceConceptName="Heart rate",TraitSourceConceptCode="250764009",TraitSourceVocabularyId="SNOMED",TraitSourceDomainId="Measurement",TraitSourceConceptClassId="Observable Entity",TraitSourceStandardConcept="NA",TraitConceptId="4239408",TraitConceptName="Heart rate",TraitConceptCode="364075005",TraitVocabularyId="SNOMED",TraitDomainId="Measurement",TraitConceptClassId="Observable Entity",TraitStandardConcept="S",TraitGroupConceptName="Measurement",TraitGroupConceptId="NA">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GCST90077560
1	930165	rs201186828	G	A	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999955:0.0502232:4.53193e-05:-0.0217538:408215:-0.006644:0.30539:0.978481
1	930204	rs148711625	G	A	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999956:0.289139:4.40944e-05:0.105264:408215:0.032148:0.30539:1.111004
1	930245	rs146327803	G	A	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999978:0.575432:2.20551e-05:0.252748:408215:0.07719:0.305397:1.287559
1	930324	rs748645146	T	C	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999993:0.0286508:7.34907e-06:-0.0315114:408215:-0.009624:0.305401:0.96898
1	935806	rs754952420	G	A	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999994:0.426929:6.12422e-06:0.38299:408215:0.116966:0.305401:1.466663
1	935833	rs144490434	C	G	.	.	.	SE:LP:AF:ES:SS:SES:SSE:OR	0.999983:3.47707:1.71478e-05:0.924088:408215:0.28222:0.305398:2.519569
```

This file contains the following data:

- Header: The header contains file-level metadata (`##field=value` lines e.g. `##GenomeBuild=GRCh38`), sumstats field descriptions (`##FORMAT=` lines), study-level metadata field descriptions (`##META=` lines), harmonised study-level data and OMOP-mapped trait data in (`##SAMPLE` lines).
- FORMAT columns should contain at least the following fields:
    * LP: -log10 p-value for effect estimate
    * ES: Effect size estimate relative to the alternative allele
    * SE: Standard error of effect size estimate
    * SS: Variant-specific number of samples/individuals with called genotypes used to test association with specified trait

    The following FORMAT fields will be available depending on the original data and the options selected:
    * AF: Alternative allele frequency in trait subset
    * NC: Variant-specific number of cases used to estimate genetic effect (binary traits only)
    * AC: Alternative allele count in the trait subset
    * OR: Odds ratio of effect (if `--coef_conversion true`)
    * SES: standardised BETA (if `--standardsise true`)
    * SSE: standardised BETA Standard Error (if `--standardise true`)
    * EZ: Z-score provided if it was used to derive the ES and SE fields
    * SI: Accuracy score of association statistics imputation
    * IS: I-squared heterogeneity for meta-GWAS

Other output:

- `results/QC_plots`: several plots and descriptive statistics files.
- `results/MultiQC`: `<study_id>_multiqc_report.html` files (one per study) as multiQC reports, with the plots available in `results/QC_plots`.
- `results/conversion`: it will keep a copy of the harmonised VCF after BETA to OR conversion and/or standardisation.
- `results/<ebi | ieu>`: it will keep a copy of the original summary statistics files downloaded from the source as well as the extracted metadata and the first harmonisation using `MungeSumstats`. Only present if `--keep_intermediate_files true`.
- `results/metadata`: mapped metadata traits files. Only present if `--keep_intermediate_files true`.
- `results/harmonised/filtered`: fully harmonised and filtered VCF.
- `results/harmonised/hail`: harmonised (and filtered if requested) GWAS summary statistics converted to `Hail` matrix format.


## 2 - <ins>Options</ins>
<a name="opt"></a>
See pipeline help (`nextflow run main.nf --help`) for a detailed explanation of the available options:

```
Usage:
The typical command for running the pipeline is as follows:
nextflow run main.nf --input input.txt [Options]

Inputs:
--input                    Input file (path). Newline delimited list of IDs or file paths.
--gwas_source              Type/source of GWAS input data. It should be one of the following supported
                           strings: 'ebi', 'ieu', 'gwas_vcf', 'gwas_tables'.
Options:
--standardise              Whether to perform BETA and SE standardisation (bool).
                           (default: false)
--coef_conversion          Whether to perform the coefficient conversion, from BETA
                           to Odds Ratio (bool).
                           (default: false)
--map_traits               Whether to map the traits in the input data to SNOMED terms using the OMOP
                           common data model (bool).
                           (default: false)
--munge_method             Method to use for munging/harmonising input data. Either use the tool
                           MungeSumstats ('MSS'), or use a simplified method ('simple').
                           (default: MSS)
--force_munge              Force munging/harmonisation of input when `--gwas_source` is
                           `gwas_vcf'. When 'false', assume GWAS VCF input is already harmonised 
                           and skip munging step. (bool).
                           (default: false)
--miss_percent_allow       Missingness proportion allowed for a column from input data. Each
                           column above this proportion will be dropped (int / float).
                           (default: 0.1)
--keep_intermediate_files  Whether to keep intermediate files (bool).
                           (default: false)
--dbsnp                    Version of dbSNP database to use when harmonising variants. Supported
                           versions: '144', '155'.
                           (default: 155)
--omop_vocabulary          Path or link to an OMOP vocabulary DB file (path).
                           (default: https://omopvocabs.s3.eu-west-1.amazonaws.com/omop-vocab-clean.zip)
--convert_to_hail          Whether to convert the harmonised VCF to Hail MatrixTable format (bool).
                           (default: false)
--take_n_studies           Take n studies from the given input file to run the pipeline. Set to '-1' to use
                           all studies in the input file. (int).
                           (default: -1)

GWAS filtering options:
--filter_beta_smaller_than  Exclude variants with BETA smaller than the specified value (float)
                            (Default: null)
--filter_beta_greater_than  Exclude variants with BETA greater than the specified value (float)
                            (Default: null)
--filter_p_greater_than     Exclude variants with a P value greater than the specified value (float)
                            (Default: null)
--filter_freq_smaller_than  Exclude variants with a alternate allele frequency smaller than the specified
                            value. (float)
                            (Default: null)
--filter_freq_greater_than  Exclude variants with a alternate allele frequency greater than the specified
                            value. (float)
                            (Default: null)
--filter_missing_info       Exclude variants with a missing value for INFO column (bool)
                            (Default: false)
--filter_info               Exclude variants with a INFO value smaller than the one specified (float).
                            (Default: null)

Resource options:
--max_cpus         Maximum number of CPUs (int)
                   (default: 2)  
--max_memory       Maximum memory (memory unit)
                   (default: 16 GB)
--max_time         Maximum time (time unit)
                   (default: 8h)
--outdir           Output directory(path)
                   (default: results)
```

## 3 - <ins>Usage</ins>  

Inside the folder [conf/](./conf) pre-curated configurations of parameters to execute the pipeline in different modes are provided.
Every configuration file can be run, from the root of the local clone of the repository, using one of the following commands.

> NOTE: Adding `-with-docker` or `-with-singularity` is required because all of the dependencies are used from the containers.

### 3.1 - <ins>Running with Docker or Singularity</ins>
## Running with Docker

```bash
nextflow run main.nf -profile <any-profile> -with-docker
```
Example:

```bash
nextflow run main.nf -profile <any-profile> -with-docker
```

## Running with Singularity

```bash
nextflow run main.nf -profile <any-profile> -with-singularity
```
Example:

```bash
nextflow run main.nf -profile <any-profile> -with-singularity
```

### 3.2 - <ins>Running with cloudos-cli</ins>
## Running with cloudos-cli to submit to CloudOS platform on AWS

### Set up of environmental variables

This step only needs to be set up once. More detailed information for the `cloudos-cli` package required and optional parameters can be found in the official documentation, at https://github.com/lifebit-ai/cloudos-cli.

```bash
# Workspace specific variables
MY_API_KEY="**"
WORKSPACE_ID="**"
CLOUDOS_URL="**"
PROJECT_NAME="**"

# Workflow specific variables
WORKFLOW_NAME="bi-gwas-sumstats-harmonisation-nf"
INSTANCE_TYPE="c3.4xlarge"

CLOUDOS_CLI_OPTIONS=" --cloudos-url $CLOUDOS_URL --apikey $MY_API_KEY --workspace-id $WORKSPACE_ID --workflow-name $WORKFLOW_NAME --project-name $PROJECT_NAME --instance-type $INSTANCE_TYPE"
CLOUDOS_CLI_OPTIONS+=" --resumable --spot"
```

### Installation of `cloudos-cli`

General installation instructions can be found in https://github.com/lifebit-ai/cloudos-cli#from-github.

#### Pre-installed cloudos-cli in the HPC

Ask from your system administrator to install the package based on the official documentation. If already available in your system, you can type the commands below to verify that the package is successfully installed:

```bash
ml load cloudos
cloudos --version
```

The output will display the version of `cloudos-cli` you are using.

```console
cloudos, version <SELECTED VERSION IN USE>
```

#### Running from pre-build docker image

> NOTE: For convenience and ensuring the commands are being executing from the same network as the jobs that will be submitted, to ensure access to Lifebit CloudOS services is successful, you can use a terminal in a JupyterLab session within the Lifebit CloudOS platform. Docker is already available, and hence you can use the `cloudos-cli` from the pre-built image that is provided in the `cloudos-cli` documentation.

```bash
docker run --rm -it quay.io/lifebitaiorg/cloudos-cli:<TAG>
```

### Execution command

After passing all the required environmental variables in $CLOUDOS_CLI_OPTIONS, the command to submit programmatically a job to the Lifebit CloudOS Platform can be condensed to:

**NOTE**: For the pipeline to successfully fetch data from GWAS catalogues, (i.e. when `--gwas_source` is set to `ebi` or `ieu`) it is necessary to use the `network_cloudos` profile in order to successfully connect to these external resources. For example, to test using the example profiles `test_ebi` and `test_ieu` use `-profile test_ebi,network_cloudos` or `-profile test_ieu,network_cloudos`.

```bash
cloudos job run --nextflow-profile <any-profile>,network_cloudos $CLOUDOS_CLI_OPTIONS
```

### 3.3 - <ins>Execution examples</ins>
   
#### Harmonising GWAS summary statistics from EBI GWAS catalogue

```bash
nextflow run main.nf -profile test_ebi -with-docker
nextflow run main.nf -profile test_ebi -with-singularity
cloudos job run --nextflow-profile test_ebi,network_cloudos $CLOUDOS_CLI_OPTIONS
```

<details>
<summary>Expected output:</summary>

```
tree -fh results
results
├── [  55]  results/conversion
│   └── [ 35M]  results/conversion/GCST90077560_harmonised_conv_sumstats.vcf
├── [ 164]  results/ebi
│   ├── [ 18M]  results/ebi/34662886-GCST90077560-EFO_0004326.h.tsv.gz
│   ├── [ 25M]  results/ebi/GCST90077560_harmonised_sumstats.vcf
│   ├── [ 716]  results/ebi/GCST90077560_metadata_ebi.tsv
│   └── [  11]  results/ebi/reported_traits.tsv
├── [  74]  results/harmonised
│   ├── [  59]  results/harmonised/filtered
│   │   └── [7.6K]  results/harmonised/filtered/GCST90077560_harmonised_filtered_sumstats.vcf
│   ├── [7.6K]  results/harmonised/GCST90077560.harmonised.gwas.vcf
│   └── [  41]  results/harmonised/hail
│       └── [ 131]  results/harmonised/hail/GCST90077560_hail_matrix.mt
│           ├── [ 149]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols
│           │   ├── [ 261]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/metadata.json.gz
│           │   ├── [ 147]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/README.txt
│           │   ├── [  72]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/rows
│           │   │   ├── [ 252]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/rows/metadata.json.gz
│           │   │   └── [  39]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/rows/parts
│           │   │       └── [  26]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/rows/parts/part-0
│           │   └── [   0]  results/harmonised/hail/GCST90077560_hail_matrix.mt/cols/_SUCCESS
│           ├── [ 149]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries
│           │   ├── [ 344]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/metadata.json.gz
│           │   ├── [ 147]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/README.txt
│           │   ├── [  72]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/rows
│           │   │   ├── [ 650]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/rows/metadata.json.gz
│           │   │   └── [ 113]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/rows/parts
│           │   │       └── [2.2K]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/rows/parts/part-0-2d9bc4c1-05d7-4b1d-9519-0e01187f4521
│           │   └── [   0]  results/harmonised/hail/GCST90077560_hail_matrix.mt/entries/_SUCCESS
│           ├── [ 164]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals
│           │   ├── [  72]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/globals
│           │   │   ├── [ 240]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/globals/metadata.json.gz
│           │   │   └── [  39]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/globals/parts
│           │   │       └── [  11]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/globals/parts/part-0
│           │   ├── [ 253]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/metadata.json.gz
│           │   ├── [ 147]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/README.txt
│           │   ├── [  72]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/rows
│           │   │   ├── [ 240]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/rows/metadata.json.gz
│           │   │   └── [  39]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/rows/parts
│           │   │       └── [  11]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/rows/parts/part-0
│           │   └── [   0]  results/harmonised/hail/GCST90077560_hail_matrix.mt/globals/_SUCCESS
│           ├── [  61]  results/harmonised/hail/GCST90077560_hail_matrix.mt/index
│           │   └── [  90]  results/harmonised/hail/GCST90077560_hail_matrix.mt/index/part-0-2d9bc4c1-05d7-4b1d-9519-0e01187f4521.idx
│           │       ├── [ 589]  results/harmonised/hail/GCST90077560_hail_matrix.mt/index/part-0-2d9bc4c1-05d7-4b1d-9519-0e01187f4521.idx/index
│           │       └── [ 199]  results/harmonised/hail/GCST90077560_hail_matrix.mt/index/part-0-2d9bc4c1-05d7-4b1d-9519-0e01187f4521.idx/metadata.json.gz
│           ├── [ 372]  results/harmonised/hail/GCST90077560_hail_matrix.mt/metadata.json.gz
│           ├── [ 147]  results/harmonised/hail/GCST90077560_hail_matrix.mt/README.txt
│           ├── [ 149]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows
│           │   ├── [ 323]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/metadata.json.gz
│           │   ├── [ 147]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/README.txt
│           │   ├── [  72]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/rows
│           │   │   ├── [ 638]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/rows/metadata.json.gz
│           │   │   └── [ 113]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/rows/parts
│           │   │       └── [ 683]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/rows/parts/part-0-2d9bc4c1-05d7-4b1d-9519-0e01187f4521
│           │   └── [   0]  results/harmonised/hail/GCST90077560_hail_matrix.mt/rows/_SUCCESS
│           └── [   0]  results/harmonised/hail/GCST90077560_hail_matrix.mt/_SUCCESS
├── [  73]  results/MultiQC
│   ├── [5.7M]  results/MultiQC/GCST90077560_multiqc_report.html
│   └── [5.7M]  results/MultiQC/multiqc_report.html
├── [ 189]  results/pipeline_info
│   ├── [2.8M]  results/pipeline_info/execution_report_2022-10-11_19-59-50.html
│   ├── [7.9K]  results/pipeline_info/execution_timeline_2022-10-11_19-59-50.html
│   ├── [1.3K]  results/pipeline_info/execution_trace_2022-10-11_20-16-36.txt
│   └── [1.5K]  results/pipeline_info/pipeline_metadata_report.tsv
└── [  26]  results/QC_plots
    └── [ 326]  results/QC_plots/GCST90077560
        ├── [ 31K]  results/QC_plots/GCST90077560/boxPlot-maf-beta.png
        ├── [ 31K]  results/QC_plots/GCST90077560/boxPlot-maf-p.png
        ├── [4.4K]  results/QC_plots/GCST90077560/descriptive_statistics.tsv
        ├── [9.4K]  results/QC_plots/GCST90077560/descriptive_statistics.txt
        ├── [ 35K]  results/QC_plots/GCST90077560/frequencyCurve-Absolute_beta.png
        ├── [ 23K]  results/QC_plots/GCST90077560/frequencyCurve-MAF.png
        ├── [171K]  results/QC_plots/GCST90077560/Manhattan.png
        ├── [ 33K]  results/QC_plots/GCST90077560/QQ.png
        ├── [ 20K]  results/QC_plots/GCST90077560/scatterPlot-beta-logp.png
        ├── [ 20K]  results/QC_plots/GCST90077560/scatterPlot-maf-beta.png
        └── [ 26K]  results/QC_plots/GCST90077560/scatterPlot-maf-p.png

26 directories, 51 files
```

</details>
<br>


<details>
<summary>Format of main outputs:</summary>

```shell
head -n 60 results/harmonised/GCST90077560.harmonised.gwas.vcf
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20221011
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Variant-specific number of samples/individuals with called genotypes used to test association with specified trait">
##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">
##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">
##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">
##phasing=unphased
##source=gwas-sumstats-harmonisation-nf
##GenomeBuild=GRCh38
##META=<ID=TotalSamples,Number=1,Type=Integer,Description="Total number of Samples in the association study">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=TraitEfos,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubTitle,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubAuthor,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubDate,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubJournal,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitReported,Number=1,Type=String,Description="Trait name used as input for trait search and mapping">
##SAMPLE=<ID=GCST90077560,TotalSamples=408215,TraitReported="heart rate",TraitEfos="http://www.ebi.ac.uk/efo/EFO_0004326",PubTitle="Exome sequencing and analysis of 454,787 UK Biobank participants.",PubAuthor="Backman JD",PubDate=2021-10-18,PubJournal="Nature",TotalVariants=405602,HarmonisedVariants=273557>
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##bcftools_filterVersion=7cd83b7+htslib-
##bcftools_filterCommand=filter '-e FORMAT/ES[*]<-0.3' temp_sumstats.vcf; Date=Tue Oct 11 20:13:29 2022
##bcftools_filterCommand=filter '-e FORMAT/ES[*]>0.3' temp_sumstats.vcf; Date=Tue Oct 11 20:13:31 2022
##bcftools_filterCommand=filter '-e FORMAT/LP[*]>0.5' temp_sumstats.vcf; Date=Tue Oct 11 20:13:32 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]<0.3' temp_sumstats.vcf; Date=Tue Oct 11 20:13:33 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]>0.8' temp_sumstats.vcf; Date=Tue Oct 11 20:13:34 2022
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  GCST90077560
1       99851033        rs2307130       A       G       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00194:0.109883:0.50213:0.000606106:408215:0.000185:0.305996:1.00061
2       70935956        rs11681642      T       C       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00026:0.0910552:0.445335:-0.000514012:408215:-0.000157:0.305483:0.999486
2       85322745        rs4240199       A       G       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00093:0.272346:0.590608:-0.00134922:408215:-0.000412:0.305687:0.998652
3       108916126       rs10933973      C       A       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       0.999904:0.04873:0.316295:-0.0003073:408215:-9.4e-05:0.305374:0.999693
5       1240642 rs7447815       C       G       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       0.996376:0.0533449:0.406968:0.00031752:408215:9.7e-05:0.304296:1.00032
5       33751349        rs1530507       A       T       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       0.99888:0.472848:0.754393:0.00237989:408215:0.000727:0.305061:1.00238
7       21543345        rs2285943       G       T       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00081:0.299042:0.502701:0.00142937:408215:0.000437:0.305651:1.00143
```

</details>

#### Harmonising GWAS summary statistics from IEU Open GWAS project

```bash
nextflow run main.nf -profile test_ieu -with-docker
nextflow run main.nf -profile test_ieu -with-singularity
cloudos job run --nextflow-profile test_ieu,network_cloudos $CLOUDOS_CLI_OPTIONS
```

<details>
<summary>Expected output:</summary>

```
tree -fh results
results
├── [  54]  results/conversion
│   └── [387M]  results/conversion/ukb-b-14043_harmonised_conv_sumstats.vcf
├── [  73]  results/harmonised
│   ├── [  58]  results/harmonised/filtered
│   │   └── [171M]  results/harmonised/filtered/ukb-b-14043_harmonised_filtered_sumstats.vcf
│   ├── [  40]  results/harmonised/hail
│   │   └── [ 131]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt
│   │       ├── [ 149]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols
│   │       │   ├── [ 261]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/metadata.json.gz
│   │       │   ├── [ 147]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/README.txt
│   │       │   ├── [  72]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows
│   │       │   │   ├── [ 252]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/metadata.json.gz
│   │       │   │   └── [  39]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/parts
│   │       │   │       └── [  25]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/parts/part-0
│   │       │   └── [   0]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/_SUCCESS
│   │       ├── [ 149]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries
│   │       │   ├── [ 339]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/metadata.json.gz
│   │       │   ├── [ 147]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/README.txt
│   │       │   ├── [  72]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows
│   │       │   │   ├── [ 635]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/metadata.json.gz
│   │       │   │   └── [ 113]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/parts
│   │       │   │       └── [ 54M]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/parts/part-0-67807801-2428-4f28-98a7-ac0a723e2d7f
│   │       │   └── [   0]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/_SUCCESS
│   │       ├── [ 164]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals
│   │       │   ├── [  72]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals
│   │       │   │   ├── [ 240]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/metadata.json.gz
│   │       │   │   └── [  39]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/parts
│   │       │   │       └── [  11]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/parts/part-0
│   │       │   ├── [ 253]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/metadata.json.gz
│   │       │   ├── [ 147]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/README.txt
│   │       │   ├── [  72]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows
│   │       │   │   ├── [ 240]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/metadata.json.gz
│   │       │   │   └── [  39]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/parts
│   │       │   │       └── [  11]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/parts/part-0
│   │       │   └── [   0]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/_SUCCESS
│   │       ├── [  61]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/index
│   │       │   └── [  90]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-67807801-2428-4f28-98a7-ac0a723e2d7f.idx
│   │       │       ├── [ 18M]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-67807801-2428-4f28-98a7-ac0a723e2d7f.idx/index
│   │       │       └── [ 206]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-67807801-2428-4f28-98a7-ac0a723e2d7f.idx/metadata.json.gz
│   │       ├── [ 366]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/metadata.json.gz
│   │       ├── [ 147]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/README.txt
│   │       ├── [ 149]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows
│   │       │   ├── [ 327]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/metadata.json.gz
│   │       │   ├── [ 147]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/README.txt
│   │       │   ├── [  72]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows
│   │       │   │   ├── [ 632]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/metadata.json.gz
│   │       │   │   └── [ 113]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/parts
│   │       │   │       └── [ 15M]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/parts/part-0-67807801-2428-4f28-98a7-ac0a723e2d7f
│   │       │   └── [   0]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/_SUCCESS
│   │       └── [   0]  results/harmonised/hail/ukb-b-14043_hail_matrix.mt/_SUCCESS
│   └── [171M]  results/harmonised/ukb-b-14043.harmonised.gwas.vcf
├── [ 138]  results/ieu
│   ├── [  52]  results/ieu/reported_traits.tsv
│   ├── [269M]  results/ieu/ukb-b-14043_harmonised_sumstats.vcf
│   ├── [ 424]  results/ieu/ukb-b-14043_metadata_ieu.tsv
│   └── [ 93M]  results/ieu/ukb-b-14043.vcf.gz
├── [  72]  results/MultiQC
│   ├── [5.6M]  results/MultiQC/multiqc_report.html
│   └── [5.6M]  results/MultiQC/ukb-b-14043_multiqc_report.html
├── [ 189]  results/pipeline_info
│   ├── [2.8M]  results/pipeline_info/execution_report_2022-10-11_17-12-33.html
│   ├── [7.9K]  results/pipeline_info/execution_timeline_2022-10-11_17-12-33.html
│   ├── [1.3K]  results/pipeline_info/execution_trace_2022-10-11_17-38-55.txt
│   └── [1.5K]  results/pipeline_info/pipeline_metadata_report.tsv
└── [  25]  results/QC_plots
    └── [ 326]  results/QC_plots/ukb-b-14043
        ├── [ 34K]  results/QC_plots/ukb-b-14043/boxPlot-maf-beta.png
        ├── [ 30K]  results/QC_plots/ukb-b-14043/boxPlot-maf-p.png
        ├── [4.5K]  results/QC_plots/ukb-b-14043/descriptive_statistics.tsv
        ├── [9.4K]  results/QC_plots/ukb-b-14043/descriptive_statistics.txt
        ├── [ 24K]  results/QC_plots/ukb-b-14043/frequencyCurve-Absolute_beta.png
        ├── [ 33K]  results/QC_plots/ukb-b-14043/frequencyCurve-MAF.png
        ├── [ 52K]  results/QC_plots/ukb-b-14043/Manhattan.png
        ├── [ 24K]  results/QC_plots/ukb-b-14043/QQ.png
        ├── [ 21K]  results/QC_plots/ukb-b-14043/scatterPlot-beta-logp.png
        ├── [ 23K]  results/QC_plots/ukb-b-14043/scatterPlot-maf-beta.png
        └── [ 19K]  results/QC_plots/ukb-b-14043/scatterPlot-maf-p.png

26 directories, 51 files
```

</details>
<br>


<details>
<summary>Format of main outputs:</summary>

```shell
head -n 60 results/harmonised/ukb-b-14043.harmonised.gwas.vcf
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20221011
##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">
##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">
##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">
##phasing=unphased
##source=gwas-sumstats-harmonisation-nf
##GenomeBuild=GRCh37
##META=<ID=Population,Number=1,Type=String,Description="Ancestral population of individuals in the association study">
##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">
##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">
##META=<ID=TotalSamples,Number=1,Type=Integer,Description="Total number of Samples in the association study">
##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=TraitTags,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitEfos,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=StudySample,Number=1,Type=String,Description="Information about study sample, e.g. participants' sex.">
##META=<ID=PubAuthor,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubDate,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=Comments,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitReported,Number=1,Type=String,Description="Trait name used as input for trait search and mapping">
##SAMPLE=<ID=ukb-b-14043,TotalSamples=361264,TotalCases=2094,TotalControls=359170,Population="European",StudySample="Males and Females",TraitReported="Illnesses of siblings: Alzheimer's disease/dementia",TraitEfos=NA,TraitTags=NA,PubAuthor="Ben Elsworth; MRC-IEU",PubDate=2018,Comments="20111#10: Output from GWAS pipeline using Phesant derived variables from UKBiobank",StudyType="CaseControl",TotalVariants=3291295,HarmonisedVariants=3275417>
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##bcftools_filterVersion=7cd83b7+htslib-
##bcftools_filterCommand=filter '-e FORMAT/ES[*]<-0.3' temp_sumstats.vcf; Date=Tue Oct 11 17:33:51 2022
##bcftools_filterCommand=filter '-e FORMAT/ES[*]>0.3' temp_sumstats.vcf; Date=Tue Oct 11 17:34:09 2022
##bcftools_filterCommand=filter '-e FORMAT/LP[*]>0.5' temp_sumstats.vcf; Date=Tue Oct 11 17:34:22 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]<0.3' temp_sumstats.vcf; Date=Tue Oct 11 17:34:33 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]>0.8' temp_sumstats.vcf; Date=Tue Oct 11 17:34:40 2022
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ukb-b-14043
1       91536   rs6702460       G       T       .       PASS    .       LP:SE:AF:ES:SES:SSE:OR  0.0809219:0.000309493:0.456754:6.81371e-05:0.000298:0.001351:1.00007
1       706368  rs12029736      A       G       .       PASS    .       LP:SE:AF:ES:SES:SSE:OR  0.00877392:0.000219453:0.515732:6.11875e-06:2.7e-05:0.000958:1.00001
1       763394  rs3115847       G       A       .       PASS    .       LP:SE:AF:ES:SES:SSE:OR  0.455932:0.000257204:0.706822:-0.000242389:-0.001058:0.001123:0.999758
```

</details>


#### Harmonising GWAS summary statistics from GWAS VCF input

```bash
nextflow run main.nf -profile test_gwas_vcf_basic -with-docker
nextflow run main.nf -profile test_gwas_vcf_basic -with-singularity
cloudos job run --nextflow-profile test_gwas_vcf_basic,network_cloudos $CLOUDOS_CLI_OPTIONS
```

<details>
<summary>Expected output:</summary>

```
tree -fh results
results
├── [  54]  results/conversion
│   └── [557K]  results/conversion/ukb-b-14043_harmonised_conv_sumstats.vcf
├── [  45]  results/harmonised
│   └── [556K]  results/harmonised/ukb-b-14043.harmonised.gwas.vcf
├── [  72]  results/MultiQC
│   ├── [5.7M]  results/MultiQC/multiqc_report.html
│   └── [5.7M]  results/MultiQC/ukb-b-14043_multiqc_report.html
├── [ 189]  results/pipeline_info
│   ├── [2.8M]  results/pipeline_info/execution_report_2022-10-11_17-11-51.html
│   ├── [6.9K]  results/pipeline_info/execution_timeline_2022-10-11_17-11-51.html
│   ├── [ 915]  results/pipeline_info/execution_trace_2022-10-11_17-19-09.txt
│   └── [1.5K]  results/pipeline_info/pipeline_metadata_report.tsv
└── [  25]  results/QC_plots
    └── [ 326]  results/QC_plots/ukb-b-14043
        ├── [ 37K]  results/QC_plots/ukb-b-14043/boxPlot-maf-beta.png
        ├── [ 32K]  results/QC_plots/ukb-b-14043/boxPlot-maf-p.png
        ├── [ 439]  results/QC_plots/ukb-b-14043/descriptive_statistics.tsv
        ├── [9.1K]  results/QC_plots/ukb-b-14043/descriptive_statistics.txt
        ├── [ 29K]  results/QC_plots/ukb-b-14043/frequencyCurve-Absolute_beta.png
        ├── [ 33K]  results/QC_plots/ukb-b-14043/frequencyCurve-MAF.png
        ├── [ 42K]  results/QC_plots/ukb-b-14043/Manhattan.png
        ├── [ 35K]  results/QC_plots/ukb-b-14043/QQ.png
        ├── [ 23K]  results/QC_plots/ukb-b-14043/scatterPlot-beta-logp.png
        ├── [ 66K]  results/QC_plots/ukb-b-14043/scatterPlot-maf-beta.png
        └── [ 68K]  results/QC_plots/ukb-b-14043/scatterPlot-maf-p.png

6 directories, 19 files
```

</details>
<br>


<details>
<summary>Format of main outputs:</summary>

```shell
head -n 60 results/harmonised/ukb-b-14043.harmonised.gwas.vcf
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20220905
##FORMAT=<ID=LP,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">
##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">
##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">
##phasing=unphased
##source=gwas-sumstats-harmonisation-nf
##GenomeBuild=GRCh37
##META=<ID=Population,Number=1,Type=String,Description="Ancestral population of individuals in the association study">
##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">
##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">
##META=<ID=TotalSamples,Number=1,Type=Integer,Description="Total number of Samples in the association study">
##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=TraitTags,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitEfos,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=StudySample,Number=1,Type=String,Description="Information about study sample, e.g. participants' sex.">
##META=<ID=PubAuthor,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=PubDate,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=Comments,Number=1,Type=String,Description="info from GWAS catalogues">
##META=<ID=TraitReported,Number=1,Type=String,Description="Trait name used as input for trait search and mapping">
##SAMPLE=<ID=ukb-b-14043,TotalSamples=361264,TotalCases=2094,TotalControls=359170,Population="European",StudySample="Males and Females",TraitReported="Illnesses of siblings: Alzheimer's disease/dementia",TraitEfos=NA,TraitTags=NA,PubAuthor="Ben Elsworth; MRC-IEU",PubDate=2018,Comments="20111#10: Output from GWAS pipeline using Phesant derived variables from UKBiobank",StudyType="CaseControl",TotalVariants=3291295,HarmonisedVariants=3275417>
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##bcftools_normVersion=7cd83b7+htslib-
##bcftools_normCommand=norm --multiallelic -any -O v -o ukb-b-14043_harmonised_sumstats.vcf mini.gwas.vcf.gz; Date=Tue Oct 11 17:16:20 2022
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ukb-b-14043
22      19033532        rs5992354       T       C       .       .       .       LP:SE:AF:ES:SES:SSE:OR  1.06048:0.000184461:0.376733:-0.000315197:-0.036857:0.021569:0.999685
22      19033768        rs1051248       G       C       .       .       .       LP:SE:AF:ES:SES:SSE:OR  0.0757207:0.000230819:0.183981:-4.69561e-05:-0.005491:0.02699:0.999953
22      19034079        rs715539        A       C,T     .       .       .       LP:SE:AF:ES:SES:SSE:OR  0.148742,0.886057:0.000230016,0.000224936:0.191493,0.204865:-8.49849e-05,-0.00033861:-0.009937,-0.039594:0.026896,0.026302:0.999915,0.999661
22      19034240        rs2012929       A       G       .       .       .       LP:SE:AF:ES:SES:SSE:OR  0.443698:0.000182256:0.402145:-0.000165982:-0.019409:0.021312:0.999834
22      19034355        rs574435982     G       C       .       .       .       LP:SE:AF:ES:SES:SSE:OR  0.244125:0.000197806:0.290243:-0.000111599:-0.01305:0.02313:0.999888
22      19034829        rs2000996       A       G       .       .       .       LP:SE:AF:ES:SES:SSE:OR  0.431798:0.000182248:0.402135:-0.000164767:-0.019267:0.021311:0.999835
22      19035072        rs17810649      A       G       .       .       .       LP:SE:AF:ES:SES:SSE:OR  1.25181:0.000226516:0.192683:-0.000432101:-0.050527:0.026487:0.999568
```
<details>
<summary>Expected output:</summary>

```
.
├── [  54]  ./conversion
│   └── [462M]  ./conversion/ukb-b-14043_harmonised_conv_sumstats.vcf
├── [  83]  ./harmonised
│   ├── [  58]  ./harmonised/filtered
│   │   └── [154M]  ./harmonised/filtered/ukb-b-14043_harmonised_filtered_sumstats.vcf
│   ├── [  40]  ./harmonised/hail
│   │   └── [ 131]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt
│   │       ├── [ 149]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols
│   │       │   ├── [ 261]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/metadata.json.gz
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/.metadata.json.gz.crc
│   │       │   ├── [ 147]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/README.txt
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/.README.txt.crc
│   │       │   ├── [  72]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows
│   │       │   │   ├── [ 252]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/metadata.json.gz
│   │       │   │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/.metadata.json.gz.crc
│   │       │   │   └── [  39]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/parts
│   │       │   │       ├── [  17]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/parts/part-0
│   │       │   │       └── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/rows/parts/.part-0.crc
│   │       │   ├── [   0]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/_SUCCESS
│   │       │   └── [   8]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/cols/._SUCCESS.crc
│   │       ├── [ 149]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries
│   │       │   ├── [ 372]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/metadata.json.gz
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/.metadata.json.gz.crc
│   │       │   ├── [ 147]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/README.txt
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/.README.txt.crc
│   │       │   ├── [  72]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows
│   │       │   │   ├── [ 680]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/metadata.json.gz
│   │       │   │   ├── [  16]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/.metadata.json.gz.crc
│   │       │   │   └── [ 113]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/parts
│   │       │   │       ├── [ 44M]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/parts/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be
│   │       │   │       └── [353K]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/rows/parts/.part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.crc
│   │       │   ├── [   0]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/_SUCCESS
│   │       │   └── [   8]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/entries/._SUCCESS.crc
│   │       ├── [ 164]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals
│   │       │   ├── [  72]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals
│   │       │   │   ├── [ 240]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/metadata.json.gz
│   │       │   │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/.metadata.json.gz.crc
│   │       │   │   └── [  39]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/parts
│   │       │   │       ├── [  11]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/parts/part-0
│   │       │   │       └── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/globals/parts/.part-0.crc
│   │       │   ├── [ 253]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/metadata.json.gz
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/.metadata.json.gz.crc
│   │       │   ├── [ 147]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/README.txt
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/.README.txt.crc
│   │       │   ├── [  72]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows
│   │       │   │   ├── [ 240]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/metadata.json.gz
│   │       │   │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/.metadata.json.gz.crc
│   │       │   │   └── [  39]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/parts
│   │       │   │       ├── [  11]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/parts/part-0
│   │       │   │       └── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/rows/parts/.part-0.crc
│   │       │   ├── [   0]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/_SUCCESS
│   │       │   └── [   8]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/globals/._SUCCESS.crc
│   │       ├── [  61]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index
│   │       │   └── [  90]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.idx
│   │       │       ├── [ 14M]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.idx/index
│   │       │       ├── [109K]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.idx/.index.crc
│   │       │       ├── [ 206]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.idx/metadata.json.gz
│   │       │       └── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/index/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.idx/.metadata.json.gz.crc
│   │       ├── [ 394]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/metadata.json.gz
│   │       ├── [ 147]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/README.txt
│   │       ├── [ 149]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows
│   │       │   ├── [ 327]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/metadata.json.gz
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/.metadata.json.gz.crc
│   │       │   ├── [ 147]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/README.txt
│   │       │   ├── [  12]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/.README.txt.crc
│   │       │   ├── [  72]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows
│   │       │   │   ├── [ 632]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/metadata.json.gz
│   │       │   │   ├── [  16]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/.metadata.json.gz.crc
│   │       │   │   └── [ 113]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/parts
│   │       │   │       ├── [6.4M]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/parts/part-0-3a863d74-e50a-408a-ac6f-a2a283df17be
│   │       │   │       └── [ 51K]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/rows/parts/.part-0-3a863d74-e50a-408a-ac6f-a2a283df17be.crc
│   │       │   ├── [   0]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/_SUCCESS
│   │       │   └── [   8]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/rows/._SUCCESS.crc
│   │       └── [   0]  ./harmonised/hail/ukb-b-14043_hail_matrix.mt/_SUCCESS
│   └── [462M]  ./harmonised/ukb-b-14043_fully_harmonised_sumstats.vcf
├── [ 146]  ./ieu
│   ├── [303M]  ./ieu/ukb-b-14043_harmonised_sumstats.vcf
│   ├── [ 387]  ./ieu/ukb-b-14043_metadata_ieu.tsv
│   ├── [   7]  ./ieu/ukb-b-14043_sample_size.txt
│   └── [ 93M]  ./ieu/ukb-b-14043.vcf.gz
├── [  88]  ./metadata
│   ├── [ 715]  ./metadata/ukb-b-14043_harmonised_metadata.txt
│   └── [1.3K]  ./metadata/ukb-b-14043_mapped_metadata.txt
├── [  72]  ./MultiQC
│   ├── [5.5M]  ./MultiQC/multiqc_report.html
│   └── [5.5M]  ./MultiQC/ukb-b-14043_multiqc_report.html
├── [ 189]  ./pipeline_info
│   ├── [2.8M]  ./pipeline_info/execution_report_2022-09-09_16-07-25.html
│   ├── [8.5K]  ./pipeline_info/execution_timeline_2022-09-09_16-07-25.html
│   ├── [1.6K]  ./pipeline_info/execution_trace_2022-09-09_16-29-32.txt
│   └── [1.5K]  ./pipeline_info/pipeline_metadata_report.tsv
└── [  25]  ./QC_plots
    └── [ 326]  ./QC_plots/ukb-b-14043
        ├── [ 34K]  ./QC_plots/ukb-b-14043/boxPlot-maf-beta.png
        ├── [ 31K]  ./QC_plots/ukb-b-14043/boxPlot-maf-p.png
        ├── [4.5K]  ./QC_plots/ukb-b-14043/descriptive_statistics.tsv
        ├── [9.4K]  ./QC_plots/ukb-b-14043/descriptive_statistics.txt
        ├── [ 24K]  ./QC_plots/ukb-b-14043/frequencyCurve-Absolute_beta.png
        ├── [ 33K]  ./QC_plots/ukb-b-14043/frequencyCurve-MAF.png
        ├── [ 47K]  ./QC_plots/ukb-b-14043/Manhattan.png
        ├── [ 25K]  ./QC_plots/ukb-b-14043/QQ.png
        ├── [ 18K]  ./QC_plots/ukb-b-14043/scatterPlot-beta-logp.png
        ├── [ 23K]  ./QC_plots/ukb-b-14043/scatterPlot-maf-beta.png
        └── [ 19K]  ./QC_plots/ukb-b-14043/scatterPlot-maf-p.png

27 directories, 77 files
```

</details>
<br>



CloudOS Example run: https://staging.lifebit.ai/app/jobs/631b644db2b00e0155866a60

</details>

#### Harmonising GWAS summary statistics from GWAS Tables input

```bash
nextflow run main.nf -profile test_gwas_tables_all -with-docker
nextflow run main.nf -profile test_gwas_tables_all -with-singularity
cloudos job run --nextflow-profile test_gwas_tables_all,network_cloudos $CLOUDOS_CLI_OPTIONS
```
</details>
<br>

CloudOS Example run: https://staging.lifebit.ai/app/jobs/631b68dab2b00e01558675c3


<details>
<summary>Expected output:</summary>

```
# Full [$ tree -fh results/] output reveals: 236 directories, 523 files
# Below is a reduced output view to compare against using [-L 3] tree  parameter
$ tree -fh -L 3 results
results
├── [4.0K]  results/conversion
│   ├── [1.2M]  results/conversion/allancs-gwas_bin-BOLT-LMM_harmonised_conv_sumstats.vcf
│   ├── [1.1M]  results/conversion/allancs-gwas_bin-fastGWA-GLMM_harmonised_conv_sumstats.vcf
│   ├── [1.0M]  results/conversion/allancs-gwas_bin-PLINK2-GLM_harmonised_conv_sumstats.vcf
│   ├── [1.1M]  results/conversion/allancs-gwas_bin-regenie_harmonised_conv_sumstats.vcf
│   ├── [1.2M]  results/conversion/allancs-gwas_bin-SAIGE_harmonised_conv_sumstats.vcf
│   ├── [1.2M]  results/conversion/allancs-gwas_qt-BOLT-LMM_harmonised_conv_sumstats.vcf
│   ├── [1.1M]  results/conversion/allancs-gwas_qt-fastGWA-GLMM_harmonised_conv_sumstats.vcf
│   ├── [1.0M]  results/conversion/allancs-gwas_qt-PLINK2-GLM_harmonised_conv_sumstats.vcf
│   ├── [1.1M]  results/conversion/allancs-gwas_qt-regenie_harmonised_conv_sumstats.vcf
│   ├── [1.2M]  results/conversion/allancs-gwas_qt-SAIGE_harmonised_conv_sumstats.vcf
│   ├── [362K]  results/conversion/allancs-notransform-Hail-GWAS_harmonised_conv_sumstats.vcf
│   └── [1.4M]  results/conversion/European-NA-METAL_harmonised_conv_sumstats.vcf
├── [4.0K]  results/gwas_vcf_input
│   ├── [892K]  results/gwas_vcf_input/allancs-gwas_bin-BOLT-LMM_harmonised_sumstats.vcf
│   ├── [772K]  results/gwas_vcf_input/allancs-gwas_bin-fastGWA-GLMM_harmonised_sumstats.vcf
│   ├── [738K]  results/gwas_vcf_input/allancs-gwas_bin-PLINK2-GLM_harmonised_sumstats.vcf
│   ├── [798K]  results/gwas_vcf_input/allancs-gwas_bin-regenie_harmonised_sumstats.vcf
│   ├── [914K]  results/gwas_vcf_input/allancs-gwas_bin-SAIGE_harmonised_sumstats.vcf
│   ├── [871K]  results/gwas_vcf_input/allancs-gwas_qt-BOLT-LMM_harmonised_sumstats.vcf
│   ├── [752K]  results/gwas_vcf_input/allancs-gwas_qt-fastGWA-GLMM_harmonised_sumstats.vcf
│   ├── [733K]  results/gwas_vcf_input/allancs-gwas_qt-PLINK2-GLM_harmonised_sumstats.vcf
│   ├── [795K]  results/gwas_vcf_input/allancs-gwas_qt-regenie_harmonised_sumstats.vcf
│   ├── [858K]  results/gwas_vcf_input/allancs-gwas_qt-SAIGE_harmonised_sumstats.vcf
│   ├── [251K]  results/gwas_vcf_input/allancs-notransform-Hail-GWAS_harmonised_sumstats.vcf
│   ├── [996K]  results/gwas_vcf_input/European-NA-METAL_harmonised_sumstats.vcf
│   └── [  17]  results/gwas_vcf_input/study_id
├── [4.0K]  results/harmonised
│   ├── [421K]  results/harmonised/allancs-gwas_bin-BOLT-LMM.harmonised.gwas.vcf
│   ├── [385K]  results/harmonised/allancs-gwas_bin-fastGWA-GLMM.harmonised.gwas.vcf
│   ├── [354K]  results/harmonised/allancs-gwas_bin-PLINK2-GLM.harmonised.gwas.vcf
│   ├── [383K]  results/harmonised/allancs-gwas_bin-regenie.harmonised.gwas.vcf
│   ├── [5.7K]  results/harmonised/allancs-gwas_bin-SAIGE.harmonised.gwas.vcf
│   ├── [235K]  results/harmonised/allancs-gwas_qt-BOLT-LMM.harmonised.gwas.vcf
│   ├── [221K]  results/harmonised/allancs-gwas_qt-fastGWA-GLMM.harmonised.gwas.vcf
│   ├── [195K]  results/harmonised/allancs-gwas_qt-PLINK2-GLM.harmonised.gwas.vcf
│   ├── [224K]  results/harmonised/allancs-gwas_qt-regenie.harmonised.gwas.vcf
│   ├── [4.3K]  results/harmonised/allancs-gwas_qt-SAIGE.harmonised.gwas.vcf
│   ├── [118K]  results/harmonised/allancs-notransform-Hail-GWAS.harmonised.gwas.vcf
│   ├── [300K]  results/harmonised/European-NA-METAL.harmonised.gwas.vcf
│   ├── [4.0K]  results/harmonised/filtered
│   │   ├── [421K]  results/harmonised/filtered/allancs-gwas_bin-BOLT-LMM_harmonised_filtered_sumstats.vcf
│   │   ├── [385K]  results/harmonised/filtered/allancs-gwas_bin-fastGWA-GLMM_harmonised_filtered_sumstats.vcf
│   │   ├── [354K]  results/harmonised/filtered/allancs-gwas_bin-PLINK2-GLM_harmonised_filtered_sumstats.vcf
│   │   ├── [383K]  results/harmonised/filtered/allancs-gwas_bin-regenie_harmonised_filtered_sumstats.vcf
│   │   ├── [5.7K]  results/harmonised/filtered/allancs-gwas_bin-SAIGE_harmonised_filtered_sumstats.vcf
│   │   ├── [235K]  results/harmonised/filtered/allancs-gwas_qt-BOLT-LMM_harmonised_filtered_sumstats.vcf
│   │   ├── [221K]  results/harmonised/filtered/allancs-gwas_qt-fastGWA-GLMM_harmonised_filtered_sumstats.vcf
│   │   ├── [195K]  results/harmonised/filtered/allancs-gwas_qt-PLINK2-GLM_harmonised_filtered_sumstats.vcf
│   │   ├── [224K]  results/harmonised/filtered/allancs-gwas_qt-regenie_harmonised_filtered_sumstats.vcf
│   │   ├── [4.3K]  results/harmonised/filtered/allancs-gwas_qt-SAIGE_harmonised_filtered_sumstats.vcf
│   │   ├── [118K]  results/harmonised/filtered/allancs-notransform-Hail-GWAS_harmonised_filtered_sumstats.vcf
│   │   └── [300K]  results/harmonised/filtered/European-NA-METAL_harmonised_filtered_sumstats.vcf
│   └── [4.0K]  results/harmonised/hail
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_bin-BOLT-LMM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_bin-fastGWA-GLMM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_bin-PLINK2-GLM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_bin-regenie_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_bin-SAIGE_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_qt-BOLT-LMM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_qt-fastGWA-GLMM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_qt-PLINK2-GLM_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_qt-regenie_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-gwas_qt-SAIGE_hail_matrix.mt
│       ├── [4.0K]  results/harmonised/hail/allancs-notransform-Hail-GWAS_hail_matrix.mt
│       └── [4.0K]  results/harmonised/hail/European-NA-METAL_hail_matrix.mt
├── [4.0K]  results/MultiQC
│   ├── [5.9M]  results/MultiQC/allancs-gwas_bin-BOLT-LMM_multiqc_report.html
│   ├── [5.9M]  results/MultiQC/allancs-gwas_bin-fastGWA-GLMM_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_bin-PLINK2-GLM_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_bin-regenie_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_bin-SAIGE_multiqc_report.html
│   ├── [5.9M]  results/MultiQC/allancs-gwas_qt-BOLT-LMM_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_qt-fastGWA-GLMM_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_qt-PLINK2-GLM_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_qt-regenie_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-gwas_qt-SAIGE_multiqc_report.html
│   ├── [5.8M]  results/MultiQC/allancs-notransform-Hail-GWAS_multiqc_report.html
│   ├── [5.7M]  results/MultiQC/European-NA-METAL_multiqc_report.html
│   └── [5.7M]  results/MultiQC/multiqc_report.html
├── [4.0K]  results/pipeline_info
│   ├── [3.0M]  results/pipeline_info/execution_report_2022-11-08_11-21-43.html
│   ├── [271K]  results/pipeline_info/execution_timeline_2022-11-08_11-21-43.html
│   ├── [ 12K]  results/pipeline_info/execution_trace_2022-11-08_11-21-43.txt
│   ├── [   3]  results/pipeline_info/execution_trace_2022-11-08_11-27-51.txt
│   └── [ 915]  results/pipeline_info/pipeline_metadata_report.tsv
└── [4.0K]  results/QC_plots
    ├── [4.0K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/boxPlot-maf-p.png
    │   ├── [ 623]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/descriptive_statistics.txt
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/frequencyCurve-MAF.png
    │   ├── [193K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/Manhattan.png
    │   ├── [ 35K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/QQ.png
    │   ├── [ 27K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/scatterPlot-beta-logp.png
    │   ├── [ 42K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/scatterPlot-maf-beta.png
    │   └── [ 53K]  results/QC_plots/allancs-gwas_bin-BOLT-LMM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/boxPlot-maf-beta.png
    │   ├── [ 33K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/boxPlot-maf-p.png
    │   ├── [ 620]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/descriptive_statistics.txt
    │   ├── [ 35K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/frequencyCurve-Absolute_beta.png
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/frequencyCurve-MAF.png
    │   ├── [207K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/Manhattan.png
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/QQ.png
    │   ├── [ 24K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/scatterPlot-beta-logp.png
    │   ├── [ 43K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/scatterPlot-maf-beta.png
    │   └── [ 40K]  results/QC_plots/allancs-gwas_bin-fastGWA-GLMM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/boxPlot-maf-beta.png
    │   ├── [ 33K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/boxPlot-maf-p.png
    │   ├── [ 621]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/descriptive_statistics.txt
    │   ├── [ 35K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/frequencyCurve-Absolute_beta.png
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/frequencyCurve-MAF.png
    │   ├── [204K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/QQ.png
    │   ├── [ 23K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/scatterPlot-beta-logp.png
    │   ├── [ 29K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/scatterPlot-maf-beta.png
    │   └── [ 41K]  results/QC_plots/allancs-gwas_bin-PLINK2-GLM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_bin-regenie
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-regenie/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_bin-regenie/boxPlot-maf-p.png
    │   ├── [ 621]  results/QC_plots/allancs-gwas_bin-regenie/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_bin-regenie/descriptive_statistics.txt
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_bin-regenie/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-regenie/frequencyCurve-MAF.png
    │   ├── [205K]  results/QC_plots/allancs-gwas_bin-regenie/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-regenie/QQ.png
    │   ├── [ 23K]  results/QC_plots/allancs-gwas_bin-regenie/scatterPlot-beta-logp.png
    │   ├── [ 45K]  results/QC_plots/allancs-gwas_bin-regenie/scatterPlot-maf-beta.png
    │   └── [ 39K]  results/QC_plots/allancs-gwas_bin-regenie/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_bin-SAIGE
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_bin-SAIGE/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_bin-SAIGE/boxPlot-maf-p.png
    │   ├── [ 621]  results/QC_plots/allancs-gwas_bin-SAIGE/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_bin-SAIGE/descriptive_statistics.txt
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-SAIGE/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-SAIGE/frequencyCurve-MAF.png
    │   ├── [202K]  results/QC_plots/allancs-gwas_bin-SAIGE/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_bin-SAIGE/QQ.png
    │   ├── [ 23K]  results/QC_plots/allancs-gwas_bin-SAIGE/scatterPlot-beta-logp.png
    │   ├── [ 41K]  results/QC_plots/allancs-gwas_bin-SAIGE/scatterPlot-maf-beta.png
    │   └── [ 41K]  results/QC_plots/allancs-gwas_bin-SAIGE/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/boxPlot-maf-p.png
    │   ├── [ 621]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/descriptive_statistics.txt
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/frequencyCurve-MAF.png
    │   ├── [195K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/Manhattan.png
    │   ├── [ 35K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/QQ.png
    │   ├── [ 26K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/scatterPlot-beta-logp.png
    │   ├── [ 38K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/scatterPlot-maf-beta.png
    │   └── [ 54K]  results/QC_plots/allancs-gwas_qt-BOLT-LMM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM
    │   ├── [ 35K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/boxPlot-maf-beta.png
    │   ├── [ 33K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/boxPlot-maf-p.png
    │   ├── [ 618]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/descriptive_statistics.txt
    │   ├── [ 38K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/frequencyCurve-Absolute_beta.png
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/frequencyCurve-MAF.png
    │   ├── [206K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/QQ.png
    │   ├── [ 22K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/scatterPlot-beta-logp.png
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/scatterPlot-maf-beta.png
    │   └── [ 38K]  results/QC_plots/allancs-gwas_qt-fastGWA-GLMM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/boxPlot-maf-beta.png
    │   ├── [ 33K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/boxPlot-maf-p.png
    │   ├── [ 624]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/descriptive_statistics.txt
    │   ├── [ 38K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/frequencyCurve-Absolute_beta.png
    │   ├── [ 37K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/frequencyCurve-MAF.png
    │   ├── [207K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/QQ.png
    │   ├── [ 25K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/scatterPlot-beta-logp.png
    │   ├── [ 30K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/scatterPlot-maf-beta.png
    │   └── [ 43K]  results/QC_plots/allancs-gwas_qt-PLINK2-GLM/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_qt-regenie
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_qt-regenie/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_qt-regenie/boxPlot-maf-p.png
    │   ├── [ 619]  results/QC_plots/allancs-gwas_qt-regenie/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_qt-regenie/descriptive_statistics.txt
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-regenie/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-regenie/frequencyCurve-MAF.png
    │   ├── [205K]  results/QC_plots/allancs-gwas_qt-regenie/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-regenie/QQ.png
    │   ├── [ 21K]  results/QC_plots/allancs-gwas_qt-regenie/scatterPlot-beta-logp.png
    │   ├── [ 31K]  results/QC_plots/allancs-gwas_qt-regenie/scatterPlot-maf-beta.png
    │   └── [ 41K]  results/QC_plots/allancs-gwas_qt-regenie/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-gwas_qt-SAIGE
    │   ├── [ 34K]  results/QC_plots/allancs-gwas_qt-SAIGE/boxPlot-maf-beta.png
    │   ├── [ 32K]  results/QC_plots/allancs-gwas_qt-SAIGE/boxPlot-maf-p.png
    │   ├── [ 619]  results/QC_plots/allancs-gwas_qt-SAIGE/descriptive_statistics.tsv
    │   ├── [9.1K]  results/QC_plots/allancs-gwas_qt-SAIGE/descriptive_statistics.txt
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-SAIGE/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-gwas_qt-SAIGE/frequencyCurve-MAF.png
    │   ├── [207K]  results/QC_plots/allancs-gwas_qt-SAIGE/Manhattan.png
    │   ├── [ 40K]  results/QC_plots/allancs-gwas_qt-SAIGE/QQ.png
    │   ├── [ 21K]  results/QC_plots/allancs-gwas_qt-SAIGE/scatterPlot-beta-logp.png
    │   ├── [ 33K]  results/QC_plots/allancs-gwas_qt-SAIGE/scatterPlot-maf-beta.png
    │   └── [ 44K]  results/QC_plots/allancs-gwas_qt-SAIGE/scatterPlot-maf-p.png
    ├── [4.0K]  results/QC_plots/allancs-notransform-Hail-GWAS
    │   ├── [ 40K]  results/QC_plots/allancs-notransform-Hail-GWAS/boxPlot-maf-beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-notransform-Hail-GWAS/boxPlot-maf-p.png
    │   ├── [4.4K]  results/QC_plots/allancs-notransform-Hail-GWAS/descriptive_statistics.tsv
    │   ├── [9.6K]  results/QC_plots/allancs-notransform-Hail-GWAS/descriptive_statistics.txt
    │   ├── [ 34K]  results/QC_plots/allancs-notransform-Hail-GWAS/frequencyCurve-Absolute_beta.png
    │   ├── [ 36K]  results/QC_plots/allancs-notransform-Hail-GWAS/frequencyCurve-MAF.png
    │   ├── [153K]  results/QC_plots/allancs-notransform-Hail-GWAS/Manhattan.png
    │   ├── [ 36K]  results/QC_plots/allancs-notransform-Hail-GWAS/QQ.png
    │   ├── [ 21K]  results/QC_plots/allancs-notransform-Hail-GWAS/scatterPlot-beta-logp.png
    │   ├── [ 41K]  results/QC_plots/allancs-notransform-Hail-GWAS/scatterPlot-maf-beta.png
    │   └── [ 47K]  results/QC_plots/allancs-notransform-Hail-GWAS/scatterPlot-maf-p.png
    └── [4.0K]  results/QC_plots/European-NA-METAL
        ├── [ 34K]  results/QC_plots/European-NA-METAL/boxPlot-maf-beta.png
        ├── [ 35K]  results/QC_plots/European-NA-METAL/boxPlot-maf-p.png
        ├── [4.2K]  results/QC_plots/European-NA-METAL/descriptive_statistics.tsv
        ├── [9.5K]  results/QC_plots/European-NA-METAL/descriptive_statistics.txt
        ├── [ 31K]  results/QC_plots/European-NA-METAL/frequencyCurve-Absolute_beta.png
        ├── [ 30K]  results/QC_plots/European-NA-METAL/frequencyCurve-MAF.png
        ├── [ 99K]  results/QC_plots/European-NA-METAL/Manhattan.png
        ├── [ 35K]  results/QC_plots/European-NA-METAL/QQ.png
        ├── [ 25K]  results/QC_plots/European-NA-METAL/scatterPlot-beta-logp.png
        ├── [ 28K]  results/QC_plots/European-NA-METAL/scatterPlot-maf-beta.png
        └── [ 46K]  results/QC_plots/European-NA-METAL/scatterPlot-maf-p.png

32 directories, 199 files
```

</details>
<br>

<details>
<summary>Format of main outputs:</summary>

```shell
head -n 60 results/harmonised/allancs-gwas_bin-SAIGE.harmonised.gwas.vcf
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20221011
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=AC,Number=A,Type=Integer,Description="Alternative allele count in the trait subset">
##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Variant-specific number of samples/individuals with called genotypes used to test association with specified trait">
##FORMAT=<ID=NC,Number=A,Type=Integer,Description="Variant-specific number of cases used to estimate genetic effect (binary traits only)">
##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">
##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">
##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">
##phasing=unphased
##source=gwas-sumstats-harmonisation-nf
##GenomeBuild=GRCh38
##META=<ID=Population,Number=1,Type=String,Description="Ancestral population of individuals in the association study">
##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">
##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">
##META=<ID=TotalSamples,Number=1,Type=Integer,Description="Total number of Samples in the association study">
##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=Method,Number=1,Type=String,Description="GWAS analysis method">
##META=<ID=GeneticModel,Number=1,Type=String,Description="Genetic model used for GWAS association test">
##META=<ID=Covariates,Number=1,Type=String,Description="Covariates used in association analysis">
##META=<ID=TraitReported,Number=1,Type=String,Description="Trait name used as input for trait search and mapping">
##SAMPLE=<ID=allancs-gwas_bin-SAIGE,TotalCases=43,TotalControls=297,Population="allancs",StudyType="CaseControl",Method="SAIGE",TotalSamples=297,GeneticModel="additive",Covariates="QUANT_TRAIT, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,",TraitReported="BIN_TRAIT",TotalVariants=17583,HarmonisedVariants=8840>
##contig=<ID=1>
##contig=<ID=2>
##bcftools_filterVersion=7cd83b7+htslib-
##bcftools_filterCommand=filter '-e FORMAT/ES[*]<-0.3' temp_sumstats.vcf; Date=Tue Oct 11 14:46:06 2022
##bcftools_filterCommand=filter '-e FORMAT/ES[*]>0.3' temp_sumstats.vcf; Date=Tue Oct 11 14:46:06 2022
##bcftools_filterCommand=filter '-e FORMAT/LP[*]>0.5' temp_sumstats.vcf; Date=Tue Oct 11 14:46:06 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]<0.3' temp_sumstats.vcf; Date=Tue Oct 11 14:46:06 2022
##bcftools_filterCommand=filter '-e FORMAT/AF[*]>0.8' temp_sumstats.vcf; Date=Tue Oct 11 14:46:06 2022
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  allancs-gwas_bin-SAIGE
1       742429  rs3094315       G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.292615:211:0.408397:0.689706:0.25128:340:43:0.012806:0.014913:1.28567
1       1011278 rs3737728       A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.358479:136:0.175753:0.8:-0.15415:340:43:-0.007856:0.018269:0.857143
1       1125105 rs9729550       A       C       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.292238:293:0.254509:0.430882:-0.171836:340:43:-0.008757:0.014893:0.842117
1       2146222 rs2460000       A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.30221:220:0.239965:0.676471:0.169234:340:43:0.008625:0.015402:1.1844
1       2172330 rs260512        G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.293569:294:0.212048:0.567647:-0.148196:340:43:-0.007553:0.014961:0.862262
1       2440254 rs2246732       A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.304581:245:0.180724:0.360294:0.13416:340:43:0.006837:0.015522:1.14358
1       2699024 rs4648356       C       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.305976:239:0.192401:0.351471:-0.142209:340:43:-0.007247:0.015593:0.86744
1       2819411 rs1869970       A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.301866:269:0.34358:0.395588:-0.226359:340:43:-0.011536:0.015384:0.797432
1       2843619 rs16823228      A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.29221:253:0.00032943:0.372059:0.000277696:340:43:1.4e-05:0.014892:1.00028
1       3099011 rs7539775       G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.305952:274:0.00240135:0.597059:0.0021144:340:43:0.000108:0.015592:1.00212
1       3248973 rs882430        G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.298285:237:0.129643:0.348529:-0.0982289:340:43:-0.005006:0.015201:0.906441
1       3273388 rs4648485       G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.313251:324:0.300929:0.476471:0.211227:340:43:0.010765:0.015964:1.23519
1       3346265 rs2483250       A       C       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.356191:144:0.10603:0.788235:0.0979248:340:43:0.004991:0.018153:1.10288
1       3562634 rs2181484       A       G       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.351177:172:0.0123686:0.747059:-0.0123606:340:43:-0.00063:0.017897:0.987715
1       3770956 rs10797348      C       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.311846:258:0.182419:0.620588:0.138468:340:43:0.007057:0.015893:1.14851
1       3939005 rs4654428       A       C       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.31862:263:0.277166:0.613235:-0.20095:340:43:-0.010241:0.016238:0.817953
1       4160278 rs6667065       G       A       .       PASS    .       SE:AC:LP:AF:ES:SS:NC:SES:SSE:OR 0.301578:216:0.0317345:0.317647:0.0266687:340:43:0.001359:0.015369:1.02703
```

</details>

## 4 - <ins>Additional information</ins>

### 4.1 - <ins>Stress testing</ins>
## Stress testing

bi-gwas-sumstats-harmonisation-nf pipeline has been tested with 2, 10, 50, 100 studies using `--gwas_source = "ieu"` and `--gwas_source = "ebi"` modes.

Output:

| Test config    | N studies | GWAS source |  Runtime\* | Instance type\*\* |
| -------------- | ----------- |----------- | ----------- |----------- |
| `conf/stress_test/test_ieu_2.config`     | 2           | IEU    |    00:52:40         |   `c3.4xlarge`          |
| `conf/stress_test/test_ieu_10.config`  | 10           |   IEU |   01:18:46        |    `c3.4xlarge`         |
| `conf/stress_test/test_ieu_50.config`     | 50           |      IEU |  01:29:52     |    `c4.4xlarge`         |
| `conf/stress_test/test_ieu_100.config`     | 100           |   IEU   |    04:18:50      |    `c4.4xlarge`         |
| `conf/stress_test/test_ebi_2.config`     | 2           |   EBI |   00:21:34       |    `c3.4xlarge`         |
| `conf/stress_test/test_ebi_10.config`     | 10           |   EBI |   00:40:52       |    `c3.4xlarge`         |
| `conf/stress_test/test_ebi_50.config`     | 50           |   EBI |   02:53:46       |     `c3.4xlarge`        |

\* The runtime does not scale linearly with the number of variants because the analyses were run in parallel, and on different timepoints more than one AWS instances were utilised. This autoscaling is managed by the Nextflow executor.

\*\* Due to constraints of instance type capacity, the instance type varies between `c3.4xlarge` and `c4.4xlarge`.
