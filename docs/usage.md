# bi-gwas-sumstats-harmonisation-nf

<!-- This README.md is the single page user documentation for this pipeline. -->

This pipeline is designed to fetch and harmonise GWAS summary statistics from [EBI GWAS catalog](https://www.ebi.ac.uk/gwas/) and [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/).

## Pipeline description

## Input

This pipeline accepts as an input a newline separated text file, with a study ID per line for EBI and IEU database retrieval, or paths to GWAS files in the case of GWAS VCF and GWAS table inputs.

Example from EBI:
```
GCST001969
GCST90077560
```

Example from IEU:
```
ukb-b-14043
ieu-a-297
```

Example for GWAS VCF:
```
s3://lifebit-featured-datasets/pipelines/downstream-omics/sumstats-harmonisation/testdata/mini.gwas.vcf.gz
```

Example for GWAS tables:
```
s3://lifebit-featured-datasets/pipelines/downstream-omics/sumstats-harmonisation/testdata/gwas_table_inputs/allancs-gwas_bin-bolt_lmm.tsv
```

## Output

The expected output is a harmonised GWAS summary statistics VCF file, with metadata from the original resources added as part of the header.

Example output (few lines):

```shell
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
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  GCST90077560
1       99851033        rs2307130       A       G       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00194:0.109883:0.50213:0.000606106:408215:0.000185:0.305996:1.00061
2       70935956        rs11681642      T       C       .       PASS    .       SE:LP:AF:ES:SS:SES:SSE:OR       1.00026:0.0910552:0.445335:-0.000514012:408215:-0.000157:0.305483:0.999486
```

## Usage

Typical command for running the pipeline for EBI studies:

```
nextflow run main.nf \
    --input "testdata/input_ebi.txt" \
    --gwas_source "ebi" \
    --standardise true \
    --coef_conversion true \
    --keep_intermediate_files true \
    -with-docker
```

Typical command for running the pipeline for IEU studies:

```
nextflow run main.nf \
    --input "testdata/input_ieu.txt" \
    --gwas_source "ieu" \
    --standardise true \
    --coef_conversion true \
    --keep_intermediate_files true \
    -with-docker
```

Typical command for running the pipeline with GWAS VCF input:

```
nextflow run main.nf \
    --input "testdata/input_gwas_vcf.txt" \
    --gwas_source "gwas_vcf" \
    --force_munge false \
    --standardise true \
    --coef_conversion true \
    --keep_intermediate_files true \
    -with-docker
```

Typical command for running the pipeline with GWAS tables input:

```
nextflow run main.nf \
    --input "testdata/input_gwas_tables.txt" \
    --gwas_source "gwas_tables" \
    --standardise true \
    --coef_conversion true \
    --keep_intermediate_files true \
    -with-docker
```

All commands will produce a very similar output to the example above.

## Options

- **`--input`**: path or link to the input file. [Required].
- **`--gwas_source`**: GWAS source from where to fetch the data. It should be one of the following supported strings: 'ebi', 'ieu'. [Required].
- **`--standardise`**: whether to perform the BETA and SE standardisation. Default : `false`.
- **`--coef_conversion`**: whether to perform the coefficient conversion, from BETA to Odds Ratio. Default : `false`.
- **`--miss_percent_allow`**: missigness percentage allowed for a column from EBI data. Each column above this percentage will be dropped. Default : `10`.
- **`--keep_intermediate_files`**: whether to keep intermediate files. Default : `false`.
- **`--omop_vocabulary`**: path or link to an OMOP vocabulary DB file. Default : `https://omopvocabs.s3.eu-west-1.amazonaws.com/omop_vocabulary.db`.
- **`--outdir`**: output directory. Default: `results`.
- **`--take_n_studies`**: the number of studies to actually process from the input file. Useful for reduce the computation time for large input files, as only the indicated first `n` studies of the input will be used. Default : `-1` (take **all** the studies from input).
- **`--convert_to_hail`**: whether to convert the harmonised or harmonised and filtered VCF to `Hail` `MatrixTable` format and store it in <outdir>/harmonised/hail. Default : `false`.
- **`--force_munge`**: Force munging/harmonisation of input when `--gwas_source` is `gwas_vcf`. When 'false', assume GWAS VCF input is already harmonised and skip munging step. Default: `false`.
- **`--map_traits`**: Whether to map the traits in the input data to SNOMED terms using the OMOP common data model. Default: `false`.
- **`--dbsnp`**: Version of dbSNP database to use when harmonising variants. Supported versions: '144', '155'. Default: `155`.

Check the pipeline help section (`nextflow main.nf --help`) for all the updated options and their default values.

<!-- For Sphinx doc, This option will be auto rendered help() section from Nextflow main.nf in the doc build -->

## QC plotting

QC plots and descriptive statisitics for a GWAS sumstats table can be generated using the the script `bin/gwas_qc_plots.py` with the following command line options.
```
options:
  -h, --help            show this help message and exit
  --out_dir <directory>
                        Location of qc plots/files
  --sumstat <tsv/csv/vcf file>
                        Location of the sumstats file
  --study_identifier identifier <String>
                        GWAS study identifier prefix. Default: "CUSTOM_GWAS".
  --n-total N <Int>     N total
  --n-cases N <Int>     N cases
  --c-chrom C_CHROM     Name of chromosome table column. Default: "chrom".
  --c-pos C_POS         Name of SNP position table column. Default: "pos".
  --c-info C_INFO       Name of imputation info table column. Default: "imp_info".
  --c-maf C_MAF         Name of MAF table column. Default: "maf".
  --c-beta C_BETA       Name of beta/effect size table column. Default: "beta".
  --c-se C_SE           Name of beta/effect size standard error table column. Default: "se".
  --c-pval C_PVAL       Name of association test Pvalue table column. Default: "pval".
```

A docker container is available which contains all needed dependencies for this script. The follwing is an example for running the script using docker with a test dataset.

```shell
aws s3 cp s3://omics-example-datasets/pipelines/gwas-sumstats-harmonisation/saige_results.csv ./
docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -w $PWD -t quay.io/lifebitai/gwas_qc_plots \
python bin/gwas_qc_plots.py --sumstat saige_results.csv \
--out_dir testdata --study_identifier GWAS_test_saige \
--n-cases 11 --n-total 158  \
--c-chrom CHR --c-pos POS \
--c-info imputationInfo --c-maf AF_Allele2 \
--c-beta BETA --c-se SE --c-pval p.value
```

## VCF filtering

The harmonised VCF can be filtered using the following command line options:

- **`--filter_beta_smaller_than`**: Exclude variants with BETA smaller than the specified value. Default: `null`.
- **`--filter_beta_greater_than`**: Exclude variants with BETA greater than the specified value. Default: `null`.
- **`--filter_p_greater_than`**: Exclude variants with a P value greater than the specified value. Default: `null`.
- **`--filter_freq_smaller_than`**: Exclude variants with a alternate allele frequency smaller than the specified value. Default: `null`.
- **`--filter_freq_greater_than`**: Exclude variants with a alternate allele frequency greater than the specified value. Default: `null`.
- **`--filter_missing_info`**: Exclude variants with a missing value for INFO column. Default: `false`.
- **`--filter_info`**: Exclude variants with a INFO value smaller than the one specified. Default: `null`.

The following example will exclude all variants with `BETA` smaller than -0.3 or greater than 0.3, with `P` value greater than 0.5 and with alternative allele frequencies smaller than 0.3 or greater than 0.8:

```shell
nextflow run main.nf -profile test_ebi -with-docker
```

<!------------------
Build of this doc in github handle by - .github/workflows/build-deploy-doc.yml

To build this doc locally follow these steps.

Needs to have installed - 
1. sphinx
2. sphinx-rtd-theme
3. nextflow

Supposing your currently in base directory of the pipeline -
```
cd docs && bash src/pre-build.sh
cp README.md src
cd src && make html 
```
index.html will be generated in `docs/src/build/html` folder
-->
