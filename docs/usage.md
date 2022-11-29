# bi-gwas-sumstats-harmonisation-nf

<!-- This README.md is the single page user documentation for this pipeline. -->

This pipeline is designed to fetch and harmonise GWAS summary statistics from [EBI GWAS catalog](https://www.ebi.ac.uk/gwas/) and [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/).

## Pipeline description

## Input

This pipeline accepts as an input a text file, with an study ID per line.

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

## Output

The expected output is a harmonised GWAS summary statistics VCF file, with metadata from the original resources added as part of the header.

Example output (few lines):

```shell
##fileformat=VCFv4.1
##fileDate=20220225
##phasing=unphased
##genome_build=GRCh38
##source=VariantAnnotation 1.40.0
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters">
##FORMAT=<ID=SNP,Number=1,Type=String,Description=".">
##FORMAT=<ID=BETA,Number=1,Type=Float,Description=".">
##FORMAT=<ID=FRQ,Number=1,Type=Float,Description=".">
##FORMAT=<ID=SE,Number=1,Type=Float,Description=".">
##FORMAT=<ID=P,Number=1,Type=Float,Description=".">
##FORMAT=<ID=STD_BETA,Number=1,Type=Float,Description="Standardised BETA">
##FORMAT=<ID=STD_BETA_SE,Number=1,Type=Float,Description="Standardised BETA SE">
##FORMAT=<ID=OR_,Number=1,Type=Float,Description="OR converted from BETA">
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
##META<pubmedId=23583979>
##META<title=Identification of heart rate-associated loci and their effects on cardiac conduction and rhythm disorders.>
##META<author=den Hoed M>
##META<publication=Nat Genet>
##META<publicationDate=2013-04-14T00:00:00Z>
##META<associationCount=32>
##META<note=NA>
##META<mr=NA>
##META<year=NA>
##META<group_name=NA>
##META<consortium=NA>
##META<sex=NA>
##META<priority=NA>
##META<population=NA>
##META<unit=NA>
##META<nsnp=NA>
##META<sample_size=NA>
##META<build=NA>
##META<id=GCST001969>
##META<subcategory=NA>
##META<category=NA>
##META<ontology=NA>
##META<trait=Heart rate>
##META<mappedLabel=heart rate>
##META<mappedUri=http://www.ebi.ac.uk/efo/EFO_0004326>
##trait=<reported="heart rate">
##trait=<source_origin="Inferred">
##trait=<source_concept_code="250764009">
##trait=<source_vocabulary="SNOMED">
##trait=<source_concept_id="40462158">
##trait=<source_concept_name="Heart rate">
##trait=<source_standard_concept="NA">
##trait=<concept_id="4239408">
##trait=<concept_code="364075005">
##trait=<concept_name="Heart rate">
##trait=<vocabulary="SNOMED">
##trait=<standard_concept="S">
##trait=<background_reported="NA">
##trait=<background_source_origin="NA">
##trait=<background_source_concept_code="NA">
##trait=<background_source_vocabulary="NA">
##trait=<background_source_concept_id="NA">
##trait=<background_source_concept_name="NA">
##trait=<background_concept_id="NA">
##trait=<background_concept_code="NA">
##trait=<background_concept_name="NA">
##trait=<background_vocabulary="NA">
##trait=<background_domain_name="NA">
##trait=<background_concept_class_name="NA">
##trait=<background_standard_concept="NA">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GWAS
1	785910	.	G	C	.	.	.	SNP:BETA:FRQ:SE:P:STD_BETA:STD_BETA_SE:OR_	rs12565286:0.1632:0.0609:0.2577:0.5266:0.002836:0.004479:1.344486
1	788439	.	T	A	.	.	.	SNP:BETA:FRQ:SE:P:STD_BETA:STD_BETA_SE:OR_	rs11804171:0.141:0.0611:0.2581:0.5848:0.002451:0.004486:1.291424
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

Both commands will produce a very similar output to the example above.

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

```
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
```

The following example will exclude all variants with `BETA` smaller than -0.3 or greater than 0.3, with `P` value greater than 0.5 and with alternative allele frequencies smaller than 0.3 or greater than 0.8:

```shell
nextflow run main.nf -profile test_ebi
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
