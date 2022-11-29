#!/usr/bin/env nextflow
/*
========================================================================================
                         iudlgnt-gwas-sumstats-utils
========================================================================================
iudlgnt-gwas-sumstats-utils
----------------------------------------------------------------------------------------
*/

// Help message

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --input input.txt [Options]
    
    Inputs:
    --input                    Input file (path). It expects a study ID per line file.

    Options:
    --gwas_source              GWAS source from where to fetch the data. It should be one of
                               the following supported strings: 'ebi', 'ieu'.
                               (default: $params.gwas_source)
    --standardise              Whether to perform the BETA and SE standardisation (bool)
                               (default: $params.standardise)
    --coef_conversion          Whether to perform the coefficient conversion, from BETA
                               to Odds Ratio (bool)
                               (default: $params.coef_conversion) 
    --miss_percent_allow       Missigness percentage allowed for a column from EBI data. Each
                               column above this percentage will be dropped (int / float)
                               (default: $params.miss_percent_allow)
    --keep_intermediate_files  Whether to keep intermediate files (bool)
                               (default: $params.keep_intermediate_files)
    --omop_vocabulary          Path or link to an OMOP vocabulary DB file (path)
                               (default: $params.omop_vocabulary)
    --convert_to_hail          Whether to convert the harmonised VCF to Hail MatrixTable format (bool)
                               (default: $params.convert_to_hail)
    --med_memory               Memory assignment for the more memory intensive tasks (string)
                               (default: 9.GB)
    --take_n_studies           Take n studies from the given input file to run the pipeline (int)
                               (default: $params.take_n_studies)

    GWAS filtering options:
    --filter_beta_smaller_than  Exclude variants with BETA smaller than the specified value (float)
                                (Default: $params.filter_beta_smaller_than)
    --filter_beta_greater_than  Exclude variants with BETA greater than the specified value (float)
                                (Default: $params.filter_beta_greater_than)
    --filter_p_greater_than     Exclude variants with a P value greater than the specified value (float)
                                (Default: $params.filter_p_greater_than)
    --filter_freq_smaller_than  Exclude variants with a alternate allele frequency smaller than the specified
                                value. (float)
                                (Default: $params.filter_freq_smaller_than)
    --filter_freq_greater_than  Exclude variants with a alternate allele frequency greater than the specified
                                value. (float)
                                (Default: $params.filter_freq_greater_than)
    --filter_missing_info       Exclude variants with a missing value for INFO column (bool)
                                (Default: $params.filter_missing_info)
    --filter_info               Exclude variants with a INFO value smaller than the one specified (float).
                                (Default: $params.filter_info)

    Resource options:
    --med_memory       Memory assignment for the more memory intensive tasks (string)
                       (default: $params.med_memory)
    --max_cpus         Maximum number of CPUs (int)
                       (default: $params.max_cpus)  
    --max_memory       Maximum memory (memory unit)
                       (default: $params.max_memory)
    --max_time         Maximum time (time unit)
                       (default: $params.max_time)
    --outdir           Output directory(path)
                       (default: $params.outdir)
    See here for more info: https://github.com/lifebit-ai/bi-gwas-sumstats-harmonisation-nf/blob/dev/README.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}



/*--------------------------------------------------------
  Defining and showing header with all params information 
----------------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Output dir']                                  = params.outdir
summary['Launch dir']                                  = workflow.launchDir
summary['Working dir']                                 = workflow.workDir
summary['Script dir']                                  = workflow.projectDir
summary['User']                                        = workflow.userName
summary['Input file']                                  = params.input
summary['GWAS source']                                 = params.gwas_source
summary['Standardise BETA and SE']                     = params.standardise
summary['Coefficient conversion']                      = params.coef_conversion
summary['EBI missigness percent allowed']              = params.miss_percent_allow
summary['Keep intermediate files']                     = params.keep_intermediate_files
summary['OMOP vocabulary DB']                          = params.omop_vocabulary
summary['Convert to Hail']                             = params.convert_to_hail
summary['Filter BETA smaller than']                    = params.filter_beta_smaller_than
summary['Filter BETA greater than']                    = params.filter_beta_greater_than
summary['Filter P value greater than']                 = params.filter_p_greater_than
summary['Filter Alt AF smaller than']                  = params.filter_freq_smaller_than
summary['Filter Alt AF greater than']                  = params.filter_freq_greater_than
summary['Filter missing INFO']                         = params.filter_missing_info
summary['Filter INFO smaller than']                    = params.filter_info

log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



/*-------------------------------------------------
  Setting up introspection variables and channels  
----------------------------------------------------*/

// Importantly, in order to successfully introspect:
// - This needs to be done first `main.nf`, before any (non-head) nodes are launched. 
// - All variables to be put into channels in order for them to be available later in `main.nf`.

ch_repository         = Channel.of(workflow.manifest.homePage)
ch_commitId           = Channel.of(workflow.commitId ?: "Not available is this execution mode. Please run 'nextflow run ${workflow.manifest.homePage} [...]' instead of 'nextflow run main.nf [...]'")
ch_revision           = Channel.of(workflow.manifest.version)

ch_scriptName         = Channel.of(workflow.scriptName)
ch_scriptFile         = Channel.of(workflow.scriptFile)
ch_projectDir         = Channel.of(workflow.projectDir)
ch_launchDir          = Channel.of(workflow.launchDir)
ch_workDir            = Channel.of(workflow.workDir)
ch_userName           = Channel.of(workflow.userName)
ch_commandLine        = Channel.of(workflow.commandLine)
ch_configFiles        = Channel.of(workflow.configFiles)
ch_profile            = Channel.of(workflow.profile)
ch_container          = Channel.of(workflow.container)
ch_containerEngine    = Channel.of(workflow.containerEngine)



/*----------------------------------------------------------------
  Setting up additional variables used for documentation purposes  
-------------------------------------------------------------------*/

Channel
    .of(params.raci_owner)
    .set { ch_raci_owner } 

Channel
    .of(params.domain_keywords)
    .set { ch_domain_keywords }



/*----------------------
  Setting up input data  
-------------------------*/

// Check for gwas_source
if (params.gwas_source != 'ebi' && params.gwas_source != 'ieu') {
    exit 1, "Please, provide one of the valid --gwas_source parameters: ebi, ieu."
}

// Check if any filter was selected
if (params.filter_beta_smaller_than || params.filter_beta_greater_than || params.filter_p_greater_than || params.filter_freq_smaller_than || params.filter_freq_greater_than || params.filter_missing_info || params.filter_info) {
    filter_activated = true
} else {
    filter_activated = false
}

// Define channels from repository files

projectDir = workflow.projectDir

// Define Channels from input
ch_input = Channel.fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find input file containing study IDs: ${params.input}" }
    .splitCsv()
    .flatten()
    .take(params.take_n_studies) //default is -1 i.e. take all files (but param is useful for testing with fewer files)

// Channel for omop vocabulary
ch_omop_vocabulary = Channel.value(file(params.omop_vocabulary))

// Channels for scripts
ch_ebi_script = Channel.fromPath("${projectDir}/bin/ebi_fetch.sh")
ch_ieu_script = Channel.fromPath("${projectDir}/bin/ieu_fetch.sh")
ch_metadata_harmonisation_script = Channel.fromPath("${projectDir}/bin/metadata_harmonisation.py")
ch_munge_ebi_script = Channel.fromPath("${projectDir}/bin/munge_ebi.R")
ch_munge_ieu_script = Channel.fromPath("${projectDir}/bin/munge_ieu.R")
ch_conv_script = Channel.fromPath("${projectDir}/bin/convert_coeff.py")
ch_map_traits_script = Channel.fromPath("${projectDir}/bin/map_traits.R")
ch_qc_plots_script = Channel.fromPath("${projectDir}/bin/gwas_qc_plots.py")
ch_vcf2hail_script = Channel.fromPath("${projectDir}/bin/vcf2hail.py")

/*-----------
  Processes  
--------------*/

// Do not delete this process
// Create introspection report

process obtain_pipeline_metadata {
  publishDir "${params.tracedir}", mode: "copy"

  input:
  val repository from ch_repository
  val commit from ch_commitId
  val revision from ch_revision
  val script_name from ch_scriptName
  val script_file from ch_scriptFile
  val project_dir from ch_projectDir
  val launch_dir from ch_launchDir
  val work_dir from ch_workDir
  val user_name from ch_userName
  val command_line from ch_commandLine
  val config_files from ch_configFiles
  val profile from ch_profile
  val container from ch_container
  val container_engine from ch_containerEngine
  val raci_owner from ch_raci_owner
  val domain_keywords from ch_domain_keywords

  output:
  file("pipeline_metadata_report.tsv") into ch_pipeline_metadata_report
  
  shell:
  '''
  echo "Repository\t!{repository}"                  > temp_report.tsv
  echo "Commit\t!{commit}"                         >> temp_report.tsv
  echo "Revision\t!{revision}"                     >> temp_report.tsv
  echo "Script name\t!{script_name}"               >> temp_report.tsv
  echo "Script file\t!{script_file}"               >> temp_report.tsv
  echo "Project directory\t!{project_dir}"         >> temp_report.tsv
  echo "Launch directory\t!{launch_dir}"           >> temp_report.tsv
  echo "Work directory\t!{work_dir}"               >> temp_report.tsv
  echo "User name\t!{user_name}"                   >> temp_report.tsv
  echo "Command line\t!{command_line}"             >> temp_report.tsv
  echo "Configuration file(s)\t!{config_files}"    >> temp_report.tsv
  echo "Profile\t!{profile}"                       >> temp_report.tsv
  echo "Container\t!{container}"                   >> temp_report.tsv
  echo "Container engine\t!{container_engine}"     >> temp_report.tsv
  echo "RACI owner\t!{raci_owner}"                 >> temp_report.tsv
  echo "Domain keywords\t!{domain_keywords}"       >> temp_report.tsv

  awk 'BEGIN{print "Metadata_variable\tValue"}{print}' OFS="\t" temp_report.tsv > pipeline_metadata_report.tsv
  '''
}

if (params.gwas_source == 'ebi') {
    process fetch_from_ebi {
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ebi/", mode: 'copy'
        }

        input:
        val(study_id) from ch_input
        each file(ebi_script) from ch_ebi_script

        output:
        set val("${study_id}"), file("${study_id}_metadata_ebi.tsv") into ch_metadata
        set val("${study_id}"), file("*${study_id}*.h.tsv.gz"), file("${study_id}_sample_size.txt") into ch_sumstats

        script:
        """
        bash $ebi_script $study_id
        """
    }

    process munge_ebi {
        label 'munge'
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ebi/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats), file(sample_size) from ch_sumstats
        each file(munge_ebi_script) from ch_munge_ebi_script

        output:
        set val(study_id), file("${study_id}_harmonised_sumstats.vcf") into ch_harmonised_sumstats

        script:
        """
        N=\$(cat ${sample_size})
        Rscript $munge_ebi_script $sumstats "${study_id}_harmonised_sumstats.vcf" ${params.miss_percent_allow} \$N
        # Get rid of repeated lines (result is already sorted)
        uniq "${study_id}_harmonised_sumstats.vcf" > temp.vcf
        mv temp.vcf "${study_id}_harmonised_sumstats.vcf"
        """
    }
}

if (params.gwas_source == 'ieu') {
    process fetch_from_ieu {
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ieu/", mode: 'copy'
        }

        input:
        val(study_id) from ch_input
        each file(ieu_script) from ch_ieu_script

        output:
        set val("${study_id}"), file("${study_id}_metadata_ieu.tsv") into ch_metadata
        set val("${study_id}"), file("${study_id}.vcf.gz"), file("${study_id}_sample_size.txt") into ch_sumstats
 
        script:
        """
        bash $ieu_script $study_id
        """
    }

    process munge_ieu {
        label 'munge'
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ieu/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats), file(sample_size) from ch_sumstats
        each file(munge_ieu_script) from ch_munge_ieu_script

        output:
        set val(study_id), file("${study_id}_harmonised_sumstats.vcf") into ch_harmonised_sumstats

        script:
        """
        N=\$(cat ${sample_size})
        Rscript $munge_ieu_script $sumstats "${study_id}_harmonised_sumstats.vcf" \$N
        # Get rid of repeated lines (result is already sorted)
        uniq "${study_id}_harmonised_sumstats.vcf" > temp.vcf
        mv temp.vcf "${study_id}_harmonised_sumstats.vcf"
        """
    }
}

process metadata_harmonisation {
    label 'python'
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/metadata/", mode: 'copy'
    }

    input:
    set val(study_id), file(metadata) from ch_metadata
    each file(metadata_harmonisation_script) from ch_metadata_harmonisation_script

    output:
    set val("${study_id}"), file("${study_id}_harmonised_metadata.txt") into ch_harmonised_metadata

    script:
    """
    python $metadata_harmonisation_script -i $metadata -o ${study_id}_harmonised_metadata.txt
    """
}

if (params.standardise && params.coef_conversion) {
    process standardisation_beta2or {
        label 'python'
        publishDir "${params.outdir}/conversion/", mode: 'copy'

        input:
        set val(study_id), file(sumstats) from ch_harmonised_sumstats
        each file(conv_script) from ch_conv_script

        output:
        set val("${study_id}"), file("${study_id}_harmonised_conv_sumstats.vcf") into ch_conv_sumstats
        
        script:
        """
        # Detect the FRQ column necessary to `--standardise`
        FRQ_COLUMN_EXISTS=`awk '/^[^#]/ {if(\$9 ~ /FRQ/) a=1; else a=0}END{print a}' $sumstats`

        # Decide and convert
        if [ \$FRQ_COLUMN_EXISTS == "1" ]; then
            python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --standardise --beta2or
            echo "[MSG] BETA and SE standardisation performed."
            echo "[MSG] BETA to OR conversion performed."
        else
            python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --beta2or
            echo "[MSG] BETA to OR conversion performed."
            echo "[WARNING] BETA and SE standardisation not performed as FRQ value not present."
        fi
        """
    }
} else if (params.standardise) {
    process standardisation {
        label 'python'
        publishDir "${params.outdir}/conversion/", mode: 'copy'

        input:
        set val(study_id), file(sumstats) from ch_harmonised_sumstats
        each file(conv_script) from ch_conv_script

        output:
        set val("${study_id}"), file("${study_id}_harmonised_conv_sumstats.vcf") into ch_conv_sumstats
        
        script:
        """
        # Detect the FRQ column necessary to `--standardise`
        FRQ_COLUMN_EXISTS=`awk '/^[^#]/ {if(\$9 ~ /FRQ/) a=1; else a=0}END{print a}' $sumstats`

        # Decide and convert
        if [ \$FRQ_COLUMN_EXISTS == "1" ]; then
            python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --standardise
            echo "[MSG] BETA and SE standardisation performed."
        else
            cp $sumstats ${study_id}_harmonised_conv_sumstats.vcf
            echo "[WARNING] BETA and SE standardisation not performed as FRQ value not present."
        fi
        """
    }
} else if (params.coef_conversion) {
    process beta2or {
        label 'python'
        publishDir "${params.outdir}/conversion/", mode: 'copy'

        input:
        set val(study_id), file(sumstats) from ch_harmonised_sumstats
        each file(conv_script) from ch_conv_script

        output:
        set val("${study_id}"), file("${study_id}_harmonised_conv_sumstats.vcf") into ch_conv_sumstats
        
        script:
        """
        # BETA value should be always present
        python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --beta2or
        """
    }
} else {
    ch_conv_sumstats = ch_harmonised_sumstats
}

process map_traits {
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/metadata/", mode: 'copy'
    }

    input:
    set val(study_id), file(metadata) from ch_harmonised_metadata
    file(omop_vocabulary) from ch_omop_vocabulary
    each file(mapTraits) from ch_map_traits_script

    output:
    set val("${study_id}"), file("${study_id}_mapped_metadata.txt") into ch_mapped_metadata

    script:
    """
    # Collect traits from metadata
    sed -e 's/##meta=<//g' -e 's/>//g' $metadata > stripped_metadata.txt
    echo "traitName" > traits.csv
    if [ ${params.gwas_source} == 'ebi' ]; then
        awk -F '=' '{if(\$1 == "mappedLabel") print \$2}' stripped_metadata.txt >> traits.csv
    elif [ ${params.gwas_source} == 'ieu' ]; then
        awk -F '=' '{if(\$1 == "trait") print \$2}' stripped_metadata.txt >> traits.csv
    else
        echo "[ERROR] GWAS source not provided"
        exit 1
    fi
    # Map
    unzip ${omop_vocabulary}
    Rscript $mapTraits --traits_csv=traits.csv --vocabulary_folder=omop-vocab-files --snomed_grouping_level=4
    # Insert traits to the metadata
    sed '1d' mapped.csv | sed 's/trait_//g' | awk -F ',' '{print "##trait=<"\$2"=\\""\$3"\\">"}' > traits_mapped.txt
    cat $metadata traits_mapped.txt > "${study_id}_mapped_metadata.txt"
    """

}

ch_retrieve_data = ch_mapped_metadata.join(ch_conv_sumstats)

process insert_metadata {
    publishDir "${params.outdir}/harmonised/", mode: 'copy'

    input:
    set val(study_id), file(metadata), file(sumstats) from ch_retrieve_data

    output:
    set val(study_id), file("${study_id}_fully_harmonised_sumstats.vcf") into (ch_full_sumstats, ch_full_sumstats_2)

    script:
    """
    # Insert metadata into sumstats
    grep '^##' $sumstats > sumstats_header.txt
    # Correct some errors detected for some of the input files
    sed -i 's/ID=BETA,Number=1,Type=String/ID=BETA,Number=1,Type=Float/' sumstats_header.txt
    sed -i 's/ID=FRQ,Number=1,Type=String/ID=FRQ,Number=1,Type=Float/' sumstats_header.txt
    sed -i 's/ID=SE,Number=1,Type=String/ID=SE,Number=1,Type=Float/' sumstats_header.txt
    cat sumstats_header.txt $metadata > "${study_id}_full_header.txt"
    grep -v '^##' $sumstats > headless_sumstats.vcf
    cat "${study_id}_full_header.txt" headless_sumstats.vcf > "${study_id}_fully_harmonised_sumstats.vcf" 
    """
}

process make_qc_plots {
    publishDir "${params.outdir}/QC_plots", mode: 'copy'
    
    input:
    set val(study_id), file("fully_harmonised_sumstats.vcf") from ch_full_sumstats
    each file('gwas_qc_plots.py') from ch_qc_plots_script

    output:
    tuple val(study_id), file("${study_id}/*") into ch_qc_plots

    script:
    """
    # Check whether the VCF contains FRQ column
    if `tail -1 fully_harmonised_sumstats.vcf | grep -q "FRQ"`; then
        python gwas_qc_plots.py \
            --sumstat fully_harmonised_sumstats.vcf \
            --c-maf FRQ \
            --c-beta BETA \
            --c-se SE \
            --c-pval P \
            --out_dir ./ \
            --study_identifier ${study_id}
    else
        mkdir ${study_id}
        touch ${study_id}/${study_id}_no_plots.txt
    fi
    """
}

if (filter_activated) {
    process filter_sumstats {
        publishDir "${params.outdir}/harmonised/filtered/", mode: 'copy'
        echo true
        
        input:
        set val(study_id), file(sumstats) from ch_full_sumstats_2

        output:
        tuple val(study_id), file("${study_id}_harmonised_filtered_sumstats.vcf") into ch_filtered_sumstats

        script:
        beta_smaller_filter = params.filter_beta_smaller_than || params.filter_beta_smaller_than == 0 ? "-e FORMAT/BETA[*]<$params.filter_beta_smaller_than" : ""
        beta_greater_filter = params.filter_beta_greater_than || params.filter_beta_greater_than == 0 ? "-e FORMAT/BETA[*]>$params.filter_beta_greater_than" : ""
        p_filter = params.filter_p_greater_than || params.filter_p_greater_than == 0 ? "-e FORMAT/P[*]>$params.filter_p_greater_than" : ""
        frq_smaller_filter = params.filter_freq_smaller_than || params.filter_freq_smaller_than == 0 ? "-e FORMAT/FRQ[*]<$params.filter_freq_smaller_than" : ""
        frq_greater_filter = params.filter_freq_greater_than || params.filter_freq_greater_than == 0 ? "-e FORMAT/FRQ[*]>$params.filter_freq_greater_than" : ""
        missing_info_filter = params.filter_missing_info ? "FILTER" : "DO_NOTHING"
        info_filter = params.filter_info || params.filter_info == 0 ? params.filter_info : "NO_FILTER"
        """
        cp $sumstats temp_sumstats.vcf
        # Check for the FRQ column
        if tail -1 temp_sumstats.vcf | grep -q "FRQ"; then
            expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter" "$frq_smaller_filter" "$frq_greater_filter")
        else
            echo "[WARNING] FRQ column not found, all allele frequency filters will be ignored"
            expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter")
        fi
        # Apply filters sequentially, otherwise bcftools would apply OR logic
        for e in "\${expression_list[@]}"; do
            if [ "\$e" != "" ]; then
                bcftools filter "\$e" temp_sumstats.vcf > temp_sumstats_new.vcf
                mv temp_sumstats_new.vcf temp_sumstats.vcf
            fi
        done
        # Apply missing info filter if required
        if [ $missing_info_filter == "FILTER" ]; then
            awk '{if ((\$1 ~ /^#/) || (\$8 != ".")) { print \$0 }}' temp_sumstats.vcf > temp_sumstats_new.vcf
            mv temp_sumstats_new.vcf temp_sumstats.vcf
        fi
        # Apply info filter if required
        if [ $info_filter != "NO_FILTER" ]; then
            awk -v var="$info_filter" '{if ((\$1 ~ /^#/) || (\$8 >= var)) { print \$0 }}' temp_sumstats.vcf > temp_sumstats_new.vcf
            mv temp_sumstats_new.vcf temp_sumstats.vcf
        fi
        mv temp_sumstats.vcf "${study_id}_harmonised_filtered_sumstats.vcf"
        """
    }
} else {
    ch_filtered_sumstats = ch_full_sumstats_2
}

if (params.convert_to_hail) {
    process vcf2hail {
        publishDir "${params.outdir}/harmonised/hail/", mode: 'copy'

        input:
        set val(study_id), file(sumstats) from ch_filtered_sumstats
        each file(vcf2hail) from ch_vcf2hail_script

        output:
        tuple val(study_id), file("${study_id}_hail_matrix.mt/*") into ch_hail_sumstats

        script:
        """
        GENOME=\$( grep -e '##genome_build=' $sumstats | cut -d "=" -f 2 )
        python $vcf2hail -v $sumstats -g \$GENOME -o "${study_id}_hail_matrix.mt"
        """
    }
}
ch_report_dir = Channel.value(file("${projectDir}/bin/report"))
process create_report {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'
  
  input:
  set val(study_id), file("${study_id}/*") from ch_qc_plots
  file(report_dir) from ch_report_dir 

  output:
  file("${study_id}_multiqc_report.html") into ch_report_outputs_all
  file("multiqc_report.html") into ch_report_outputs1

  script:
  """
  cp -r ${report_dir}/* .
  cp ${study_id}/* .
  # Skip when no plots were generated due to lack of FRQ column
  if [ -f "${study_id}_no_plots.txt" ]; then
    touch multiqc_report.html
    touch ${study_id}_multiqc_report.html
  else
    for i in `ls *.png`; do name=`basename \$i .png | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
    for i in `ls *.tsv`; do name=`basename \$i .tsv | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
    cat file_list.txt | tr "\n" "," | sed 's/,\$//g' > file_list1.txt
    # Generates the report
    Rscript -e "rmarkdown::render('report.Rmd', params = list(`cat file_list1.txt`))"
    cp report.html multiqc_report.html
    mv report.html ${study_id}_multiqc_report.html
  fi
  """
}

// When the pipeline is not run locally
// Ensure trace report is output in the pipeline results (in 'pipeline_info' folder)

userName = workflow.userName

if ( userName == "ubuntu" || userName == "ec2-user") {
  workflow.onComplete {

  def trace_timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')

  traceReport = file("/home/${userName}/nf-out/trace.txt")
  traceReport.copyTo("results/pipeline_info/execution_trace_${trace_timestamp}.txt")
  }
}
