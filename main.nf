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
    --input                    Input file (path). Newline delimited list of IDs or file paths, or a 
                               single supported sumstats file.
    --input_type               Whether the file passed to --input should be treated as a file containing 
                               a newline delimited list of paths to sumstats files to process ('list') or
                               a single sumstats file to process ('single').
                               (default: $params.input_type) 
    --gwas_source              Type/source of GWAS input data. It should be one of the following supported
                               strings: 'ebi', 'ieu', 'gwas_vcf', 'gwas_table'.
    Options:
    --standardise              Whether to perform BETA and SE standardisation (bool).
                               (default: $params.standardise)
    --coef_conversion          Whether to perform the coefficient conversion, from BETA
                               to Odds Ratio (bool).
                               (default: $params.coef_conversion)
    --map_traits               Whether to map the traits in the input data to SNOMED terms using the OMOP
                               common data model (bool).
                               (default: $params.map_traits)
    --munge_method             Method to use for munging/harmonising input data. Either use the tool
                               MungeSumstats ('MSS'), or use a simplified method ('simple').
                               (default: $params.munge_method)
    --force_munge              Force munging/harmonisation of input when `--gwas_source` is
                               `gwas_vcf'. When 'false', assume GWAS VCF input is already harmonised 
                               and skip munging step. (bool).
                               (default: $params.force_munge)
    --miss_percent_allow       Missingness proportion allowed for a column from input data. Each
                               column above this proportion will be dropped (int / float).
                               (default: $params.miss_percent_allow)
    --keep_intermediate_files  Whether to keep intermediate files (bool).
                               (default: $params.keep_intermediate_files)
    --dbsnp                    Version of dbSNP database to use when harmonising variants. Supported
                               versions: '144', '155'.
                               (default: $params.dbsnp)
    --omop_vocabulary          Path or link to an OMOP vocabulary DB file (path).
                               (default: $params.omop_vocabulary)
    --convert_to_hail          Whether to convert the harmonised VCF to Hail MatrixTable format (bool).
                               (default: $params.convert_to_hail)
    --take_n_studies           Take n studies from the given input file to run the pipeline. Set to '-1' to use
                               all studies in the input file. (int).
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
    --filter_info               Exclude variants with a INFO/SI value smaller than the one specified (float).
                                (Default: $params.filter_info)

    Resource options:
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
summary['EBI missingness percent allowed']             = params.miss_percent_allow
summary['Keep intermediate files']                     = params.keep_intermediate_files
summary['OMOP vocabulary DB']                          = params.omop_vocabulary
summary['Convert to Hail']                             = params.convert_to_hail
summary['Filter BETA smaller than']                    = params.filter_beta_smaller_than
summary['Filter BETA greater than']                    = params.filter_beta_greater_than
summary['Filter P value greater than']                 = params.filter_p_greater_than
summary['Filter Alt AF smaller than']                  = params.filter_freq_smaller_than
summary['Filter Alt AF greater than']                  = params.filter_freq_greater_than
summary['Filter missing SI']                           = params.filter_missing_info
summary['Filter INFO/SI smaller than']                 = params.filter_info

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
if (!(params.gwas_source in ['ebi', 'ieu', 'gwas_vcf', 'gwas_table'])) {
    exit 1, "Please, provide one of the valid --gwas_source parameters: ebi, ieu, gwas_vcf, gwas_table."
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
if (params.gwas_source in ['ebi', 'ieu']) {
    ch_input = Channel.fromPath(params.input)
        .ifEmpty { exit 1, "Cannot find input file containing study IDs: ${params.input}" }
        .splitCsv()
        .flatten()
        .take(params.take_n_studies) //default is -1 i.e. take all files (but param is useful for testing with fewer files)
} else if (params.gwas_source in ['gwas_vcf'] && params.input_type == 'list' ) {
    ch_sumstats = Channel.fromPath(params.input)
        .ifEmpty { exit 1, "Cannot find input file containing GWAS VCF paths: ${params.input}" }
        .splitCsv()
        .flatten()
        .map { vcf -> file(vcf) }
        .take(params.take_n_studies) //default is -1 i.e. take all files (but param is useful for testing with fewer files)
} else if (params.gwas_source in ['gwas_vcf'] && params.input_type == 'single') {
    ch_sumstats = Channel.fromPath(params.input)
} else if (params.gwas_source in ['gwas_table'] && params.input_type == 'list') {
    ch_gwas_tables = Channel.fromPath(params.input)
        .ifEmpty { exit 1, "Cannot find input file containing GWAS Table paths: ${params.input}" }
        .splitCsv()
        .flatten()
        .map { table -> file(table) }
        .take(params.take_n_studies) //default is -1 i.e. take all files (but param is useful for testing with fewer files)
} else if (params.gwas_source in ['gwas_table'] && params.input_type == 'single') {
    ch_gwas_tables = Channel.fromPath(params.input)
} else if (params.gwas_source in ['gwas_table'] && params.input_type == 'regenie_folder_for_lb_connect') {
    gwas_results_dir = params.input.split(',').collect()
    ch_gwas_results_dir = Channel.fromPath(gwas_results_dir)
    //ch_gwas_results_dir.view()
    //ch_gwas_tables
}

if (params.input_type == 'regenie_folder_for_lb_connect'){
    // not to have same name collision for "results" directory sufix
    process stageResults {
        input:
        file(gwas_results_dir) from ch_gwas_results_dir
	    //echo true
	
        output:
        file("${uuid}_notransform.regenei") into ch_gwas_tables

        script:
        uuid = UUID.randomUUID().toString()
        """
        mv $gwas_results_dir/allancs/notransform/regenie/*.regenie ${uuid}_notransform.regenei
        """
    }
}

// Channel for omop vocabulary
ch_omop_vocabulary = Channel.value(file(params.omop_vocabulary))

// Channels for SNPlocs
ch_snplocs_grch38 =  Channel.value(file(params.snplocs_grch38))
ch_snplocs_grch37 =  Channel.value(file(params.snplocs_grch37))

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
        set val("${study_id}"), file("*${study_id}*.h.tsv.gz"), file("${study_id}_metadata_ebi.tsv") into ch_sumstats
        file("reported_traits.tsv") into ch_reported_traits

        script:
        """
        bash $ebi_script $study_id
        """
    }

    process munge_ebi {
        label 'munge'
        label "high_memory"
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ebi/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats), file(metadata) from ch_sumstats
        each file('munge_ebi.R') from Channel.fromPath("${projectDir}/bin/munge_ebi.R")
        each file('munge_funcs.R') from Channel.fromPath("${projectDir}/bin/munge_funcs.R")
        each file('harmonise_metadata_funcs.R') from Channel.fromPath("${projectDir}/bin/harmonise_metadata_funcs.R")
        each file('field_descriptions.tsv') from Channel.fromPath("${projectDir}/assets/field_descriptions.tsv")
        each file(snplocs_grch37) from ch_snplocs_grch37
        each file(snplocs_grch38) from ch_snplocs_grch38

        output:
        set val(study_id), file("${study_id}_harmonised_sumstats.vcf") into ch_harmonised_sumstats

        script:
        """
        R CMD INSTALL $snplocs_grch37
        R CMD INSTALL $snplocs_grch38

        Rscript munge_ebi.R \
            --ss_table $sumstats \
            --meta_table $metadata \
            --method ${params.munge_method} \
            --study_id $study_id \
            --output "${study_id}_harmonised_sumstats.vcf" \
            --col_miss ${params.miss_percent_allow} \
            --dbsnp ${params.dbsnp} \
            --field_descriptions field_descriptions.tsv \
            --ncpus ${task.cpus}
        
        gunzip -S .bgz ${study_id}_harmonised_sumstats.vcf.bgz
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
        set val("${study_id}"), file("${study_id}.vcf.gz"), file("${study_id}_metadata_ieu.tsv") into ch_sumstats
        file("reported_traits.tsv") into ch_reported_traits
 
        script:
        """
        bash $ieu_script $study_id
        """
    }

    process munge_ieu {
        label 'munge'
        label "high_memory"
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/ieu/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats), file(metadata) from ch_sumstats
        each file('munge_ieu.R') from Channel.fromPath("${projectDir}/bin/munge_ieu.R")
        each file('munge_funcs.R') from Channel.fromPath("${projectDir}/bin/munge_funcs.R")
        each file('harmonise_metadata_funcs.R') from Channel.fromPath("${projectDir}/bin/harmonise_metadata_funcs.R")
        each file('field_descriptions.tsv') from Channel.fromPath("${projectDir}/assets/field_descriptions.tsv")
        each file(snplocs_grch37) from ch_snplocs_grch37
        each file(snplocs_grch38) from ch_snplocs_grch38

        output:
        set val(study_id), file("${study_id}_harmonised_sumstats.vcf") into ch_harmonised_sumstats

        script:
        """
        R CMD INSTALL $snplocs_grch37
        R CMD INSTALL $snplocs_grch38

        Rscript munge_ieu.R \
            --gwas_vcf $sumstats \
            --meta_table $metadata \
            --method ${params.munge_method} \
            --study_id $study_id \
            --output "${study_id}_harmonised_sumstats.vcf" \
            --col_miss ${params.miss_percent_allow} \
            --dbsnp ${params.dbsnp} \
            --field_descriptions field_descriptions.tsv \
            --ncpus ${task.cpus}
        
        gunzip -S .bgz ${study_id}_harmonised_sumstats.vcf.bgz
        """
    }
}

if (params.gwas_source == 'gwas_vcf') {
    if (params.force_munge) {
        process munge_gwas_vcf {
            label 'munge'
            label "high_memory"
            if (params.keep_intermediate_files) {
                publishDir "${params.outdir}/gwas_vcf_input/", mode: 'copy'
            }

            input:
            file(sumstats) from ch_sumstats
            each file('munge_gwasvcf.R') from Channel.fromPath("${projectDir}/bin/munge_gwasvcf.R")
            each file('munge_funcs.R') from Channel.fromPath("${projectDir}/bin/munge_funcs.R")
            each file('field_descriptions.tsv') from Channel.fromPath("${projectDir}/assets/field_descriptions.tsv")
            each file(snplocs_grch37) from ch_snplocs_grch37
            each file(snplocs_grch38) from ch_snplocs_grch38

            output:
            set file('study_id'), file("*_harmonised_sumstats.vcf") into ch_expanded_gwas_vcf
            file("reported_traits.tsv") into ch_reported_traits

            script:
            """
            R CMD INSTALL $snplocs_grch37
            R CMD INSTALL $snplocs_grch38

            if [[ $sumstats == *gz ]]; then 
                study_id=\$(zcat $sumstats | awk '\$1=="#CHROM"{print \$10; quit}')
            else
                study_id=\$(awk '\$1=="#CHROM"{print \$10; quit}' $sumstats)
            fi
            echo -n \$study_id > study_id

            Rscript munge_gwasvcf.R \
                --gwas_vcf $sumstats \
                --method ${params.munge_method} \
                --output "\${study_id}_harmonised_sumstats.vcf" \
                --dbsnp ${params.dbsnp} \
                --field_descriptions field_descriptions.tsv \
                --ncpus ${task.cpus}
            
            gunzip -S .bgz \${study_id}_harmonised_sumstats.vcf.bgz
            """
        }
    }
    
    if (!params.force_munge) {
        process expand_gwas_vcf {
            label 'bcftools'

            input:
            file(sumstats) from ch_sumstats

            output:
            set file("study_id"), file("*_harmonised_sumstats.vcf") into ch_expanded_gwas_vcf
            file("reported_traits.tsv") into ch_reported_traits

            script:
            """
            if [[ $sumstats == *gz ]]; then 
                study_id=\$(zcat $sumstats | awk '\$1=="#CHROM"{print \$10; quit}')
            else
                study_id=\$(awk '\$1=="#CHROM"{print \$10; quit}' $sumstats)
            fi
            echo -n \$study_id > study_id

            bcftools norm --multiallelic -any $sumstats -O v -o \${study_id}_harmonised_sumstats.vcf

            sed -n -E '/##SAMPLE=/ s/.*TraitReported=([^>=]+)[>=].*/\\1/p' \${study_id}_harmonised_sumstats.vcf \
            | sed -E 's/,[^,]*//' \
            > reported_traits.tsv
            """
        }
    }
    ch_harmonised_sumstats = ch_expanded_gwas_vcf.map{ [it[0].text] + [it[1]] }
}

if (params.gwas_source == 'gwas_table') {
    process munge_gwas_table {
        label 'munge'
        label "high_memory"
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/gwas_vcf_input/", mode: 'copy'
        }

        input:
        file(gwas_table) from ch_gwas_tables
        each file('munge_gwastable.R') from Channel.fromPath("${projectDir}/bin/munge_gwastable.R")
        each file('transform_for_munge.sh') from Channel.fromPath("${projectDir}/bin/transform_for_munge.sh")
        each file('harmonise_metadata_funcs.R') from Channel.fromPath("${projectDir}/bin/harmonise_metadata_funcs.R")
        each file('munge_funcs.R') from Channel.fromPath("${projectDir}/bin/munge_funcs.R")
        each file('field_descriptions.tsv') from Channel.fromPath("${projectDir}/assets/field_descriptions.tsv")
        each file(snplocs_grch37) from ch_snplocs_grch37
        each file(snplocs_grch38) from ch_snplocs_grch38

        output:
        set file(study_id), file("*_harmonised_sumstats_${uuid}.vcf") into ch_harmonised_table_vcf

        script:
        // not to have same name collision for "results" directory sufix
        uuid = UUID.randomUUID().toString()
        """
        # Generate study_id
        population=\$(grep -m1 -h "#Population" $gwas_table | awk -F "\\t" '{print \$2}')
        studytag=\$(grep -m1 -h "#StudyTag" $gwas_table | awk -F "\\t" '{print \$2}')
        if [ -z \$studytag ]; then
            studytag="NA"
        fi
        method=\$(grep -m1 -h "#Method" $gwas_table | awk -F "\\t" '{print \$2}')
        study_identifier="\${population}-\${studytag}-\${method}"

        echo -n \$study_identifier > study_id

        # Get genome_build from table if available
        genome_build=\$(grep -m1 -h "#GenomeBuild" $gwas_table | awk -F "\\t" '{print \$2}')
        if [[ -z \$genome_build ]]; then
            genome_build=\$(
                grep -ho 'genome_build=GRCh[0-9][0-9]' $gwas_table | awk -F '=' '{print \$2}'
            )
        fi
        metagwas="FALSE"
        if [[ \$method == "METAL" ]]; then
            metagwas="TRUE"
        fi

        # Get total samples from table
        totalcases=\$(grep -m1 -h "#TotalCases" $gwas_table | awk -F "\\t" '{print \$2}')
        totalcontrols=\$(grep -m1 -h "#TotalControls" $gwas_table | awk -F "\\t" '{print \$2}')
        totalsamples=\$((totalcases + totalcontrols))

        # Get available metadata for munge
        sed -n -e '/^[^#]/ q' -e 's/^##// p' $gwas_table |
        awk 'BEGIN {FS=OFS="\\t"} { for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") \$i; } \
        END{ for (i=1; i<=NF; i++) print RtoC[i] }' | sed -e 's/^\\t//g' > meta_table.tsv

        # Run transformation script to make table mungeable based on gwas method
        bash transform_for_munge.sh $gwas_table > sumstat.tmp

        # Install only necessary ref genomes for munge script if possible
        if [ -z \$genome_build ]; then
            R CMD INSTALL $snplocs_grch38
            R CMD INSTALL $snplocs_grch37
        elif [[ \$genome_build == "GRCh38" ]]; then
            R CMD INSTALL $snplocs_grch38
        elif [[ \$genome_build == "GRCh37" ]]; then
            R CMD INSTALL $snplocs_grch37
        fi

        # Run the gwas_table munging script
        Rscript munge_gwastable.R \
            --gwas_table=sumstat.tmp \
            --meta_table=meta_table.tsv \
            --outfile=\${study_identifier}_harmonised_sumstats_${uuid}.vcf \
            --genome_build="\$genome_build" \
            --total_samples=\$totalsamples \
            --dbsnp=${params.dbsnp} \
            --study_id=\${study_identifier} \
            --metagwas=\$metagwas \
            --ncpus=${task.cpus} 2> \${study_identifier}.munge.log

        gunzip -S .bgz \${study_identifier}_harmonised_sumstats_${uuid}.vcf.bgz
        """
    }
    ch_harmonised_table_vcf.map{ [it[0].text] + [it[1]] }
    .into {ch_harmonised_sumstats; ch_harmonised_sumstats_2}
    //ch_harmonised_sumstats_2.view()
}

if (params.standardise && params.coef_conversion) {
    process standardisation_beta2or {
        label 'python'
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/conversion/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats) from ch_harmonised_sumstats
        each file(conv_script) from ch_conv_script

        output:
        set val("${study_id}"), file("${study_id}_harmonised_conv_sumstats_${uuid}.vcf") into ch_conv_sumstats
        
        script:
        // not to have same name collision for "results" directory sufix
        uuid = UUID.randomUUID().toString()
        """
        # Detect the FRQ column necessary to `--standardise`
        FRQ_COLUMN_EXISTS=\$(awk '/^[^#]/ {if(\$9 ~ /AF/) a=1; else a=0; exit}END{print a}' $sumstats)

        # Decide and convert
        if [ \$FRQ_COLUMN_EXISTS == "1" ]; then
            python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats_${uuid}.vcf --standardise --beta2or
            echo "[MSG] BETA and SE standardisation performed."
            echo "[MSG] BETA to OR conversion performed."
        else
            python $conv_script -v $sumstats -o ${study_id}_harmonised_conv_sumstats_${uuid}.vcf --beta2or
            echo "[MSG] BETA to OR conversion performed."
            echo "[WARNING] BETA and SE standardisation not performed as FRQ value not present."
        fi
        """
    }
} else if (params.standardise) {
    process standardisation {
        label 'python'
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/conversion/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats) from ch_harmonised_sumstats
        each file(conv_script) from ch_conv_script

        output:
        set val("${study_id}"), file("${study_id}_harmonised_conv_sumstats.vcf") into ch_conv_sumstats
        
        script:
        """
        # Detect the FRQ column necessary to `--standardise`
        FRQ_COLUMN_EXISTS=\$(awk '/^[^#]/ {if(\$9 ~ /AF/) a=1; else a=0; exit}END{print a}' $sumstats)

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
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/conversion/", mode: 'copy'
        }

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


if (params.map_traits) {

    // combine all traits into a single csv for improved efficiency
    //(loading omop vocabulary is the majority of the processing here)

    ch_collected_traits = ch_reported_traits
        .map{it.text}
        .collect()
        .map{it.join('\n')}

    process map_traits {
        label "mid_memory"
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/metadata/", mode: 'copy'
        }

        input:
        file('reported_traits.list') from ch_collected_traits
        file(omop_vocabulary) from ch_omop_vocabulary
        each file("map_traits.R") from Channel.fromPath("${projectDir}/bin/map_traits.R")

        output:
        file("mapped.csv") into ch_mapped_traits

        script:
        """
        unzip ${omop_vocabulary}

        echo "traitName" > traits.csv
        cat reported_traits.list >> traits.csv

        Rscript map_traits.R --traits_csv=traits.csv --vocabulary_folder=omop-vocab-files --snomed_grouping_level=4
        """
    }

    process insert_traits {
        if (params.keep_intermediate_files) {
            publishDir "${params.outdir}/metadata/", mode: 'copy'
        }

        input:
        set val(study_id), file(sumstats_vcf) from ch_conv_sumstats
        each file("mapped.csv") from ch_mapped_traits
        each file("add_trait_metadata.R") from Channel.fromPath("${projectDir}/bin/add_trait_metadata.R")
        each file('field_descriptions.tsv') from Channel.fromPath("${projectDir}/assets/field_descriptions.tsv")

        output:
        set val(study_id), file("${study_id}_harmonised_conv_omop_sumstats.vcf") into ch_omop_gwas_vcf

        script:
        """
        sed '/^#CHROM/ q' $sumstats_vcf > header.txt

        Rscript add_trait_metadata.R \
            --gwas_vcf_head header.txt \
            --trait_table mapped.csv \
            --field_descriptions field_descriptions.tsv

        sed '/^#/ d' $sumstats_vcf  \
        | cat updated_header.txt - \
        > ${study_id}_harmonised_conv_omop_sumstats.vcf
        """
    }

} else {
    ch_omop_gwas_vcf = ch_conv_sumstats
}


ch_omop_gwas_vcf.into{
    ch_full_sumstats;
    ch_full_sumstats_2;
}


process make_qc_plots {
    label "mid_memory"
    publishDir "${params.outdir}/QC_plots", mode: 'copy'
    
    input:
    set val(study_id), file(fully_harmonised_sumstats_vcf) from ch_full_sumstats
    each file('gwas_qc_plots.py') from ch_qc_plots_script

    output:
    tuple val(study_id), file("${study_id}/*") into ch_qc_plots

    script:
    """
    # Check whether the VCF contains FRQ column
    if `tail -1 fully_harmonised_sumstats.vcf | grep -q "AF"`; then
        python gwas_qc_plots.py \
            --sumstat $fully_harmonised_sumstats_vcf \
            --c-maf AF \
            --c-beta ES \
            --c-se SE \
            --c-pval LP \
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
        label 'bcftools'
        publishDir "${params.outdir}/harmonised/filtered/", mode: 'copy'
        echo true
        
        input:
        set val(study_id), file(sumstats) from ch_full_sumstats_2

        output:
        tuple val(study_id), file("${study_id}_harmonised_filtered_sumstats_${uuid}.vcf") into ch_filtered_sumstats

        script:
        // not to have same name collision for "results" directory sufix
        uuid = UUID.randomUUID().toString()
        beta_smaller_filter = params.filter_beta_smaller_than || params.filter_beta_smaller_than == 0 ? "-e FORMAT/ES[*]<$params.filter_beta_smaller_than" : ""
        beta_greater_filter = params.filter_beta_greater_than || params.filter_beta_greater_than == 0 ? "-e FORMAT/ES[*]>$params.filter_beta_greater_than" : ""
        p_filter = params.filter_p_greater_than || params.filter_p_greater_than == 0 ? "-e FORMAT/LP[*]>$params.filter_p_greater_than" : ""
        frq_smaller_filter = params.filter_freq_smaller_than || params.filter_freq_smaller_than == 0 ? "-e FORMAT/AF[*]<$params.filter_freq_smaller_than" : ""
        frq_greater_filter = params.filter_freq_greater_than || params.filter_freq_greater_than == 0 ? "-e FORMAT/AF[*]>$params.filter_freq_greater_than" : ""
        info_filter = params.filter_info ? "-e FORMAT/SI[*]<$params.filter_info" : ""
        missing_info_filter = params.filter_missing_info || params.filter_missing_info == 0 ? "-e FORMAT/SI[*]='.'" : ""
        """
        cp $sumstats temp_sumstats.vcf
        # Filter expressions to pass into bcftools filter
        # Check for the AF and SI columns
        if tail -1 temp_sumstats.vcf | grep -q "AF"; then
            expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter" "$frq_smaller_filter" "$frq_greater_filter")
        else
            echo "[WARNING] AF column not found, all allele frequency filters will be ignored"
            expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter")
        fi
        if tail -1 temp_sumstats.vcf | grep -q "SI"; then
            expression_list+=("$info_filter" "$missing_info_filter")
        else
            echo "[WARNING] SI column not found, all imputation accuracy filters will be ignored"
        fi

        # Apply filters sequentially, otherwise bcftools would apply OR logic
        for e in "\${expression_list[@]}"; do
            if [ "\$e" != "" ]; then
                bcftools filter "\$e" temp_sumstats.vcf > temp_sumstats_new.vcf
                mv temp_sumstats_new.vcf temp_sumstats.vcf
            fi
        done
        mv temp_sumstats.vcf "${study_id}_harmonised_filtered_sumstats_${uuid}.vcf"
        """
    }
} else {
    ch_filtered_sumstats = ch_full_sumstats_2
}

ch_filtered_sumstats.into{
    ch_filtered_sumstats_for_hail;
    ch_filtered_sumstats_for_collapse;
    ch_filtered_sumstats_for_metagwas_metal
}

process convert_vcf_metal {
    label 'vcf'
    publishDir "${params.outdir}/converted_vcf", mode: "copy"

    input:
    set val(study_name), file(study) from ch_filtered_sumstats_for_metagwas_metal

    output:
    file("${study.baseName}_${uuid}.converted.csv") into all_input_studies_ch

    script:
    // not to have same name collision for "results" directory sufix
    uuid = UUID.randomUUID().toString()
    """
    # Process metadata from VCF file
    grep -m1 "^##SAMPLE" $study | \
    sed -e "s/^##SAMPLE=<//" -e "s/>\$//" -e "s/\\"//g" -e "s/=/\\t/g" | \
    sed -E "s/([^ ]),([A-Z])/\\1\\n\\2/g" | \
    while read line; do echo "##\$line"; done > "${study.baseName}_${uuid}.converted.csv"

    echo "${params.chr_col},${params.pos_col},${params.rsid_col},${params.a1_col},${params.a2_col},${params.metal_freq_col},${params.beta_col},${params.se_col},${params.pval_col}" >> "${study.baseName}_${uuid}.converted.csv"
    bcftools norm --multiallelic -any -Ou $study | \
    bcftools query -f'%CHROM,%POS,[%ID],%REF,%ALT,[%AF],[%ES],[%SE],[%LP]\n' | \
        awk 'BEGIN{FS=OFS=","} {
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,sprintf("%.10f",10^-\$9);
        }' >> "${study.baseName}_${uuid}.converted.csv"
    """
}

/*-----------------------------
  Setting up extra METAL flags
-------------------------------*/

// Initialise variable to store optional parameters
extra_flags = ""

// 1 - METAL options for describing input files

if ( params.flip ) { extra_flags += "FLIP \n" }

// 2 - METAL options for filtering input files

if ( params.addfilter ) { extra_flags += "ADDFILTER ${params.addfilter}\n" }
if ( params.removefilters ) { extra_flags += "REMOVEFILTERS  \n" }

// 3 - METAL options for sample size weighted meta-analysis

if ( params.weightlabel ) { extra_flags += "WEIGHTLABEL ${params.weightlabel}\n" }
if ( params.defaultweight ) { extra_flags += "DEFAULTWEIGHT ${params.defaultweight}\n" }
if ( params.minweight ) { extra_flags += "MINWEIGHT ${params.minweight}\n" }

// 4 - METAL options for inverse variance weighted meta-analysis

if ( params.stderrlabel ) { extra_flags += "STDERRLABEL ${params.stderrlabel}\n" }
if ( params.scheme ) { extra_flags += "SCHEME ${params.scheme}\n" }

// 5 - METAL options to enable tracking of allele frequencies

if ( params.averagefreq ) { extra_flags += "AVERAGEFREQ ${params.averagefreq}\n" }
if ( params.minmaxfreq ) { extra_flags += "MINMAXFREQ ${params.minmaxfreq}\n" }
if ( params.freqlabel ) { extra_flags += "FREQLABEL ${params.freqlabel}\n" }

// 6 - METAL options to enable tracking of user defined variables

if ( params.customvariable ) { extra_flags += "CUSTOMVARIABLE ${params.customvariable}\n" }
if ( params.label ) { extra_flags += "LABEL ${params.label}\n" }

// 7 - METAL options to enable explicit strand information

if ( params.usestrand ) { extra_flags += "USESTRAND ${params.usestrand}\n" }
if ( params.strandlabel ) { extra_flags += "STRANDLABEL ${params.strandlabel}\n"  }

// 8 - METAL options for automatic genomic control correction of input statistics

if ( params.genomiccontrol ) { extra_flags += "GENOMICCONTROL ${params.genomiccontrol}\n" }

// 9 - METAL options for general analysis control  

if ( params.outfile ) { extra_flags += "OUTFILE ${params.outfile}\n"}
if ( params.maxwarnings ) { extra_flags += "MAXWARNINGS ${params.maxwarnings}\n" }
if ( params.verbose ) { extra_flags += "VERBOSE ${params.verbose}\n"}
if ( params.logpvalue ) { extra_flags += "LOGPVALUE ${params.logpvalue}\n" }

// 10 - METAL options for general run control not available (pipeline is not currently developed to handle this)

/*------------------------------
  Running METAL (meta-analysis)
--------------------------------*/

// NB: this process must be "padded to the wall" to allow for extra flags to be properly inserted

process run_metal {
publishDir "${params.outdir}/METAL", mode: "copy"

input:
file(study) from all_input_studies_ch.collect()
file('fill_chr_pos_metal_output.R') from Channel.fromPath("${projectDir}/bin/fill_chr_pos_metal_output.R", followLinks: false)

output:
file("METAANALYSIS*") into results_ch
file("METAANALYSIS*.with_pstns.tsv") into ch_metal_output_for_plots

shell:
'''
# 1 - Dynamically obtain files to process
touch process_commands.txt
touch metadata_tmp.txt
for csv in $(ls *.csv)
do 
echo "PROCESS $csv" >> process_commands.txt
if ( head -1 $csv | grep -q "^##" ); then
    sed -n '/^##/!q;p' $csv >> metadata_tmp.txt
    sed -i '/^##/d' $csv
fi
done
echo "##Method\tMETAL" > metadata.txt
cat metadata_tmp.txt | sort -u | \
awk '{a[$1]=a[$1] FS $2} END{for(i in a) print i a[i]}' | \
sed 's/ /'$'\t''/' | \
sed '/^##Method/d' >> metadata.txt
process_commands=$(cat process_commands.txt)
# 2 - Make METAL script 
cat > metal_command.txt <<EOF
MARKER !{params.rsid_col}
ALLELE !{params.a1_col} !{params.a2_col}
EFFECT !{params.beta_col}
PVALUE !{params.pval_col} 
SEPARATOR COMMA
AVERAGEFREQ ON
MINMAXFREQ ON
FREQLABEL !{params.metal_freq_col}
!{extra_flags}
$process_commands
ANALYZE !{params.heterogeneity ? "HETEROGENEITY" : ""}
QUIT
EOF
# 3 - Run METAL
echo "Running METAL"
metal metal_command.txt
Rscript fill_chr_pos_metal_output.R \
--metal_output METAANALYSIS1.TBL \
--study_file_list $(echo !{study} | tr ' ' ',') \
--study_file_cols '!{params.rsid_col},!{params.chr_col},!{params.pos_col}'
# 4 - Add metadata
cat metadata.txt METAANALYSIS1.TBL \
> metadata_new.txt && mv metadata_new.txt METAANALYSIS1.TBL
cat metadata.txt METAANALYSIS1.TBL.with_pstns.tsv \
> metadata_new.txt && mv metadata_new.txt METAANALYSIS1.TBL.with_pstns.tsv
'''
}

process metal_report {
    label 'report'
    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    file("results.tsv") from ch_metal_output_for_plots
    file(report_files) from Channel.fromPath("${projectDir}/bin/report/*").collect()
    file('genes.gff3') from Channel.fromPath("${params.gff3_for_locuszoom}")

    output:
    file("multiqc_report.html")

    script:
    """
    # Check for empty results.tsv where filtering was rather strict
    if ( tail -1 results.tsv | grep -q "MarkerName" ); then
      cat > multiqc_report.html <<- EOM
      <!DOCTYPE html> 
      <html lang=???en???> 
	     <head></head> 
	      <body> 
		      <p><strong><em>There is no data available to generate report</em></strong></p> 
	      </body> 
      </html>
    EOM
      exit 0
    fi
    Rscript -e "wd=getwd();rmarkdown::render('metal_report.Rmd', \
      params=list(
        input_table='results.tsv',
        analyse_het=${params.heterogeneity ? "T" : "F"},
        analyse_beta=${params.scheme == "STDERR" ? "T" : "F"},
        p_val_thresh=${params.p_val_threshold},
        gene_gff='genes.gff3',
        ensembl_pop='${params.ld_ensembl_pop}',
        top_n=${params.plot_top_n}), \
      clean=F, intermediates_dir=wd, output_dir=wd, knit_root_dir=wd)"
    
    mv metal_report.html multiqc_report.html
    """
  }

process collapse_vcf {
    label 'python'
    publishDir "${params.outdir}/harmonised/", mode: 'copy'

    input:
    set val(study_id), file(gwas_vcf) from ch_filtered_sumstats_for_collapse
    each file('collapse_vcf.py') from Channel.fromPath("${projectDir}/bin/collapse_vcf.py")

    output:
    file("${study_id}.harmonised.gwas.vcf")
    
    script:
    """
    python3 collapse_vcf.py $gwas_vcf -o ${study_id}.harmonised.gwas.vcf
    """
}

if (params.convert_to_hail) {
    process vcf2hail {
        publishDir "${params.outdir}/harmonised/hail/", mode: 'copy'

        input:
        set val(study_id), file(sumstats) from ch_filtered_sumstats_for_hail
        each file(vcf2hail) from ch_vcf2hail_script

        output:
        tuple val(study_id), file("${study_id}_hail_matrix.mt/*") into ch_hail_sumstats

        script:
        """
        GENOME=\$( grep -e '##GenomeBuild=' $sumstats | cut -d "=" -f 2 )
        python $vcf2hail -v $sumstats -g \$GENOME -o "${study_id}_hail_matrix.mt"
        """
    }
}

// ch_report_dir = Channel.value(file("${projectDir}/bin/report"))
// process create_report {
//   publishDir "${params.outdir}/MultiQC", mode: 'copy'
  
//   input:
//   set val(study_id), file("${study_id}/*") from ch_qc_plots
//   file(report_dir) from ch_report_dir 

//   output:
//   file("${study_id}_multiqc_report.html") into ch_report_outputs_all
//   file("multiqc_report.html") into ch_report_outputs1

//   script:
//   """
//   cp -r ${report_dir}/* .
//   cp ${study_id}/* .
//   # Skip when no plots were generated due to lack of FRQ column
//   if [ -f "${study_id}_no_plots.txt" ]; then
//     touch multiqc_report.html
//     touch ${study_id}_multiqc_report.html
//   else
//     for i in `ls *.png`; do name=`basename \$i .png | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
//     for i in `ls *.tsv`; do name=`basename \$i .tsv | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
//     cat file_list.txt | tr "\n" "," | sed 's/,\$//g' > file_list1.txt
//     # Generates the report
//     Rscript -e "rmarkdown::render('report.Rmd', params = list(`cat file_list1.txt`))"
//     cp report.html multiqc_report.html
//     mv report.html ${study_id}_multiqc_report.html
//   fi
//   """
// }
