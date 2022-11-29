// See https://www.opendevstack.org/ods-documentation/ for usage and customization.

@Library('ods-jenkins-shared-library@4.x') _

node {
  dockerRegistry = env.DOCKER_REGISTRY
}

/*
Parameters related to CloudOS run
*/

// [REQUIRED] Pipeline specific variables - must be updated
CLOUDOS_WORKFLOW_NAME = "iudlgnt-gwas-sumstats-utils"

// Test profiles
CLOUDOS_PROFILE_1="network_cloudos,test_ebi"
CLOUDOS_PROFILE_2="network_cloudos,test_ebi_dbsnp155"
CLOUDOS_PROFILE_3="network_cloudos,test_ieu"
CLOUDOS_PROFILE_4="network_cloudos,test_ieu_dbsnp155"
CLOUDOS_PROFILE_5="network_cloudos,test_gwas_vcf_basic"
CLOUDOS_PROFILE_6="network_cloudos,test_gwas_tables_all"

CLOUDOS_PROFILE_LIST="$CLOUDOS_PROFILE_1 $CLOUDOS_PROFILE_2 $CLOUDOS_PROFILE_3 $CLOUDOS_PROFILE_4 $CLOUDOS_PROFILE_5 $CLOUDOS_PROFILE_6"

// [REQUIRED] Environment Specific Variables - must be updated
CLOUDOS_URL = "https://prod.cloudos.aws.boehringer.com"
CLOUDOS_WORKSPACE_ID = "620fc193a487941c70e8a2b0"
//CLOUD_OS_TOKEN = "xxxx" // don't define! This is stored as a secret env in OpenShift

// [OPTIONAL] Variables with sensible defaults - could be optionally parameters
CLOUDOS_PROJECT_NAME = "JenkinsCI" // the project must already exist in the workspace
CLOUDOS_INSTANCE_TYPE="r4.2xlarge"

// [CLOUDOS CMD]
CLOUDOS_RUN_CMD = "cloudos job run"
CLOUDOS_RUN_CMD += " -k \${CLOUD_OS_TOKEN} -c ${CLOUDOS_URL} --workspace-id ${CLOUDOS_WORKSPACE_ID} --project-name ${CLOUDOS_PROJECT_NAME}"
CLOUDOS_RUN_CMD += " --workflow-name ${CLOUDOS_WORKFLOW_NAME}"
CLOUDOS_RUN_CMD += " --instance-type ${CLOUDOS_INSTANCE_TYPE} --resumable --spot --wait-completion"


odsComponentPipeline(
  podContainers: [// <editor-fold desc="podContainers>
                  containerTemplate(
                    name: 'jnlp', // do not change, see https://github.com/jenkinsci/kubernetes-plugin#constraints
                    image: "${dockerRegistry}/ods/jenkins-agent-base:4.x",
                    workingDir: '/tmp',
                    resourceRequestCpu: '2',
                    resourceLimitCpu: '10',
                    resourceRequestMemory: '2Gi',
                    resourceLimitMemory: '16Gi',
                    alwaysPullImage: true,
                    /* groovylint-disable-next-line GStringExpressionWithinString */
                    args: '${computer.jnlpmac} ${computer.name}'
                  ),
                  containerTemplate(
                    name: 'cloud-os',
                    image: "quay.io/lifebitaiorg/cloudos-py:v0.0.8bi",
                    workingDir: '/tmp',
                    alwaysPullImage: true,
                    ttyEnabled: true,
                    command: 'sleep',
                    args: 'inf'
                  )
  ], // </editor-fold desc="podContainers>
  branchToEnvironmentMapping: [
    'feature/': 'dev',
    'bugfix/': 'dev',

    'master': 'dev',
    'staging': 'test',
    'production': 'prod',
  ]
) { context ->


// run this in the "normal" ODS closures (logically):
  // build in ODS (OpenShift)
  // have a known image tag available
  // run test with that image
  // add another tag to successfully tested image

interrogateEnvironment(context)


// context.targetProject == iudlgnt-dev
// context.projectId == iudlgnt
// context.environment == dev

    // General concept to pre generate
    // def foo = map[string]Closure
    // foo.one = {echo 'bar '}
    // foo.two = {echo 'bar '}
    // foo.three = {echo 'bar '}
    // parallel( foo )

    runCloudOSJob(context, "${context.projectId}-cd-cloudos-token") // devSecret | test.. | prod
    // parallel(
    //   firstLabel: {
    //     runCloudOSJob(context, "${context.projectId}-cd-cloudos-token") // devSecret | test.. | prod
    //   },
    //   secondLabel: {
    //     runCloudOSJob(context, "${context.projectId}-cd-cloudos-token") // devSecret | test.. | prod
    //   },
    // )
}

private void interrogateEnvironment(def context){
  container('cloud-os'){
    sh 'pwd'
    sh 'ls -lah'
    sh 'python --version'
  }
}


private void runCloudOSJob(def context, String secretName){

  withEnv(["HTTP_PROXY=${env.HTTP_PROXY}", "HTTPS_PROXY=${env.HTTPS_PROXY}", "NO_PROXY=${env.NO_PROXY}",]) {
    withCredentials([
      string(credentialsId: secretName, variable: 'CLOUD_OS_TOKEN')
    ]){
      container('cloud-os'){
        output = sh(
          returnStdout: true,
          script:"""
            export HTTP_PROXY=http://appaccess-zscaler.boehringer.com:80 && \
            export HTTPS_PROXY=http://appaccess-zscaler.boehringer.com:80 && \
            export NO_PROXY=localhost,.boehringer.com,*.boehringer.com,10.,10.*,172.20.,172.20.*,0,1,2,3,4,5,6,7,8,9 && \
            GIT_COMMIT=`git log -n 1 --pretty=format:'%H'`
            for CLOUDOS_PROFILE in ${CLOUDOS_PROFILE_LIST}; do
                ${CLOUDOS_RUN_CMD} --nextflow-profile \${CLOUDOS_PROFILE} --job-name \${CLOUDOS_PROFILE} --git-commit \${GIT_COMMIT}
            done
            """).trim()

        }
        echo "XXXXXXXXXXXXXXX"
        echo "${output}"
        echo "XXXXXXXXXXXXXXX"
     }
  }
}
