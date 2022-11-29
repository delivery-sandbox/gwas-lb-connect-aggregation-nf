> Note: None of the template header sections should be removed and each should be addressed

## Overview

Does this

## Purpose

To achieve that

## Acceptance Criteria/Jira Ticket

[DEL-XXXX]()

## Changes

- Implements X
- Refactors Y
- Adds/Removes Z

## Client CloudOS Tests

- [ ] [-profile test_ebi_dbsnp155]()
- [ ] [-profile test_ebi]()
- [ ] [-profile test_gwas_tables_all]()
- [ ] [-profile test_gwas_vcf_list]()
- [ ] [-profile test_gwas_vcf_single_file]()
- [ ] [-profile test_ieu_dbsnp155]()
- [ ] [-profile test_ieu]()

## PR Checklist

#### Author to check:

- [ ] The PR title is informative and begins with appropriate `[Fix:|Feat:|Doc:]` based on PR type
- [ ] This PR is raised to `dev` branch for DSL1 or `DSL2-dev` for DSL2 translations
- [ ] The description contains the acceptance criteria or link to the associated Jira ticket
- [ ] The description contains links to CloudOS tests covering all test cases in client environment

#### Reviewer to check:

- [ ] A new test profile covering new functions for acceptance criteria is created if it does not already exist
- [ ] Nextflow configs have associated profiles named to match the config basename
- [ ] The [Jenkinsfile](https://github.com/lifebit-ai/gwas-sumstats-harmonisation-nf/blob/dev/Jenkinsfile) has been updated with current Nextflow profiles
- [ ] The [internal_lifebit_cloudos_ci.yml](https://github.com/lifebit-ai/gwas-sumstats-harmonisation-nf/blob/dev/.github/workflows/internal_lifebit_cloudos_ci.yml) has been updated to run all current Nextflow profiles
- [ ] Redundant tests have been removed
- [ ] The [README.md](https://github.com/lifebit-ai/gwas-sumstats-harmonisation-nf/blob/dev/docs/README.md) has been updated with parameter changes and current `results/` and conforms to the README.md guidelines
- [ ] The code is logical, readable, reasonably optimal and meets the acceptance criteria
