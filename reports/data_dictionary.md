# Data Dictionary

This dictionary describes the columns present in the processed dataset.

| Column Name                     | Source File         | Data Type   | % Missing   | Description                                                                      |
|:--------------------------------|:--------------------|:------------|:------------|:---------------------------------------------------------------------------------|
| patient_id                      | patients_master.csv | str         | 0.00%       | Unique identifier for the patient (TCGA case submitter ID).                      |
| age_at_diagnosis                | patients_master.csv | float64     | 1.46%       | Age of the patient in days at the time of primary diagnosis.                     |
| gender                          | patients_master.csv | str         | 0.09%       | Sex of the patient at birth.                                                     |
| race                            | patients_master.csv | str         | 8.74%       | Race of the patient.                                                             |
| ethnicity                       | patients_master.csv | str         | 15.94%      | Ethnicity of the patient.                                                        |
| ajcc_pathologic_stage           | patients_master.csv | str         | 1.09%       | AJCC pathologic stage of the tumor.                                              |
| ajcc_pathologic_t               | patients_master.csv | str         | 0.09%       | AJCC pathologic T (tumor) stage.                                                 |
| ajcc_pathologic_n               | patients_master.csv | str         | 0.09%       | AJCC pathologic N (node) stage.                                                  |
| ajcc_pathologic_m               | patients_master.csv | str         | 0.09%       | AJCC pathologic M (metastasis) stage.                                            |
| primary_diagnosis               | patients_master.csv | str         | 0.09%       | Text description of the primary diagnosis.                                       |
| morphology                      | patients_master.csv | str         | 0.09%       | ICD-O-3 morphology code.                                                         |
| vital_status                    | patients_master.csv | str         | 0.09%       | Survival state of the patient (Alive/Dead).                                      |
| days_to_death                   | patients_master.csv | float64     | 86.25%      | Number of days from diagnosis to death.                                          |
| days_to_last_follow_up_clinical | patients_master.csv | float64     | 99.91%      | Number of days to last follow-up recorded in clinical data.                      |
| days_to_last_follow_up          | patients_master.csv | float64     | 0.18%       | Number of days to latest follow-up recorded in follow_up data.                   |
| patient_id                      | survival.csv        | str         | 0.00%       | Unique identifier for the patient (TCGA case submitter ID).                      |
| OS_time                         | survival.csv        | float64     | 0.00%       | Overall survival time in days (days to death or last follow-up).                 |
| OS_event                        | survival.csv        | float64     | 0.00%       | Overall survival event indicator (1 = Dead, 0 = Alive).                          |
| patient_id                      | samples.csv         | str         | 0.00%       | Unique identifier for the patient (TCGA case submitter ID).                      |
| sample_barcode                  | samples.csv         | str         | 0.00%       | Unique identifier for the biological sample.                                     |
| sample_type_code                | samples.csv         | int64       | 0.00%       | 2-digit sample type code extracted from the barcode.                             |
| sample_type                     | samples.csv         | str         | 0.00%       | Text description of the sample type (e.g., Primary Tumor, Blood Derived Normal). |