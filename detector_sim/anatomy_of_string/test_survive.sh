#!/bin/bash

for i in {0..23}
do
    hadd ../single_muon_jobs/job_${i}/copy.root ../single_muon_jobs/job_${i}/data_hadd.root
    ./analysis ${i}
done