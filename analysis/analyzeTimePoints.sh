#!/bin/bash

export DATA="../data"
export RESULT="../results"
export PLOT="../plots"

## Compare mutant to wild type

Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E14.5.csv --testCountFile=${DATA}/mut_E14.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E14.5_mut
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E15.5.csv --testCountFile=${DATA}/mut_E15.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E15.5_mut
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E16.5.csv --testCountFile=${DATA}/mut_E16.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E16.5_mut
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E17.5.csv --testCountFile=${DATA}/mut_E17.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E17.5_mut
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_NB.csv --testCountFile=${DATA}/mut_NB.csv --results=${RESULT} --plots=${PLOT} --experimentName=NB_mut

## Compare heterozygotes to wild type

Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E14.5.csv --testCountFile=${DATA}/het_E14.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E14.5_het
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E15.5.csv --testCountFile=${DATA}/het_E15.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E15.5_het
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E16.5.csv --testCountFile=${DATA}/het_E16.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E16.5_het
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_E17.5.csv --testCountFile=${DATA}/het_E17.5.csv --results=${RESULT} --plots=${PLOT} --experimentName=E17.5_het
Rscript differentialExpression.R --controlCountFile=${DATA}/wt_NB.csv --testCountFile=${DATA}/het_NB.csv --results=${RESULT} --plots=${PLOT} --experimentName=NB_het

## Combine all results, remove the header
cd ${RESULT}
cat *_results.csv  > all_results.csv
head -n 1 all_results.csv > result_header.txt
cat all_results.csv | grep -v logFC > a.txt; mv a.txt all_results.csv
cat result_header.txt all_results.csv > a.txt; mv a.txt all_results.csv 

rm result_header.txt ## Clean-up

## Plot fold-changes in gene expression versus time
cd "../analysis"
Rscript genesVsTime.R -i ../results/all_results.csv -p ../plots/ -l g17090,g23871