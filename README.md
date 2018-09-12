# Deconvolution analysis for neuronal-specific miRNA expression

This reprosition is for the deconvolution analysis for Chiara's project

### To-do
  - [x] Understand the difference between the WT and henn-1 mutant, because there is a good correlation between the log2FC in both genotypes. The starting point
    will be to compare the henn1-mutant without promoter, henn1-mutant with ubiquitous promoter and N2 without promoter (with ubiquitous promoter data ?? not sure exits)
    - the conclusion is not very clear; because the highly enriched guys are not very highly expressed in the untreated samples, suggesting they are not simply due to the unefficient removal for unmethylated ones in the oxidation step
    - But it is clear that in the WT background, since the oxidation is shared by more species, e.g. siRNAs, piRNAs; consequently, the oxidation efficiency is lower for miRNAs. So the log2FC in WT is always lower than the log2FC in the henn1.mutant. But there are are some particular guys, which can be highly methylated (probably) already in the background. The reason could be the Chiara's assumption.
    
  - [ ] Double check the batch correction works
    - The batch correction works more or less; but the the batch correction for rab-3 is much improved but still not perfect, in spit of limma or combat batch correction.   
    - The reason requires further understanding of batch correction methods.
    - double check the data before and after batch correction for example lsy-6, mir-791 and mir-790 
    - By comparing the log and linear scaling with examples lsy-6, mir-791, mir-790, the observation is that the total amount of specific microRNA defined by pan-neurons is much large than other samples, consequently, the linear model must assign the rest to some neurons, which makes false positives; in contrast, in log-scale the solution will be much sparse and resonable.
    
  
  - [ ] Integrate the sample qualities into the linear model
  
  - [ ] Optimize the fitting and parameter to have final results

### Deconvolution analysis
`deconvolution_analysis.R`

- the binary fraction (or proportion) matrix works better than the actual cell number matrix in which the big groups are 
  hardly penelized (not clear about the reason)
- We stil have the trouble for pan-neuron sample, one of the replicates have too high gene expression in lsy-6 at least, and also higher in general.
  For the moment, we try to remove this sample to test 

### time series analsyis
`miRNA_time_series_spikeIns.R`
- versions:   
  - "_miRNAs_timeSeries_spikeIn_R5922_R6016_20180620"  
    -correct the spike-in normalization    
  - "_miRNAs_timeSeries_spikeIn_R5922_R6016_20180719"  
    - correct the concentration for spike ins and the units in the table
    - make heatmaps for treated samples across stages, trying to address the main question  
      which miRNAs are for newly-born neurons and which ones are deleted.
    - reorganize the code : 
      - data and design process; 
      - spike-in normalization
      - QC for cpm and spike-in normalization
      - enrichment analysis
      - some real analysis
  - "_miRNAs_timeSeries_spikeIn_R5922_R6016_20180802"        
    consider the cell number c(230, 293, 301, 302) for late L1, L2, L3 and L4
    and meanwhile scale the expression into the range (0, 1)
    

## update
  - 20180822: 
    new rab-3 in WT background with spike-in cames; to integrate the new data into the existing analysis, I just need to run the
    `miRNA_enrichmentAnalysis_for_neuronClass.R` which compares the WT and henn1-mutant and also enrichment analysis for both background;
    `prepare_ExpressionMatrix.R` to compare the spike-in normalization and also the piRNA normalization; as long as the new rab-3 sample was added,
    the matrix can be used directly by `deconvolution_analysis.R`
    
    This will be the same for the last sample to come.
    
  - 20180831: 
    - the last sample of unc-3 is sequenced; so if all samples look good, the data should be complete. 
    - IMPORTANT NOTE: after discussion with Thomas, there is a way to extract read counts for each piRNAs; the count information is in the tables like
      `72762_ATCACG_CCRBHANXX_6_20180822B_20180822.trimmed.virus.bam.tailor.bam.seqCnt.txt.gz`; An easy solution suggested by Thomas is to modify the 
      one of R codes, by specifiying the piRNAs together with miRNAs; but anyway this needs to be double checked. 
    
    
    