# Deconvolution analysis for neuronal-specific miRNA expression

This reprosition is for the deconvolution analysis for Chiara's project

### Deconvolution analysis
`deconvolution_analysis.R`

- the binary fraction (or proportion) matrix works better than the actual cell number matrix in which the big groups are 
  hardly penelized (not clear about the reason)


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
  - 20180822: new rab-3 in WT background with spike-in cames; to integrate the new data into the existing analysis, I just need to run the
    `miRNA_enrichmentAnalysis_for_neuronClass.R` which compares the WT and henn1-mutant and also enrichment analysis for both background;
    `prepare_ExpressionMatrix.R` to compare the spike-in normalization and also the piRNA normalization; as long as the new rab-3 sample was added,
    the matrix can be used directly by `deconvolution_analysis.R`
    
    This will be the same for the last sample to come.