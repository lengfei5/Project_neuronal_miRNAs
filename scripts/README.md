# Deconvolution analysis for neuronal-specific miRNA expression

### `miRNA_time_series_spikeIns.R`
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

