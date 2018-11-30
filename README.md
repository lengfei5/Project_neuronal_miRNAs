# Deconvolution analysis for neuronal-specific miRNA expression

This reprosition is for the deconvolution analysis for Chiara's project

### To-do and Done
  - [x] Understand the difference between the WT and henn-1 mutant, because there is a good correlation between the log2FC in both genotypes. The starting point
    will be to compare the henn1-mutant without promoter, henn1-mutant with ubiquitous promoter and N2 without promoter (with ubiquitous promoter data ?? not sure exits)
    - the conclusion is not very clear; because the highly enriched guys are not very highly expressed in the untreated samples, suggesting they are not simply due to the unefficient removal for unmethylated ones in the oxidation step
    - But it is clear that in the WT background, since the oxidation is shared by more species, e.g. siRNAs, piRNAs; consequently, the oxidation efficiency is lower for miRNAs. So the log2FC in WT is always lower than the log2FC in the henn1.mutant. But there are are some particular guys, which can be highly methylated (probably) already in the background. The reason could be the Chiara's assumption.
    
  - [x] Select tuning parameter for elastic-net penalty   
    - [x] Test a global parameter lambda2 for ridge and gene-specific parameter lambda1 for lasso with gcdnet 
    - [x] Test gene-specifc parameter alpha, find the method to select the optimal alpha for each gene
      - test adaptive elastic-net and extract results from the `msaenet` function
      - use the cv to select alpha, specifically, prefix a vector of alpha, and then run cv.glment and find the minimum of MSE; finally, the alpha yielding the lowest MSE will be selected  
  - [x] After testing, one conclusion drawn is that the global alpha in glment works better than other choices, namely, global lambda2 in gcdnet and gene-specific alpha in glmnet until now;
      - “one-standard-error” rule turns out to be more efficient way to have sparse model, even compared with BIC, at least for mir7-91 example
      - so at the end, the global alpha and one-standard-error method is chosen to fit the data
      - [x] try BIC, AIC, EBIC with propre degree of freedom for gcdnet with global lambda2
      - [x] try BIC, AIC, EBIC for glmnet with global alpha
      
 - [x] Improve the piRNA normalization by using individual piRNA read count
    - To count the reads for miRNA (and piRNAs), there are just one R function to do it in Thomas' pipleline.  
      - functions.miRNASummarize.R (most comprehensive counting, considering different 5' or isoforms, count only the most abundance 5' or isoform and ignore 
        others )
      - functions.miRNASummarize.no.fixation.R (pool different 5' or isoforms) 
      - total different R function in nextflow pipeline (suppose to give same results as functions.miRNASummarize.no.fixation.R)
    - at the end my own function was coded and the conclusion for this step is that the size.factors calcluated from piRNA count table is the very similar to piRNA library size.   
    
  - [x] Double check the batch correction works 
    - The batch correction works more or less; but the the batch correction for rab-3 is much improved but still not perfect, in spit of limma or combat batch correction.   
    - The reason requires further understanding of batch correction methods.    
    - [x] double check the data before and after batch correction for example lsy-6, mir-791 and mir-790 
    - [x] By comparing the log and linear scaling with examples lsy-6, mir-791, mir-790, the observation is that the total amount of specific microRNA defined by pan-neurons is much large than other samples, consequently, the linear model must assign the rest to some neurons, which makes false positives; in contrast, in log-scale the solution will be much sparse and resonable, whereas the interpretation will be non-trivial
    - Conclusion: the combat batch correction works quite well and the log ratio (treated/untreated) does not solve the problem for rab-3 sample
  
  - [ ] Improve the model fitting, especially for known examples, lsy-6 
    - [ ] recheck the linear model, especailly the intercept for the background; perhaps we can consider promoter-specific background, which take into count the           miRNAs in the background compositions and also the promoter mythelation efficiencies.
    
    - [ ] Integrate the sample qualities or sample variance into into the linear model ??
    - [ ] post-filtering (i.e. removing non-significant ones)
  
  - [ ] Optimize the fitting and parameter
  - [ ] save the final tables and plots


### code structure
- `miRNA_enrichAnalysis_for_neuronClassess.R`
  the first code to run after raw data processing. It prepare the count table and design matrix for all samples
  Then it runs the enrichment analysis for each tissue and promoter (as in the Nature Methods paper)
- `prepare_ExpressionMatrix.R`
  normalize the miRNA expression with piRNA or test spike-in 
  prepare the expression matrix for deconvolution analysis

- `deconvolution_analysis.R`
  The main function for deconvolution analysis

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
    

## Update
  - 20181120:
    - R7011 for the new pan-neuron smaples is done, in which the pan-neuron with rab-3 promoter and with another promoter rgef-1
    
    
    
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
    
    
    