library(edgeR)
library(limma)

# DE analysis using edgeR function

edgeR_analysis = function(samples_info, count_values, pair, out_filename = FALSE) {
  
  # input: 
  # samples_info - table of samples including columns-'name','treatment'
  # count_values - table of count-values from Salmon, for each gene and sample
  #     NOTE: samples_info and count_values - for a specific species and control+stimuation
  # pair = c('control','stimuation')
  
  # create DGE objects
  de_obj = DGEList(counts=count_values, group=samples_info$treatment)
  de_obj = calcNormFactors(de_obj, method = "TMM")
  de_obj = estimateDisp(de_obj)
  
  # run exact test
  et = exactTest(de_obj, pair = pair)
  
  # add QValue
  et$table$QValue = p.adjust(et$table$PValue, "BH")
  
  # order by FDR
  et$table = et$table[order(et$table$QValue),]
  
  # write output
  if (out_filename == FALSE) {
    return(et$table)
  }
  write.table(et$table, file = out_filename, sep = "\t")
}
