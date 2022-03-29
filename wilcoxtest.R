
dataset = c('sc_10x_5cl')

# Load the single-cell label
label_dir = paste('/label/',dataset[1],'_labelName.txt',sep='')
labelName = read.table(label_dir,sep = "\t",header = FALSE)
print(dim(labelName))
cellLine = unique(labelName$V1)
cellLine_num = length(cellLine)  

# Load the single cell dataset
originalSC_dir = paste('/processed_data/scranNorm/',dataset[1],'_genebycell.txt', sep='')
originalData = read.table(originalSC_dir,sep = "\t",header = TRUE)
print(dim(originalData))
originalData <- as.matrix(originalData)

####test all cellline in the dataset
for (i in 1:(cellLine_num-1)) {
  for (ii in (i+1):cellLine_num){
    print(paste(cellLine[i],'vs',cellLine[ii]))
    ## wilcoxon
    pval = sapply(1:nrow(originalData), function(i1) {
      cellLine_one = originalData[i1,labelName$V1 %in% cellLine[i]]
      cellLine_two = originalData[i1,labelName$V1 %in% cellLine[ii]]
      names(cellLine_one) <- names(cellLine_two) <- NULL
      if (identical(cellLine_one,cellLine_two)) {
        1
      }else{
        res = wilcox.test(cellLine_one,cellLine_two)
        res$p.value
      }
    })
    ##output the number of diffgene that p.adjust less than 0.05
    fdr = p.adjust(pval, method='fdr')
    fdr_list=list(fdr)
    diffgene_raw=which(fdr_list[[1]] < 0.05)
    print(length(diffgene_raw))
  }
}






