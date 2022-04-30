
dataset = c('sc_10x_5cl','FiveEncode_GSE81861')
method = c('bayes','deepimpute','SAVER','scimpute','VIPER','bayNorm')
bulk_dataset = c('CellBench_GSE86337_averagelog2TPM','FiveEncode_averagelog2TPM')
dataset_number=length(dataset)
method_number=length(method)


pearsonCor_imputed=list()
pearsonCor_raw=list()
pearsonCor_pseudo=list()

for (ii in 1:dataset_number){
  print(dataset[ii])
  #####load original data, load imputed data, load bulk data
  bulk_dir = paste('/home/csq/processed_data/bulk/',bulk_dataset[ii],'.txt', sep='')
  label_dir = paste('/home/csq/processed_data/label/',dataset[ii],'_labelName.txt',sep='')
  originalSC_dir = paste('/home/csq/processed_data/scranNorm/',dataset[ii],'_genebycell.txt', sep='')
  # Load the single-cell label
  labelName = read.table(label_dir,sep = "\t",header = FALSE)
  # Load the single-cell dataset
  originalData = read.table(originalSC_dir,sep = "\t",header = TRUE)
  # Load the bulk dataset
  if (dataset[ii] %in% c('sc_10x','sc_celseq2','sc_dropseq')){
    averagebulkData = read.table(bulk_dir,sep = "\t",header = TRUE)
    averagebulkData = averagebulkData[,c('HCC827','H1975','H2228')]
  }else if (dataset[ii] == '293T_Jurkat'){
    bulkData = read.table(bulk_dir,sep = "\t",header = TRUE)
    print(dim(bulkData))
    colnames(bulkData)=c('293T','293T','Jurkat','Jurkat')
    bulkname = unique(colnames(bulkData))
    print(bulkname)
    averagebulkData = sapply(bulkname,function(i) rowMeans(bulkData[, colnames(bulkData) %in% i]))
  }else{
    averagebulkData = read.table(bulk_dir,sep = "\t",header = TRUE)
  }
  ##output dataset dim
  print(dim(labelName))
  print(dim(originalData))
  print(dim(averagebulkData))
  
  
 
  for (j in 1:method_number) {
    print(method[j])
    # Load the imputed dataset
    if (method[j] == 'scimpute'){
      imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[ii],"_genebycell_count.txtscimpute_count.txt", sep='')
      imputedData = read.table(imputedSC_dir,sep = " ",header = TRUE,row.names = 1)
      libsize = colSums(imputedData)
      normal_data <- sweep(imputedData,2,libsize,'/')
      imputedData <- log2(normal_data + 1)
    }else if (method[j] == 'DCA'){
      imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[ii],"_genebycell.txt", sep='')
      imputedData = read.table(imputedSC_dir,sep = "\t",header = TRUE,row.names=1)
      libsize = colSums(imputedData) 
      normal_data <- sweep(imputedData,2,libsize,'/')
      imputedData <- log2(normal_data + 1)
    }else if (method[j] == 'bayNorm'){
      imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[ii],"_genebycell.txt", sep='')
      imputedData = read.table(imputedSC_dir,sep = "\t",header = TRUE,row.names=1)
      imputedData <- log2(imputedData + 1)
    }else{
      imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[ii],"_genebycell.txt", sep='')
      imputedData = read.table(imputedSC_dir,sep = "\t",header = TRUE,row.names = 1)
    }
    
   
    matchID = intersect(rownames(originalData), intersect(rownames(imputedData),rownames(averagebulkData)))
    matchoriginalData = originalData[matchID, ]
    matchbulkData = averagebulkData[matchID, ]
    matchimputedData = imputedData[matchID, ]
    cellLine = unique(labelName$V1)  
    
    
    matchimputedData = as.data.frame(matchimputedData)
    spearmanCor=list()
    Cor = vector()
    spearmanCor_mean = vector()
    k1 = 0
    for(kid in cellLine){
      #print(kid)
      k1=k1+1
      for(i in 1:dim(matchimputedData[, labelName$V1 %in% kid])[2]){
        Cor[i] = cor(matchimputedData[, labelName$V1 %in% kid][i],matchbulkData[,(colnames(matchbulkData)) %in% kid],method='pearson')
      }
      #print(mean(CorBYS))
      spearmanCor[[kid]] = Cor
      spearmanCor_mean[kid] = mean(Cor)
      Cor= vector()
    }
    pearsonCor_imputed[[method[j]]]<-append(pearsonCor_imputed[[method[j]]],spearmanCor_mean)
    print(spearmanCor_mean)
  }
  
  
  spearmanCor_raw=list()
  Cor_raw = vector()
  spearmanCor_raw_mean = vector()
  k = 0
  for(kid in cellLine){
    k=k+1
    for(i in 1:dim(matchoriginalData[, labelName$V1 %in% kid])[2]){
      Cor_raw[i] = cor(matchoriginalData[, labelName$V1 %in% kid][i],matchbulkData[,(colnames(matchbulkData)) %in% kid],method='pearson')
    }
    spearmanCor_raw[[kid]] = Cor_raw
    spearmanCor_raw_mean[kid] = mean(Cor_raw)
    Cor_raw = vector()
    pearsonCor_raw[[kid]]<-append(pearsonCor_raw[[kid]],spearmanCor_raw_mean[kid])
  }
  print(spearmanCor_raw_mean)
  
  
  
  pseudoCor=vector()
  pseudobulk = sapply(cellLine,function(i) rowMeans(matchoriginalData[, labelName$V1 %in% i]))
  colnames(pseudobulk) = cellLine
  for(kid in cellLine){
    pseudoCor[kid] = cor(pseudobulk[,colnames(pseudobulk)==kid],matchbulkData[,(colnames(matchbulkData)) %in% kid],method='pearson')
    pearsonCor_pseudo[[kid]]<-append(pearsonCor_pseudo[[kid]],pseudoCor[kid])
  }
  print(pseudoCor)
}




