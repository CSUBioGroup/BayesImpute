dataset = c('FiveEncode_GSE81861_2cl')
method = c('Bayes','scimpute','DrImpute','ALRA','scRMD','MAGIC','DCA','deepimpute','SAVER')
method_number = length(method)

# Load the single-cell label
label_dir = paste('/label/',dataset[1],'_labelName.txt',sep='')
labelName = read.table(label_dir,sep = "\t",header = FALSE)
print(dim(labelName))
cellLine = unique(labelName$V1)
cellLine_num = length(cellLine) 
print(cellLine)

cellLine[1]='A549'
cellLine[2]='IMR90'


print(cellLine)
# Load the original single cell dataset
originalSC_dir = paste('/home/csq/processed_data/scranNorm/',dataset[1],'_genebycell.txt', sep='')
originalData = read.table(originalSC_dir,sep = "\t",header = TRUE)
print(dim(originalData))
originalData <- as.matrix(originalData)

#bulkdata
bulk_dir = '/home/csq/processed_data/bulk/FiveEncode_rmMTrawTPM.txt'
# Load the bulk dataset
bulkData = read.table(bulk_dir,sep = "\t",header = TRUE,row.names = 1)
colnames(bulkData)=c('A549','A549','IMR90','IMR90')
print(dim(bulkData))
print(colnames(bulkData))
cellLine_bulk = unique(colnames(bulkData))
print(cellLine_bulk)
print(cellLine==cellLine_bulk)


unoin_gene=intersect(rownames(bulkData),rownames(originalData))
print(length(unoin_gene))


diffgene_raw_num = list()
diffgene_bulk_num = list()
overlap_raw_bulk = list()
overlap_imputed_bulk = list()
diffgene_fdr_imputed_ascend = list()
diffgene_fdr_ascend = list()
diffgene_fdr_bulk_ascend = list()


k=5000
for (i in 1:(cellLine_num-1)) {
  for (ii in (i+1):cellLine_num){
    print(paste(cellLine[i],'vs',cellLine[ii]))
    
    for (j in 1:method_number) {
      print(method[j])
      # Load the imputed dataset
      if (method[j] == 'scimpute'){
        imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[1],"_genebycell_count.txtscimpute_count.txt", sep='')
        imputedData = read.table(imputedSC_dir,sep = " ",header = TRUE,row.names = 1)
        libsize = colSums(imputedData)
        normal_data <- sweep(imputedData,2,libsize,'/')
        imputedData <- log2(normal_data + 1)
        imputedData = imputedData[unoin_gene, ]
      }else if (method[j] == 'DCA'){
        imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[1],"_genebycell.txt", sep='')
        imputedData = read.table(imputedSC_dir,sep = "\t",header = TRUE,row.names=1)
        libsize = colSums(imputedData)
        normal_data <- sweep(imputedData,2,libsize,'/')
        imputedData <- log2(normal_data + 1)
        imputedData = imputedData[unoin_gene, ]
      }else{
        imputedSC_dir =paste('/home/csq/Bayesmodel/imputed_data/',method[j],'/',dataset[1],"_genebycell.txt", sep='')
        imputedData = read.table(imputedSC_dir,sep = "\t",header = TRUE,row.names = 1)
        imputedData = imputedData[unoin_gene, ]
      }
      print(dim(imputedData))
      ##input imputed dataset,the raw is gene,the column is cell
      imputedData <- as.matrix(imputedData)
      
      ####test  cellline in the dataset
      ## wilcoxon-imputedData
      pval_imputed = sapply(1:nrow(imputedData), function(i1) {
        cellLine1_imputed = imputedData[i1,labelName$V1 %in% cellLine[i]]
        cellLine2_imputed = imputedData[i1,labelName$V1 %in% cellLine[ii]]
        names(cellLine1_imputed) <- names(cellLine2_imputed) <- NULL
        if (identical(cellLine1_imputed,cellLine2_imputed)) {
          1
        }else{
          res = wilcox.test(cellLine1_imputed,cellLine2_imputed)
          res$p.value
        }
      })
      ##output the number of diffgene less than 0.05
      fdr_imputed = p.adjust(pval_imputed, method='fdr')
      fdr_list_imputed=list(fdr_imputed)

      fdr_imputed_ascend=order(fdr_list_imputed[[1]])
      diffgene_fdr_imputed_ascend[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_fdr_imputed_ascend[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_imputed_ascend)
      fdr_imputed_top=fdr_imputed_ascend[1:k]
      diffgene_imputed_num[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_imputed_num[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_imputed_top)
      
      ## wilcoxon-originalData
      originalData = originalData[unoin_gene, ]
      pval = sapply(1:nrow(originalData), function(i1) {
        cellLine1_original = originalData[i1,labelName$V1 %in% cellLine[i]]
        cellLine2_original = originalData[i1,labelName$V1 %in% cellLine[ii]]
        names(cellLine1_original) <- names(cellLine2_original) <- NULL
        if (identical(cellLine1_original,cellLine2_original)) {
          1
        }else{
          res = wilcox.test(cellLine1_original,cellLine2_original)
          res$p.value
        }
      })
      ##output the number of diffgene that p.adjust less than 0.05
      fdr = p.adjust(pval, method='fdr')
      fdr_list=list(fdr)
      # diffgene_raw=which(fdr_list[[1]] < 0.05)
      # print(length(diffgene_raw))
   
      fdr_ascend=order(fdr_list[[1]])
      diffgene_fdr_ascend[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_fdr_ascend[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_ascend)
      fdr_top=fdr_ascend[1:k]
      diffgene_raw_num[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_raw_num[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_top)
      
      
      ##DEseq_bulkdata
      bulkData1 <- bulkData[unoin_gene,colnames(bulkData) %in% c(cellLine_bulk[i],cellLine_bulk[ii])]
      a1=sum(str_detect(colnames(bulkData1),cellLine_bulk[i]))
      a2=sum(str_detect(colnames(bulkData1),cellLine_bulk[ii]))     
      condition <- factor(c(rep(cellLine_bulk[i],a1),rep(cellLine_bulk[ii],a2)))
      coldata <- data.frame(row.names = colnames(bulkData1), condition)
      dds <- DESeqDataSetFromMatrix(countData=bulkData1, colData=coldata, design=~condition)
      dds <- DESeq(dds) 
      res <- results(dds)
      
      
      alpha=0.05
      padj_list=list(res$padj) 
      fdr_bulk_ascend=order(padj_list[[1]])
      diffgene_fdr_bulk_ascend[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_fdr_bulk_ascend[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_bulk_ascend)
      fdr_bulk_top=fdr_bulk_ascend[1:k]
      diffgene_bulk_num[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(diffgene_bulk_num[[paste(cellLine[i],'vs',cellLine[ii])]],fdr_bulk_top)

      bulk_diffgene_top1000=rownames(bulkData)[fdr_bulk_top]
      data_top1000=rownames(originalData)[fdr_top]
      data_imputed_top1000=rownames(imputedData)[fdr_imputed_top]
      overlap_raw_top1000=length(intersect(bulk_diffgene_top1000,data_top1000))
      overlap_imputed_top1000=length(intersect(bulk_diffgene_top1000,data_imputed_top1000))
      
      overlap_raw_bulk[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(overlap_raw_bulk[[paste(cellLine[i],'vs',cellLine[ii])]],overlap_raw_top1000)
      overlap_imputed_bulk[[paste(cellLine[i],'vs',cellLine[ii])]]<-append(overlap_imputed_bulk[[paste(cellLine[i],'vs',cellLine[ii])]],overlap_imputed_top1000)
      print(length(bulk_diffgene_top1000))
      print(overlap_raw_top1000)
      print(overlap_imputed_top1000)
    }
  }
}

