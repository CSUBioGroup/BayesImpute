
normalizeData <- function(x, y = x) {
  sf <- colSums(y)/1000000
  return(sweep(x, 2, sf, "/"))
}

#calculate gene_wise correlation
get.cor.gene <- function(X, Y) {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ]))
}

#calculate cell_wise correlation
get.cor.cell <- function(X, Y) {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i]))
}



###############################################################################

method = c('bayes')
dataset = c('baron1')
method_number=length(method)
dataset_number=length(dataset)
gene_cor_list=list()
cell_cor_list=list()
gene_cmd_list=list()
cell_cmd_list=list()
for (i in 1:method_number){
  print(method[i])
  for (ii in 1:dataset_number){
    print(dataset[ii])
    if (method[i] == 'downsample_lognorm'){
     
      data_imp_dir = paste('/home/csq/Bayesmodel/',method[i],'/',dataset[ii],'.csv',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1)
    }
    else{
      data_imp_dir = paste('/home/csq/Bayesmodel/imputed_data/',method[i],'/',dataset[ii],'_imputed.csv',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1)
    
    }
    data_ref_dir=paste('/home/csq/Bayesmodel/reference_data/','baron','_ref.csv',sep='')
    data_ref=read.csv(data_ref_dir, header = T, row.names = 1)
    data_ref<-log(normalizeData(data_ref)+1)
    data_imp=as.matrix(data_imp)
    data_ref=as.matrix(data_ref)
    
    #calculate gene_wise correlation
    gene_cor<-get.cor.gene(data_ref, data_imp)
    gene_cor_mean<-mean(gene_cor)
    #gene_cor_list[[method[i]]]<-append(gene_cor_list[[method[i]]],gene_cor_list)
    #paste("gene_cor:",gene_cor)
    print(gene_cor_mean)
    
    #calculate cell_wise correlation
    cell_cor<-get.cor.cell(data_ref, data_imp)
    cell_cor_mean<-mean(cell_cor)
    #cell_cor_list[[method[i]]]<-append(cell_cor_list[[method[i]]],cell_cor)
    #paste("cell_cor:",cell_cor)
    print(cell_cor_mean)
    
  }
}
  
