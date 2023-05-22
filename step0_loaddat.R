library(dplyr)
library(stringr)
rm(list = ls())  
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 
## batch reading data
### Set the data path and sample name
#### data path
file_list = list.files('../AM/filtered/')


#### sample name
# split function
samples = str_split(file_list,"_",simplify = T)[,2]
setwd("../AM/filtered/")
lapply(unique(samples),function(x){
  y = file_list[grepl(x,file_list)]
  folder=paste(str_split(y[1],'_',simplify = T)[,2],collapse = '') 
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"genes.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})
setwd("/home/data/t070421/adenomyosis/scrna/")
