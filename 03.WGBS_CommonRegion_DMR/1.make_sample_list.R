#/work1/home/joweihsieh/20190731_HDA6/WGBS/raw_data/rawData/bam

library(magrittr)
Filelist=list.files(pattern=".CGmap.gz")

list2 = Filelist %>% sapply(function(x)strsplit(x,'\\.')[[1]][1])%>% sapply(function(x)strsplit(x,'-')[[1]][2]) 


datalist=list()
for (ff in seq(1:length(Filelist))){
        File1=Filelist[ff]
        code=paste(list2[ff],"",File1)
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)
write.table(dataset,"sample_list.txt",quote=F,sep="\t",row.names=F,col.names=F)

