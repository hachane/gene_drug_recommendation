
getGene<-function(a){
  pstar = regexpr("star", a)
  substr(a, 1, pstar - 1)
}

dataOrig=read.table("CYPs_hap_frq_output_2020oct", header=FALSE,stringsAsFactors=FALSE)
v2new=gsub("\\*","star",dataOrig$V2) # to avoid special treatment to '*'
data=data.frame("V1"=dataOrig$V1,"V2"=v2new,"V3"=dataOrig$V3,stringsAsFactors=FALSE)

#pop=c("YRI","CEU","SAS","JPT","CHB","CHS","KHV")
#popNSamples=c(186,183,661,105,108,171,206)

pop=c("CEU","CHB","CHS","JPT","SAS","YRI","KHV")
#popNSamples=c(183,108,171,105,661,186,206)
popNSamples=c(99,103,105,104,489,108,206)

geneObj=unique(getGene(data$V2))
gene=c()
for(i in 1:length(geneObj)){
  gene[i] = toString(geneObj[i])
}

gene=gene[order(substr(gene, 1, 6))]
nGene=length(gene)
nPop=length(pop)

pValueTbl=data.frame(Gene=character(),Hap=character(),
                     CEU=double(),
                     CHB=double(),
                     CHS=double(),
                     JPT=double(),
                     SAS=double(),
                     YRI=double(),                     
                     stringsAsFactors=FALSE)

for(i in 1:length(gene)){
  print(paste("Analyzing ", gene[i], "...", sep=""))
  hapObj=unique(data$V2[grepl(paste(gene[i],"star",sep=""),data$V2)])
  nHap=length(hapObj)
  
  # compute COUNTS matrix haplotype (row) by population (col)
  popCount=matrix(nrow=nHap, ncol=nPop) 
  dimnames(popCount) = list(hapObj,pop)
  for(j in 1:length(pop)){
    p=pop[j]
    prows=data[data$V1 == p,]
    for(h in hapObj){
      pf = prows[prows$V2 == h,]
      if(nrow(pf) == 1){
        popCount[h, p] = trunc(pf$V3 * 2 * popNSamples[j] / 100)
      } else popCount[h, p] = 0
    }
  }

  write.table(popCount, file=paste(gene[i],'.pop.hap.2020oct.tsv'), quote=FALSE, sep='\t', row.names=TRUE)
  
  print(popCount[1,])
  
  for(hid in 1:nrow(popCount)){
#    test=pairwise.prop.test(unlist(popCount[hid,]), colSums(popCount), p.adjust.method = "bonferroni")
    test=pairwise.prop.test(unlist(popCount[hid,]), 2 * popNSamples)
    print(test)
#    test=pairwise.prop.test(unlist(popCount[hid,]), popNSamples)
    pValueTbl[nrow(pValueTbl) + 1,] = c(gene[i],hapObj[hid],
                                        test[[3]][6,1],
                                        test[[3]][6,2],
                                        test[[3]][6,3],
                                        test[[3]][6,4],
                                        test[[3]][6,5],
                                        test[[3]][6,6])
  }
}

print(pValueTbl)
write.table(pValueTbl, file='pairwise.prop.test.redo3.adjust.2020oct.tsv', quote=FALSE, sep='\t', row.names=FALSE)
