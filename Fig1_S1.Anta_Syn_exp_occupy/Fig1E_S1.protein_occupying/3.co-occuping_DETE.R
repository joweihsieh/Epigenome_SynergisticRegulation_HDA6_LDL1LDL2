hda6 = read.table("DHMG_TE_hda6.txt",header=T)
ldl1 = read.table("DHMG_TE_ldl1.txt",header=T)
dete = read.table("DETE_L12a.txt",header=T)
dete$hda6 = hda6[match(dete$Row.names,hda6$Row.names),'DHMG']
dete$ldl1 = ldl1[match(dete$Row.names,ldl1$Row.names),'DHMG']
dete$cobind = ifelse(dete$hda6 == 'Yes' & dete$ldl1 == 'Yes','Yes','No')
dete$DETE2 = ifelse(dete$DETE == 'Up', 'Up', 'No')

bed = read.table("tair10_TE.bed")
TE_HDA6_binding = bed[bed$V4 %in% hda6[hda6$DHMG=='Yes','Row.names'],]
TE_LDL1_binding = bed[bed$V4 %in% ldl1[ldl1$DHMG=='Yes','Row.names'],]
TE_HDA6LDL1_cobinding = TE_HDA6_binding[TE_HDA6_binding$V4 %in% TE_LDL1_binding$V4,]
write.table(TE_HDA6_binding,file = "TE_HDA6_bound.bed",sep = '\t',quote=F, row.names=F, col.names=F)
write.table(TE_LDL1_binding,file = "TE_LDL1_bound.bed",sep = '\t',quote=F, row.names=F, col.names=F)
write.table(TE_HDA6LDL1_cobinding,file = "TE_HDA6LDL1_cobound.bed",sep = '\t',quote=F, row.names=F, col.names=F)



### hypergeometric test
q = nrow(dete[dete$DETE=='Up' & dete$cobind =='Yes',])
m = nrow(dete[dete$DETE=='Up',])
n = nrow(dete[dete$DETE!='Up',])
k = nrow(dete[dete$cobind =='Yes',])
p = 1 - phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)


table(dete$DETE2, dete$cobind)
#chisq.test(table(dete$DETE2, dete$cobind))


