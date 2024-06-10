
##CG
CG =  read.table("Common_regions_CG_200bp_deltaM_pvalues.txt", sep = "\t", header = T)

##CHG
CHG =  read.table("Common_regions_CHG_200bp_deltaM_pvalues.txt", sep = "\t", header = T)

##CHH
CHH =  read.table("Common_regions_CHH_200bp_deltaM_pvalues.txt", sep = "\t", header = T)

#######
#genebody and TEbody
system('bash -c "bedtools intersect -wa -wb -a Common_regions_CG_200bp_deltaM_pvalues.txt -b ../tair10_genebody_20210514.bed > Common_regions_200bp_CG_GeneandTEbody.txt"')

system('bash -c "bedtools intersect -wa -wb -a Common_regions_CHG_200bp_deltaM_pvalues.txt -b ../tair10_genebody_20210514.bed > Common_regions_200bp_CHG_GeneandTEbody.txt"')

system('bash -c "bedtools intersect -wa -wb -a Common_regions_CHH_200bp_deltaM_pvalues.txt -b ../tair10_genebody_20210514.bed > Common_regions_200bp_CHH_GeneandTEbody.txt"')


system('bash -c "bedtools intersect -wa -wb -a Common_regions_CG_200bp_deltaM_pvalues.txt -b ../tair10_promoter_20210514.bed > Common_regions_200bp_CG_promoter.txt"')

system('bash -c "bedtools intersect -wa -wb -a Common_regions_CHG_200bp_deltaM_pvalues.txt -b ../tair10_promoter_20210514.bed > Common_regions_200bp_CHG_promoter.txt"')

system('bash -c "bedtools intersect -wa -wb -a Common_regions_CHH_200bp_deltaM_pvalues.txt -b ../tair10_promoter_20210514.bed > Common_regions_200bp_CHH_promoter.txt"')


CG_b = read.table("Common_regions_200bp_CG_GeneandTEbody.txt", sep = "\t")
CHG_b = read.table("Common_regions_200bp_CHG_GeneandTEbody.txt", sep = "\t")
CHH_b = read.table("Common_regions_200bp_CHH_GeneandTEbody.txt", sep = "\t")

colnames(CG_b)[1:length(CG)] = colnames(CG)
colnames(CHG_b)[1:length(CHG)] = colnames(CHG)
colnames(CHH_b)[1:length(CHH)] = colnames(CHH)

CG_b$ID = paste0(CG_b$Chr, "_", CG_b$RegionStart, "_", CG_b$RegionEnd)
CHG_b$ID = paste0(CHG_b$Chr, "_", CHG_b$RegionStart, "_", CHG_b$RegionEnd)
CHH_b$ID = paste0(CHH_b$Chr, "_", CHH_b$RegionStart, "_", CHH_b$RegionEnd)


colnames(CG_b)[33] = "gene_id"
colnames(CHG_b)[33] = "gene_id"
colnames(CHH_b)[33] = "gene_id"

# TE or gene?

#TE = read.table("/work1/home/joweihsieh/20201103_hda6jmj28/2replciates/tair10_TEbody_20210514.bed", sep = "\t", header = F)

type = read.table("/work1/home/joweihsieh/20201103_hda6jmj28/2replciates/gene_type_20210514.txt", sep = "\t", header = F)


colnames(type) = c("gene_id", "type")
CG_b = merge(CG_b, type, by = "gene_id", all.x = T)
CHG_b = merge(CHG_b, type, by = "gene_id", all.x = T)
CHH_b = merge(CHH_b, type, by = "gene_id", all.x = T)



# promoter
CG_p = read.table("Common_regions_200bp_CG_promoter.txt", sep = "\t")
CHG_p = read.table("Common_regions_200bp_CHG_promoter.txt", sep = "\t")
CHH_p = read.table("Common_regions_200bp_CHH_promoter.txt", sep = "\t")

colnames(CG_p)[1:length(CG)] = colnames(CG)
colnames(CHG_p)[1:length(CHG)] = colnames(CHG)
colnames(CHH_p)[1:length(CHH)] = colnames(CHH)

CG_p$ID = paste0(CG_p$Chr, "_", CG_p$RegionStart, "_", CG_p$RegionEnd)
CHG_p$ID = paste0(CHG_p$Chr, "_", CHG_p$RegionStart, "_", CHG_p$RegionEnd)
CHH_p$ID = paste0(CHH_p$Chr, "_", CHH_p$RegionStart, "_", CHH_p$RegionEnd)


colnames(CG_p)[33] = "gene_id"
colnames(CHG_p)[33] = "gene_id"
colnames(CHH_p)[33] = "gene_id"

CG_p = merge(CG_p, type, by = "gene_id", all.x = T)
CHG_p = merge(CHG_p, type, by = "gene_id", all.x = T)
CHH_p = merge(CHH_p, type, by = "gene_id", all.x = T)


###### promoter or body?

CG$ID = paste0(CG$Chr, "_", CG$RegionStart, "_", CG$RegionEnd)
CHG$ID = paste0(CHG$Chr, "_", CHG$RegionStart, "_", CHG$RegionEnd)
CHH$ID = paste0(CHH$Chr, "_", CHH$RegionStart, "_", CHH$RegionEnd)

for (i in 1:nrow(CG)){
		if (CG[i, "ID"]%in%CG_b$ID){
					CG[i, "body"] = "Yes"
	} else{CG[i, "body"] = "No"}
}

for (i in 1:nrow(CHG)){
		if (CHG[i, "ID"]%in%CHG_b$ID){
					CHG[i, "body"] = "Yes"
	} else{CHG[i, "body"] = "No"}
}

for (i in 1:nrow(CHH)){
		if (CHH[i, "ID"]%in%CHH_b$ID){
					CHH[i, "body"] = "Yes"
	} else{CHH[i, "body"] = "No"}
}



for (i in 1:nrow(CG)){
		if (CG[i, "ID"]%in%CG_p$ID){
					CG[i, "promoter"] = "Yes"
	} else{CG[i, "promoter"] = "No"}
}

for (i in 1:nrow(CHG)){
		if (CHG[i, "ID"]%in%CHG_p$ID){
					CHG[i, "promoter"] = "Yes"
	} else{CHG[i, "promoter"] = "No"}
}

for (i in 1:nrow(CHH)){
		if (CHH[i, "ID"]%in% $ID){
					CHH[i, "promoter"] = "Yes"
	} else{CHH[i, "promoter"] = "No"}
}



CG_1 = merge(CG, CG_b[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
CG_2 = merge(CG_1, CG_p[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
colnames(CG_2)[c(33:36)] = c("gene_id_body", "type_body", "gene_id_promoter", "type_promoter")

CHG_1 = merge(CHG, CHG_b[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
CHG_2 = merge(CHG_1, CHG_p[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
colnames(CHG_2)[c(33:36)] = c("gene_id_body", "type_body", "gene_id_promoter", "type_promoter")

CHH_1 = merge(CHH, CHH_b[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
CHH_2 = merge(CHH_1, CHH_p[, c("gene_id", "ID", "type")], by = "ID", all.x = T)
colnames(CHH_2)[c(33:36)] = c("gene_id_body", "type_body", "gene_id_promoter", "type_promoter")


write.table(CG_2, "Common_regions_CG_200bp_deltaM_pvalues_genomic_type.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(CHG_2, "Common_regions_CHG_200bp_deltaM_pvalues_genomic_type.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(CHH_2, "Common_regions_CHH_200bp_deltaM_pvalues_genomic_type.txt", sep = "\t", col.names = T, row.names = F, quote = F)


####################
# DMR (we did not use DMR in this study)
####################
treatments = 4
replicates = 2
cutoff = 0.2
# delta: 26:29
# p-value: 17:20


DMR_hda6_Col_CG = CG_2[(abs(CG_2[, 27]) >= cutoff&(CG_2[, 18]< 0.05)), ]
DMR_ldl12_Col_CG = CG_2[(abs(CG_2[, 28]) >= cutoff&(CG_2[, 19]< 0.05)), ]
DMR_hda6ldl12_Col_CG = CG_2[(abs(CG_2[, 29]) >= cutoff&(CG_2[, 20]< 0.05)), ]


write.table(DMR_hda6_Col_CG, "CG_DMR_200bp_02_hda6.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_ldl12_Col_CG, "CG_DMR_200bp_02_ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_hda6ldl12_Col_CG, "CG_DMR_200bp_02_hda6ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)



DMR_hda6_Col_CHG = CHG_2[(abs(CHG_2[, 27]) >= cutoff&(CHG_2[, 18]< 0.05)), ]
DMR_ldl12_Col_CHG = CHG_2[(abs(CHG_2[, 28]) >= cutoff&(CHG_2[, 19]< 0.05)), ]
DMR_hda6ldl12_Col_CHG = CHG_2[(abs(CHG_2[, 29]) >= cutoff&(CHG_2[, 20]< 0.05)), ]


write.table(DMR_hda6_Col_CHG, "CHG_DMR_200bp_02_hda6.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_ldl12_Col_CHG, "CHG_DMR_200bp_02_ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_hda6ldl12_Col_CHG, "CHG_DMR_200bp_02_hda6ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)

cutoff2 = 0.1
DMR_hda6_Col_CHH = CHH_2[(abs(CHH_2[, 27]) >= cutoff2&(CHH_2[, 18]< 0.05)), ]
DMR_ldl12_Col_CHH = CHH_2[(abs(CHH_2[, 28]) >= cutoff2&(CHH_2[, 19]< 0.05)), ]
DMR_hda6ldl12_Col_CHH = CHH_2[(abs(CHH_2[, 29]) >= cutoff2&(CHH_2[, 20]< 0.05)), ]


write.table(DMR_hda6_Col_CHH, "CHH_DMR_200bp_01_hda6.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_ldl12_Col_CHH, "CHH_DMR_200bp_01_ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DMR_hda6ldl12_Col_CHH, "CHH_DMR_200bp_01_hda6ldl12.txt", sep = "\t", col.names = T, row.names = F, quote = F)


