WT1 = read.table("ATAC_WT_1_body.tab")
WT2 = read.table("ATAC_WT_2_body.tab")
hda61 = read.table("ATAC_hda6_1_body.tab")
hda62 = read.table("ATAC_hda6_2_body.tab")
hda6ldl121 = read.table("ATAC_hda6ldl12_1_body.tab")
hda6ldl122 = read.table("ATAC_hda6ldl12_2_body.tab")
ldl121 = read.table("ATAC_ldl12_1_body.tab")
ldl122 = read.table("ATAC_ldl12_2_body.tab")

RNA = read.table("Exprssed_TE.txt",header=T)
RNA$RNA_WT = (RNA$C_6 + RNA$C_7 + RNA$C_7_2) /3
RNA$RNA_hda6 = (RNA$a_6 + RNA$a_7) /2
RNA$RNA_ldl12 = (RNA$L12_6 + RNA$L12_7) / 2
RNA$RNA_hda6ldl12 = (RNA$L12a_6 + RNA$L12a_7) /2


df = data.frame(row.names = WT1$V1, WT = (WT1$V6 + WT2$V6) /2, hda6 = (hda61$V6 + hda62$V6) /2, ldl12 = (ldl121$V6 + ldl122$V6) /2, hda6ldl12 = (hda6ldl121$V6 + hda6ldl122$V6) /2, WT1=WT1$V6, WT2=WT2$V6, hda61=hda61$V6, hda62=hda62$V6, ldl121 = ldl121$V6, ldl122 = ldl122$V6, hda6ldl121 = hda6ldl121$V6, hda6ldl122 = hda6ldl122$V6)
df$WT_max = apply(df[,c('WT1','WT2')],1,max)
df$hda6_min = apply(df[,c('hda61','hda62')],1,min)
df$ldl12_min = apply(df[,c('hda61','ldl122')],1,min)

#df$hda6_opened = ifelse(df$WT < df$hda6, 'hda6_opened','closed')
#df$ldl12_opened = ifelse(df$WT < df$ldl12, 'ldl12_opened','closed')

df$hda6_opened = ifelse(df$WT_max < df$hda6_min, 'hda6_opened','closed')
df$ldl12_opened = ifelse(df$WT_max < df$ldl12_min, 'ldl12_opened','closed')

df2 = merge(df, RNA[,c(1,24:27)],by.x = 'row.names',by.y = 'Row.names')
row.names(df2) = df2$Row.names
df = df2[,2:ncol(df2)]

table(df$hda6_opened)
table(df$ldl12_opened)
table(df$hda6_opened, df$ldl12_opened)


# double-locked
df_001 = df[row.names(df) %in% TE[TE$types == '0_0_1','Row.names'],]
table(df_001$hda6_opened)
table(df_001$ldl12_opened)
table(df_001$hda6_opened, df_001$ldl12_opened)

# HDA6-predominant
df_101 = df[row.names(df) %in% TE[TE$types == '1_0_1','Row.names'],]
table(df_101$hda6_opened)
table(df_101$ldl12_opened)
table(df_101$hda6_opened, df_101$ldl12_opened)


