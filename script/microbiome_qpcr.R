#######################
# set working directory

setwd('/Users/hou/Desktop/tci/microbiome_eve/')

#########################
# read mirobiome qpcr raw 

microbiome = read.csv("191028/20191028-G2190910.csv")

################################
# table of sample and each genes
table_micro = matrix(NA, nrow = 14*3*3, ncol = 3)

# add column and row name
rownames(table_micro) = microbiome$Sample.Name[1:126]

colnames(table_micro) = unique(microbiome$Target.Name)

# add 16s CT
table_micro[, 1] = microbiome$CT[1:126]

# add Lacto CT
table_micro[, 2] = microbiome$CT[127:252]

# add LPI CT
table_micro[, 3] = microbiome$CT[253:378]

write.csv(table_micro, '191028/result/20191028-G2190910_ct_table.csv')


##############
# table of dCT

# mean of 3 replicate in 16S
row_num = seq(1, dim(table_micro)[1], by = 3)

mean_16s = matrix(NA, nrow = 42)

for(i in 1:length(row_num)){
  
  mean_16s[i] = mean(table_micro[row_num[i]:(row_num[i]+2),1])

}


sample_row = seq(1, 126, 3)
rownames(mean_16s) = microbiome$Sample.Name[sample_row]
colnames(mean_16s) = c('mean_16s')

# table of dct
dct = matrix(NA, nrow = 126, ncol = 2)

rownames(dct) = rownames(table_micro)
colnames(dct) = colnames(table_micro)[2:3]

for (row in 1:length(mean_16s)){
  
  for (col in 1:2){
  
    dct[(row*3-2):(row*3),col] = table_micro[(row*3-2):(row*3), col+1] - mean_16s[row]
  
  }
}

write.csv(dct, '191028/result/20191028-G2190910_dct.csv')

###############
# table of ddCT

# mean of W0 dCT

mean_w0_dct = matrix(NA, nrow = 14 , ncol = 2)

sample_row2 = seq(1, 42, 3)
rownames(mean_w0_dct) = rownames(dct[1:42,])[sample_row2]
colnames(mean_w0_dct) = colnames(dct)


for (row in 1:14){
  
  for (col in 1:2){
     
    mean_w0_dct[row, col] = mean(dct[(row*3-2):(row*3), col])
  
  }
}

# table ddCT

dct_W2 = dct[43:84,]
dct_W4 = dct[85:126,]


ddct_W2 = matrix(NA, nrow = 42, ncol = 2)
rownames(ddct_W2) = rownames(dct_W2)
  
ddct_W4 = matrix(NA, nrow = 42, ncol = 2)
rownames(ddct_W4) = rownames(dct_W4)


for (row in 1:dim(mean_w0_dct)[1]){
  
  for (col in 1:2){
    
    ddct_W2[(row*3-2):(row*3),col] = dct_W2[(row*3-2):(row*3), col] - mean_w0_dct[row,col]
    ddct_W4[(row*3-2):(row*3),col] = dct_W4[(row*3-2):(row*3), col] - mean_w0_dct[row,col]
  
  }
}


ddct = rbind(ddct_W2, ddct_W4)
colnames(ddct) = colnames(dct)

write.csv(ddct, '191028/result/20191028-G2190910_ddct.csv')

#############
# fold change

fold_change = 2^(-ddct)

fold_change_sci = format(fold_change, scientific = FALSE)

write.csv(fold_change_sci, '191028/result/20191028-G2190910_fold_change.csv')


########
# result

#Avg
Avg = matrix(NA, nrow = 3, ncol = 2)
rownames(Avg) = c('W0', 'W2', 'W4')
colnames(Avg) = c('Lacto', 'LPI')

Avg[1,] = 1

Avg[2, 1] = mean(fold_change[1:42, 1])
Avg[3, 1] = mean(fold_change[43:84, 1])

Avg[2, 2] = mean(fold_change[1:42, 2])
Avg[3, 2] = mean(fold_change[43:84, 2])

# Std

std = matrix(NA, nrow = 3, ncol = 2)
rownames(Avg) = c('W0', 'W2', 'W4')
colnames(Avg) = c('Lacto', 'LPI')

std[1,] = 0

std[2, 1] = sd(fold_change[1:42, 1])
std[3, 1] = sd(fold_change[43:84, 1])

std[2, 2] = sd(fold_change[1:42, 2])
std[3, 2] = sd(fold_change[43:84, 2])

# SEM

sem = std/sqrt(3)



####################
# plot group barplot
library(ggplot2)

# create the dataframe for plots
bacteria = c(rep("Lacto" , 3) , rep("LPI" , 3))
week = rep(c("W0" , "W2", "W4") , 2)
value = c(Avg)

data = data.frame(bacteria, week,value)

# group barplot

ggplot(data, aes(fill = week, y = value, x = bacteria)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab('Bacteria') +
  ylab('Fold change') +
  scale_fill_manual(values=c("#AED6F1", "#3498DB", "#21618C"))

ggplot(data, aes(fill = week, y = value, x = bacteria)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab('Bacteria') +
  ylab('Fold change') 

