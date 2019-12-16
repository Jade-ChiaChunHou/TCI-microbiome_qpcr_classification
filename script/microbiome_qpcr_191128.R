library(varhandle)
library(stringr)
#######################   
# set working directory

setwd('/home/tcigene/jade/microbiome_eve/raw/EVE_G20190325')

result_path = "../../result/191128/"
sample_name = "EVE-G20190325"
#########################
# read mirobiome qpcr raw 

microbiome = read.csv("EVE-G20190325.csv")
microbiome = unfactor(microbiome)

for (row in 1:dim(microbiome)[1]){
  
  if(microbiome[row,3] == "Undetermined"){
    
      microbiome[row,3] = "40"
    
  }
}


##############################
# order by bacteria in columns

bac_name = unique(microbiome$Target.Name)

bac_table = matrix(NA, nrow = (dim(microbiome)[1] / length(bac_name)), ncol = length(bac_name))

rownames(bac_table) = microbiome[which(microbiome$Target.Name == "16S"),c("Sample.Name")]
colnames(bac_table) = bac_name


for (col in 1:length(bac_name)){
  
  bac_table[,col] = microbiome[which(microbiome$Target.Name == bac_name[col]),c("CT")]
  
}

#bac1 = microbiome[which(microbiome$Target.Name == bac_name[1]),c("Sample.Name", "CT")]

path = paste(result_path, sample_name, "_ct_table.csv", sep = "")
write.csv(bac_table, path)


##############
# table of dCT

# mean of 3 replicate in 16S
row_num = seq(1, dim(bac_table)[1], by = 3)

mean_16s = matrix(NA, nrow = (dim(bac_table)[1]/3))


for(i in 1:length(row_num)){
  
  print(bac_table[row_num[i]:(row_num[i]+2),7])
  mean_16s[i] = mean(as.numeric(bac_table[row_num[i]:(row_num[i]+2), which(bac_name == "16S")]))

}


sample_row = seq(1, dim(bac_table)[1], 3)
rownames(mean_16s) = rownames(bac_table)[sample_row]
colnames(mean_16s) = c('mean_16s')

# table of dct
dct = matrix(NA, nrow = dim(bac_table)[1], ncol = (dim(bac_table)[2] - 1))

rownames(dct) = rownames(bac_table)
colnames(dct) = colnames(bac_table)[which(bac_name != "16S")]

for (row in 1:length(mean_16s)){
  
  for (col in 1:dim(dct)[2]){
  
    dct[(row*3-2):(row*3),col] = as.numeric(bac_table[(row*3-2):(row*3), col]) - mean_16s[row]
  
  }
}

path = paste(result_path, sample_name, "_dct.csv", sep = "")
write.csv(dct, path)

###############
# table of ddCT

# timepoint 1: get dct first timepoint
timepoint1_dct = dct[which((str_sub(rownames(dct), start= -1)) == "1"), ]

# mean of W0 dCT

mean_timepoint1_dct = matrix(NA, nrow = (dim(timepoint1_dct)[1] / 3) , ncol = dim(timepoint1_dct)[2])

sample_row2 = seq(1, dim(timepoint1_dct)[1], 3)
rownames(mean_timepoint1_dct) = rownames(timepoint1_dct)[sample_row2]
colnames(mean_timepoint1_dct) = colnames(timepoint1_dct)


for (row in 1:dim(mean_timepoint1_dct)[1]){
  
  for (col in 1:dim(mean_timepoint1_dct)[2]){
     
    mean_timepoint1_dct[row, col] = mean(timepoint1_dct[(row*3-2):(row*3), col])
  
  }
}


############
# table ddCT

# timepoint 2: get dct second timepoint
timepoint2_dct = dct[which((str_sub(rownames(dct), start= -1)) == "2"), ]

# timepoint 3: get dct second timepoint
timepoint3_dct = dct[which((str_sub(rownames(dct), start= -1)) == "3"), ]


ddct_timepoint2 = matrix(NA, nrow = dim(timepoint2_dct)[1], ncol = dim(timepoint2_dct)[2])
rownames(ddct_timepoint2) = rownames(timepoint2_dct)
colnames(ddct_timepoint2) = colnames(timepoint2_dct)

ddct_timepoint3 = matrix(NA, nrow = dim(timepoint3_dct)[1], ncol = dim(timepoint3_dct)[2])
rownames(ddct_timepoint3) = rownames(timepoint3_dct)
colnames(ddct_timepoint3) = colnames(timepoint3_dct)


for (row in 1:dim(mean_timepoint1_dct)[1]){
  
  for (col in 1:dim(mean_timepoint1_dct)[2]){
    
    ddct_timepoint2[(row*3-2):(row*3),col] = timepoint2_dct[(row*3-2):(row*3), col] - mean_timepoint1_dct[row,col]
    ddct_timepoint3[(row*3-2):(row*3),col] = timepoint3_dct[(row*3-2):(row*3), col] - mean_timepoint1_dct[row,col]
  
  }
}


ddct = rbind(ddct_timepoint2, ddct_timepoint3)

path = paste(result_path, sample_name, "_ddct.csv", sep = "")
write.csv(ddct, path)


#############
# fold change

fold_change = 2^(-ddct)

fold_change_sci = format(fold_change, scientific = FALSE)

path = paste(result_path, sample_name, "_fold_change.csv", sep = "")
write.csv(fold_change_sci, path)


########
# result

#Avg

avg = matrix(NA, nrow = (dim(fold_change)[1] / 3), ncol = dim(fold_change)[2])
colnames(avg) = colnames(fold_change)

sample_row3 = seq(1, dim(fold_change)[1], 3)
rownames(avg) = rownames(fold_change)[sample_row3]


for (row in 1:dim(avg)[1]){
  
  for (col in 1:dim(avg)[2]){
    
    avg[row, col] = mean(fold_change[(row*3-2):(row*3), col])
    
  }
}

avg_sci = format(avg, scientific = FALSE)

path = paste(result_path, sample_name, "_avg_sci.csv", sep = "")
write.csv(avg_sci, path)


# SEM = Std / sqrt(sample_size)

std = matrix(NA, nrow = (dim(fold_change)[1] / 3), ncol = dim(fold_change)[2])
colnames(std) = colnames(fold_change)
rownames(std) = rownames(fold_change)[sample_row3]


for (row in 1:dim(std)[1]){
  
  for (col in 1:dim(std)[2]){
    
    std[row, col] = sd(fold_change[(row*3-2):(row*3), col]) / sqrt(3)
    
  }
}

std_sci = format(std, scientific = FALSE)

path = paste(result_path, sample_name, "_std_sci.csv", sep = "")
write.csv(std_sci, path)



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













#######################
# For Selected ct table
#######################

#bac_table2 = read.csv('EVE-G20190325_ct_table_v2.csv')
bac_table2 = read.csv('EVE-G20190325_ct_table_v3.csv')
#rownames(bac_table2) = microbiome[which(microbiome$Target.Name == "16S"),c("Sample.Name")]
rowname = bac_table2$X
bac_table2 = bac_table2[,2:8]

# the repeat number after selection
sample_table = data.frame(table(rowname))

path = paste(result_path, sample_name, "_sample_table.csv", sep = "")
write.csv(sample_table, path)

##############
# table of dCT

# mean of 3 replicate in 16S
row_num = seq(1, dim(bac_table2)[1], by = 3)

mean_16s = matrix(NA, nrow = (dim(bac_table2)[1]/3))


for(i in 1:length(row_num)){
  
  print(bac_table2[row_num[i]:(row_num[i]+2),7])
  mean_16s[i] = mean(as.numeric(bac_table2[row_num[i]:(row_num[i]+2), which(bac_name == "16S")]), na.rm = T)
  
}


sample_row = seq(1, dim(bac_table2)[1], 3)
rownames(mean_16s) = rowname[sample_row]
colnames(mean_16s) = c('mean_16s')

# table of dct
dct = matrix(NA, nrow = dim(bac_table2)[1], ncol = (dim(bac_table2)[2] - 1))

rownames(dct) = rowname
colnames(dct) = colnames(bac_table2)[which(bac_name != "16S")]

for (row in 1:length(mean_16s)){
  
  for (col in 1:dim(dct)[2]){
    
    dct[(row*3-2):(row*3),col] = as.numeric(bac_table2[(row*3-2):(row*3), col]) - mean_16s[row]
    
  }
}

path = paste(result_path, sample_name, "_dct_select_v3.csv", sep = "")
write.csv(dct, path)

###############
# table of ddCT

# timepoint 1: get dct first timepoint
timepoint1_dct = dct[which((str_sub(rownames(dct), start= -1)) == "1"), ]

# mean of W0 dCT

mean_timepoint1_dct = matrix(NA, nrow = (dim(timepoint1_dct)[1] / 3) , ncol = dim(timepoint1_dct)[2])

sample_row2 = seq(1, dim(timepoint1_dct)[1], 3)
rownames(mean_timepoint1_dct) = rownames(timepoint1_dct)[sample_row2]
colnames(mean_timepoint1_dct) = colnames(timepoint1_dct)


for (row in 1:dim(mean_timepoint1_dct)[1]){
  
  for (col in 1:dim(mean_timepoint1_dct)[2]){
    
    mean_timepoint1_dct[row, col] = mean(timepoint1_dct[(row*3-2):(row*3), col], na.rm = T)
    
  }
}


############
# table ddCT

# timepoint 2: get dct second timepoint
timepoint2_dct = dct[which((str_sub(rownames(dct), start= -1)) == "2"), ]

# timepoint 3: get dct second timepoint
timepoint3_dct = dct[which((str_sub(rownames(dct), start= -1)) == "3"), ]


ddct_timepoint2 = matrix(NA, nrow = dim(timepoint2_dct)[1], ncol = dim(timepoint2_dct)[2])
rownames(ddct_timepoint2) = rownames(timepoint2_dct)
colnames(ddct_timepoint2) = colnames(timepoint2_dct)

ddct_timepoint3 = matrix(NA, nrow = dim(timepoint3_dct)[1], ncol = dim(timepoint3_dct)[2])
rownames(ddct_timepoint3) = rownames(timepoint3_dct)
colnames(ddct_timepoint3) = colnames(timepoint3_dct)


for (row in 1:dim(mean_timepoint1_dct)[1]){
  
  for (col in 1:dim(mean_timepoint1_dct)[2]){
    
    ddct_timepoint2[(row*3-2):(row*3),col] = timepoint2_dct[(row*3-2):(row*3), col] - mean_timepoint1_dct[row,col]
    ddct_timepoint3[(row*3-2):(row*3),col] = timepoint3_dct[(row*3-2):(row*3), col] - mean_timepoint1_dct[row,col]
    
  }
}


ddct = rbind(ddct_timepoint2, ddct_timepoint3)

path = paste(result_path, sample_name, "_ddct_select_v3.csv", sep = "")
write.csv(ddct, path)


#############
# fold change

fold_change = 2^(-ddct)

fold_change_sci = format(fold_change, scientific = FALSE)

path = paste(result_path, sample_name, "_fold_change_select_v3.csv", sep = "")
write.csv(fold_change_sci, path)


########
# result

#Avg for each sample

avg = matrix(NA, nrow = (dim(fold_change)[1] / 3), ncol = dim(fold_change)[2])
colnames(avg) = colnames(fold_change)

sample_row3 = seq(1, dim(fold_change)[1], 3)
rownames(avg) = rownames(fold_change)[sample_row3]


for (row in 1:dim(avg)[1]){
  
  for (col in 1:dim(avg)[2]){
    
    avg[row, col] = mean(fold_change[(row*3-2):(row*3), col], na.rm = T)
    
  }
}

avg_sci = format(avg, scientific = FALSE)

path = paste(result_path, sample_name, "_avg_select_v3.csv", sep = "")
write.csv(avg_sci, path)

##########################################
# AVG for A, B, C team with timepoint 2, 3

# avg A & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "A"))
avg_a_2 = avg_sci[row, ]

# avg A & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "A"))
avg_a_3 = avg_sci[row, ]

# avg B & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "B"))
avg_b_2 = avg_sci[row, ]

# avg B & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "B"))
avg_b_3 = avg_sci[row, ]

# avg C & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "C"))
avg_c_2 = avg_sci[row, ]

# avg C & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "C"))
avg_c_3 = avg_sci[row, ]


# mean of each team in different timepoint

mean_a_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_a_2)[2]){
  
  mean_a_2[1, col] = mean(as.numeric(avg_a_2[, col]), na.rm = T)
  
}


mean_a_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_a_3)[2]){
  
  mean_a_3[1, col] = mean(as.numeric(avg_a_3[, col]), na.rm = T)
  
}


mean_b_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_b_2)[2]){
  
  mean_b_2[1, col] = mean(as.numeric(avg_b_2[, col]), na.rm = T)
  
}


mean_b_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_b_3)[2]){
  
  mean_b_3[1, col] = mean(as.numeric(avg_b_3[, col]), na.rm = T)
  
}


mean_c_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_c_2)[2]){
  
  mean_c_2[1, col] = mean(as.numeric(avg_c_2[, col]), na.rm = T)
  
}


mean_c_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(avg_c_3)[2]){
  
  mean_c_3[1, col] = mean(as.numeric(avg_c_3[, col]), na.rm = T)
  
}


mean = rbind(mean_a_2, mean_a_3, mean_b_2, mean_b_3, mean_c_2, mean_c_3)
colnames(mean) = colnames(fold_change)
rownames(mean) = c("A_2", "A_3", "B_2", "B_3", "C_2", "C_3")


path = paste(result_path, sample_name, "_mean_v3.csv", sep = "")
write.csv(mean, path)



###############################
# SEM = Std / sqrt(sample_size)

# std A & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "A"))
std_a_2 = avg_sci[row, ]

# std A & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "A"))
std_a_3 = avg_sci[row, ]

# std B & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "B"))
std_b_2 = avg_sci[row, ]

# std B & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "B"))
std_b_3 = avg_sci[row, ]

# std C & timepoint 2: get avg from team A & timepoint 2
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "2"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "C"))
std_c_2 = avg_sci[row, ]

# std C & timepoint 3: get avg from team A & timepoint 3
row = intersect(which((str_sub(rownames(avg_sci), start= -1)) == "3"), which((str_sub(rownames(avg_sci), start= 1, end = 1)) == "C"))
std_c_3 = avg_sci[row, ]

# delete the na row
std_c_3 = std_c_3[-8,]


# sem

sem_a_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_a_2)[2]){
  
  sem_a_2[1, col] = sd(std_a_2[, col]) / sqrt(dim(std_a_2)[1])
  
}


sem_a_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_a_3)[2]){
  
  sem_a_3[1, col] = sd(std_a_3[, col]) / sqrt(dim(std_a_3)[1])
  
}


sem_b_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_b_2)[2]){
  
  sem_b_2[1, col] = sd(std_b_2[, col]) / sqrt(dim(std_b_2)[1])
  
}


sem_b_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_b_3)[2]){
  
  sem_b_3[1, col] = sd(std_b_3[, col]) / sqrt(dim(std_b_3)[1])
  
}


sem_c_2 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_c_2)[2]){
  
  sem_c_2[1, col] = sd(std_c_2[, col]) / sqrt(dim(std_c_2)[1])
  
}


sem_c_3 = matrix(NA, nrow = 1, ncol = dim(fold_change)[2])

for (col in 1:dim(std_c_3)[2]){
  
  sem_c_3[1, col] = sd(std_c_3[, col]) / sqrt(dim(std_c_3)[1])
  
}

  
sem = rbind(sem_a_2, sem_a_3, sem_b_2, sem_b_3, sem_c_2, sem_c_3)
colnames(sem) = colnames(fold_change)
rownames(sem) = c("A_2", "A_3", "B_2", "B_3", "C_2", "C_3")


path = paste(result_path, sample_name, "_sem_v3.csv", sep = "")
write.csv(sem, path)

