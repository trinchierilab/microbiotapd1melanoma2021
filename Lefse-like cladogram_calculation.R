# Create file for LEFSE-like cladogram
# Calculate stats for all levels

# get level names
library(tidyverse)
library(tidyr)
merged_df<- left_join(featuretable, LKT_counts)
merged_df<- left_join(merged_df, LKT_HZ_double_q)
merged_df<- merged_df %>% drop_na()
levels<- colnames(featuretable)[-c(1:3)]
for (i in 1:7){
  
  aggregate(temp$PPM, by=list(Category=temp$Phylum), FUN=sum)
  
  temp<- merged_df[c(levels[i],12)]
  assign()
}

### Stats for cladogram

length(which((subset(merged_df, merged_df$Class=="c__Bacteroidia")$HR) > 1.5))
# Result 35
length(which((subset(merged_df, merged_df$Class=="c__Bacteroidia")$HR) < 0.66666))
# Result 3
obs <- c(35,3)
exp <- c(0.5,0.5)
chisq.test(obs, p=exp)
# Result:
#Chi-squared test for given probabilities
#data:  obs
#X-squared = 26.947, df = 1, p-value = 2.091e-07

merged_df<- subset(merged_df, merged_df$Kingdom== "k__Bacteria")
> View(merged_df)

# template 
cbind(paste(merged_df$Kingdom[2],merged_df$Phylum[2], sep = "."), merged_df$HR[2], merged_df$Second_q[2])

stats<- cbind(merged_df$HR, merged_df$Second_q, merged_df$PPM)

phylum<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum, sep = ".")),stats))
class<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class, sep = ".")),stats))
order<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order, sep = ".")),stats))
family<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family, sep = ".")),stats))
genus<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family,merged_df$Genus, sep = ".")), stats))
species<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family,merged_df$Genus,merged_df$Species, sep = ".")),stats))

colnames(phylum)<- c("Phylum", "HR", "double_q", "PPM")
colnames(class)<- c("Class", "HR", "double_q", "PPM")
colnames(order)<- c("Order", "HR", "double_q", "PPM")
colnames(family)<- c("Family", "HR", "double_q", "PPM")
colnames(genus)<- c("Genus", "HR", "double_q", "PPM")
colnames(species)<- c("Species", "HR", "double_q", "PPM")

phylum.names<- unique(phylum[,1])
class.names<- unique(class[,1])
order.names<- unique(order[,1])
family.names<- unique(family[,1])
genus.names<- unique(genus[,1])



length(which((subset(phylum, phylum$Phylum==(as.character(phylum.names[2])[,2]) > 1.5))))

subset(phylum, phylum$Phylum==(phylum.names[2]))[,2]

length(which(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[2]))[,2])) > 1.5))
length(which(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[2]))[,2])) > 1.5))


# Result 35
length(which((subset(merged_df, merged_df$Class=="c__Bacteroidia")$HR) < 0.66666))
# Result 3
obs <- c(35,3)
exp <- c(0.5,0.5)
chisq.test(obs, p=exp)


tax_levels<- colnames(merged_df)[c(4:10)]
unique(as.character(merged_df[,tax_levels[2]]))


#Sum levels example
aggregate(as.numeric(as.character(phylum$PPM)), by=list(Category=phylum$phylum), FUN=sum)

# metaanalysis of p-values for each taxon
pchisq(105.0389, 70)

# total number of taxa for metaanalysis
n_taxa<- length(phylum$Phylum)

z<- as.numeric(as.character(phylum[c(which(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[5]))[,2])) > 1.5)),][,3]))


# filter(df, fct %in% vc)

count.taxa<- c()
for (i in 1:n_taxa){
  subset.taxa<- subset(phylum, phylum$Phylum==(phylum.names[i]))
  taxa.qount<- cbind(as.character(phylum.names[i]), length(subset.taxa$Phylum))
  count.taxa<- rbind(count.taxa, taxa.qount)
}
count.taxa<- as.data.frame(count.taxa)
count.taxa[,2]<- as.numeric(as.character(count.taxa$V2))
# Select taxa which has more than "w" taxa
# set taxa number cutoff
w<- 3
# subset taxa that pass the cutoff
selected.taxa<- as.character((subset(count.taxa, count.taxa$V2 > (w-1)))[,1])

library(dplyr)
Phylum<- filter(phylum, Phylum %in% selected.taxa)


###########################################





phylum<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum, sep = ".")),stats))
class<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class, sep = ".")),stats))
order<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order, sep = ".")),stats))
family<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family, sep = ".")),stats))
genus<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family,merged_df$Genus, sep = ".")), stats))
species<- as.data.frame(cbind((paste(merged_df$Kingdom,merged_df$Phylum,merged_df$Class,merged_df$Order,merged_df$Family,merged_df$Genus,merged_df$Species, sep = ".")),stats))

colnames(phylum)<- c("Phylum", "HR", "double_q", "PPM")
colnames(class)<- c("Class", "HR", "double_q", "PPM")
colnames(order)<- c("Order", "HR", "double_q", "PPM")
colnames(family)<- c("Family", "HR", "double_q", "PPM")
colnames(genus)<- c("Genus", "HR", "double_q", "PPM")
colnames(species)<- c("Species", "HR", "double_q", "PPM")

phylum.names<- unique(phylum[,1])
class.names<- unique(class[,1])
order.names<- unique(order[,1])
family.names<- unique(family[,1])
genus.names<- unique(genus[,1])


### Phylum calculation
phylum.cladogram.input<- data.frame()

#p-values

n_taxa<- length(unique(phylum$Phylum))

for (i in 1:n_taxa){
  vector.HR<- as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[i]))[,2]))
  protective_vector<- vector.HR[c(which(vector.HR<0.6666))]
  detrimental_vector<- vector.HR[c(which(vector.HR>1.5))]
  protective<- length(protective_vector)
  detrimental<- length(detrimental_vector)
  obs <- c((0.01+protective),(0.01+detrimental))
  exp <- c(0.5,0.5)
  x<- chisq.test(obs, p=exp)
  discordance_p<- x$p.value
  vector.PPM<- as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  if (protective>detrimental){
    p<- as.numeric(as.character(phylum[c(which(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[i]))[,2])) < 0.66666)),][,3])) 
    HR_mean<- mean(protective_vector)
  } else {
    p<- as.numeric(as.character(phylum[c(which(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[i]))[,2])) > 1.5)),][,3]))
    HR_mean<- mean(detrimental_vector)
  }
  ### For responder status in LEFse cladogram
  if (protective>detrimental){
    taxa_status<- "Protective"
  } else {
    taxa_status<- "Detrimental"
  }
  # Calculate p from chidist - Fisher's method from the group that goes same direction
  ### nonsign<- n_taxa - length(p)
  #### nonsign_vector<-  c(rep(0.99,nonsign))
  
  chidist<- sum(-2*(sapply(p, log)))
  chi_p<- 1 - pchisq(chidist, (2*length(p)))
  
  ### Calculate sum of the reads per taxa
  sum.reads.taxa<- sum(as.numeric(as.character(subset(phylum, phylum$Phylum==(phylum.names[i]))[,4])))
  # Extract taxa name into tax_name file
  taxa_name<- unique(subset(phylum, phylum$Phylum==(phylum.names[i]))[1])
  if (chi_p<0.01 && discordance_p<0.05) {
    temp_value<- cbind(taxa_name, log.PPM, taxa_status, as.factor(HR_mean), as.factor(chi_p), chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  } else {
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  
  phylum.cladogram.input<- rbind(phylum.cladogram.input, temp_value)
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
### Class calculation ################################################################################################

class.cladogram.input<- data.frame()

#p-values

n_taxa<- length(unique(class$Class))

for (i in 1:n_taxa){
  vector.HR<- as.numeric(as.character(subset(class, class$Class==(class.names[i]))[,2]))
  protective_vector<- vector.HR[c(which(vector.HR<0.6666))]
  detrimental_vector<- vector.HR[c(which(vector.HR>1.5))]
  protective<- length(protective_vector)
  detrimental<- length(detrimental_vector)
  obs <- c((0.01+protective),(0.01+detrimental))
  exp <- c(0.5,0.5)
  x<- chisq.test(obs, p=exp)
  discordance_p<- x$p.value
  vector.PPM<- as.numeric(as.character(subset(class, class$Class==(class.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  if (protective>detrimental){
    p<- as.numeric(as.character(class[c(which(as.numeric(as.character(subset(class, class$Class==(class.names[i]))[,2])) < 0.66666)),][,3])) 
    HR_mean<- mean(protective_vector)
  } else {
    p<- as.numeric(as.character(class[c(which(as.numeric(as.character(subset(class, class$Class==(class.names[i]))[,2])) > 1.5)),][,3]))
    HR_mean<- mean(detrimental_vector)
  }
  ### For responder status in LEFse cladogram
  if (protective>detrimental){
    taxa_status<- "Protective"
  } else {
    taxa_status<- "Detrimental"
  }
  # Calculate p from chidist - Fisher's method from the group that goes same direction
  ### nonsign<- n_taxa - length(p)
  #### nonsign_vector<-  c(rep(0.99,nonsign))
  
  chidist<- sum(-2*(sapply(p, log)))
  chi_p<- 1 - pchisq(chidist, (2*length(p)))
  
  ### Calculate sum of the reads per taxa
  sum.reads.taxa<- sum(as.numeric(as.character(subset(class, class$Class==(class.names[i]))[,4])))
  # Extract taxa name into tax_name file
  taxa_name<- unique(subset(class, class$Class==(class.names[i]))[1])
  if (chi_p<0.01 && discordance_p<0.05) {
    temp_value<- cbind(taxa_name, log.PPM, taxa_status, as.factor(HR_mean), as.factor(chi_p), chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  } else {
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  
  class.cladogram.input<- rbind(class.cladogram.input, temp_value)
}


######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
### Order calculation ################################################################################################

order.cladogram.input<- data.frame()

#p-values

n_taxa<- length(unique(order$Order))

for (i in 1:n_taxa){
  vector.HR<- as.numeric(as.character(subset(order, order$Order==(order.names[i]))[,2]))
  protective_vector<- vector.HR[c(which(vector.HR<0.6666))]
  detrimental_vector<- vector.HR[c(which(vector.HR>1.5))]
  protective<- length(protective_vector)
  detrimental<- length(detrimental_vector)
  obs <- c((0.01+protective),(0.01+detrimental))
  exp <- c(0.5,0.5)
  x<- chisq.test(obs, p=exp)
  discordance_p<- x$p.value
  vector.PPM<- as.numeric(as.character(subset(order, order$Order==(order.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  if (protective>detrimental){
    p<- as.numeric(as.character(order[c(which(as.numeric(as.character(subset(order, order$Order==(order.names[i]))[,2])) < 0.66666)),][,3])) 
    HR_mean<- mean(protective_vector)
  } else {
    p<- as.numeric(as.character(order[c(which(as.numeric(as.character(subset(order, order$Order==(order.names[i]))[,2])) > 1.5)),][,3]))
    HR_mean<- mean(detrimental_vector)
  }
  ### For responder status in LEFse cladogram
  if (protective>detrimental){
    taxa_status<- "Protective"
  } else {
    taxa_status<- "Detrimental"
  }
  # Calculate p from chidist - Fisher's method from the group that goes same direction
  ### nonsign<- n_taxa - length(p)
  #### nonsign_vector<-  c(rep(0.99,nonsign))
  
  chidist<- sum(-2*(sapply(p, log)))
  chi_p<- 1 - pchisq(chidist, (2*length(p)))
  
  ### Calculate sum of the reads per taxa
  sum.reads.taxa<- sum(as.numeric(as.character(subset(order, order$Order==(order.names[i]))[,4])))
  # Extract taxa name into tax_name file
  taxa_name<- unique(subset(order, order$Order==(order.names[i]))[1])
  if (chi_p<0.01 && discordance_p<0.05) {
    temp_value<- cbind(taxa_name, log.PPM, taxa_status, as.factor(HR_mean), as.factor(chi_p), chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  } else {
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  
  order.cladogram.input<- rbind(order.cladogram.input, temp_value)
}




######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
### Family calculation ###############################################################################################

family.cladogram.input<- data.frame()

#p-values

n_taxa<- length(unique(family$Family))

for (i in 1:n_taxa){
  vector.HR<- as.numeric(as.character(subset(family, family$Family==(family.names[i]))[,2]))
  protective_vector<- vector.HR[c(which(vector.HR<0.6666))]
  detrimental_vector<- vector.HR[c(which(vector.HR>1.5))]
  protective<- length(protective_vector)
  detrimental<- length(detrimental_vector)
  obs <- c((0.01+protective),(0.01+detrimental))
  exp <- c(0.5,0.5)
  x<- chisq.test(obs, p=exp)
  discordance_p<- x$p.value
  vector.PPM<- as.numeric(as.character(subset(family, family$Family==(family.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  if (protective>detrimental){
    p<- as.numeric(as.character(family[c(which(as.numeric(as.character(subset(family, family$Family==(family.names[i]))[,2])) < 0.66666)),][,3])) 
    HR_mean<- mean(protective_vector)
  } else {
    p<- as.numeric(as.character(family[c(which(as.numeric(as.character(subset(family, family$Family==(family.names[i]))[,2])) > 1.5)),][,3]))
    HR_mean<- mean(detrimental_vector)
  }
  ### For responder status in LEFse cladogram
  if (protective>detrimental){
    taxa_status<- "Protective"
  } else {
    taxa_status<- "Detrimental"
  }
  # Calculate p from chidist - Fisher's method from the group that goes same direction
  ### nonsign<- n_taxa - length(p)
  #### nonsign_vector<-  c(rep(0.99,nonsign))
  
  chidist<- sum(-2*(sapply(p, log)))
  chi_p<- 1 - pchisq(chidist, (2*length(p)))
  
  ### Calculate sum of the reads per taxa
  sum.reads.taxa<- sum(as.numeric(as.character(subset(family, family$Family==(family.names[i]))[,4])))
  # Extract taxa name into tax_name file
  taxa_name<- unique(subset(family, family$Family==(family.names[i]))[1])
  if (chi_p<0.01 && discordance_p<0.05) {
    temp_value<- cbind(taxa_name, log.PPM, taxa_status, as.factor(HR_mean), as.factor(chi_p), chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  } else {
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  
  family.cladogram.input<- rbind(family.cladogram.input, temp_value)
}





######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
### Genus calculation ###############################################################################################



genus.cladogram.input<- data.frame()

#p-values

n_taxa<- length(unique(genus$Genus))

for (i in 1:n_taxa){
  vector.HR<- as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,2]))
  protective_vector<- vector.HR[c(which(vector.HR<0.6666))]
  detrimental_vector<- vector.HR[c(which(vector.HR>1.5))]
  protective<- length(protective_vector)
  detrimental<- length(detrimental_vector)
  taxa_name<- unique(subset(genus, genus$Genus==(genus.names[i]))[1])
  vector.PPM<- as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  
  if (sum(protective+detrimental)>2){
  
  obs <- c((0.01+protective),(0.01+detrimental))
  exp <- c(0.5,0.5)
  x<- chisq.test(obs, p=exp)
  discordance_p<- x$p.value
  vector.PPM<- as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,4]))
  log.PPM<- log(sum(vector.PPM),2)
  
  if (protective>detrimental){
    p<- as.numeric(as.character(genus[c(which(as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,2])) < 0.66666)),][,3])) 
    HR_mean<- mean(protective_vector)
  } else {
    p<- as.numeric(as.character(genus[c(which(as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,2])) > 1.5)),][,3]))
    HR_mean<- mean(detrimental_vector)
  }
  ### For responder status in LEFse cladogram
  if (protective>detrimental){
    taxa_status<- "Protective"
  } else {
    taxa_status<- "Detrimental"
  }
  # Calculate p from chidist - Fisher's method from the group that goes same direction
  ### nonsign<- n_taxa - length(p)
  #### nonsign_vector<-  c(rep(0.99,nonsign))
  
  chidist<- sum(-2*(sapply(p, log)))
  chi_p<- 1 - pchisq(chidist, (2*length(p)))
  
  ### Calculate sum of the reads per taxa
  sum.reads.taxa<- sum(as.numeric(as.character(subset(genus, genus$Genus==(genus.names[i]))[,4])))
  # Extract taxa name into tax_name file
  taxa_name<- unique(subset(genus, genus$Genus==(genus.names[i]))[1])
  if (chi_p<0.01 && discordance_p<0.05) {
    temp_value<- cbind(taxa_name, log.PPM, taxa_status, as.factor(HR_mean), as.factor(chi_p), chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  } else {
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  }else{
    chi_p<- 1
    discordance_p<- 1
    temp_value<- cbind(taxa_name, log.PPM, "", "", "-", chi_p, discordance_p)
    colnames(temp_value)<- c("taxa", "log.PPM", "taxa_status", "HR_mean", "chi_p", "chi_p_val", "discordance_p")
  }
  
  genus.cladogram.input<- rbind(genus.cladogram.input, temp_value)
  }
  





