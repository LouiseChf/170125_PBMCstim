#09.01.2025 - Louise CHORFI 
# Graph exploratoires, expression de CD49a par les lymphos CD8 PBL et TIl après stim ou non,
# addition de nouvelles conditions de stim : aCD3, IL-2,IL-4, IL-15, TGFb

#### set up env. and data ####

#clean environment and load packages:
rm(list=ls())
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd(dir = "C:/Users/l_chorfi/Documents/TRM - CD49a/stims/")
exportpath <- "C:/Users/l_chorfi/Documents/TRM - CD49a/stims/R_export"

#import data from excel
TIL <- read.table(file = "C:/Users/l_chorfi/Documents/TRM - CD49a/stims/TIL_080125_allconditions.csv",
                  header = TRUE, sep =";")
PBMC <- read.table(file = "C:/Users/l_chorfi/Documents/TRM - CD49a/stims/PBMC_080125_newpanel2.csv",
                 header = TRUE, sep =";")

str(TIL)
str(PBMC) # !! small values in 10E3 make the observation characters

#clean data
#PBMC <- separate(PBMC, X, into= c(NA, "celltype", "sample","treatment", NA), sep = "_",)
TIL <- separate(TIL, X, into= c("celltype",NA, NA, "sample", "treatment", NA), sep = "_",)

TIL$sample <- gsub("new panel 2-", "", TIL$sample)
TIL$sample <- gsub("new panel-TOX-APC-", "", TIL$sample)
TIL$treatment <- gsub("NS", "NT", TIL$treatment)
TIL$treatment <- gsub("CD3", "aCD3", TIL$treatment)
TIL$treatment <- gsub("S", "aCD3", TIL$treatment)
TIL$treatment <- gsub("aa", "a", TIL$treatment)

PBMC$sample_ID <- gsub("new panel 2-", "", PBMC$sample_ID)
PBMC$sample_ID <- gsub("new panel-TOX-APC-", "", PBMC$sample_ID)
PBMC$treatment <- gsub("NS", "NT", PBMC$treatment)
PBMC$treatment <- gsub("CD3", "aCD3", PBMC$treatment)
PBMC$treatment <- gsub("S", "aCD3", PBMC$treatment)
PBMC$treatment <- gsub("NT+","", PBMC$treatment, fixed = TRUE)
PBMC$treatment <- gsub("^IL4$","NT+IL4", PBMC$treatment)
PBMC$treatment <- gsub("^IL12$","NT+IL12", PBMC$treatment)
PBMC$treatment <- gsub("^IL15$","NT+IL15", PBMC$treatment)
PBMC$treatment <- gsub("^TGFb$","NT+TGFb", PBMC$treatment)
PBMC$freq_CD103pos_CD49apos.in.CD4 <- as.numeric(PBMC$freq_CD103pos_CD49apos.in.CD4) 
PBMC$freq_CD103pos_CD49apos.in.CD4 [is.na(PBMC$freq_CD103pos_CD49apos.in.CD4)] <- 0

colnames(TIL)
colnames(PBMC)
aov_summ <- data.frame()
tukey_df <- data.frame()

#### TIL : all plots in one script 8) ####
data <- TIL
for (i in 4:15){
  i_name <- colnames(data[i])
  
  #paired data
  temppaired <- data %>%
    select(celltype, sample, treatment, i)
  temppaired <-  pivot_wider(temppaired, names_from = treatment, values_from = 4)
  assign( x= paste0("paired_", i_name), value = temppaired)
  
  
  #graph
  binwidth <- diff(range(data[i])) / 40 
  set.seed(2)
  plot <- ggplot(data, aes(x = treatment, y = data[i][[1]] , fill = treatment))+
    geom_boxplot(alpha = 0.3) + 
    geom_dotplot(aes(), binaxis = "y", stackdir = "center", binwidth = binwidth)+
    scale_x_discrete(limits = c("NT", "aCD3", "TGFb","aCD3+TGFb"))+
    theme_classic() + 
    labs(x="treatment", y= i_name, 
         title = paste0( "TIL stimulations ", i_name))
  #  geom_segment( data = paired_CD8freq, aes(x = 1, xend = 2, y = NT , yend = aCD3),
  #                color = "gray50", alpha = 0.5)
  
  
  plot
  ggsave(filename = paste0(i_name, ".png"), 
         path = exportpath)

  #stat
  temp_anova <- aov(data[i][[1]] ~ treatment, data = data)
  temp_aovsumm <- data.frame( param = i_name,(summary(temp_anova)[[1]]))
  aov_summ <- rbind(aov_summ,temp_aovsumm)
  
  if(temp_aovsumm$Pr..F.[1] < 0.05) {
    tukey <- TukeyHSD(temp_anova)
    tukey_temp <- data.frame( param = i_name, tukey[[1]])
    tukey_temp$comp <- rownames(tukey_temp) 
    tukey_df <- rbind(tukey_df, tukey_temp)
  }  else {
    
  }
  
  #display stat on plot and save plot 
#  signif <-  tukey_temp$p.adj[] < 0.05
 
# for (i in signif){
#   if( i == TRUE ){
#    plot <- plot + geom_text( data = tukey_temp$p.adj[] ,
#                                    aes(x = 1.5, y = 110, label = paste0("p = ", signif(pvalue, 3))),
#                                   inherit.aes = FALSE,  
#                                  size = 3) +
#        geom_segment( aes(x = 1, xend = 2, y = 105 , yend = 105),
#                      color = "black", alpha = 1)
#       } else {
#    plot <- plot 

#      }
#  }

    
  ggsave(filename = paste0(i_name, ".png"), 
         path = exportpath)
  
  #clean environment
  rm(tukey_temp)
  rm( list = paste0("paired_", i_name) )
  
}

write.table(aov_summ, "/R_export/TIL_aov_summ.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tukey_df, "/R_export/TIL_Tukey_tests.csv", sep = "\t", row.names = FALSE, quote = FALSE)

#### PBMC : all plots in one script !####
data <- PBMC
data$treatment <- factor(data$treatment)

## I- stats 

for (i in 4:ncol(data)) {
  i_name <- colnames(data[i])
  treatmentcond <- unique(data$treatment )
  
  #stat
  temp_anova <- aov(data[i][[1]] ~ treatment, data = data)
  temp_aovsumm <- data.frame( param = i_name,(summary(temp_anova)[[1]]))
  aov_summ <- rbind(aov_summ,temp_aovsumm)
  
  if(temp_aovsumm$Pr..F.[1] < 0.05) {
    tukey <- TukeyHSD(temp_anova)
    tukey_temp <- data.frame( param = i_name, tukey[[1]])
    tukey_temp$comp <- rownames(tukey_temp) 
    tukey_df <- rbind(tukey_df, tukey_temp)
    rm(tukey_temp)
  }  else {
    tukey_temp <- data.frame( param = i_name, diff = NA,
                              lwr = NA, upr = NA, p.adj = "Anova NS",
                              comp = NA)
    tukey_df <- rbind(tukey_df, tukey_temp)
    rm(tukey_temp)
  }
  rm(temp_anova, temp_aovsumm)
}

write.table(aov_summ, "R_export/PBMC_aov_summ.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tukey_df, "R_export/PBMC_Tukey_tests.csv", sep = "\t", row.names = FALSE, quote = FALSE)


## II - Basic graph 
for (i in 4:ncol(data)){
  i_name <- colnames(data[i])
  treatmentcond <- unique(data$treatment)
  
  temp_anova <- compare_means(data[i][[1]] ~ treatment, data = data, method = "anova")
  
  #paired data
  temppaired <- data %>%
    select(celltype, sample_ID, treatment, i)
  temppaired <-  pivot_wider(temppaired, names_from = treatment, values_from = 4)
  assign( x= paste0("paired_", i_name), value = temppaired)
  
  #graph
  binwidth <- diff(range(data[i])) / 50 
  set.seed(2)
  plot <- ggplot(data, aes(x = treatment, y = data[i][[1]] , fill = treatment))+
    geom_boxplot(alpha = 0.3 ) + 
    stat_compare_means( method = "anova") +
   geom_dotplot(aes(), binaxis = "y", stackdir = "center", binwidth = binwidth)+
   theme_classic() + 
    labs(x="treatment", y= i_name, 
         title = paste0( "PBMC stimulations ", i_name))
 


  
  plot
  ggsave(filename = paste0(i_name, ".png"), 
         path = exportpath)
  
  #clean environment
  rm( list = paste0("paired_", i_name) )
  
}

## III- Advanced plot : display stat 

signif <-  tukey_df$p.adj[] < 0.05

for (i in signif){
  for (j in seq_along(signif)){
    if( i == TRUE ){
      plot <- plot + geom_text( data = tukey_temp ,
                                aes(x = j+0.5, y = 110, label = paste0("p = ", 
                                                                       signif(p.adj[j], 3))),
                                inherit.aes = FALSE,  
                                size = 3) +
        geom_segment( aes(x = j, xend = j+1, y = 105 , yend = 105),
                      color = "black", alpha = 1)
      plot
    } else {
      plot <- plot + geom_text(aes(x = j+0.5 , y = 110, label = paste0("NS")),
                               inherit.aes = FALSE,  
                               size = 3) +
        geom_segment( aes(x = j, xend = j+1, y = 105 , yend = 105),
                      color = "black", alpha = 1)
      plot
      
      ggsave(filename = paste0(i_name, "_tuckey.png"), 
             path = exportpath)
    }
  }
}

