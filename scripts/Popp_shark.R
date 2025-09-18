########## Brian, Rita, Carl Meyer shark and reefs ##################

#load packages and data
# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load all the packages you are missing!
pacman::p_load("knitr", "ggord", "MASS", "tidyverse", "plyr", "dplyr", "gridExtra", 
               "cowplot", "ggplot2", "R2jags", "MixSIAR", "remotes")

# load data
df<-read.csv("data/Popp_Trevl_Gshark.csv")

# parse data
raw.df<- df %>% dplyr::select(SampleID, Group, Thr, Val, Leu, Phe, Lys)
raw.df<-na.omit(raw.df)

norm.df<- df %>% dplyr::select(SampleID, Group, Thrn, Valn, Leun, Phen, Lysn)
norm.df<-na.omit(norm.df)
########## 
########## LDA

raw.df$SampleID<-as.factor(raw.df$SampleID)
raw.df$Group<-as.factor(raw.df$Group)

# remove all samples, keep source
prod <- raw.df[!(raw.df$Group == "GT" |raw.df$Group == "GS"),]
prod$Group<-droplevels(prod$Group) # drop sample levels

#Run an LDA with a jackknifing model fit to look at error rate
#'CV = TRUE' makes the LDA run a jackknifed - leave one out - model fit
All.lda <- lda(Group ~Thr + Val + Leu + Phe + Lys, data = prod, CV = TRUE)

#Create a table which compares the classification from the LDA model to the actual classification
All.reclass <- table(prod$Group, All.lda$class) #100% success

#Total percent of samples correctly classified is the sum of the diagonal of this table
sum(diag(prop.table(All.reclass))) #100% !

mean(All.lda$class == prod$Group) # mean accuracy

#Percent of each producer Group correctly reclassified: 100% for all groups
diag(prop.table(All.reclass, 1))

#Create a training LDA function from the library data
#Note - you can't use the 'All.lda' object above because the 'CV = TRUE' command was used to create it, and for some reason this won't work with the predict() function
All.train <- lda(Group~ Thr + Val + Leu + Phe + Lys, data = prod)
All.train

# save as txt file
sink("output/Popp_LDA_Train_Info.txt")
print(All.train)
sink()

#plot LD1 and LD2 to show the sources
prod$Group <-factor(prod$Group, levels=c("PelagicPlankton", "ReefPlankton", "Coral"))

Source.LDA <- ggord(All.train, prod$Group, arrow=0,
                    txt = NULL, veclsz = 0, vectyp="blank", # remove vectors and arrows
                    cols = c("dodgerblue", "springgreen4","coral"), grp_title = "Producers",
                    ellipse_pro = 0.95, xlim = c(-10,14), ylim = c(-8,6), alpha=0.5) + theme_classic() + theme(legend.position = 'top')

#view plot, not used in manuscript
Source.LDA
dev.copy(pdf, "figures/Popp_LDA.pdf", height=6, width=8)
dev.off()

#Write a csv file for the coefficients of LD1/2
lda.info <- as.data.frame(All.train$scaling)
write.csv(lda.info, "output/Popp_LDA_loadings.csv", row.names=FALSE)

#Create a data frame with these LDA coordinates
AllProdPredict <- data.frame(SampleID = prod$SampleID, Group = prod$Group, LD.class = "Source", predict(All.train)$x)

#Write a csv file for the LDA coordinates
write.csv(AllProdPredict, "output/Popp_LDA_coords.csv", row.names = FALSE)


######## ######## ######## ######## #####
######## Classify Trevaly and Sharks ####
######## ######## ######## ######## #####
# remove all samples, keep source
shk.df <- raw.df[(raw.df$Group == "GT" |raw.df$Group == "GS"),]
shk.df$Group<-droplevels(shk.df$Group) # drop sample levels

shk.predict <- predict(object = All.train, newdata = shk.df)

shk.predict.data <- data.frame(SampleID = shk.df$SampleID, Group = shk.df$Group, LD.class = shk.predict$class, shk.predict$x)

All.reclass <- table(shk.df$Group, shk.predict$class) 
#  GS: 1 as PelagicPlank, 3 as ReefPlankton
#  GT: 7 as PelagicPlank, 6 as ReefPlankton

###### ###### write and compile data

# write a csv file for the zooplankton LDA info
write.csv(shk.predict.data, "output/Popp_LDA_info.csv", row.names=FALSE)

# combine the shark data and source LDA #
LDA.df<-rbind(shk.predict.data, AllProdPredict)

# combine with raw and mean-normalized data and rearrange
raw.norm.LD.data <- merge(df, LDA.df, by = "SampleID", all.x = TRUE)

# remove extra column and rename the other
raw.norm.LD.data<- raw.norm.LD.data %>% 
  dplyr::select(-Group.y)
raw.norm.LD.data <- raw.norm.LD.data %>% 
  dplyr::rename("Group" = "Group.x")

###### ###### ######
# export all data
write.csv(raw.norm.LD.data, "output/Popp_LDA.data.compiled.csv")
###### ###### ######


#### plot it
raw.norm.LD.data$Group<-factor(raw.norm.LD.data$Group, 
                      levels=c("PelagicPlankton", "ReefPlankton", "Coral", 
                               "GT", "GS"))

lda.5.colors <- c("dodgerblue", "springgreen4","coral", "gold", "darkgoldenrod")


LDA.scat.Group <- ggplot(raw.norm.LD.data, aes(x = LD1, y = LD2))+
  geom_point(size = 3.5, color="black",
             aes(fill = Group, color = Group, shape = Group), alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.90, linewidth = 0.5, aes(lty = Group, color = Group))+
  scale_color_manual(values = lda.5.colors)+
  scale_fill_manual(values = lda.5.colors)+
  scale_shape_manual(values = c(24, 24, 24, 21, 21))+
  scale_linetype_manual(values = c(2,2,2,0,0))+
  xlab("LD1 (79.1%)") + ylab("LD2 (20.9%)") +
  ggtitle("LDA: Measured d13C-EAA")+
  theme_classic() + guides(lty = "none")

LDA.scat.Group
dev.copy(pdf, "figures/Popp_LDA_sorc_shark.pdf", width=7, height=6)
dev.off()


################################################
## MIXSIAR #####################################
################################################

# create a new data frame with desired columns for the zooplankton consumers (mixtures)
Mix.df.EAA <- raw.norm.LD.data %>%
  dplyr::select("SampleID", "Group", "Thrn", "Valn", "Leun", "Phen", "Lysn")

# remove sources
Mix.df.EAA.shk <- Mix.df.EAA[!(Mix.df.EAA$Group == "ReefPlankton" | Mix.df.EAA$Group == "PelagicPlankton" |
                                  Mix.df.EAA$Group == "Coral"),]
Mix.df.EAA.shk$Group <- droplevels(Mix.df.EAA.shk$Group)
Mix.df.EAA.shk<-na.omit(Mix.df.EAA.shk)

# write and export a csv of the zooplankton data frame
write.csv(Mix.df.EAA.shk, "output/Popp_Mix.df.EAAn.shark.csv", row.names = FALSE)

###########
# SOURCES #
###########

# NORMALIZED DATA
# Reef and pelagic plankton, coral

# create a new data frame with desired columns for the primary producers (sources)
Source.df.EAA <-  Mix.df.EAA[(Mix.df.EAA$Group == "ReefPlankton" | Mix.df.EAA$Group == "PelagicPlankton" |
                               Mix.df.EAA$Group == "Coral"),]

Source.df.EAA <- Source.df.EAA %>%
  dplyr::select("Group", "Thrn", "Valn", "Leun", "Phen", "Lysn")

##### source mean and SD
### Thr mean and SD
# calculate mean Thr values for the sources
mean.Thr.df.EAA <- aggregate(Thrn ~ Group, data = Source.df.EAA, FUN = mean)
colnames(mean.Thr.df.EAA) <- c("Group", "Mean.Thrn")

# calculate standard deviations of Thr values for the sources
sd.Thr.df.EAA <- aggregate(Thrn ~ Group, data = Source.df.EAA, FUN = sd)
colnames(sd.Thr.df.EAA) <- c("Group", "SD.Thrn")

### Val mean and SD
# calculate mean Val values for the sources
mean.Val.df.EAA <- aggregate(Valn ~ Group, data = Source.df.EAA, FUN = mean)
colnames(mean.Val.df.EAA) <- c("Group", "Mean.Valn")

# calculate standard deviations of Val values for the sources
sd.Val.df.EAA <- aggregate(Valn ~ Group, data = Source.df.EAA, FUN = sd)
colnames(sd.Val.df.EAA) <- c("Group", "SD.Valn")

### Leu mean and SD
# calculate mean Leu values for the sources
mean.Leu.df.EAA <- aggregate(Leun ~ Group, data = Source.df.EAA, FUN = mean)
colnames(mean.Leu.df.EAA) <- c("Group", "Mean.Leun")

# calculate standard deviations of Leu values for the sources
sd.Leu.df.EAA <- aggregate(Leun ~ Group, data = Source.df.EAA, FUN = sd)
colnames(sd.Leu.df.EAA) <- c("Group", "SD.Leun")

### Phe mean and SD
# calculate mean Phe values for the sources
mean.Phe.df.EAA <- aggregate(Phen ~ Group, data = Source.df.EAA, FUN = mean)
colnames(mean.Phe.df.EAA) <- c("Group", "Mean.Phen")

# calculate standard deviations of Phe values for the sources
sd.Phe.df.EAA <- aggregate(Phen ~ Group, data = Source.df.EAA, FUN = sd)
colnames(sd.Phe.df.EAA) <- c("Group", "SD.Phen")

### Lys mean and SD
# calculate mean Phe values for the sources
mean.Lys.df.EAA <- aggregate(Lysn ~ Group, data = Source.df.EAA, FUN = mean)
colnames(mean.Lys.df.EAA) <- c("Group", "Mean.Lysn")

# calculate standard deviations of Phe values for the sources
sd.Lys.df.EAA <- aggregate(Lysn ~ Group, data = Source.df.EAA, FUN = sd)
colnames(sd.Lys.df.EAA) <- c("Group", "SD.Lysn")

#### sample size
# calculate the number of samples for each of the sources
df.size.EAA <- aggregate(Lysn ~ Group, data = Source.df.EAA, FUN = length)
colnames(df.size.EAA) <- c("Group", "Lys.size")

# bind the values calculated above into one data frame
source.agg.df.EAA <- cbind(df.size.EAA, 
                           mean.Thr.df.EAA[2], sd.Thr.df.EAA[2], 
                           mean.Val.df.EAA[2], sd.Val.df.EAA[2], 
                           mean.Leu.df.EAA[2], sd.Leu.df.EAA[2], 
                           mean.Phe.df.EAA[2], sd.Phe.df.EAA[2], 
                           mean.Lys.df.EAA[2], sd.Lys.df.EAA[2])
colnames(source.agg.df.EAA) <- c("", "n", "MeanThrn", "SDThrn", "MeanValn", "SDValn", 
                                 "MeanLeun", "SDLeun", "MeanPhen", "SDPhen", "MeanLysn", "SDLysn")

# write and export a csv of the source data frame
write.csv(source.agg.df.EAA, "output/EAAn_Mix/Popp_source.agg.df.EAAn.csv", row.names = FALSE, na = "")


# load zooplankton (consumer/mixture) data and assign factors
# SampleID is a fixed factor
cons.EAAn <- load_mix_data(filename = "output/EAAn_Mix/Popp_Mix.df.EAAn.shark.csv",
                          iso_names = c("Thrn", "Valn", "Leun", "Phen", "Lysn"),
                          factors = "SampleID",
                          fac_random = FALSE,
                          fac_nested = FALSE,
                          cont_effects = NULL)

# load source data
source.EAAn <- load_source_data(filename = "output/EAAn_Mix/Popp_source.agg.df.EAAn.csv",
                               source_factors = NULL,
                               conc_dep = FALSE,
                               data_type = "means",
                               cons.EAAn)

#Load TDF data
discr.EAAn <- load_discr_data(filename = "output/EAAn_Mix/MixSIAR_EAAn_TDFs.csv", cons.EAAn)


################## ################## ################# 
################## RUN MIXSIAR MODEL BY SAMPLEID #######

# define the model and error structure and write JAGS model file
# PROCESS ERROR ONLY - HAVE TO DO THIS WHEN RUNNING MIXTURE POINTS INDIVIDUALLY
model_filename.EAA <- "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_Model_SampleID_fixed.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename.EAA, resid_err, process_err, cons.EAAn, source.EAAn)

################################################
############# run the JAGS model

set.seed(12)
zoop.jags.mod.EAAn <- run_model(run = "normal", cons.EAAn, source.EAAn, discr.EAAn, model_filename.EAA, alpha.prior = 1, resid_err, process_err)

#save RDS object of JAGS model
saveRDS(zoop.jags.mod.EAAn, file = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_Model_SampleID_rjags.rds")

# code below will read back in the RDS object
zoop.jags.mod.EAAn <- readRDS(file = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_Model_SampleID_rjags.rds")

# process diagnostics, summary stats, and posterior plots
output_options_EAAn <- list(
  summary_save = TRUE,                   
  summary_name = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_summary_stats.EAAn",
  sup_post = FALSE,                       
  plot_post_save_pdf = FALSE,              
  plot_post_name = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_posterior_density.EAAn",   
  sup_pairs = TRUE,                       
  plot_pairs_save_pdf = FALSE,             
  plot_pairs_name = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_pairs_plot.EAAn",         
  sup_xy = FALSE,                          
  plot_xy_save_pdf = TRUE,                
  plot_xy_name = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_xy_plot.EAAn",             
  gelman = TRUE,                         
  heidel = FALSE,                         
  geweke = TRUE,                         
  diag_save = TRUE,                       
  diag_name = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_diagnostics.EAAn",
  indiv_effect = FALSE,
  plot_post_save_png = FALSE,            
  plot_pairs_save_png = FALSE,            
  plot_xy_save_png = FALSE,
  graphics.off())

# set max to get all data
options(max.print = 3000)
out.JAGS <- output_JAGS(zoop.jags.mod.EAAn, cons.EAAn, source.EAAn, output_options_EAAn)

# get the output stats for the model
# output_stats(zoop.jags.mod.EAnA, cons.EAAn, source.EAAn, output_options_EAAn)

######### a posteriori aggregation of sources
# combine benthic algae and POM into one autochthonous source
shark.in.out.combined.EAAn <- combine_sources(zoop.jags.mod.EAAn, cons.EAAn, source.EAAn, alpha.prior = 1,
                                          groups = list(Inshore = c("ReefPlankton", "Coral"),
                                                        Offshore = "PelagicPlankton"))

# print and save the interval plot
plot_intervals(shark.in.out.combined.EAAn, toplot = "fac1")
dev.copy(pdf, "output/MixSIAR_models/EAAn_SampleID_model/combines_output_intervals.EAAn.pdf", width=10, height=10)
dev.off()

# summary statistics for a posteriori aggregation of sources (autochthonous vs allochthonous)
summary_stat(shark.in.out.combined.EAAn, toprint = "fac1", meanSD = TRUE, savetxt = TRUE,
             quantiles = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
             filename = "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_summary_stats_a_posteriori.EAAn")


################################################
### CLEAN MIXSIAR MODEL BY SAMPLEID RESULTS ###
################################################

SAMP.mixSIAR.EAAn.output.3 <- read.table(
  "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_summary_stats.EAAn.txt",
  sep = '\t', header = FALSE, skip = 6)

SAMP.mixSIAR.EAAn.output.3 <- as.data.frame(SAMP.mixSIAR.EAAn.output.3)
colnames(SAMP.mixSIAR.EAAn.output.3) <- "full.data"

# remove all extra spaces using "squish", make them 1 space only
library(stringr)
SAMP.mixSIAR.EAAn.output.3$full.data <- str_squish(SAMP.mixSIAR.EAAn.output.3$full.data)

# rename columns
library(plyr)
samp.mix.EAAn.out.3 <- separate(SAMP.mixSIAR.EAAn.output.3, full.data,  
                               into = c("factor", "Mean", "SD", "2.5.perc", "5.perc", "25.perc", "50.perc",
                                        "75.perc", "95.perc", "97.5.perc"), sep=" ")

# parse the factors here into 3 columns
new.cols.EAAn.3 <- stringr::str_split_fixed(samp.mix.EAAn.out.3$factor, "\\.", 3) %>%
  as.data.frame() %>%
  setNames(c("proportion", "SampleID", "Source"))

new.cols.EAAn.3 <- new.cols.EAAn.3 %>%
  dplyr::select("SampleID", "Source")

# combine the new, separated columns with the data, verify the levels match, then remove "factor"
Samp.out.cleaned.EAAn.3 <- cbind(new.cols.EAAn.3, samp.mix.EAAn.out.3[1:10])
Samp.out.cleaned.EAAn.3 <- Samp.out.cleaned.EAAn.3 %>%
  dplyr::select(-factor)

Samp.out.cleaned.EAAn.3<-Samp.out.cleaned.EAAn.3[-1,]

write.csv(Samp.out.cleaned.EAAn.3, "output/EAAn_Mix/MixSIAR_SampleID.df.EAAn.3sources.csv", row.names = FALSE)

################## ################## ################## 
################## ready to plot ################## 

# bring in some metadata to make sense of it
Mix.samp.EAAn.3 <- merge.data.frame(raw.norm.LD.data, Samp.out.cleaned.EAAn.3, by="SampleID")
Mix.samp.EAAn.3$Source <- factor(Mix.samp.EAAn.3$Source, 
                                 levels = c("PelagicPlankton", "ReefPlankton", "Coral"))
Mix.samp.EAAn.3$Mean<-as.numeric(Mix.samp.EAAn.3$Mean) 

# plot it: some general formatting here for the plot...
Fig.formatting <- 
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.key.size = unit(0.6, "cm"))


## plot the 3 sources and shark consumers
Mix.samp.plot.EAAn.3 <- ggplot(data = Mix.samp.EAAn.3, aes(x = Group, y = Mean, fill = Source)) +
  geom_boxplot(alpha = 0.9, color = "black")+
  geom_point(aes(fill = Source), pch = 21, alpha = 0.5, color = "black", 
             position = position_dodge(0.75))+
  scale_fill_manual(values = c("dodgerblue", "springgreen4","coral")) +
  ylab("Proportional Contribution") +
  xlab(NULL) +
  ylim(0,1) +
  guides(fill = guide_legend(title = expression(paste(AA[ESS], " ", Source)))) +
  theme_classic() + Fig.formatting + guides(color = "none")

#print and save
Mix.samp.plot.EAAn.3
dev.copy(pdf, "figures/Popp_Fig.mixsiar.3sour.EAAn.pdf", width = 6, height = 7)
dev.off()


################################ #################### ##########
############ a posteriori pooling: 2 sources ###################

### CLEAN MIXSIAR MODEL BY SAMPLEID RESULTS ###

SAMP.mixSIAR.EAAn.output <- read.table(
  "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_summary_stats_a_posteriori.EAAn.txt",
  sep = '\t', header = FALSE, skip = 6)

SAMP.mixSIAR.EAAn.output <- as.data.frame(SAMP.mixSIAR.EAAn.output)
colnames(SAMP.mixSIAR.EAAn.output) <- "full.data"

# remove all extra spaces using "squish", make them 1 space only
library(stringr)
SAMP.mixSIAR.EAAn.output$full.data <- str_squish(SAMP.mixSIAR.EAAn.output$full.data)

# rename columns
library(plyr)
samp.mix.EAAn.out <- separate(SAMP.mixSIAR.EAAn.output, full.data,  
                             into = c("factor", "Mean", "SD", "2.5.perc", "5.perc", "25.perc", "50.perc",
                                      "75.perc", "95.perc", "97.5.perc"), sep=" ")

# parse the factors here into 3 columns
new.cols.EAAn <- stringr::str_split_fixed(samp.mix.EAAn.out$factor, "\\.", 3) %>%
  as.data.frame() %>%
  setNames(c("proportion", "Source", "SampleID"))

new.cols.EAAn <- new.cols.EAAn %>%
  dplyr::select("SampleID", "Source")

# combine the new, separated columns with the data, verify the levels match, then remove "factor"
Samp.out.cleaned.EAAn <- cbind(new.cols.EAAn, samp.mix.EAAn.out[1:10])
Samp.out.cleaned.EAAn <- Samp.out.cleaned.EAAn %>%
  dplyr::select(-factor)

write.csv(Samp.out.cleaned.EAAn, "output/MixSIAR_models/EAAn_SampleID_model/MixSIAR_SampleID_a_posteriori.df.EAAn.csv", row.names = FALSE)

samp.mix.df<-Samp.out.cleaned.EAAn
samp.mix.reduced <- samp.mix.df %>% 
  dplyr::select(SampleID, Source, Mean, SD)

################# ################### ################ ################

# bring in some metadata to make sense of it
Mix.samp.EAAn.2 <- merge.data.frame(raw.norm.LD.data, samp.mix.reduced, by="SampleID")
Mix.samp.EAAn.2$Source <- factor(Mix.samp.EAAn.2$Source, 
                                 levels = c("Offshore", "Inshore"))
Mix.samp.EAAn.2$Mean<-as.numeric(Mix.samp.EAAn.2$Mean) 

write.csv(Mix.samp.EAAn.2, "output/Popp_Mix.samp.EAAn.2.csv")

########### plot
Mix.samp.plot.EAAn.InOff <- ggplot(data = Mix.samp.EAAn.2, aes(x = Group, y = Mean, fill = Source)) +
  geom_boxplot(alpha = 0.9, color = "black")+
  geom_point(aes(fill = Source), pch = 21, alpha = 0.5, color = "black", position = position_dodge(0.75))+
  scale_fill_manual(values = c("dodgerblue","springgreen4")) +
  ylab("Proportional Contribution") +
  xlab(NULL) +
  ylim(0,1) +
  guides(fill = guide_legend(title = expression(paste(AA[ESS], " ", Source)))) +
  theme_classic() + Fig.formatting + guides(color = "none")

#print and save grouped with 3 categories
Mix.samp.plot.EAAn.InOff
dev.copy(pdf, "figures/Popp_Fig.mixsiar.2sour.EAAn.pdf", width = 4, height = 7)
dev.off()

##########

write.csv(Mix.samp.EAAn.2, "output/Popp_Mix.samp.EAAn.2.csv")

Inshore.df<-Mix.samp.EAAn.2[(Mix.samp.EAAn.2$Source=="Inshore"),]

# group mean model
Inshore.mod<-lm(Mean~Group, data=Inshore.df) #p=0.974
anova(Inshore.mod)

# residency by mean model
Res.Inshore.mod<-lm(Mean~Residency.index, data=Inshore.df) #p=0.0145
summary(Res.Inshore.mod)
anova(Res.Inshore.mod)

Inshore.prop.plant<-ggplot(Inshore.df, aes(x=Residency.index, y=Mean)) +
  geom_smooth(method = lm, color="seagreen", alpha=0.2, fill="mediumseagreen")+
  geom_point(color="seagreen", size=2.5, aes(shape=Group))+
  ylab("Proportion Inshore-EAA")+
  xlab("Residency Index")+
  annotate("text", x = 0.08, y = 0.68, label = "Adj.R2=0.29, p=0.015", size = 3) +
  theme_classic()

Inshore.prop.plant
ggsave("figures/Popp_Inshore.prop.plantEAAn.pdf", width=5, height=7)


