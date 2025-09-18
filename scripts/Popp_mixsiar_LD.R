################################################
## MIXSIAR #####################################
################################################

# create a new data frame with desired columns for the zooplankton consumers (mixtures)
Mix.df.LD <- raw.norm.LD.data %>%
  dplyr::select("SampleID", "Group", "LD1", "LD2")

# remove sources
Mix.df.LD.shk <- Mix.df.LD[!(Mix.df.LD$Group == "ReefPlankton" | Mix.df.LD$Group == "PelagicPlankton" |
                               Mix.df.LD$Group == "Coral"),]
Mix.df.LD.shk$Group <- droplevels(Mix.df.LD.shk$Group)
Mix.df.LD.shk<-na.omit(Mix.df.LD.shk)

# write and export a csv of the zoyoplankton data frame
write.csv(Mix.df.LD.shk, "output/Popp_Mix.df.LD.shark.csv", row.names = FALSE)

###########
# SOURCES #
###########

# NORMALIZED DATA
# Reef and pelagic plankton, coral

# create a new data frame with desired columns for the primary producers (sources)
Source.df.LD <-  Mix.df.LD[(Mix.df.LD$Group == "ReefPlankton" | Mix.df.LD$Group == "PelagicPlankton" |
                                Mix.df.LD$Group == "Coral"),]

Source.df.LD <- Source.df.LD %>%
  dplyr::select("Group", "LD1", "LD2")

##### source mean and SD
# LD1
# calculate mean Thr values for the sources
mean.LD1.df <- aggregate(LD1 ~ Group, data = Source.df.LD, FUN = mean)
colnames(mean.LD1.df) <- c("Group", "Mean.LD1")

# calculate standard deviations of Thr values for the sources
sd.LD1.df <- aggregate(LD1 ~ Group, data = Source.df.LD, FUN = sd)
colnames(sd.LD1.df) <- c("Group", "SD.LD1")

# calculate mean Thr values for the sources
mean.LD2.df <- aggregate(LD2 ~ Group, data = Source.df.LD, FUN = mean)
colnames(mean.LD2.df) <- c("Group", "Mean.LD2")

# calculate standard deviations of Thr values for the sources
sd.LD2.df <- aggregate(LD2 ~ Group, data = Source.df.LD, FUN = sd)
colnames(sd.LD2.df) <- c("Group", "SD.LD2")

#### sample size
# calculate the number of samples for each of the sources
df.size.LD <- aggregate(LD2 ~ Group, data = Source.df.LD, FUN = length)
colnames(df.size.LD) <- c("Group", "LD2.size")

# bind the values calculated above into one data frame
Mix.df.source.agg.df.LD <- cbind(df.size.LD, 
                           mean.LD1.df[2], sd.LD1.df[2], 
                           mean.LD2.df[2], sd.LD2.df[2])
colnames(Mix.df.source.agg.df.LD) <- c("", "n", "MeanLD1", "SDLD1", "MeanLD2", "SDLD2")

# write and export a csv of the source data frame
write.csv(Mix.df.source.agg.df.LD, "output/LD_MIX/Popp_Mix.df.source.agg.df.LD.csv", row.names = FALSE, na = "")


# load zooplankton (consumer/mixture) data and assign factors
# SampleID is a fixed factor
cons.LD <- load_mix_data(filename = "output/LD_MIX/Popp_Mix.df.LD.shark.csv",
                           iso_names = c("LD1", "LD2"),
                           factors = "SampleID",
                           fac_random = FALSE,
                           fac_nested = FALSE,
                           cont_effects = NULL)

# load source data
source.LD <- load_source_data(filename = "output/LD_MIX/Popp_Mix.df.source.agg.df.LD.csv",
                                source_factors = NULL,
                                conc_dep = FALSE,
                                data_type = "means",
                                cons.LD)

#Load TDF data
discr.LD <- load_discr_data(filename = "output/LD_MIX/MixSIAR_LD_TDFs.csv", cons.LD)


################## ################## ################## 
################## RUN MIXSIAR MODEL BY SAMPLEID #######

# define the model and error structure and write JAGS model file
# PROCESS ERROR ONLY - HAVE TO DO THIS WHEN RUNNING MIXTURE POINTS INDIVIDUALLY
model_filename.LD <- "output/MixSIAR_models/LD_SampleID_model/MixSIAR_Model_SampleID_fixed.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename.LD, resid_err, process_err, cons.LD, source.LD)

################################################
############# run the JAGS model

set.seed(12)
GTGS.jags.mod.LD <- run_model(run = "normal", cons.LD, source.LD, discr.LD, model_filename.LD, alpha.prior = 1, resid_err, process_err)

#save RDS object of JAGS model
saveRDS(GTGS.jags.mod.LD, file = "output/MixSIAR_models/LD_SampleID_model/GTGS.jags.mod.LD.rds")

# code below will read back in the RDS object
GTGS.jags.mod.LD <- readRDS(file = "output/MixSIAR_models/LD_SampleID_model/GTGS.jags.mod.LD.rds")

# process diagnostics, summary stats, and posterior plots
output_options_LD <- list(
  summary_save = TRUE,                   
  summary_name = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_summary_stats.LD",
  sup_post = FALSE,                       
  plot_post_save_pdf = FALSE,              
  plot_post_name = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_posterior_density.LD",   
  sup_pairs = FALSE,                       
  plot_pairs_save_pdf = TRUE,             
  plot_pairs_name = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_pairs_plot.LD",         
  sup_xy = FALSE,                          
  plot_xy_save_pdf = FALSE,                
  plot_xy_name = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_xy_plot.LD",             
  gelman = TRUE,                         
  heidel = FALSE,                         
  geweke = TRUE,                         
  diag_save = TRUE,                       
  diag_name = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_diagnostics.LD",
  indiv_effect = FALSE,
  plot_post_save_png = FALSE,            
  plot_pairs_save_png = FALSE,            
  plot_xy_save_png = FALSE,
  graphics.off())

# set max to get all data
options(max.print = 3000)
out.JAGS.LD <- output_JAGS(GTGS.jags.mod.LD, cons.LD, source.LD, output_options_LD)

# get the output stats for the model
# output_stats(zoop.jags.mod.EAnA, cons.LD, source.LD, output_options_LD)

######### a posteriori aggregation of sources
# combine benthic algae and POM into one autochthonous source
GTGS.in.out.combined.LD <- combine_sources(GTGS.jags.mod.LD, cons.LD, source.LD, alpha.prior = 1,
                                              groups = list(Inshore = c("ReefPlankton", "Coral"),
                                                            Offshore = "PelagicPlankton"))

# print and save the interval plot
plot_intervals(GTGS.in.out.combined.LD, toplot = "fac1")
dev.copy(pdf, "output/MixSIAR_models/LD_SampleID_model/combines_output_intervals.LD.pdf", width=10, height=10)
dev.off()

# summary statistics for a posteriori aggregation of sources (inshore vs offshore)
summary_stat(GTGS.in.out.combined.LD, toprint = "fac1", meanSD = TRUE, savetxt = TRUE,
             quantiles = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
             filename = "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_summary_stats_a_posteriori.LD")


################################################
### CLEAN MIXSIAR MODEL BY SAMPLEID RESULTS ###
################################################

SAMP.mixSIAR.LD.output.3 <- read.table(
  "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_summary_stats.LD.txt",
  sep = '\t', header = FALSE, skip = 6)

SAMP.mixSIAR.LD.output.3 <- as.data.frame(SAMP.mixSIAR.LD.output.3)
colnames(SAMP.mixSIAR.LD.output.3) <- "full.data"

# remove all extra spaces using "squish", make them 1 space only
library(stringr)
SAMP.mixSIAR.LD.output.3$full.data <- str_squish(SAMP.mixSIAR.LD.output.3$full.data)

# rename columns
library(plyr)
samp.mix.LD.out.3 <- separate(SAMP.mixSIAR.LD.output.3, full.data,  
                                into = c("factor", "Mean", "SD", "2.5.perc", "5.perc", "25.perc", "50.perc",
                                         "75.perc", "95.perc", "97.5.perc"), sep=" ")

# parse the factors here into 3 columns
new.cols.LD.3 <- stringr::str_split_fixed(samp.mix.LD.out.3$factor, "\\.", 3) %>%
  as.data.frame() %>%
  setNames(c("proportion", "SampleID", "Source"))

new.cols.LD.3 <- new.cols.LD.3 %>%
  dplyr::select("SampleID", "Source")

# combine the new, separated columns with the data, verify the levels match, then remove "factor"
Samp.out.cleaned.LD.3 <- cbind(new.cols.LD.3, samp.mix.LD.out.3[1:10])
Samp.out.cleaned.LD.3 <- Samp.out.cleaned.LD.3 %>%
  dplyr::select(-factor)

Samp.out.cleaned.LD.3<-Samp.out.cleaned.LD.3[-1,]

write.csv(Samp.out.cleaned.LD.3, "output/MixSIAR_SampleID.df.LD.3sources.csv", row.names = FALSE)

################## ################## ################## 
################## ready to plot ################## 

# bring in some metadata to make sense of it
Mix.samp.LD.3 <- merge.data.frame(raw.norm.LD.data, Samp.out.cleaned.LD.3, by="SampleID")
Mix.samp.LD.3$Source <- factor(Mix.samp.LD.3$Source, 
                                 levels = c("PelagicPlankton", "ReefPlankton", "Coral"))
Mix.samp.LD.3$Mean<-as.numeric(Mix.samp.LD.3$Mean) 

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
Mix.samp.plot.LD.3 <- ggplot(data = Mix.samp.LD.3, aes(x = Group, y = Mean, fill = Source)) +
  geom_boxplot(alpha = 0.9, color = "black")+
  geom_point(aes(fill = Source), pch = 21, alpha = 0.5, color = "black", 
             position = position_dodge(0.75))+
  scale_fill_manual(values = c("dodgerblue", "springgreen4","coral")) +
  ylab("Proportional Contribution") +
  xlab(NULL) +
  ylim(0,1) +
  ggtitle("3 Sources on LD1 and LD2") +
  guides(fill = guide_legend(title = expression(paste(AA[ESS], " ", Source)))) +
  theme_classic() + Fig.formatting + guides(color = "none")

#print and save
Mix.samp.plot.LD.3
dev.copy(pdf, "figures/Popp_Fig.mixsiar.3sour.LD.pdf", width = 6, height = 7)
dev.off()


################################ #################### ##########
############ a posteriori pooling: 2 sources ###################

### CLEAN MIXSIAR MODEL BY SAMPLEID RESULTS ###

SAMP.mixSIAR.LD.output <- read.table(
  "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_summary_stats_a_posteriori.LD.txt",
  sep = '\t', header = FALSE, skip = 6)

SAMP.mixSIAR.LD.output <- as.data.frame(SAMP.mixSIAR.LD.output)
colnames(SAMP.mixSIAR.LD.output) <- "full.data"

# remove all extra spaces using "squish", make them 1 space only
library(stringr)
SAMP.mixSIAR.LD.output$full.data <- str_squish(SAMP.mixSIAR.LD.output$full.data)

# rename columns
library(plyr)
samp.mix.LD.out <- separate(SAMP.mixSIAR.LD.output, full.data,  
                              into = c("factor", "Mean", "SD", "2.5.perc", "5.perc", "25.perc", "50.perc",
                                       "75.perc", "95.perc", "97.5.perc"), sep=" ")

# parse the factors here into 3 columns
new.cols.LD <- stringr::str_split_fixed(samp.mix.LD.out$factor, "\\.", 3) %>%
  as.data.frame() %>%
  setNames(c("proportion", "Source", "SampleID"))

new.cols.LD <- new.cols.LD %>%
  dplyr::select("SampleID", "Source")

# combine the new, separated columns with the data, verify the levels match, then remove "factor"
Samp.out.cleaned.LD <- cbind(new.cols.LD, samp.mix.LD.out[1:10])
Samp.out.cleaned.LD <- Samp.out.cleaned.LD %>%
  dplyr::select(-factor)

write.csv(Samp.out.cleaned.LD, "output/MixSIAR_models/LD_SampleID_model/MixSIAR_SampleID_a_posteriori.df.LD.csv", row.names = FALSE)

samp.mix.df<-Samp.out.cleaned.LD
samp.mix.reduced <- samp.mix.df %>% 
  dplyr::select(SampleID, Source, Mean, SD)

################# ################### ################ ################

# bring in some metadata to make sense of it
Mix.samp.LD.2 <- merge.data.frame(raw.norm.LD.data, samp.mix.reduced, by="SampleID")
Mix.samp.LD.2$Source <- factor(Mix.samp.LD.2$Source, 
                                 levels = c("Offshore", "Inshore"))
Mix.samp.LD.2$Mean<-as.numeric(Mix.samp.LD.2$Mean) 

write.csv(Mix.samp.LD.2, "output/Popp_Mix.samp.LD.2.csv")

########### plot
Mix.samp.plot.LD.InOff <- ggplot(data = Mix.samp.LD.2, aes(x = Group, y = Mean, fill = Source)) +
  geom_boxplot(alpha = 0.9, color = "black")+
  geom_point(aes(fill = Source), pch = 21, alpha = 0.5, color = "black", position = position_dodge(0.75))+
  scale_fill_manual(values = c("dodgerblue","springgreen4")) +
  ylab("Proportional Contribution") +
  xlab(NULL) +
  ylim(0,1) +
  guides(fill = guide_legend(title = expression(paste(AA[ESS], " ", Source)))) +
  theme_classic() + Fig.formatting + guides(color = "none")

#print and save grouped with 3 categories
Mix.samp.plot.LD.InOff
dev.copy(pdf, "figures/Popp_Fig.mixsiar.2sour.LD.pdf", width = 5, height = 7)
dev.off()

##########

write.csv(Mix.samp.LD.2, "output/LD_MIX/Popp_Mix.samp.LD.2.csv")

Inshore.LD.df<-Mix.samp.LD.2[(Mix.samp.LD.2$Source=="Inshore"),]

# group mean model
Inshore.LD.mod<-lm(Mean~Group, data=Inshore.LD.df) #p=0.0.457
anova(Inshore.LD.mod)

# residency by mean model
Res.Inshore.LD.mod<-lm(Mean~Residency.index, data=Inshore.LD.df) #p=0.0145
summary(Res.Inshore.LD.mod)
anova(Res.Inshore.LD.mod)

Inshore.prop.plant.LD<-ggplot(Inshore.LD.df, aes(x=Residency.index, y=Mean)) +
  geom_smooth(method = lm, color="seagreen", alpha=0.2, fill="mediumseagreen")+
  geom_point(color="seagreen", size=2.5, aes(shape=Group))+
  ylab("Proportion Inshore-EAA")+
  xlab("Residency Index")+
  annotate("text", x = 0.08, y = 0.88, label = "Adj.R2=0.20, p=0.043", size = 3) +
  theme_classic()

Inshore.prop.plant.LD
ggsave("figures/Inshore.prop.plant.LD.pdf", width=7, height=5)


