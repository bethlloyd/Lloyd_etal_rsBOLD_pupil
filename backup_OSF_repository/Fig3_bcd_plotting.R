# Author: Beth Lloyd
rm(list= ls(all.names = TRUE))
gc()

# libraries
library(ggplot2) # makes graphs
library(dplyr)
library(psych)
library(plyr)
library(ggpubr)
library(psycho)
library(rlist)
library(tidyverse)
library(data.table)
library(heemod)
library(phonTools)
library(ggsci)
library(apa)
library(rstatix)
library(Hmisc)
# ggplot theme
mytheme <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(size=1, colour = "black"),
                 text = element_text(size=30,colour = "black"),
                 strip.background = element_blank())
                 #axis.title.x=element_blank(), 
                 #axis.text.x=element_blank(), 
                 #axis.ticks.x=element_blank(),
                 #legend.title = element_blank())
# -----------------------------------------------------------------------------------------------------

# SET PATHS
group.colors =c(LC='#00C1AA', VTA='#ED8141', SN='#9590FF', DR='#FF62BC', MR='#5BB300', BF='goldenrod1', ACC='tan3', OCC='grey', Pons='#a2bffe')
#home comp
datapath = "D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure3-SourceData"
save_folder = 'D:\\NYU_RS_LC\\stats\\OSF_fig_output'

# Make Figure 3 B
# -------------------------

# load in data
rsHRF_data  <- (paste(datapath, 'Fig3_B_rsHRF_extracted_HRFs_both_days.csv', sep = '\\'))
rsHRF_data <- read.csv(rsHRF_data, header = TRUE, sep = ',')
rsHRF_data <- as.data.frame(rsHRF_data)

# Get info
labels = names(rsHRF_data)
# Select variables
var_names = c("subj",
              "time",
              "LC",
              "VTA",
              "SN",
              "DR",
              "MR",
              "BF", 
              "ACC",
              "OCC",
              "Pons")
dat = rsHRF_data[,match(var_names,labels)]

# Make data frame
df_dat = data.frame(dat)

# Make long data format [ csp_preTlapse:csm_stim_mean ]
long_dat = tidyr::gather(df_dat, value=Val, key=Names, LC:Pons, factor_key=TRUE)

# Add variables and change names 
names(long_dat) = c("subj","time", "ROI","signal")
#have to do this:
long_dat$time <- as.character(long_dat$time)
long_dat$time <- as.numeric(long_dat$time)
long_dat$time <- as.factor(long_dat$time)
long_dat$ROI <- as.factor(long_dat$ROI)

long_dat$time <- as.character(long_dat$time)
long_dat$time <- as.numeric(long_dat$time)
long_dat$subj <- substr(long_dat$subj, nchar(long_dat$subj) - 3 + 1, nchar(long_dat$subj))
long_dat$subj <- as.factor(long_dat$subj)

## Plot timcourse - facet wrap 
g3 <-long_dat %>%
  mutate(ROI = fct_relevel(ROI, "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", 'Pons')) %>%
  ggplot(aes(x=time,y= signal, color = subj)) +
  geom_line(aes(color=subj), size = 0.1) + 
  #geom_ribbon(alpha=0.2) + 
  mytheme +
  ylab("Normalised BOLD response") +
  xlab("time (s)")  +
  #scale_color_brewer(palette="Dark2") +
  #scale_linetype_manual(values=c("solid", "solid"))+
  #scale_y_continuous(limits=c(-0.2,1.2))+
  scale_x_continuous(limits=c(0,20))+
  theme(legend.title=element_blank()) +
  labs(title = "") + 
  facet_wrap(~ROI) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")

g3a <- g3 + stat_summary(aes(group=ROI), fun.y=mean, geom="line", colour="black", size=1)
#ggsave(paste(save_folder,"indiv_subjs_overlayed_HRF_both_days.eps", sep=""), g3a, width = 12, height = 12)



# Make Figure 3 C
# -------------------------

# load in the number of events 
rsHRF_events_path  <- (paste(datapath, 'Fig3_C_rsHRF_roi_event_num_both_days.csv', sep = '\\'))
rsHRF_events <- read.csv(rsHRF_events_path, header = TRUE, sep = ',')
rsHRF_events <- as.data.frame(rsHRF_events)


# Get info
labels = names(rsHRF_events)
# Select variables
var_names = c("subj",
              "LC",
              "VTA",
              "SN",
              "DR",
              "MR",
              "BF", 
              "ACC",
              "OCC",
              "Pons")
dat = rsHRF_events[,match(var_names,labels)]

# Make data frame
df_dat = data.frame(dat)

# Make long data format [ csp_preTlapse:csm_stim_mean ]
long_dat = tidyr::gather(df_dat, value=Val, key=Names, LC:Pons, factor_key=TRUE)

# Add variables and change names 
names(long_dat) = c("subj", "ROI","events")
g_events <- long_dat %>%
  mutate(ROI = fct_relevel(ROI, "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", "Pons")) %>%
  ggplot(aes(x = ROI, y = events, fill = ROI)) +
  #geom_boxplot(lwd=0.2) + 
  #geom_point(aes(color=ROI), size=2, alpha=.4) +
  geom_jitter(width = 0.1, aes(color=ROI), size=1.5, alpha=1) + 
  #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
  #scale_y_continuous(limits=c(3, 8.3)) +
  #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
  stat_summary(fun="mean",colour = "black", size = 1.5,geom = "point") +
  stat_summary(fun.data =mean_se,colour = "black", width=0,size = 1,geom = "errorbar") +
  scale_color_manual(values=group.colors) +
  mytheme                 +                                 
  labs(title = "", x = "", y = "Number of psuedo events")   +
  theme(legend.position = "none") 
#geom_violin(alpha = 0, aes(color=ROI, fill = ROI))

#ggsave(paste(save_folder,"\\ROI_rsHRf_event_number.eps", sep=""), g_events, width = 6, height = 4)



# Make Figure 3 D
# -------------------------
#input data
TTP_data  <- (paste(datapath, 'Fig3_D_rsHRF_roi_TTP_both_days.csv', sep = '\\'))
TTP_data <- read.csv(TTP_data, header = TRUE, sep = ',')
TTP_data <- as.data.frame(TTP_data)


# Get info
labels = names(TTP_data)
# Select variables
var_names = c("subj",
              "LC",
              "VTA",
              "SN",
              "DR",
              "MR",
              "BF", 
              "ACC",
              "OCC",
              "Pons")
dat = TTP_data[,match(var_names,labels)]

# Make data frame
df_dat = data.frame(dat)

# Make long data format [ csp_preTlapse:csm_stim_mean ]
long_dat = tidyr::gather(df_dat, value=Val, key=Names, LC:Pons, factor_key=TRUE)

# Add variables and change names 
names(long_dat) = c("subj", "ROI","time")



# Run statistics on TTP differences 

# descriptive statistics:
long_dat %>%
  group_by(ROI) %>%
  summarise_at(vars(time), list(name = mean))

#------------- ANOVAS 
# check for main effect of model on T-stat
m1 = afex::aov_ez("subj", "time", 
                  long_dat, 
                  within = c("ROI"))
summary(m1)



#-------------- follow up t-test 

# Pairwise comparisons between models
# Paired t-test is used because we have repeated measures by time
stat.test <- long_dat %>%
  pairwise_t_test(
    time ~ ROI, paired = TRUE,
    p.adjust.method = "fdr"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test


g4 <- long_dat %>%
  mutate(ROI = fct_relevel(ROI, "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC")) %>%
  ggplot(aes(x = ROI, y = time, fill = ROI)) +
  #geom_point(aes(color=ROI), size=2, alpha=.4) +
  geom_jitter(width = 0.1, aes(color=ROI), size=1.5, alpha=1) + 
  #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
  scale_y_continuous(limits=c(3, 8.3)) +
  #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
  stat_summary(fun="mean",colour = "black", size = 1.5,geom = "point") +
  stat_summary(fun.data =mean_se,colour = "black", width=0,size = 1,geom = "errorbar") +
  scale_color_manual(values=group.colors) +
  mytheme                 +                                 
  labs(title = "", x = "", y = "time to peak (s)")   +
  theme(legend.position = "none") 
  #geom_violin(alpha = 0, aes(color=ROI, fill = ROI))



#geom_hline(yintercept=0, linetype="dashed", color = "black")
#ggsave(paste(save_folder,"ROI_TTP_", bf, "_", day_dat, ".png", sep=""), g4, width = 8, height = 8)
#ggsave(paste(save_folder,"ROI_TTP_both_days.eps", sep=""), g4, width = 6, height = 5)
  






