rm(list= ls(all.names = TRUE))
gc()

library(ggplot2)
library(dplyr)
library(Hmisc)
library(tidyverse)
library(RColorBrewer)
library(psych)
library(rstatix)
library(apa)
library(RColorBrewer)
library(ggcorrplot)

group.colors =c(LC='#00C1AA', VTA='#ED8141', SN='#9590FF', DR='#FF62BC', MR='#5BB300', BF='goldenrod1', ACC='tan3', OCC='grey', Pons='#a2bffe')

# ggplot theme
mytheme <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(size=1, colour = "black"),
                 text = element_text(size=22,colour = "black"),
                 strip.background = element_blank(),
                 axis.title.x=element_blank(), 
                 axis.text.x=element_blank(), 
                 axis.ticks.x=element_blank(),
                 legend.title = element_blank())


Fig2_path <- paste("D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure2-SourceData")

## Create Figure 2 B
## --------------------------

output_path <- "D:\\NYU_RS_LC\\stats\\OSF_fig_output"
day = c("day1", 'day2')

for (c_day in day){
  # load in data and rename column
  TSNR_dat <- read.delim(paste(Fig2_path, paste('Fig2_B_tSNR_',c_day,'.csv',sep=""),
                                    sep = '\\'), header = TRUE, sep = ',')

  # get data into long format 
  labels = names(TSNR_dat)
  # Select variables
  var_names = c("subj",
                "LC",
                'DR',
                "MR", 
                "VTA", 
                "SN",
                "BF",
                "OCC", 
                "ACC",
                "Pons")
  
  dat = TSNR_dat[,match(var_names,labels)]
  
  # Make data frame
  df_dat = data.frame(dat)
  
  # Make long data format [ csp_preTlapse:csm_stim_mean ]
  long_dat = tidyr::gather(df_dat, value=Val, key=Names, LC:Pons, factor_key=TRUE)
  
  names(long_dat)[names(long_dat) == 'Names'] <- 'ROI'
  names(long_dat)[names(long_dat) == 'Val'] <- 'tSNR'
  
  long_dat$ROI <- as.factor(long_dat$ROI)
  
  long_dat_summary <- describeBy(long_dat$tSNR, 
                                 group = list(long_dat$ROI), mat = TRUE)
  
  plot1 <-subset(long_dat, !is.na(tSNR)) %>%
    mutate(ROI = fct_relevel(ROI, "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", "Pons")) %>%
    ggplot(aes(x = ROI, y = tSNR, fill = ROI)) +
    #geom_point(aes(color=ROI), size=2, alpha=.4) +
    geom_jitter(width = 0.1, aes(color=ROI), size=2, alpha=1) + 
    #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
    scale_y_continuous(limits=c(0, 120)) +
    #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
    stat_summary(fun="mean",colour = "black", size = 4,geom = "point", alpha=1) +
    stat_summary(fun.data="mean_se",colour = "black", width=0,size = 1,geom = "errorbar", alpha=1) +
    scale_color_manual(values=group.colors) +
    #scale_color_brewer(palette="Dark2") +
    mytheme                 +                                 
    labs(title = "", x = "", y = "tSNR")   +
    theme(legend.position = "none") 
  
  ggsave(paste(output_path, paste("tSNR_", c_day, '.eps',  sep = ''), sep='\\'),  height = 6, width = 4,plot1)
  
  
  stat.test <- long_dat %>%
    pairwise_t_test(
      tSNR ~ ROI, paired = TRUE,
      p.adjust.method = "fdr"
    ) %>%
    select(-df, -statistic, -p) # Remove details
  stat.test
  
  
}


## Create Figure 2 C
## --------------------------

options(scipen=999)

p_val = list(numeric())
uber_collect_p <-list()
uber_collect_p <- append(uber_collect_p, list(p_val))

r_val = list(numeric())
uber_collect_r <-list()
uber_collect_r <- append(uber_collect_r, list(r_val))

# load in data and rename column
BS_corr_dat_day1 <- read.delim(paste(Fig2_path, paste('Fig2_C_day1_BS_partial_PONS_correlations.csv',sep=""),
                                     sep = '\\'), header = TRUE, sep = ',')
BS_corr_dat_day2 <- read.delim(paste(Fig2_path, paste('Fig2_C_day2_BS_partial_PONS_correlations.csv',sep=""),
                                     sep = '\\'), header = TRUE, sep = ',')


# get data into long format 
labels = names(BS_corr_dat_day1)
# Select variables
var_names = c("subj",
              "LC_VTA",
              'LC_SN',
              "LC_DR", 
              "LC_MR", 
              "VTA_SN", 
              "VTA_DR", 
              "VTA_MR",
              "SN_DR",
              "SN_MR",
              "DR_MR",
              "LC_BF",
              "VTA_BF",
              "SN_BF",
              "DR_BF",
              "MR_BF")

dat_day1 = BS_corr_dat_day1[,match(var_names,labels)]
dat_day2 = BS_corr_dat_day2[,match(var_names,labels)]

# Make data frame
df_dat_day1 = data.frame(dat_day1)
df_dat_day2 = data.frame(dat_day2)


# Make long data format [ csp_preTlapse:csm_stim_mean ]
long_dat_day1 = tidyr::gather(df_dat_day1, value=Val, key=Names, LC_VTA:MR_BF, factor_key=TRUE)
names(long_dat_day1)[names(long_dat_day1) == 'Names'] <- 'BS_ROIs'
names(long_dat_day1)[names(long_dat_day1) == 'Val'] <- 'r_value_day1'
long_dat_day2 = tidyr::gather(df_dat_day2, value=Val, key=Names, LC_VTA:MR_BF, factor_key=TRUE)
names(long_dat_day2)[names(long_dat_day2) == 'Names'] <- 'BS_ROIs'
names(long_dat_day2)[names(long_dat_day2) == 'Val'] <- 'r_value_day2'

day_both_days <- left_join(long_dat_day1, long_dat_day2, by = c('subj', 'BS_ROIs'))

day_both_days$r_values_average <- rowMeans(day_both_days[3:4], na.rm=TRUE)
average_corr <- describeBy(day_both_days$r_values_average, group = list(day_both_days$BS_ROIs), mat = TRUE, na.rm=TRUE)

#write.csv(average_corr, paste(Fig2_path, 'average_corr.csv'))

# perform fisher r to z transform
day_both_days$fisher_z <- fisherz(day_both_days$r_values_average)
# group by correlation- 
long_dat_summary <- describeBy(r_values_average ~ BS_ROIs, mat = TRUE, data= day_both_days)


# run ANOVA
m1 = afex::aov_ez("subj", "fisher_z", 
                  day_both_days, 
                  within = c("BS_ROIs"))
anova_apa(m1)

# T-tests

P_mods = c("LC_VTA",
           'LC_SN',
           "LC_DR", 
           "LC_MR", 
           "VTA_SN", 
           "VTA_DR", 
           "VTA_MR",
           "SN_DR",
           "SN_MR",
           "DR_MR",
           "LC_BF",
           "VTA_BF",
           "SN_BF",
           "DR_BF",
           "MR_BF")
i = 1

# run indiv t-tests on each roi
for (c_mod in P_mods) { 
  
  dat_model = subset(day_both_days, BS_ROIs == c_mod)
  
  # one-way t-tests
  stat.test_p1 <- 
    t.test(dat_model$fisher_z,
           mu = 0) # Remove details
  
  #t_apa(stat.test_p1)
  
  print(paste('p value for t-test correlation', c_mod, 'is', t_apa(stat.test_p1), sep =" "))
  p_val$p_value <- as.numeric(as.character(unlist(stat.test_p1[3][1])))
  uber_collect_p <- append(uber_collect_p, list(p_val))
  
} 


# bind p-vals and correct for number of tests (FDR-correction) (per roi)
df_pval <- as.data.frame(unlist(uber_collect_p))
df_pval_corr <- p.adjust(df_pval$`unlist(uber_collect_p)`, method = "fdr", n = length(df_pval$`unlist(uber_collect_p)`))

# plot the correlations
day_both_days$BS_ROIs <- as.factor(day_both_days$BS_ROIs)

plot1 <-subset(day_both_days, !is.na(fisher_z)) %>%
  ggplot(aes(x = BS_ROIs, y = fisher_z, fill = BS_ROIs)) +
  #geom_point(aes(color=ROI), size=2, alpha=.4) +
  geom_jitter(width = 0.1, aes(color=BS_ROIs), size=2, alpha=.4) + 
  #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
  scale_y_continuous(limits=c(-1, 1)) +
  #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
  stat_summary(fun="mean",colour = "black", size = 2,geom = "point") +
  stat_summary(fun.data="mean_se",colour = "black", width=0,size = 1,geom = "errorbar") +
  #scale_color_brewer(palette="Dark2") +
  mytheme                 +                                 
  labs(title = "", x = "", y = "Fisher Z correlation coefficient")   +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")

#ggsave(paste(dat_path, dat_type, 'plots', paste("BS_correlations_both_days.png",  sep = ''), sep='\\'),  height = 6, width = 12, plot1)


# make correlation matrix 
mat_arr <- array(dim=c(6,6))
mat_arr[1,1] = 1 # LC-LC
mat_arr[2,1] = long_dat_summary$mean[1] # VTA-LC
mat_arr[1,2] = long_dat_summary$mean[1] # VTA-LC

mat_arr[3,1] = long_dat_summary$mean[2] # SN-LC
mat_arr[1,3] = long_dat_summary$mean[2] # SN-LC

mat_arr[4,1] = long_dat_summary$mean[3] # DR-LC
mat_arr[1,4] = long_dat_summary$mean[3] # DR-LC

mat_arr[5,1] = long_dat_summary$mean[4] # MR-LC
mat_arr[1,5] = long_dat_summary$mean[4] # MR-LC

mat_arr[2,2] = 1 # VTA-VTA
mat_arr[3,2] = long_dat_summary$mean[5] # SN-VTA
mat_arr[2,3] = long_dat_summary$mean[5] # SN-VTA

mat_arr[4,2] = long_dat_summary$mean[6] # DR-VTA
mat_arr[2,4] = long_dat_summary$mean[6] # DR-VTA

mat_arr[5,2] = long_dat_summary$mean[7] # MR-VTA
mat_arr[2,5] = long_dat_summary$mean[7] # MR-VTA

mat_arr[3,3] = 1 # SN-SN
mat_arr[4,3] = long_dat_summary$mean[8] # DR-SN
mat_arr[3,4] = long_dat_summary$mean[8] # DR-SN

mat_arr[5,3] = long_dat_summary$mean[9] # MR-SN
mat_arr[3,5] = long_dat_summary$mean[9] # MR-SN

mat_arr[4,4] = 1 # DR-DR
mat_arr[5,4] = long_dat_summary$mean[10] # MR-DR
mat_arr[4,5] = long_dat_summary$mean[10] # MR-DR

mat_arr[5,5] = 1 # MR-MR
mat_arr[6,1] = long_dat_summary$mean[11] # LC-BF
mat_arr[1,6] = long_dat_summary$mean[11] # LC-BF

mat_arr[6,2] = long_dat_summary$mean[12] # VTA_BF
mat_arr[2,6] = long_dat_summary$mean[12] # VTA_BF

mat_arr[6,3] = long_dat_summary$mean[13] # SN_BF
mat_arr[3,6] = long_dat_summary$mean[13] # SN_BF

mat_arr[6,4] = long_dat_summary$mean[14] # DR_BF
mat_arr[4,6] = long_dat_summary$mean[14] # DR_BF

mat_arr[6,5] = long_dat_summary$mean[15] # MR_BF
mat_arr[5,6] = long_dat_summary$mean[15] # MR_BF
mat_arr[6,6] = 1 # BF-BF



#COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)



rownames(mat_arr) <- c("LC", "VTA", "SN", "DR", "MR", "BF")
colnames(mat_arr) <- c("LC", "VTA", "SN", "DR", "MR", "BF")

mat_arr<-as.matrix(mat_arr)
diag(mat_arr) <- NA  
#png(height=1800, width=1800, file=paste(dat_path, 'figures', 'BS_correlations', paste("BS_correlations_matrix_circles.png"), sep="\\"), type = "cairo")
# make basic correlation plot - add p-values manually
p = ggcorrplot(mat_arr,type = "lower")  

p + scale_fill_gradient2(limit = c(-0.4,0.4), low = "blue", high =  "red", mid = "white", midpoint = 0)





## Create Figure 2 D
## --------------------------

# load in data
LC_dat <- read.delim(paste(Fig2_path, 'Fig2_D_LC_intensities.csv',
                             sep = '\\'), header = TRUE, sep = ',')

save_fig <- 'D:\\NYU_RS_LC\\stats\\figures\\LC_intensity'

LC_dat$ID <- seq.int(nrow(LC_dat))
# get data into long format 
labels = names(LC_dat)
# Select variables
var_names = c("ID", "LC_cnr_intensity")

dat = LC_dat[,match(var_names,labels)]

# Make data frame
df_dat = data.frame(dat)

# Make long data format [ csp_preTlapse:csm_stim_mean ]
long_dat = tidyr::gather(df_dat, value=Val, key=Names, LC_cnr_intensity, factor_key=TRUE)


plot1 <-long_dat %>%
  ggplot(aes(x = Names, y = Val)) +
  #geom_point(size=2, alpha=.4) +
  geom_jitter(width = 0.02, size=1, alpha=1) + 
  #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
  scale_y_continuous(limits=c(0, 0.35)) +
  #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
  stat_summary(fun="mean",colour = "grey", size = 1.75,geom = "point", alpha=1) +
  stat_summary(fun.data="mean_se",colour = "grey", width=0,size = 0.5,geom = "errorbar", alpha=1) +
  #scale_color_brewer(palette="Dark2") +
  mytheme                 +                                 
  labs(title = "", x = "", y = "Locus Coeruleus CNR (a.u.)")   +
  theme(legend.position = "none")
ggsave(paste(save_fig, paste('LC_CNR_plot.eps',  sep = ''), sep='\\'),  height = 3, width = 2,plot1)