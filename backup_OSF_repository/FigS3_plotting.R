rm(list= ls(all.names = TRUE))
gc()

library(ggplot2)
library(dplyr)
library(Hmisc)
library(ggplot2) # makes graphs
library(dplyr)
library(psych)
library(plyr)
library(ggpubr)
library(psycho)
library(rlist)
library(tidyverse)
library(data.table)
library(lme4)
library(apa)

# ggplot theme
mytheme <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(size=1, colour = "black"),
                 text = element_text(size=26,colour = "black"),
                 strip.background = element_blank(),
                 #axis.title.x=element_blank(), 
                 #axis.text.x=element_blank(), 
                 axis.ticks.x=element_blank(),
                 legend.title = element_blank())


data_path <- "D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure5-S3-SourceData"

# list ROIs to loop through
rois = c("LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", "Pons")

# prepare lists
p_val = list(numeric())
uber_collect_p <-list()
uber_collect_p <- append(uber_collect_p, list(p_val))

#make output dir
output_plots_dir = 'D:\\NYU_RS_LC\\stats\\OSF_fig_output'
dir.create(file.path(output_plots_dir), showWarnings = FALSE)

# make plot fig
setEPS()                                             # Set postscript arguments
postscript(paste(output_plots_dir, "XcorrPLOTS.eps",  sep='\\'),  width=4, height=5) 

par(mfrow=c(3,3),
    oma = c(5,4,0,0) + 0.7,
    mar = c(0,0,1,1) + 0.7)
i=0

session = c('sess1', 'sess2')

for (sess in session) {
  print(paste('running script on', sess))
  
  # loop through ROIs
  for (c_roi in rois){
    i + 1
    print(paste('running stats on ', c_roi))
    
    # load in data (Z val for pupil + BOLD signal at +4 to -4 time bins)
    day1 <- read.delim(paste(data_path, 'ses-day-1', paste("XCorr_stat_", c_roi, '_pup_size.csv',  sep = ''), 
                             sep = '\\'), header = TRUE, sep = ',')
    names(day1)[names(day1) == 'CC'] <- 'CC_D1'
    day2 <- read.delim(paste(data_path, 'ses-day-2', paste("XCorr_stat_", c_roi, '_pup_size.csv',  sep = ''), 
                             sep = '\\'), header = TRUE, sep = ',')
    names(day2)[names(day2) == 'CC'] <- 'CC_D2'
    # settings
    day1$lag <-as.factor(day1$lag)
    day2$lag <-as.factor(day2$lag)
    
    # join data from both days together 
    dat_all <- left_join(day1, day2, by = c('subj', 'lag'))
    # calculate average of both days
    dat_all$CC_ave =  rowMeans(dat_all[3:4], na.rm = TRUE)
    
    # settings
    dat_all$ROI <- as.character(c_roi)
    dat_all$subj <- as.character(dat_all$subj)
    
    # get data into long format 
    labels = names(dat_all)
    var_names = c("subj",
                  "lag",
                  "CC_D1", 
                  "CC_D2", 
                  "CC_ave", 
                  "ROI")
    dat = dat_all[,match(var_names,labels)]
    long_dat = data.frame(dat)
    # correct column names
    names(long_dat)[names(long_dat) == 'CC_D1'] <- 'CC_day1'
    names(long_dat)[names(long_dat) == 'CC_D2'] <- 'CC_day2'
    names(long_dat)[names(long_dat) == 'CC_ave'] <- 'CC_both'
    
    # loop over each time bin and run t-test (correct for multiple comparisons at end of script)
    P_mods = c("-4",
               '-3',
               "-2", 
               "-1", 
               "0", 
               "1", 
               "2",
               "3",
               "4")
    
    for (c_mod in P_mods) { 
      
      dat_model = subset(long_dat, lag == c_mod)
      
      if (sess == 'sess1') {
        # one-way t-tests
        stat.test_p1 <- 
          t.test(dat_model$CC_day1,
                 mu = 0) # Remove details
        print(paste('is the uncorrected p-value for Z corr at lag ', c_mod, 'for day 1', t_apa(stat.test_p1), sep =" "))
        p_val$p_value <- as.numeric(as.character(unlist(stat.test_p1[3][1])))
        uber_collect_p <- append(uber_collect_p, list(p_val))
      } else if (sess == 'sess2') {
        # one-way t-tests
        stat.test_p1 <- 
          t.test(dat_model$CC_day2,
                 mu = 0) # Remove details
        print(paste('is the uncorrected p-value for Z corr at lag ', c_mod, 'for day 2', t_apa(stat.test_p1), sep =" "))
        p_val$p_value <- as.numeric(as.character(unlist(stat.test_p1[3][1])))
        uber_collect_p <- append(uber_collect_p, list(p_val))
      } 
    }
    # some people excluded from session 1 or session 2 (noted as NaN) - removed before plotting
    long_dat_fin <- subset(long_dat, !is.na(CC_both))
    long_dat_fin_D1 <- subset(long_dat, !is.na(CC_day1))
    long_dat_fin_D2 <- subset(long_dat, !is.na(CC_day2))
    
    if (c_roi == 'LC') {
      colour_plot = '#00C1AA'
    } else if (c_roi == 'VTA') {
      colour_plot = '#ED8141' 
    } else if (c_roi == 'SN') {
      colour_plot = '#9590FF' 
    } else if (c_roi == 'DR') {
      colour_plot = '#FF62BC' 
    } else if (c_roi == 'MR') {
      colour_plot = '#5BB300'
    } else if (c_roi == 'BF') {
      colour_plot = 'goldenrod1'
    } else if (c_roi == 'ACC') {
      colour_plot = 'tan3'
    } else if (c_roi == 'OCC') {
      colour_plot = 'grey'
    } else if (c_roi == 'Pons') {
      colour_plot = '#a2bffe'
    }
    
    
    # compute averages 
    if (sess == 'sess1') {
      AVE_long_dat_fin <- describeBy(long_dat_fin$CC_day1, group = list(long_dat_fin$lag), mat = TRUE)
    } else if (sess == 'sess2') {  
      AVE_long_dat_fin <- describeBy(long_dat_fin$CC_day2, group = list(long_dat_fin$lag), mat = TRUE)
    }
    names(AVE_long_dat_fin)[c(2)] <- c("lag")
    # compute SEM vals 
    AVE_long_dat_fin$corr_max <- AVE_long_dat_fin$mean + AVE_long_dat_fin$se
    AVE_long_dat_fin$corr_min <- AVE_long_dat_fin$mean - AVE_long_dat_fin$se
    
    # Make cross corr plots
    ses <- AVE_long_dat_fin$mean + outer(AVE_long_dat_fin$se, c(1,-1))
    x1  = factor(AVE_long_dat_fin$lag, levels=c(-4,-3,-2,-1, 0, 1, 2, 3, 4))
    
    if (i == 4 | 8) {
      with(AVE_long_dat_fin, 
           graphics::plot(
             lag, mean, type="b",
             ylim = c(-0.10, 0.10), # different y-lim for cortical vs subcortical ROIs
             xlab = "",
             ylab = "",
             family = "A",
             main=c_roi,
             panel.first=polygon(c(lag,rev(lag)), c(ses[,1],rev(ses[,2])), 
                                 border=NA, 
                                 grid(NULL, NULL),
                                 col=colour_plot)))
      lines(AVE_long_dat_fin$lag, AVE_long_dat_fin$mean , col = "black")
    } else {
      with(AVE_long_dat_fin, 
           graphics::plot(
             lag, mean, type="b",
             ylim = c(-0.04, 0.04),
             xlab = "",
             ylab = "",
             cex.main=1.5,
             family = "A",
             main=c_roi,
             panel.first=polygon(c(lag,rev(lag)), c(ses[,1],rev(ses[,2])), 
                                 border=NA, 
                                 grid(NULL, NULL),
                                 col=colour_plot)))
      lines(AVE_long_dat_fin$lag, AVE_long_dat_fin$mean , col = "black")
    }
    
    
    abline(a=NULL, b=NULL, h=NULL, v=NULL)
    abline(h=0, col="black",lwd=1, lty=2)
    abline(v=0, col="black", lwd=1, lty=2)
  }    
  title(xlab = "lag (s)",
        ylab = "Cross-correlation (normalized)",
        outer = TRUE, line = 3)
  
  dev.off()

  
  # correct for multple comparisons (FDR-correction) on t-tests for each time bin within each day
  options(scipen=999)
  df_pval <- as.data.frame(unlist(uber_collect_p))
  df_pval_corr <- p.adjust(df_pval$`unlist(uber_collect_p)`, method = "fdr", n = length(df_pval$`unlist(uber_collect_p)`))

}