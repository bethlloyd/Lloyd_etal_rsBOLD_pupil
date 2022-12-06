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
library(ts.extend)
library(pracma)
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

XCorr_path <- paste("D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure6-SourceData")

rois = c("LC", "VTA", "SN", "ACC", "DR", "MR", "BF", "OCC")

ana_type = 'Pxy'

#create data frame for sub and cort
df_sub <- data.frame(matrix(ncol = 11, nrow = 2170))
colnames(df_sub) <- c("subj", "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", "Pons", "freq")

total_ROIs = list(subj=character(),
                  roi=character(),
                  lag=factor(),
                  CC_both = numeric())
uber_tota_ROIs = list()
uber_tota_ROIs <- append(uber_tota_ROIs, list(total_ROIs))

setEPS()                                             # Set postscript arguments
postscript('D:\\NYU_RS_LC\\stats\\OSF_fig_output\\Pxy_cpsd_stat.eps',  width=4, height=3) 


par(mfrow=c(2,4),
    oma = c(5,4,0,0) + 0.7,
    mar = c(0,0,1,1) + 0.7)
i=1
for (c_roi in rois){
  i = i + 1
  print(c_roi)

  # load in data and rename column
  day1 <- read.delim(paste(XCorr_path, paste('Pxy_cpsd_stat_ses-day1_', c_roi, '.csv',  sep = ''), 
                           sep = '\\'), header = TRUE, sep = ',')
  names(day1)[c(1:3)] <- c("subj", "freq", "power")
  day2 <- read.delim(paste(XCorr_path, paste('Pxy_cpsd_stat_ses-day2_', c_roi, '.csv',  sep = ''), 
                           sep = '\\'), header = TRUE, sep = ',')
  names(day2)[c(1:3)] <- c("subj", "freq", "power")
  
  day1$freq <-as.factor(day1$freq)
  day2$freq <-as.factor(day2$freq)
  
  
  # join data from both days together 
  dat_all <- left_join(day1, day2, by = c('subj', 'freq'))
  dat_all$power_ave =  rowMeans(dat_all[3:4], na.rm = TRUE)
  dat_all$ROI <- as.character(c_roi)
  dat_all$subj <- as.character(dat_all$subj)
  
  
  # get data into long format 
  labels = names(dat_all)
  # Select variables
  var_names = c("subj",
                "freq",
                "power.x", 
                "power.y", 
                "power_ave", 
                "ROI")
  
  dat = dat_all[,match(var_names,labels)]
  
  # Make data frame
  long_dat = data.frame(dat)
  
  names(long_dat)[names(long_dat) == 'power.x'] <- 'power_day1'
  names(long_dat)[names(long_dat) == 'power.y'] <- 'power_day2'
  names(long_dat)[names(long_dat) == 'power_ave'] <- 'power_ave'
  
  # colour settings
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
 
  # average across participants 
  long_dat <- subset(long_dat, !is.na(power_ave))
  AVE_long_dat_fin <- describeBy(long_dat$power_ave, group = list(long_dat$freq), mat = TRUE)
  names(AVE_long_dat_fin)[c(2)] <- c("freq")
  
  df_sub[, i] <- long_dat$power_ave
  
  # Create plots 
  
  ses <- AVE_long_dat_fin$mean + outer(AVE_long_dat_fin$se, c(1,-1))
  print(i)
  if (i == 5 | i == 9) {
    ylimits = c(-0.1, 0.9)
  } else {
    ylimits = c(-0.1, 0.6)
  }
  
  with(AVE_long_dat_fin, 
       plot(
         freq, mean, type="l",
         ylim = ylimits,
         xlab = "",
         ylab = "",
         main=c_roi,
         panel.first=polygon(c(freq,rev(freq)), c(ses[,1],rev(ses[,2])), 
                             border=NA, 
                             grid(NULL, NULL),
                             col=colour_plot)))
  axis(side = 1, at = seq(0,0.25,0.05),labels = T)
  if (i == 5 | i == 9) {
    axis(2, at = c(0, 0.4, 0.8))
  } else {
    axis(2, at = c(0, 0.3, 0.6))
  }
  
  
  # compute peak frequency
  print(paste('peak frequencey for ', c_roi, ' is: ', AVE_long_dat_fin$freq[which(AVE_long_dat_fin$mean == max(AVE_long_dat_fin$mean))]))
  
  # abline(a=NULL, b=NULL, h=NULL, v=NULL)
  # abline(h=0, col="black",lwd=1, lty=2)
  # abline(v=0, col="black", lwd=1, lty=2)
}    
title(xlab = "Frequency (Hz)",
      ylab = "Cross spectral power density (A.U.)",
      outer = TRUE, line = 3)

dev.off()


