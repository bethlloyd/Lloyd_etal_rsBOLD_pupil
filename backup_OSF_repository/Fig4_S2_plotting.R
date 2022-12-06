rm(list= ls(all.names = TRUE))
gc()

library(ggplot2)
library(dplyr)
library(Hmisc)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(apa)
# ggplot theme
mytheme <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(size=0.1, colour = "black"),
                 text = element_text(size=12,colour = "black"),
                 strip.background = element_blank(),
                 axis.title.x=element_blank(), 
                 #axis.text.x=element_blank(), 
                 axis.ticks.x=element_blank(),
                 legend.title = element_blank())
group.colors =c(LC='#00C1AA', VTA='#ED8141', SN='#9590FF', DR='#FF62BC', MR='#5BB300', BF='goldenrod1', ACC='tan3', OCC='grey', Pons='#a2bffe')

# to plot figure S2, change this path
data_path <- 'D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure4-SourceData\\Fig4_BC'  # ''D:\\NYU_RS_LC\\stats\\OSF_repository\\FigureS2-SourceData\\FigS2_A''

# define ROIs to loop over
roi = c("DR", "MR", "VTA", "LC", "SN", "ACC", "OCC", "Pons", "BF")

# prepare lists 
p_val = list(numeric())
uber_collect_p <-list()
uber_collect_p <- append(uber_collect_p, list(p_val))
uber_collect_anov_p <- list()
uber_collect_anov_p <- append(uber_collect_anov_p, list(p_val))
mat = matrix(ncol = 0, nrow = 0)
df_uber=data.frame(mat)


# loop through ROIs, run statistics, and plot data
for (c_roi in roi){
  print(paste('running stats on ', c_roi))
  
  #make output dir
  output_plots_dir = 'D:\\NYU_RS_LC\\stats\\OSF_fig_output'
  
  # load in data for each TTP model and rename column
  conmap_mod1 <- read.delim(paste(data_path, paste("con_stat_rsHRF_P1_", c_roi,'.csv',  sep = ''), 
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod1)[names(conmap_mod1) == 'spmT_0001.nii'] <- 'P1'
  
  conmap_mod2 <- read.delim(paste(data_path, paste("con_stat_rsHRF_P2_", c_roi, '.csv',  sep = ''), 
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod2)[names(conmap_mod2) == 'spmT_0001.nii'] <- 'P2'
  
  conmap_mod3 <- read.delim(paste(data_path, paste("con_stat_rsHRF_p3_", c_roi,  '.csv',  sep = ''),  
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod3)[names(conmap_mod3) == 'spmT_0001.nii'] <- 'P3'
  
  conmap_mod4 <- read.delim(paste(data_path,  paste("con_stat_rsHRF_P4_", c_roi, '.csv',  sep = ''),
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod4)[names(conmap_mod4) == 'spmT_0001.nii'] <- 'P4'
  
  conmap_mod5 <- read.delim(paste(data_path, paste("con_stat_rsHRF_P5_", c_roi,  '.csv',  sep = ''),
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod5)[names(conmap_mod5) == 'spmT_0001.nii'] <- 'P5'
  
  conmap_mod6 <- read.delim(paste(data_path, paste("con_stat_rsHRF_P6_", c_roi, '.csv',  sep = ''),
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod6)[names(conmap_mod6) == 'spmT_0001.nii'] <- 'P6'
  
  # combine all TTP models into 1 dataframe 
  dat_all <- left_join(conmap_mod1, conmap_mod2, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod3, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod4, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod5, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod6, by = 'subj')
  
  # get data into long format 
  labels = names(dat_all)
  # Select variables
  var_names = c("subj",
                "P1",
                "P2", 
                "P3", 
                "P4", 
                "P5",
                "P6")
  dat = dat_all[,match(var_names,labels)]
  df_dat = data.frame(dat)
  long_dat = tidyr::gather(df_dat, value=Val, key=Names, P1:P6, factor_key=TRUE)
  # rename columns
  names(long_dat)[names(long_dat) == 'Names'] <- 'model'
  names(long_dat)[names(long_dat) == 'Val'] <- 't_stat'
  # rename elements
  long_dat_fin <- long_dat %>% 
    mutate(TTPmod = case_when(
      model == 'P1' ~ 1,
      model == 'P2' ~ 2,
      model == 'P3' ~ 3,
      model == 'P4' ~ 4,
      model == 'P5' ~ 5,
      model == 'P6' ~ 6))
  
  # prep data for appending all ROIs (to run big ANOVA on TTP + ROI effects) 
  long_dat$ROI <- replicate(length(long_dat$subj), c_roi)
  df_uber <- rbind(df_uber, long_dat)
  
  
  # statistics: ANOVAS 
  # check for main effect of ttp model on T-stat
  m1 = afex::aov_ez("subj", "t_stat", 
                    long_dat_fin, 
                    within = c("model"))
  paste('anova for model FOR ', c_roi, ' is:')
  m1$anova_table
  #summary(m1)
  #anova_apa(m1)
  anov_p = m1$anova_table$`Pr(>F)`
  # append all p-values in order to correct for multple comparisons (at end of script)
  uber_collect_anov_p <- append(uber_collect_anov_p, list(anov_p))
  
  #-------------- follow up t-test 
  #settings
  P1_dat <- subset(long_dat_fin, model == 'P1')
  P2_dat <- subset(long_dat_fin, model == 'P2')
  P3_dat <- subset(long_dat_fin, model == 'P3')
  P4_dat <- subset(long_dat_fin, model == 'P4')
  P5_dat <- subset(long_dat_fin, model == 'P5')
  P6_dat <- subset(long_dat_fin, model == 'P6')
  
  # run r-test on each ttp model separately (apply FDR correction later in script)
  P_mods = c('P1', 'P2', 'P3', 'P4', 'P5', 'P6')
  i = 1
  print(paste(c_roi))
  for (c_mod in P_mods) { 
    
    dat_model = subset(long_dat_fin, model == c_mod)
    
   # one-way t-tests
    stat.test_p1 <- 
      t.test(dat_model$t_stat,
        mu = 0) # Remove details
    
    #t_apa(stat.test_p1)

    #print(paste('p value for t-test model P', i, 'is', t_apa(stat.test_p1), sep =" "))
    p_val$p_value <- as.numeric(as.character(unlist(stat.test_p1[3][1])))
    uber_collect_p <- append(uber_collect_p, list(p_val))
    
    i <- i + 1
  } 
  
  
# ---- plot data: create Figure 4 B & C
  
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
  long_dat_fin$TTPmod <- as.factor(long_dat_fin$TTPmod)
  
  # change mod names 
  long_dat_fin <- long_dat_fin %>% 
    mutate(model = case_when(
      model == 'P1'  ~ '1',
      model == 'P2' ~ '2',
      model == 'P3' ~ '3',
      model == 'P4' ~ '4',
      model == 'P5' ~ '5',
      model == 'P6' ~ '6'))
  
  plot1 <- ggplot(data=subset(long_dat_fin, !is.na(t_stat)), aes(x = model, y = t_stat)) +
    #geom_point(aes(color=ROI), size=2, alpha=.4) +
    geom_jitter(width = 0.05, size=0.0002, alpha=1, color = colour_plot) + 
    #geom_line(data = fortify(quad_mod), aes(x = long_dat_fin$TTPmod, y = .fitted)) +
    #geom_line() + 
    #geom_line(aes(color=session, group=subj), size=.3, color="black", alpha=.4) +
    scale_y_continuous(limits=c(-3, 3)) +
    #scale_x_discrete()(limits  = c(1, 6)) +
    #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
    stat_summary(fun="mean",colour = "black", size = 0.6,geom = "point") +
    stat_summary(fun.data="mean_se",colour = "black", width=0,size = 0.6,geom = "errorbar") +
    scale_colour_manual(values = colour_plot) +
    mytheme                 +                                 
    labs(title = "", x = "", y = "")   +
    theme(legend.position = "none") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.1) +
    geom_violin(alpha = 0, color='black', size = 0.2) 
  

  print(plot1)
  #ggsave(paste(output_plots_dir, paste("roi_sig_extract_", c_roi, '.eps',  sep = ''), sep='\\'), height = 2, width = 1.5, plot1)
  
  
}    


# bind p-vals from all ttp models and correct for number of tests
df_pval <- as.data.frame(unlist(uber_collect_p))
df_pval_corr <- p.adjust(df_pval$`unlist(uber_collect_p)`, method = "fdr", n = length(df_pval$`unlist(uber_collect_p)`))

# df_pval_corr = list of corrected p-values for each ttp model for each ROI 


# correct p-values from running 1 anova per ROI: (FDR-correction for 9 rois)
anov_p <- as.data.frame(unlist(uber_collect_anov_p))
df_anov_pval_corr <- p.adjust(anov_p$`unlist(uber_collect_anov_p)`, method = "fdr", n = length(anov_p$`unlist(uber_collect_anov_p)`))

# df_anov_pval_corr = list of corrected p-values for each anova per ROI 


# run 1 large ANOVA with main effects: ROI & Time 
m2 = afex::aov_ez("subj", "t_stat", 
                  df_uber, 
                  within = c("ROI", "model"))
summary(m2)

  



