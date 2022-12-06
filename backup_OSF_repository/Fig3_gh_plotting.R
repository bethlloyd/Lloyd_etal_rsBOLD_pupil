rm(list= ls(all.names = TRUE))
gc()

library(ggplot2)
library(dplyr)
library(Hmisc)
library(tidyverse)
install.packages( c( "dplyr", "ggplot2", "plyr" ) )
# ggplot theme
mytheme <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(size=1, colour = "black"),
                 text = element_text(size=34,colour = "black"),
                 strip.background = element_blank(),
                 axis.title.x=element_blank(), 
                 #axis.text.x=element_blank(), 
                 axis.ticks.x=element_blank(),
                 legend.title = element_blank())

dat_path <- paste("D:\\NYU_RS_LC\\stats\\OSF_repository\\Figure3-SourceData")

group.colors =c(LC='#00C1AA', VTA='#ED8141', SN='#9590FF', DR='#FF62BC', MR='#5BB300', BF='goldenrod1', ACC='tan3', OCC='grey', Pons='#a2bffe')

#Fig3_G = pupil size, Fig3_H = pupil derivative 
pup_type = c('Fig3_G', 'Fig3_H')

# for p vals 
p_val = list(numeric())
uber_collect_p <-list()
uber_collect_p <- append(uber_collect_p, list(p_val))


for (c_pup in pup_type){
  

  #make output dir
  output_plots_dir = paste("D:\\NYU_RS_LC\\stats", "OSF_fig_output", sep='\\')
  
  # load in data and rename column
  conmap_mod1 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_DR.csv', 
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod1)[names(conmap_mod1) == 'spmT_0001.nii'] <- 'DR'
  
  conmap_mod2 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_MR.csv', 
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod2)[names(conmap_mod2) == 'spmT_0001.nii'] <- 'MR'
  
  conmap_mod3 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_VTA.csv',  
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod3)[names(conmap_mod3) == 'spmT_0001.nii'] <- 'VTA'
  
  conmap_mod4 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_LC.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod4)[names(conmap_mod4) == 'spmT_0001.nii'] <- 'LC'
  
  conmap_mod5 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_SN.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod5)[names(conmap_mod5) == 'spmT_0001.nii'] <- 'SN'
  
  conmap_mod6 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_ACC.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod6)[names(conmap_mod6) == 'spmT_0001.nii'] <- 'ACC'
  
  conmap_mod7 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_OCC.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod7)[names(conmap_mod7) == 'spmT_0001.nii'] <- 'OCC'
  
  conmap_mod8 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_BF.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod8)[names(conmap_mod8) == 'spmT_0001.nii'] <- 'BF'
  
  conmap_mod9 <- read.delim(paste(dat_path, c_pup, 'con_stat_rsHRF_ana_PONS.csv',
                                  sep = '\\'), header = TRUE, sep = ',')
  names(conmap_mod9)[names(conmap_mod9) == 'spmT_0001.nii'] <- 'Pons'
  
  dat_all <- left_join(conmap_mod1, conmap_mod2, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod3, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod4, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod5, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod6, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod7, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod8, by = 'subj')
  dat_all <- left_join(dat_all, conmap_mod9, by = 'subj')
  
  
  # get data into long format 
  labels = names(dat_all)
  # Select variables
  var_names = c("subj",
                "DR",
                "MR", 
                "VTA", 
                "LC", 
                "SN", 
                "ACC",
                "OCC",
                "BF",
                "Pons")
  
  dat = dat_all[,match(var_names,labels)]
  
  # Make data frame
  df_dat = data.frame(dat)
  
  
  
  # Make long data format [ csp_preTlapse:csm_stim_mean ]
  long_dat = tidyr::gather(df_dat, value=Val, key=Names, DR:Pons, factor_key=TRUE)
  names(long_dat)[names(long_dat) == 'Names'] <- 'ROI'
  names(long_dat)[names(long_dat) == 'Val'] <- 't_stat'
  
  long_dat_fin <- subset(long_dat, !is.na(t_stat))

  # run t-tests on each roi (then correct for multiple comparisons)
  i = 1
  print(paste(c_pup))
  for (c_mod in tail(var_names, -1)) { 
    
    dat_model = subset(long_dat_fin, ROI == c_mod)
    
    stat.test_p1 <- 
      t.test(dat_model$t_stat,
             mu = 0) # Remove details
    
    #t_apa(stat.test_p1)
    
    print(paste('p value for t-test model roi ', c_mod, 'is', t_apa(stat.test_p1), sep =" "))
    p_val$p_value <- as.numeric(as.character(unlist(stat.test_p1[3][1])))
    print(unlist(stat.test_p1[3][1]))
    
    uber_collect_p <- append(uber_collect_p, unlist(stat.test_p1[3][1]))
    i <- i + 1
  } 
  
  # plot data for each ROI 
  plot1 <- long_dat_fin %>%
    mutate(ROI = fct_relevel(ROI, "LC", "VTA", "SN", "DR", "MR", "BF", "ACC", "OCC", "Pons")) %>%
    ggplot(aes(x = ROI, y = t_stat, fill = ROI)) +
    #geom_point(aes(color=ROI), size=2, alpha=.4) +
    geom_jitter(width = 0.1, aes(color=ROI), size=0.8) + 
    #geom_line(aes(color=session, group=subject), size=.3, color="black", alpha=.4) +
    scale_y_continuous(limits=c(-3.5,3.5)) +
    #scale_y_continuous(limits=c((min(long_dat_fin$con_stat)*1.1), (max(long_dat_fin$con_stat)*1.1))) +
    stat_summary(fun="mean",colour = "black", size = 2,geom = "point") +
    stat_summary(fun.data="mean_se",colour = "black", width=0,size = 1,geom = "errorbar") +
    scale_color_manual(values=group.colors) +
    mytheme                 +                                 
    labs(title = "", x = "", y = "t-values")   +
    theme(legend.position = "none") +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_violin(alpha = 0, aes(color='#000000', fill = ROI), size = 0.2)
  
  print(plot1)
  #ggsave(paste(output_plots_dir, paste("roi_pup_deriv_sig_extract.eps",  sep = ''), sep='\\'),width = 9, height = 8, plot1)
   
  
}

# correct for multple comparisons - FDR correction
var_names
options(scipen=999) 
corr_p_val  = unlist(uber_collect_p)
p.adjust(corr_p_val, method = "fdr", n = length(corr_p_val))
