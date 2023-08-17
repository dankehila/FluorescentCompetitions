setwd(dir = "cleaned_datasets/")
library(tidyverse);library(broom);library(segmented)


dat = read.csv('monocultures_cleaned_all_glucose.csv',row.names = 1)

cfp_mat = 0.29; ofp_mat=0.63;





# Segmented dynamic calibration ---------------------------------------------------

segreg <- function(df){
  exp(segmented(lm(OD600~RFU, data=df), seg.Z = ~RFU, npsi=1)$psi[2]) # getting the breakpoint from the segmented model
}



calib = dat %>% filter(Hours <15) %>% 
  filter(glucose > 0.05) %>% 
  group_by(Colour,Well,day) %>% 
  arrange(Colour,Well,day,Hours) %>%
  mutate(RFU =  case_when(
    Colour == "T" & Hours > first(Hours) ~ lead(RFU,n=1)*(1-cfp_mat) + lead(RFU,n=2)*cfp_mat,
    Colour == "O" &  Hours > first(Hours)  ~ lead(RFU,n=1)*(1-ofp_mat) + lead(RFU,n=2)*ofp_mat,
    #Colour == "O" &  Hours > first(Hours)  ~ RFU*(1-ofp_mat) + lead(RFU,n=1)*ofp_mat,
    T ~ as.numeric(RFU))) %>%
  filter(OD600>0 & RFU>0) %>% # about 500 datapoints below 0 due to background subtraction 
  # removed so we can log-log regress, which is favorable because it constrains to positive values
  # and also reflects well the appearance of spread which is constant on the log-log scale
  # this can be seen in the visualization section above
  mutate(OD600=log(OD600),RFU=log(RFU))%>% 
  group_by(Colour)  %>% # removing well / day information and fitting all days and wells as a single regression
  dplyr::select(c(RFU,OD600)) %>% nest %>% 
  mutate(breakpoint = data %>% map_dbl(segreg), # piecewise (segmented) regression using the function above
         filtered_data = map2(data,breakpoint,
                             ~filter(.x,RFU > as.numeric(log(.y)))), # subsetting filtered data for linear model
         lm = map(filtered_data,~lm(.x$OD600~.x$RFU)), # calibration fit for filtered data,
         sigma = map_dbl(lm,function(x) summary(x)$sigma),
         lm = map(lm,tidy)
         ) %>%  # obtaining broom::tidy summary of linear model
  unnest(lm) %>% subset(select = c(Colour,breakpoint,term,estimate,sigma)) %>% # keeping only slope and intercept info
  mutate(term = recode(term, '(Intercept)' = 'Intercept', '.x$RFU' = 'Slope'),
         Colour = paste(ifelse(Colour=='T','C',Colour),'FP',sep="")) %>% 
  pivot_wider(names_from=term,values_from=estimate) %>% as.data.frame


calib_unshifted = dat %>% filter(Hours <15) %>% 
  filter(glucose > 0.05) %>% 
  group_by(Colour,Well,day) %>% 
  arrange(Colour,Well,day,Hours) %>%
  filter(OD600>0 & RFU>0) %>% # about 500 datapoints below 0 due to background subtraction 
  # removed so we can log-log regress, which is favorable because it constrains to positive values
  # and also reflects well the appearance of spread which is constant on the log-log scale
  # this can be seen in the visualization section above
  mutate(OD600=log(OD600),RFU=log(RFU))%>% 
  group_by(Colour)  %>% # removing well / day information and fitting all days and wells as a single regression
  dplyr::select(c(RFU,OD600)) %>% nest %>% 
  mutate(breakpoint = data %>% map_dbl(segreg), # piecewise (segmented) regression using the function above
         filtered_data = map2(data,breakpoint,
                              ~filter(.x,RFU > as.numeric(log(.y)))), # subsetting filtered data for linear model
         lm = map(filtered_data,~lm(.x$OD600~.x$RFU)), # calibration fit for filtered data,
         sigma = map_dbl(lm,function(x) summary(x)$sigma),
         lm = map(lm,tidy)) %>% 
  unnest(lm) %>% subset(select = c(Colour,breakpoint,term,estimate,sigma)) %>% # keeping only slope and intercept info
  mutate(term = recode(term, '(Intercept)' = 'Intercept', '.x$RFU' = 'Slope')) %>% 
  pivot_wider(names_from=term,values_from=estimate) %>% as.data.frame

# generating predicted values from log-log calibration to plot a line over the points
calib_preds = tibble(Colour = c('GFP','OFP','CFP'),
                     max=dat %>% group_by(Colour) %>% dplyr::summarize(max=max(RFU)) %>% pull(max)) %>% 
  full_join(calib) %>% 
  mutate(RFU = map2(breakpoint,max,~seq(.x,.y,length.out=1000))) %>% unnest(RFU) %>% 
  mutate(OD600 = exp(Intercept)*RFU^Slope)

# showing the results of calibration
dat %>% #filter(Hours <15) %>% 
  group_by(Colour,Well,day) %>% #filter(glucose == 1) %>% 
  arrange(Colour,Well,day,Hours) %>%
  mutate(RFU =  case_when(
    Colour == "T" & Hours > first(Hours) ~ lead(RFU,n=1)*(1-cfp_mat) + lead(RFU,n=2)*cfp_mat,
    Colour == "O" &  Hours > first(Hours)  ~ lead(RFU,n=1)*(1-ofp_mat) + lead(RFU,n=2)*ofp_mat,
    T ~ as.numeric(RFU)),
    Colour = paste(ifelse(Colour=='T','C',Colour),'FP',sep="")) %>%   filter(OD600>0 & RFU>0) %>%
  ggplot(aes(RFU,OD600,color=as.factor(Genotype))) +facet_wrap(~Colour,scales="free") + theme_classic()+
  labs(y=expression("OD"[600]),color='Plasmid',shape="")+
  geom_point(aes(shape = day),alpha=1) + scale_x_log10() + scale_y_log10()+
  geom_vline(data=calib,
             mapping=aes(xintercept=breakpoint),linetype='dashed',color='darkred',linewidth=1) +
  geom_line(data=calib_preds,color='grey15',linewidth=1) +
  scale_colour_grey(start=0.3,end=0.8) +guides(shape = "none") 


# detection limit using residual standard error and background values --------

# even though all the RFU and OD values are background subtracted
# We still need to use background values in order to compute a new limit of detection
# if extrapolating further below the OD limits


max_bg = rbind(
  read.csv("../make_datasets/220209_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% 
  filter(Colour == "blank") %>% subset(select = c(Read,value)),
  read.csv("../make_datasets/220216_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% 
  filter(Colour == "blank") %>% subset(select = c(Read,value)),
  read.csv("../make_datasets/220223_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(grepl("([1][7-9]|[2][0-4])",Well))  %>% subset(select = -Well) %>% gather(key='Read'),
  read.csv("../make_datasets/220224_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(grepl("([1][7-9]|[2][0-4])",Well)) %>% subset(select = -Well) %>% gather(key='Read'),
  read.csv("../make_datasets/220301_compete.csv",header=T) %>%   pivot_longer(-c(Read,Time),names_to = "Well") %>% filter(grepl("A8|A9|A10|A11|A12|A13",Well))  %>% 
  filter(Time ==min(Time)) %>% subset(select = c(Read,value))) %>% group_by(Read) %>% dplyr::summarize(bg=mean(value)) %>% 
  filter(Read != 'OD600') %>% rename('Colour' = 'Read')

calib = calib %>% full_join(max_bg) %>% 
  mutate(detection.lim = bg/exp(-sigma),detection.lim.odeq = exp(Intercept)*detection.lim^Slope)
calib_unshifted = calib_unshifted %>% full_join(max_bg) %>% 
  mutate(detection.lim = bg/exp(-sigma),detection.lim.odeq = exp(Intercept)*detection.lim^Slope)

calib_preds_extrapolated =  tibble(Colour = c('GFP','OFP','CFP'),
                       max=dat %>% group_by(Colour) %>% dplyr::summarize(max=max(RFU)) %>% pull(max)) %>% 
  full_join(calib) %>% 
  mutate(RFU = map2(detection.lim,max,~seq(.x,.y,length.out=1000))) %>% unnest(RFU) %>% 
  mutate(OD600 = exp(Intercept)*RFU^Slope)

# showing the results of calibration with sigma-derived detection limit in orange
dat %>% 
  group_by(Colour,Well,day) %>% filter(glucose > 0.05) %>% filter(Hours<15) %>% 
  arrange(Colour,Well,day,Hours) %>%
  mutate(RFU =  case_when(
    Colour == "T" & Hours > first(Hours) ~ lead(RFU,n=1)*(1-cfp_mat) + lead(RFU,n=2)*cfp_mat,
    Colour == "O" &  Hours > first(Hours)  ~ lead(RFU,n=1)*(1-ofp_mat) + lead(RFU,n=2)*ofp_mat,
    T ~ as.numeric(RFU)),
    Colour = paste(ifelse(Colour=='T','C',Colour),'FP',sep="")) %>%   filter(OD600>0 & RFU>0) %>% 
  ggplot(aes(RFU,OD600)) +facet_wrap(~Colour) + 
  geom_point(aes(color=as.factor(Genotype), shape = day)) + 
 scale_x_log10() + scale_y_log10(limits = c(0.02,1.2),
                                 breaks = c(1,0.5,0.1,0.05))+
  scale_colour_grey(start=0.3,end=0.8) +guides(shape = "none") +theme_classic()+
  geom_line(data=calib_preds_extrapolated,color='#004D40',linewidth=1.5,linetype='11') +
  geom_line(data=calib_preds,color='#004D40',linewidth=1.5) +
  labs(y=expression("OD"[600]),color='Plasmid')+
  geom_vline(data=calib,mapping=aes(xintercept=breakpoint),linetype='dashed',color='darkred',linewidth=1) +
  geom_vline(data=calib,mapping=aes(xintercept=detection.lim),color='darkorange1',linewidth=1,linetype='longdash') 

write.csv(calib,'../cleaned_datasets/log_log_calib_all_glucose_shifted.csv')

# Comparing calibrations on competition data --------------------------------------------------
competition_data = read.csv('raw_competition_dat_methodspaper.csv',row.names = 1)


shifted_competedat = competition_data %>% separate(pair,into=c('col1','col2'),remove=F,sep='v') %>% 
  rename('gen1'='genA','gen2'='genB') %>% 
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
  mutate(RFU1 = case_when(
    pair %in% c('CvG','CvO') ~ (1-0.29)*lead(RFUA,n=1)+0.29*lead(RFUA,n=1),
    pair == 'OvG'  ~ (1-0.63)*(lead(RFUA,n=1) + 0.63*lead(RFUA,n=2)),
    T ~ RFUA),
    RFU2 = case_when(
      pair == 'CvO'  ~ (1-0.63)*(lead(RFUB,n=1) + 0.63*lead(RFUB,n=2)),
      T ~ RFUB)) %>%
  na.omit %>% 
  mutate(ODeq1 = case_when(
    col1=='C'~ exp(calib[3,'Intercept'])*RFU1^calib[3,'Slope'], 
    col1=='O'~ exp(calib[2,'Intercept'])*RFU1^calib[2,'Slope']),
    ODeq2 = case_when(
      col2=='G'~exp(calib[1,'Intercept'])*RFU2^calib[1,'Slope']
      ,col2=='O'~exp(calib[2,'Intercept'])*RFU2^calib[2,'Slope'])) %>% ungroup


calib_static = read.csv('static_calibration_curve.csv',row.names = 1)


shifted_competedat_static = competition_data %>% separate(pair,into=c('col1','col2'),remove=F,sep='v') %>% 
  rename('gen1'='genA','gen2'='genB') %>% 
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
  mutate(RFU1 = case_when(
    pair %in% c('CvG','CvO') ~ lead(RFUA,n=1),
    pair == 'OvG'  ~ 0.5*(lead(RFUA,n=1) + lead(RFUA,n=2)),
    T ~ RFUA),
    RFU2 = case_when(
      pair == 'CvO'  ~ 0.5*(lead(RFUB,n=1) + lead(RFUB,n=2)),
      T ~ RFUB)) %>%
  na.omit %>% 
  mutate(ODeq1 = case_when(
    col1=='C'~ exp(calib_static[1,'Intercept'])*RFU1^calib_static[1,'Slope'], 
    col1=='O'~ exp(calib_static[3,'Intercept'])*RFU1^calib_static[3,'Slope']),
    ODeq2 = case_when(
      col2=='G'~exp(calib_static[2,'Intercept'])*RFU2^calib_static[2,'Slope']
      ,col2=='O'~exp(calib_static[3,'Intercept'])*RFU2^calib_static[3,'Slope'])) %>% ungroup


## supplementary figure
full_join(shifted_competedat %>% mutate(calib='dynamic'),
          shifted_competedat_static %>% mutate(calib='static')) %>% 
  group_by(Well,day) %>%  arrange(Well,day,Hours) %>%
  mutate(sumeq = ODeq1+ODeq2) %>% mutate(dilution=paste(dilution,'X',sep='')) %>% 
  # removing three outlier wells for presentation
  filter(Hours<16 & 
           !(day == "Jan27" & Well == "B13") & !(day == 'Jan19' & Well %in%  c("P1","P9"))) %>% 
  ggplot(aes(x=sumeq,y=OD600,color=calib)) + #facet_grid(fct_rev(dilution)~.) +
  geom_point(alpha=0.5) +  theme_classic() +
  scale_color_manual(values = c("#DC990A","#a9a9a9")) +
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey20",size=1) +
  labs(y="OD600", x ="Sum ODeq",color=NULL) +
  theme(legend.key = element_rect(fill = NA, color = NA),
        text=element_text(size=rel(5),family="Helvetica",color="grey20"),
        legend.position="bottom", 
        legend.text = element_text(size=rel(4),family="Helvetica",color="grey20"),
        strip.text.x=element_text(size=rel(4),family="Helvetica",color="grey20"))+
  # scale_x_log10() + scale_y_log10()+
  guides(color=guide_legend(override.aes=list(size=5)))

# calculating mean absolute error

full_join(shifted_competedat %>% mutate(calib='dynamic'),
          shifted_competedat_static %>% mutate(calib='static')) %>% 
  group_by(Well,day) %>%  arrange(Well,day,Hours) %>%
  mutate(sumeq = ODeq1+ODeq2) %>% group_by(calib) %>% na.omit %>% 
  # removing the three outliers for calculating mean absolute error does not affect results much
  # (dynamic is still better)
  filter(!(day == "Jan27" & Well == "B13") & !(day == 'Jan19' & Well %in%  c("P1","P9"))) %>% 
  # Optional but can also see how MAE changes when looking after hour 2 and before hour 12
  # filtering for above the optical density detection limit but keeping the timeframe within growth 
  # (after 12 hours they're mostly saturated)
  # again this does not affect results qualitatively
  filter(Hours < 16) %>% 
  dplyr::summarise(MAE = mean(abs(sumeq-OD600)))



# Showing maturation correction effect ------------------------------------

unshifted_competedat = competition_data %>% separate(pair,into=c('col1','col2'),remove=F,sep='v') %>%
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% mutate(RFU1=RFUA,RFU2=RFUB) %>%
  na.omit %>% 
  mutate(ODeq1 = case_when(
    col1=='C'~ exp(calib_unshifted[3,'Intercept'])*RFU1^calib_unshifted[3,'Slope'], 
    col1=='O'~ exp(calib_unshifted[2,'Intercept'])*RFU1^calib_unshifted[2,'Slope']),
    ODeq2 = case_when(
      col2=='G'~exp(calib_unshifted[1,'Intercept'])*RFU2^calib_unshifted[1,'Slope'],
      col2=='O'~exp(calib_unshifted[2,'Intercept'])*RFU2^calib_unshifted[2,'Slope'])) %>% ungroup


# For figure 3
unshifted_competedat %>% rename('gen1'='genA','gen2'='genB') %>% 
  filter(pair == 'OvG' & gen1 == gen2) %>% filter(Hours<15) %>% 
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
  mutate(ODeq1=ODeq1/first(ODeq1),ODeq2=ODeq2/first(ODeq2),freq=log(ODeq1/ODeq2)) %>% 
  ggplot(aes(Hours,freq,group=paste(Well,day))) + geom_smooth(se=F) + facet_grid(dilution~day)+
  geom_smooth(data=shifted_competedat %>% rename('gen1'='genA','gen2'='genB') %>% 
                filter(pair == 'OvG' & gen1 == gen2) %>% filter(Hours<15) %>% 
                group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
                mutate(ODeq1=ODeq1/first(ODeq1),ODeq2=ODeq2/first(ODeq2),freq=log(ODeq1/ODeq2)),
              color='pink3',se=F) +
  geom_point(aes(y=0.5*log(OD600)),color='black') +
  geom_point(size=1,color='blue3') +
  geom_point(data=shifted_competedat %>% rename('gen1'='genA','gen2'='genB') %>% 
               filter(pair == 'OvG' & gen1 == gen2) %>% filter(Hours<15) %>% 
               group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
               mutate(ODeq1=ODeq1/first(ODeq1),ODeq2=ODeq2/first(ODeq2),freq=log(ODeq1/ODeq2)),color='pink3',size=1)

fig3_unshifted = unshifted_competedat %>% rename('gen1'='genA','gen2'='genB') %>% 
  filter(pair == 'OvG' & gen1 == gen2) %>% filter(Hours<16) %>% 
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
  # normalizing by first timepoint so that the freq curve starts at 0,0 for visualization ease
  # slope doesn't change due as log division becomes subtraction so division by first point only affects intercept
  mutate(ODeq1=ODeq1/first(ODeq1),ODeq2=ODeq2/first(ODeq2),freq=log(ODeq1/ODeq2)) 
fig3_shifted =  shifted_competedat %>% rename('gen1'='genA','gen2'='genB') %>% 
  filter(pair == 'OvG' & gen1 == gen2) %>% filter(Hours<16) %>% 
  group_by(Well,day) %>% arrange(Well,day,Hours) %>% 
  mutate(ODeq1=ODeq1/first(ODeq1),ODeq2=ODeq2/first(ODeq2),freq=log(ODeq1/ODeq2)) 

rbind(fig3_unshifted %>% mutate(shifting = 'unshifted'),
      fig3_shifted %>% mutate(shifting = 'maturation - shifted')) %>% 
  ggplot(aes(Hours,freq,group=paste(Well,day))) + 
  stat_smooth(se=F,method='loess',span=0.5,method.args = list(degree=1)) + facet_grid(~fct_rev(shifting))+
  geom_point(aes(y=0.5*log(OD600)-2),color='black') +
  geom_point(size=1,color='blue3')

rbind(fig3_unshifted %>% mutate(shifting = 'unshifted'),
      fig3_shifted %>% mutate(shifting = 'maturation - shifted')) %>% 
  ggplot(aes(Hours,freq,group=paste(Well,day))) + geom_line() + facet_grid(~fct_rev(shifting))+
  geom_point(aes(y=0.5*log(OD600)-2),color='black') +
  geom_point(size=1,color='blue3')





write.csv(shifted_competedat,'competition_data_calibrated_shifted.csv')
