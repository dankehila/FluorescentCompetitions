
# Load packages and raw data ----------------------------------------------
library(cmdstanr);library(posterior);library(tidyverse);library(bayesplot);library(ggdist);library(deSolve)
options(mc.cores = 8)
# competition data that has been corrected for differential maturation
competition_data_raw = read.csv('cleaned_datasets/competition_data_calibrated_shifted.csv',row.names = 1)


# seems like 0.001 is a good replacement for when ODeq1 or ODeq2 are zero -- necessary for lognormal fitting (positive support only)
# competition_data_raw %>% group_by(col1) %>% filter(ODeq1 > 0) %>% filter(ODeq1 == min(ODeq1))
# competition_data_raw %>% group_by(col2) %>% filter(ODeq2 > 0) %>% filter(ODeq2 == min(ODeq2))
# 
# competition_data_raw %>%  replace(is.na(.),0) %>% group_by(Well,day) %>% filter(any(ODeq1==0)) %>% mutate(ODeq1f = ifelse(ODeq1==0,0.001,ODeq1)) %>% filter(Hours<4) %>% 
#   ggplot(aes(Hours)) + geom_point(aes(y=ODeq1),color='grey30') + geom_point(aes(y=ODeq1f),color='darkred')
# competition_data_raw %>%  replace(is.na(.),0) %>% group_by(Well,day) %>% filter(any(ODeq2==0)) %>% mutate(ODeq2f = ifelse(ODeq2==0,0.001,ODeq2)) %>% filter(Hours<4) %>% 
#   ggplot(aes(Hours)) + geom_line(aes(y=ODeq2,group=paste(Well,day)),color='grey30') + geom_line(aes(y=ODeq2f,group=paste(Well,day)),color='darkred')

# replacing zeroes with 0.001 and filtering to only the first 16 timepoints -- cultures saturate around 10-12 hours.
competition_data =   competition_data_raw %>%  mutate(ID = as.numeric(as.factor(paste(Well,day)))) %>%
  # for lognormal fits: can't have OD at 0 so switching to small value justified above
  replace(is.na(.),0.001) %>% 
  # mutate(ODeq1=ifelse(ODeq1<=0,0.001,ODeq1),ODeq2=ifelse(ODeq2<=0,0.001,ODeq2)) %>% 
  filter(Hours < 16) %>% arrange(ID,Hours) #%>% filter(ID %in% sample(1200,30)) %>% mutate(ID=as.numeric(as.factor(ID)))

# high y0 --> some contaminating source of fluorescence that muddles the signal
high_y0_ID = competition_data %>% group_by(Well,day) %>% filter(Hours==min(Hours)) %>% 
  pivot_longer(c(ODeq1,ODeq2)) %>% filter(value>0.2) %>% pull(ID) %>% unique

# most of the time it's green except once orange fluorescence wih this issue
# competition_data %>% filter(ID %in% high_y0_ID) %>% ggplot(aes(Hours)) +
#   geom_point(aes(y=ODeq1,color=col1)) + geom_point(aes(y=ODeq2,color=col2)) +
#   facet_wrap(~Well+day) + geom_line(aes(y=OD600,),color='grey30') + scale_y_log10() + scale_color_manual(values = c('darkturquoise','darkgreen','darkorange')) + theme_classic()
# 


competition_data = competition_data %>% filter(!(ID %in% high_y0_ID)) %>% mutate(ID = as.numeric(as.factor(paste(Well,day)))) 

# 
# ## all growth curves visualized (supplementary figure)
# competition_data %>% ggplot(aes(Hours)) + 
#   geom_line(aes(y=ODeq1,group=paste(Well,day),color='competitor 1'),alpha=0.1) +
#   geom_line(aes(y=ODeq2,group=paste(Well,day),color='competitor 2'),alpha=0.1) +
#   scale_y_log10() + theme_classic() + scale_color_manual(values = c("#336CA1","#E2B89B")) +
#   labs(y='OD equivalents',color=NULL) + guides(color=guide_legend(override.aes=list(linewidth=5,alpha=1))) +
#   theme(legend.position = 'bottom')



# Generate dataset --------------------------------------------------------

# for y0 list expect 4 * 3 levels = 12

# for y0 list expect 4 days * 3 colors * 10 genotypes  = 120 levels
ID_y0_levels = expand.grid(day = 1:4,col = 1:3,gen = 1:10) %>% mutate(ID_y0 = (paste(day,col,gen))) %>% 
  pull(ID_y0)

design = competition_data %>% 
  subset(select=c(Well,day,col1,col2,gen1,gen2,ID,dilution)) %>% unique %>% 
  mutate(
    gen1=factor(gen1,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                "VIMR6.4","VIMR18.1","NDMR18")),
    gen2=factor(gen2,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                "VIMR6.4","VIMR18.1","NDMR18")),
    col1=factor(col1, levels = c("C","G","O")),col2=factor(col2, levels = c("C","G","O")),
    dilution = factor(dilution, levels = c(500,1250)),
    day = factor(day, levels = c('Jan11','Jan19','Jan26','Jan27')),
    ID_y10 = as.numeric(factor(paste(as.numeric(day),as.numeric(col1),as.numeric(gen1)),levels = ID_y0_levels)), 
    ID_y20 = as.numeric(factor(paste(as.numeric(day),as.numeric(col2),as.numeric(gen2)),levels = ID_y0_levels)),
    selcoef = rnorm(n())) # dummy variable for constructing model matrices only




modmat1 = model.matrix(selcoef ~ col1*gen1,data=design)
modmat2 = model.matrix(selcoef ~ col2*gen2,data=design)

#checking model matrix column names are the same
identical(gsub('gen1','gen2',(gsub('col1','col2',colnames(modmat1)))),colnames(modmat2))

#model matrix for predicting individual genotype growth rates/selective advantages
newdata =expand.grid(gen1= c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                             "VIMR6.4","VIMR18.1","NDMR18"),
                     col1 = c("C","G","O")) %>% 
  mutate(
    gen1=factor(gen1,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                "VIMR6.4","VIMR18.1","NDMR18")),
    col1=factor(col1, levels = c("C","G","O")),
    estimate_corr=rnorm(n()))

newdata_modmat = model.matrix(lm(estimate_corr ~ col1*gen1,
                                 data=newdata
))
#no epistasis model matrix
noepimodmat = newdata_modmat
# for no epistasis model matrix, setting all epistatic terms (13:30) to zero so it can be used for posterior prediction
noepimodmat[,-c(1:12)] = 0



#llg stands for lognormal, logistic growth. This model has only one K
llg_mod = cmdstan_model('stan models/one_K_model.stan')
#in case you want to run LOO, you may have to run these with force_recompile = T
llg_mod_multK = cmdstan_model('stan models/multilevel_K_model.stan')



# initial values for y0 -- rough guesses help sampler
med_y0 = competition_data %>% full_join(design %>% mutate(day = as.character(day)) %>% subset(select=c(Well,day,ID_y10,ID_y20))) %>% 
  pivot_longer(c(ODeq1,ODeq2)) %>% group_by(Well,day) %>% filter(Hours==min(Hours)) %>% 
  mutate(ID_y0 = ifelse(name == 'ODeq1',ID_y10,ID_y20)) %>% group_by(ID_y0) %>% dplyr::summarize(med=median(log(value))) %>% 
  pull(med)

#initial values for selective differences
slope_s_init = competition_data %>%filter(Hours<6) %>% group_by(Well,day) %>% nest %>% 
  mutate(s1 = map_dbl(data,~coef(lm(log(ODeq1)~Hours,data=.x))[2]),
         s2 = map_dbl(data,~coef(lm(log(ODeq2)~Hours,data=.x))[2])) %>% ungroup %>% 
  mutate(means = mean(s1+s2)/2,
         s1 = s1-means,s2=s2-means) %>% 
  subset(select=c(s1,s2)) %>% as.matrix
# initial value for differential costs
dc_init = (replace_na(coef(lm(slope_s_init[,1] ~0+ modmat1)),0) +
             replace_na(coef(lm(slope_s_init[,2] ~0+ modmat2)),0))/2
# initial value for K_array
K_array_init = competition_data %>% mutate(eqsum = ODeq1+ODeq2) %>% group_by(Well,day,ID) %>% filter(eqsum>quantile(eqsum,0.99)) %>% 
  dplyr::summarize(meanK=median(eqsum)) %>% pull(meanK) %>% log

# rough values of r and sd to inform priors
competition_data %>% filter(Hours<8) %>% group_by(Well,day) %>% 
  nest %>% mutate(tid=map(data,~broom::tidy(lm(log(OD600)~Hours,data=.x)))) %>% ungroup %>% 
  unnest(tid) %>% filter(term == 'Hours') %>% dplyr::summarize(r = mean(estimate),rsd = sd(estimate))

#same for K
competition_data %>% mutate(eqsum = ODeq1+ODeq2) %>% filter(eqsum>quantile(eqsum,0.99)) %>% 
  dplyr::summarize(mean_K = log(mean(eqsum)),sdK = (sd(log(eqsum))))




stan_data_y0 = list(
  N = nrow(competition_data),
  N_ID = nrow(design),
  IDX= competition_data$ID,
  N_y0 = length(ID_y0_levels),
  ID_y10 = design$ID_y10,
  ID_y20 = design$ID_y20,
  t0=0,
  yobs = as.matrix(competition_data[,c("ODeq1","ODeq2")]),
  ts = competition_data$Hours,
  culturing_time = 16,
  N_X = ncol(modmat1),
  X1 = modmat1,
  X2 = modmat2,
  N_pred = nrow(newdata_modmat),
  X_pred = newdata_modmat,
  X_pred_no_epistasis = noepimodmat,
  prior_r = c(0.5,0.25),
  prior_K = c(0.3,0.2),
  prior_y0 = c(median(med_y0),1.5),
  epsilon = 1e-4,
  num_steps = 1e6
  
)


# 
init_fun= function() list(log_y0 = med_y0,
                          K=0.3,sigma=0.1,sigma_dc=0.1,
                          #rbar = 0.5,
                          r = 0.5+slope_s_init,
                          dc = dc_init,
                          nu=6,L=matrix(c(1,0.5,0.5,1),ncol=2))



llg_fit =  llg_mod$sample(
  data = stan_data_y0,chains=4,iter_sampling = 1000,iter_warmup = 1000,
  init = init_fun,
)



init_fun_multK= function() list(log_y0 = med_y0,
                                K=0.3,sigma=0.1,sigma_dc=0.1,
                                r = 0.5+slope_s_init,
                                dc = dc_init,nu=6,L=matrix(c(1,0.5,0.5,1),ncol=2),
                                K_array = K_array_init,sigma_K = 2*sd(K_array_init)
)



llg_fit_multK = llg_mod_multK$sample(
  data = stan_data_y0,chains=4,iter_warmup = 1000,iter_sampling = 1000,
  init = init_fun_multK
)

#llg_fit$save_output_files(dir = 'full_models/full model outputs/20230712 antibiotic final/',basename = '230712_onek')
#llg_fit_multK$save_output_files(dir = 'full_models/full model outputs/20230712 antibiotic final/', basename = '230712_multk')


dir = 'full_models/full model outputs/20230712 antibiotic final/'
llg_fit = as_cmdstan_fit(files = paste(dir, list.files(dir,pattern = '230712_onek.*'),sep=""))
llg_fit_multK = as_cmdstan_fit(files = paste(dir, list.files(dir,pattern = '230712_multk.*'),sep=""))

# LOO CV ------------------------------------------------------------------

onek_loo = llg_fit$loo(variables = 'log_lik_ODE',cores=6,save_psis=T)
multk_loo = llg_fit_multK$loo(variables = 'log_lik_ODE',cores=6,save_psis=T)

print(loo_compare(onek_loo,multk_loo),simplify=F)



competition_data %>% 
  mutate(onek_outliers_1 = onek_loo$pointwise[1:nrow(competition_data),"influence_pareto_k"],
         onek_outliers_2 = onek_loo$pointwise[nrow(competition_data) + (1:nrow(competition_data)),"influence_pareto_k"],
         multk_outliers_1 = multk_loo$pointwise[1:nrow(competition_data),"influence_pareto_k"],
         multk_outliers_2 = multk_loo$pointwise[nrow(competition_data) + (1:nrow(competition_data)),"influence_pareto_k"]
  ) %>% 
  group_by(ID) %>% 
  filter(any(
    onek_outliers_1 > 0.7 |
      onek_outliers_2 > 0.7 |
      multk_outliers_1 > 0.7 |
      multk_outliers_2 > 0.7
  )) %>% 
  mutate(
    gen1=as.numeric(factor(gen1,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                           "VIMR6.4","VIMR18.1","NDMR18"))),
    gen2=as.numeric(factor(gen2,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                           "VIMR6.4","VIMR18.1","NDMR18")))
  ) %>% 
  ggplot(aes(Hours)) + geom_line(aes(y=ODeq1,color=paste(col1,"FP",sep=""))) + geom_line(aes(y=ODeq2,color=paste(col2,"FP",sep=""))) + 
  scale_color_manual(values = c('darkturquoise','darkgreen','darkorange')) +
  geom_point(data=. %>% ungroup %>% filter(onek_outliers_1>0.7),aes(y=ODeq1,shape='one K'),size=5) +
  geom_point(data=. %>% ungroup %>% filter(onek_outliers_2>0.7),aes(y=ODeq2,shape='one K'),size=5) +
  geom_point(data=. %>% ungroup %>% filter(multk_outliers_1>0.7),aes(y=ODeq1,shape='multiple K'),size=5) +
  geom_point(data=. %>% ungroup %>% filter(multk_outliers_2>0.7),aes(y=ODeq2,shape = 'multiple K'),size=5) +
  facet_wrap(Well~paste(paste(gen1,col1,"FP",sep=""),
                        paste(gen2,col2,"FP",sep=""),
                        sep=" vs. ")) + 
  scale_y_log10() + scale_shape_manual(values = c(5,8)) + theme_classic() +
  labs(y='ODeq',color = 'FP',shape=NULL)


outlier_points = competition_data %>% 
  mutate(onek_outliers_1 = onek_loo$pointwise[1:nrow(competition_data),"influence_pareto_k"],
         onek_outliers_2 = onek_loo$pointwise[nrow(competition_data) + (1:nrow(competition_data)),"influence_pareto_k"],
         multk_outliers_1 = multk_loo$pointwise[1:nrow(competition_data),"influence_pareto_k"],
         multk_outliers_2 = multk_loo$pointwise[nrow(competition_data) + (1:nrow(competition_data)),"influence_pareto_k"]
  ) %>% 
  filter(
    onek_outliers_1 > 0.7 |
      onek_outliers_2 > 0.7 |
      multk_outliers_1 > 0.7 |
      multk_outliers_2 > 0.7
  ) %>% mutate(outlier_point = paste(ID,Hours)) %>% pull(outlier_point)



# Posterior predictive checks ---------------------------------------------------------------------
#generating posterior predictive distributions
yrep_onek = t(apply(log(llg_fit$draws('mu',format='matrix')[1:50,]),1,
                    function(x) rlnorm(length(x),meanlog=x,
                                       sdlog =llg_fit$draws('sigma')[sample(1:4000,length(x),replace=T)])))
yrep_multk = t(apply(log(llg_fit_multK$draws('mu',format='matrix')[1:50,]),1,
                     function(x) rlnorm(length(x),meanlog=x,
                                        sdlog =llg_fit_multK$draws('sigma')[sample(1:4000,length(x),replace=T)])))

# general fit on the log scale of all data
#one K
ppc_dens_overlay(y = c(competition_data$ODeq1,competition_data$ODeq2),yrep = (yrep_onek)) + scale_x_log10()
#multilevel K
ppc_dens_overlay(  y = c(competition_data$ODeq1,competition_data$ODeq2),yrep = (yrep_multk)) + scale_x_log10()


relerror_modcompare_against_hours = 
  gather(as_tibble(sweep(yrep_onek,2,c(competition_data$ODeq1,competition_data$ODeq2))/yrep_onek),key = 'ypoint1',value='onek') %>% 
  cbind(
    gather(as_tibble(sweep(yrep_multk,2,c(competition_data$ODeq1,competition_data$ODeq2))/yrep_multk),key = 'ypoint2',value='multk')) %>% 
  mutate(Hours = rep(rep(competition_data$Hours,2),each=50),
         rep = rep(1:50,2*nrow(competition_data)),
         ypoint = parse_number(ypoint1),
         bin = cut(ypoint,floor(sqrt(nrow(competition_data))))) # binning number formula borrowed from bayesplot source code

yrep_against_hours = 
  gather(as_tibble(yrep_onek[,1:nrow(competition_data)]+yrep_onek[,nrow(competition_data)+(1:nrow(competition_data))]),key = 'ypoint1',value='onek') %>% 
  cbind(
    gather(as_tibble(yrep_multk[,1:nrow(competition_data)]+yrep_onek[,nrow(competition_data)+(1:nrow(competition_data))]),key = 'ypoint2',value='multk')) %>% 
  mutate(Hours = rep(competition_data$Hours,each=50),
         rep = rep(1:50,nrow(competition_data)),
         ypoint = parse_number(ypoint1),
         bin = cut(ypoint,floor(sqrt(nrow(competition_data))))) %>% # binning number formula borrowed from bayesplot source code
  group_by(Hours,rep,bin) %>%
  dplyr::summarize(onek = mean(onek),multk = mean(multk))


# Figure 4b
relerror_modcompare_against_hours %>% mutate(Hours=ceiling(Hours)) %>%  group_by(Hours,rep,bin) %>%
  dplyr::summarize(onek = mean(-onek),multk = mean(-multk)) %>% 
  ggplot(aes(Hours)) +# geom_point(aes(y=onek,group=paste(rep,bin),fill='one k'),alpha=0.1)
  stat_summary(aes(y=onek, color = 'one K'),fill='transparent',geom="ribbon", #linewidth=1,
               fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975),show.legend = F) +
  stat_summary(aes(y=multk, color = 'multilevel K'),fill='transparent',geom="ribbon", #linewidth=1,
               fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975),show.legend = F) +
  stat_ribbon(aes(y=onek,fill='one K'),alpha=0.2,.width=seq(0.05,0.95,.05)) +
  stat_ribbon(aes(y=multk,fill='multilevel K'),alpha=0.2,.width=seq(0.05,0.95,.05)) + 
  # stat_ribbon(aes(y=onek,fill='one K'),alpha=0.2,.width=.95) +
  # stat_ribbon(aes(y=multk,fill='multilevel K'),alpha=0.2,.width=.95) +
  scale_color_manual(values = c('#D75857','#1E88E5')) +scale_fill_manual(values = c('#D75857','#1E88E5')) +
  labs(y= expression(paste("E[",frac(y-yrep,y),"]")),fill=NULL) + geom_hline(yintercept = 0,linetype='24',color='grey30')+
  scale_y_continuous(breaks = c(-0.5,0,0.5),
                     sec.axis = sec_axis(trans = ~exp((.+0.3)*8),name = "Total ODeq",breaks = c(0,0.01,0.1,1))) + 
  geom_line(data=competition_data,aes(x=Hours,y=-0.3+log((ODeq1+ODeq2))/8,group=ID),alpha=0.01,color='grey30') +
  # geom_point(data = yrep_against_hours,aes(y=-0.5+log(onek)/8,color='one K'),alpha=0.003) +
  # geom_point(data = yrep_against_hours,aes(y=-0.5+log(multk)/8,color='multilevel K'),alpha=0.003) +
  theme_classic() + theme(text=element_text(family='Helvetica'), legend.position = 'bottom') +
  guides(color=guide_legend(override.aes=list(size=5,alpha=1))) 







# Refitting without outier points -----------------------------------------

#refitting without outliers will introduce uneven number of timepoints per competition, so have to use new model code
onek_unequal_timepoints_mod = cmdstan_model('stan models/one_K_unequal_times.stan')

#this requires outlier_points above at the LOO section, which in turn requires LOO
outlier_removed_compdata = competition_data %>% mutate(points = paste(ID,Hours)) %>% 
  filter(!(points %in% outlier_points))

start_end_times = outlier_removed_compdata %>% rownames_to_column(var='N_growthcurves') %>% 
  group_by(ID) %>% 
  dplyr::summarize(
    start_time = first(N_growthcurves),
    end_time = last(N_growthcurves))

stan_data_y0_unequal = list(
  N = nrow(outlier_removed_compdata),
  N_ID = nrow(design),
  IDX= outlier_removed_compdata$ID,
  N_y0 = length(ID_y0_levels),
  ID_y10 = design$ID_y10,
  ID_y20 = design$ID_y20,
  t0=0,
  yobs = as.matrix(outlier_removed_compdata[,c("ODeq1","ODeq2")]),
  ts = outlier_removed_compdata$Hours,
  start_time = as.numeric(start_end_times$start_time),
  end_time = as.numeric(start_end_times$end_time),
  N_X = ncol(modmat1),
  X1 = modmat1,
  X2 = modmat2,
  N_pred = nrow(newdata_modmat),
  X_pred = newdata_modmat,
  X_pred_no_epistasis = noepimodmat,
  prior_r = c(0.5,0.25),
  prior_K = c(0.3,0.2),
  prior_y0 = c(median(med_y0),1.5),
  epsilon = 1e-4,
  num_steps = 1e6
  
)


onek_refit =  onek_unequal_timepoints_mod$sample(
  data = stan_data_y0_unequal,chains=4,iter_sampling = 1000,iter_warmup = 1000,
  init = init_fun
)


w_outlier_global_params = summarize_draws(llg_fit$draws(c('dc','sigma','sigma_dc','rho[1,2]','nu','K')),default_summary_measures())
no_outlier_global_params = summarize_draws(onek_refit$draws(c('dc','sigma','sigma_dc','rho[1,2]','nu','K')),default_summary_measures())



w_outlier_global_params %>% filter(grepl('dc\\[',variable)) %>% 
  ggplot(aes(variable)) + 
  geom_pointinterval(aes(y=mean,ymin=q5,ymax=q95,color='with outliers')) +
  geom_pointinterval(data = no_outlier_global_params %>%filter(grepl('dc\\[',variable)),
                     aes(y=mean,ymin=q5,ymax=q95,color='no outliers')) +
  labs(x=NULL,y="differential costs posterior mean and 95% HPD",color=NULL) + theme_classic() +
  theme(legend.position = 'none')

w_outlier_global_params %>% filter(!grepl('dc\\[',variable)) %>% 
  mutate(variable = recode(variable, 'rho[1,2]' = 'correlation rho', 'nu' = 'degrees of freedom nu', 'sigma' = 'sigma_mu', 'sigma_dc' = 'sigma_c')) %>% 
  ggplot(aes(variable)) + facet_wrap(~variable,scales="free")+
  geom_pointinterval(aes(y=mean,ymin=q5,ymax=q95,color='with outliers')) +
  geom_pointinterval(data = no_outlier_global_params %>%filter(!grepl('dc',variable)) %>% 
                       mutate(variable = recode(variable, 'rho[1,2]' = 'correlation rho', 'nu' = 'degrees of freedom nu', 'sigma' = 'sigma_mu', 'sigma_dc' = 'sigma_c')),
                     aes(y=mean,ymin=q5,ymax=q95,color='no outliers')) + 
  labs(x=NULL,y="posterior mean and 95% HPD",color=NULL) + theme_classic() + theme(  strip.background = element_blank(),
                                                                                     strip.text.x = element_blank())


w_outlier_global_params %>% rename_with(~paste("outlier",.x,sep="_"),) %>% 
  cbind(no_outlier_global_params) %>% 
  mutate(
    sigdiff = ifelse(
      (outlier_mean>mean & outlier_q5 > q95) | 
        (mean>outlier_mean & q5 > outlier_q95),
      'significant','overlapping'
    )
  ) %>% filter(sigdiff == 'significant')


w_outlier_global_params %>% rename_with(~paste("outlier",.x,sep="_"),) %>% 
  cbind(no_outlier_global_params) %>% 
  mutate(
    sigdiff = ifelse(
      (outlier_mean>mean & outlier_mean-outlier_sd > mean+sd) | 
        (mean>outlier_mean & (mean-sd) > outlier_mean+sd),
      'significant','overlapping'
    )
  ) %>% filter(sigdiff == 'significant') 

w_outlier_local_params = summarize_draws(llg_fit$draws(c('log_y0','r')),default_summary_measures())
no_outlier_local_params = summarize_draws(onek_refit$draws(c('log_y0','r')),default_summary_measures())

w_outlier_local_params %>% rename_with(~paste("outlier",.x,sep="_"),) %>% 
  cbind(no_outlier_local_params) %>% 
  mutate(
    sigdiff = ifelse(
      (outlier_mean>mean & outlier_q5 > q95) | 
        (mean>outlier_mean & q5 > outlier_q95),
      'significant','overlapping'
    )
  ) %>% filter(sigdiff == 'significant')

w_outlier_local_params %>% rename_with(~paste("outlier",.x,sep="_"),) %>% 
  cbind(no_outlier_local_params) %>% 
  mutate(
    sigdiff = ifelse(
      (outlier_mean>mean & outlier_mean-outlier_sd > mean+sd) | 
        (mean>outlier_mean & (mean-sd) > outlier_mean+sd),
      'significant','overlapping'
    )
  ) %>% filter(sigdiff == 'significant') %>%  mutate(across(!c(outlier_variable, variable,sigdiff), exp)) 





# ODE system using deSolve for posterior prediction -----------------------


# ## ODE system -----------------------------------------------------------


#'
#'
#' This function generates the growth curves for two populations A and B
#' That grow according to the Baranyi-Roberts differential equation with shared carrying-capacity K
#' B always grows with a cost that is equivalent to the selection coefficient
#' @param A,B population size of A and B respectively
#' @param mu average growth rate of both populations
#' @param s1,s2 cost as described in paper
#' @param K carrying capacity of the overall culture containing both populations
#' @return the list containing the derivatives of both ODEs to be used by package deSolve for numerical integration
ode_density_logis <- function (time, init, parms, ...) {
  with(as.list(c(parms, init)), {
    dA <- (mu+s1)*A*(1-(A+B)/K)
    dB <- (mu+s2)*B*(1-(A+B)/K)
    list(c(dA,dB))
  })
}
grow_density_logis <- function(time, parms, ...) {
  init <- parms[c("A0","B0")] # initial values
  names(init) <- c("A", "B") # force names of state variables
  odeparms <- parms[c("K","mu","s1","s2")] # the parms of the ODE model
  out <- ode(init, time, ode_density_logis, parms = odeparms)
  out}


# check fit growh curve using mu -------------------------------------------------------
nwells = 4 # number of wells (competition experiments) for which to generate predictions
test_wells = design %>% mutate(day = as.character(day)) %>% subset(select=c(Well,day,ID_y10,ID_y20)) %>% 
  slice_sample(n=nwells) %>% inner_join(competition_data %>% mutate(rowname=1:n())) %>% 
  subset(select = c(ODeq1,ODeq2,Hours,Well,day,ID,gen1,gen2,col1,col2,ID_y10,ID_y20,rowname)) %>% 
  mutate(col1 = factor(paste(col1,'FP',sep=""),levels=c('CFP','OFP','GFP')),
         col2 = factor(paste(col2,'FP',sep=""),levels=c('CFP','OFP','GFP')))

#number of posterior draws for the actual predicted growth curve (mu) and cost-based prediction
n_sim_mu = 20;
n_sim_costs = 20;

#predict growth curves for one K model
growthcurve_preds_onek = as_tibble(
  llg_fit$draws(c(
    unique(paste('mu[',test_wells$rowname,',1]',sep="")),
    unique(paste('mu[',test_wells$rowname,',2]',sep=""))),format='matrix')[sample(4000,n_sim_mu),],
  .name_repair = ~vctrs::vec_as_names(...,repair = 'universal',quiet=T)) %>% 
  rownames_to_column(var='draw') %>% pivot_longer(-draw) %>%
  separate(name,into=c('na','rowname','competitor'),sep="\\.",convert=T) %>% subset(select=-na) %>% 
  mutate(competitor=paste('pred',competitor,sep="")) %>% pivot_wider(names_from=competitor) %>% 
  inner_join(test_wells) %>% 
  mutate(sig = as.numeric(llg_fit_multK$draws('sigma',format='matrix')[sample(4000,1),]),
         pred1 = rlnorm(n(),meanlog=log(pred1),sdlog=sig),
         pred2 = rlnorm(n(),meanlog=log(pred2),sdlog=sig),model = 'one K')

#cost based predicted growth curves for one K model
cost_predict_growthcurves_onek = llg_fit$draws(c('rbar','K',#'log_y0',
                                                 unique(paste('log_y0[',test_wells$ID_y10,']',sep="")),
                                                 unique(paste('log_y0[',test_wells$ID_y20,']',sep="")),
                                                 unique(paste('y_rep[',test_wells$ID,',1]',sep="")),
                                                 unique(paste('y_rep[',test_wells$ID,',2]',sep=""))),format='matrix')[sample(4000,n_sim_costs),] %>% 
  as_tibble(.name_repair = ~vctrs::vec_as_names(...,repair = 'universal',quiet=T)) %>% rownames_to_column(var='draw') %>% 
  pivot_longer(cols = starts_with('y_rep')) %>% 
  pivot_longer(cols = starts_with('log'),names_to='logy0ID',values_to='logy0') %>% 
  mutate(name = gsub('y_rep.','',name)) %>% 
  #  rename_with(.cols = starts_with("`log_y0"),.fn = ~gsub(`log_y0,'',.x))
  separate(name,into=c('ID','competitor',"na"),sep="\\.",convert=T) %>%subset(select=-na) %>%# mutate(ID=as.numeric(parse_number(ID))) %>% 
  mutate(competitor = ifelse(competitor == 1,'s1','s2')) %>% pivot_wider(names_from=competitor) %>% 
  inner_join(test_wells,relationship = "many-to-many") %>% 
  mutate(logy0ID = parse_number(gsub('log_y0.','',logy0ID)),
         y10 = ifelse(ID_y10 == logy0ID,exp(logy0),NA),y20 = ifelse(ID_y20 == logy0ID,exp(logy0),NA)) %>% 
  pivot_longer(c(y10,y20)) %>%  
  filter(logy0ID == ID_y10 | logy0ID == ID_y20) %>% 
  na.omit %>% 
  subset(select=-c(logy0ID,logy0)) %>% unique %>% pivot_wider() %>% 
  group_by(Well,draw,day,col1,col2,ID,rbar,K,y10,y20,s1,s2) %>% nest %>% 
  mutate(
    logis = pmap(list(data,y10,y20,rbar,K,s1,s2),
                 ~as.data.frame(grow_density_logis(
                   time= 0:16,
                   parms = c(A0 = ..2,B0=..3,mu=..4,K=exp(..5),s1=..6,s2=..7))))
  ) %>% unnest(logis) %>% rename('Hours' = 'time') %>% 
  pivot_longer(c(A,B),values_to='pred') %>% mutate(col = ifelse(name=='A',as.character(col1),as.character(col2)),model = 'one K')

#same as above but for multiple K model
growthcurve_preds_multk = as_tibble(
  llg_fit_multK$draws(c(
    unique(paste('mu[',test_wells$rowname,',1]',sep="")),
    unique(paste('mu[',test_wells$rowname,',2]',sep=""))),format='matrix')[sample(4000,n_sim_mu),],
  .name_repair = ~vctrs::vec_as_names(...,repair = 'universal',quiet=T)) %>% 
  rownames_to_column(var='draw') %>% pivot_longer(-draw) %>%
  separate(name,into=c('na','rowname','competitor'),sep="\\.",convert=T) %>% subset(select=-na) %>% 
  mutate(competitor=paste('pred',competitor,sep="")) %>% pivot_wider(names_from=competitor) %>% 
  inner_join(test_wells) %>% 
  mutate(sig = as.numeric(llg_fit_multK$draws('sigma',format='matrix')[sample(4000,1),]),
         pred1 = rlnorm(n(),meanlog=log(pred1),sdlog=sig),
         pred2 = rlnorm(n(),meanlog=log(pred2),sdlog=sig),model = 'multilevel K')

cost_predict_growthcurves_multk = llg_fit_multK$draws(c('rbar',
                                                        unique(paste('K_array[',test_wells$ID,']',sep="")),
                                                        unique(paste('log_y0[',test_wells$ID_y10,']',sep="")),
                                                        unique(paste('log_y0[',test_wells$ID_y20,']',sep="")),
                                                        unique(paste('y_rep[',test_wells$ID,',1]',sep="")),
                                                        unique(paste('y_rep[',test_wells$ID,',2]',sep=""))),format='matrix')[sample(4000,n_sim_costs),] %>% 
  as_tibble(.name_repair = ~vctrs::vec_as_names(...,repair = 'universal',quiet=T)) %>% rownames_to_column(var='draw') %>% 
  pivot_longer(cols = starts_with('y_rep'),names_transform = list(name = ~gsub('y_rep.','',.x))) %>% 
  pivot_longer(cols = starts_with('K_array'),names_to = 'K_ID',values_to = 'log_K',names_transform = list(K_ID = ~parse_number(gsub('K_array.','',.x)))) %>% 
  pivot_longer(cols = starts_with('log_y0'),names_to='logy0ID',values_to='logy0',names_transform = list(logy0ID = ~parse_number(gsub('log_y0.','',.x)))) %>% 
  #  rename_with(.cols = starts_with("`log_y0"),.fn = ~gsub(`log_y0,'',.x))
  separate(name,into=c('ID','competitor',"na"),sep="\\.",convert=T) %>%subset(select=-na) %>%# mutate(ID=as.numeric(parse_number(ID))) %>% 
  filter(K_ID == ID) %>% 
  mutate(competitor = ifelse(competitor == 1,'s1','s2')) %>% pivot_wider(names_from=competitor) %>% 
  inner_join(test_wells,relationship = "many-to-many") %>% 
  mutate(
    y10 = ifelse(ID_y10 == logy0ID,exp(logy0),NA),y20 = ifelse(ID_y20 == logy0ID,exp(logy0),NA)) %>% 
  pivot_longer(c(y10,y20)) %>%  
  filter(logy0ID == ID_y10 | logy0ID == ID_y20) %>% 
  na.omit %>% 
  subset(select=-c(logy0ID,logy0)) %>% unique %>% pivot_wider() %>% 
  group_by(Well,draw,col1,col2,day,ID,rbar,log_K,y10,y20,s1,s2) %>% nest %>% 
  mutate(
    logis = pmap(list(data,y10,y20,rbar,log_K,s1,s2),
                 ~as.data.frame(grow_density_logis(
                   time= 0:16,
                   parms = c(A0 = ..2,B0=..3,mu=..4,K=exp(..5),s1=..6,s2=..7))))
  ) %>% unnest(logis) %>% rename('Hours' = 'time') %>% 
  pivot_longer(c(A,B),values_to='pred') %>% mutate(col = ifelse(name=='A',as.character(col1),as.character(col2)),model = 'multilevel K')



# Figure 4a
test_wells %>%
  ggplot(aes(Hours,value)) + facet_grid(model~ID) + scale_color_manual(values = c('darkturquoise','darkgreen','darkorange'))+
  scale_fill_manual(values = c('darkturquoise','darkgreen','darkorange'))+
  stat_ribbon(data=growthcurve_preds_onek,aes(y=pred1,fill=col1,fill_ramp=after_stat(.width)),alpha=0.5,.width=ppoints(50))+
  stat_ribbon(data=growthcurve_preds_onek,aes(y=pred2,fill=col2,fill_ramp=after_stat(.width)),alpha=0.5,.width=ppoints(50))+
  stat_ribbon(data=growthcurve_preds_multk,aes(y=pred1,fill=col1,fill_ramp=after_stat(.width),alpha=0.5),alpha=0.25,.width=ppoints(50))+
  stat_ribbon(data=growthcurve_preds_multk,aes(y=pred2,fill=col2,fill_ramp=after_stat(.width)),alpha=0.5,.width=ppoints(50))+
  scale_fill_ramp_continuous(range = c(0.75,0.25))+
  geom_line(data=cost_predict_growthcurves_onek,
            aes(y=pred,group=paste(draw,ID,name)),color='grey20',alpha=0.3,linewidth=0.25,linetype='dashed') +
  geom_line(data=cost_predict_growthcurves_multk,
            aes(y=pred,group=paste(draw,ID,name)),color='grey20',alpha=0.3,linewidth=0.25,linetype='dashed') +
  geom_point(aes(y=ODeq1,color=col1),size=0.6) + geom_point(aes(y=ODeq2,color=col2),size=0.6) +
  coord_trans(y = "log10")+ scale_y_continuous(breaks = c(0.01,0.1,1)) +
  theme_classic() +  theme(text = element_text(family='Helvetica',size=rel(4)),legend.position='bottom',strip.text.x = element_blank(),
                           legend.text = element_text(size=rel(3))) + 
  labs(color=NULL) + 
  guides(color=guide_legend(override.aes=list(size=5))) + labs(y='ODeq')



# rhat plot ---------------------------------------------------------------

rhats = llg_fit_multK$summary(c('r','log_y0','K','s','dc','sigma','sigma_dc','nu','rho[1,2]'))$rhat
bayesplot::mcmc_rhat(rhats) + scale_x_continuous(limits = c(0.999,max(rhats))) + theme(legend.position = 'none')

# rhat /  neff plots --------------------------------------------------------------

convergence_stats_onek = summarise_draws(
  llg_fit$draws(c('r','log_y0','K','s','dc','sigma','sigma_dc','nu','rho[1,2]')),
  default_convergence_measures())

convergence_stats = summarise_draws(
  llg_fit_multK$draws(c('r','log_y0','K','K_array','s','dc','sigma','sigma_dc','sigma_K','nu','rho[1,2]')),
  default_convergence_measures())

bayesplot::mcmc_rhat(convergence_stats_onek$rhat) + scale_x_continuous(limits = c(0.999,max(convergence_stats_onek$rhat))) + theme(legend.position = 'none')
bayesplot::mcmc_rhat(convergence_stats$rhat) + scale_x_continuous(limits = c(0.999,max(convergence_stats$rhat))) + theme(legend.position = 'none')

mcmc_neff(convergence_stats_onek$ess_tail) + theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0,'cm'))
mcmc_neff(convergence_stats$ess_tail) + theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0,'cm'))

mcmc_neff(convergence_stats_onek$ess_bulk) + theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0,'cm'))
mcmc_neff(convergence_stats$ess_bulk) + theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0,'cm'))

# variance and covariance parameter  densities  -------------------------------------
#supp figure
llg_fit_multK$draws(c('sigma','sigma_dc','nu','sigma_K'),format='matrix') %>% as_tibble %>% 
  rename('σ μ' = 'sigma', 'σ c' = 'sigma_dc', 'nu' = 'nu','σ K' = 'sigma_K') %>% gather %>%#filter(key == 'nu') %>% 
  ggplot(aes(x=value)) + facet_grid(~key,scales="free_x") + 
  stat_halfeye(aes(thickness=stat(pdf)),normalize='panels')+
  theme_classic() +
  theme(text=element_text(family = 'Helvetica',size=rel(4)),axis.text = element_text(size=rel(1)),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),strip.text.x = element_text(size=rel(6))
  )+ labs(x=NULL,y='posterior density')

#Figure 4c inset
tibble(rvar = as_draws_rvars(llg_fit_multK$draws('rho[1,2]',format='matrix'))$rho[1,2]) %>% 
  ggplot(aes(xdist=rvar)) + stat_pointinterval()+
  #stat_halfeye() + #geom_vline(xintercept=0,show.legend=F) +
  theme_classic() + labs(x='Correlation coefficient',y=NULL) + scale_y_continuous(breaks = NULL)



# Predicted and residuals for s and ppcheck -------------------------------------------------------------

preds = as_tibble(llg_fit_multK$draws('s',format='draws_list')) %>% 
  mutate(ID=rep(1:nrow(design),2),which=rep(c('s1','s2'),each=nrow(design))) %>% 
  group_by(ID,which) %>% nest %>% 
  mutate(alldat = map(data,~unlist(.x)),mean=map(alldat,mean),upr=map(alldat,~mean(.x)+sd(.x)),lwr=map(alldat,~mean(.x)-sd(.x))) %>% 
  unnest(c(mean,lwr,upr)) %>% #mutate(rvar = rvar/0.65) %>% 
  subset(select = -c(alldat,data)) #%>% pivot_wider(names_from=which,values_from=c(mean,lwr,upr)) 

resids = as_tibble(llg_fit_multK$draws('resids',format='draws_list')) %>% 
  mutate(ID=rep(1:nrow(design),2),which=rep(c('s1_resid','s2_resid'),each=nrow(design))) %>% 
  group_by(ID,which) %>% nest %>% 
  mutate(alldat = map(data,~unlist(.x)),mean=map(alldat,mean),upr=map(alldat,~mean(.x)+sd(.x)),lwr=map(alldat,~mean(.x)-sd(.x))) %>% unnest(c(mean,lwr,upr)) %>% #mutate(rvar = rvar/0.65) %>% 
  subset(select = -c(alldat,data)) %>% pivot_wider(names_from=which,values_from=c(mean,lwr,upr)) 

#Figure 4c
resids %>%
  ggplot() + 
  geom_pointinterval(aes(x=mean_s1_resid,y=mean_s2_resid,ymin=lwr_s2_resid,ymax=upr_s2_resid),color='grey30',linewidth=0.01,alpha=0.75) + 
  geom_pointinterval(aes(x=mean_s1_resid,y=mean_s2_resid,xmin=lwr_s1_resid,xmax=upr_s1_resid),color='grey30',linewidth=0.01,alpha=0.75) + 
  theme_classic() + 
  geom_vline(xintercept=0) + geom_hline(yintercept=0) + labs(x='Residuals 1 posterior mean', y='Residuals 2 posterior mean')


#Figure 4D
preds%>% full_join(design) %>% 
  mutate(rep = ifelse(day %in% c('Jan11','Jan26'),'rep1','rep2')) %>% 
  subset(select = c(which,mean,upr,lwr,dilution,rep,gen1,gen2,col1,col2)) %>% 
  pivot_wider(names_from=rep,values_from=c(mean,upr,lwr)) %>% mutate(dilution = paste(dilution,'X',sep="")) %>% na.omit %>% 
  mutate(outlier = ifelse(abs(mean_rep1 - mean_rep2) > 0.15,1,0)) %>% 
  ggplot(aes(color=dilution,size=outlier)) + 
  geom_pointinterval(aes(x=mean_rep1,y=mean_rep2,xmin=lwr_rep1,xmax=upr_rep1),alpha=0.75,linewidth=0.01) +
  geom_pointinterval(aes(x=mean_rep1,y=mean_rep2,ymin=lwr_rep2,ymax=upr_rep2),alpha=0.75,linewidth=0.01) +
  geom_abline(slope=1,intercept=0) + theme_classic() + scale_color_manual(values = c('#B77F08','#1975C5'),guide = guide_legend(reverse = TRUE))+
  labs(x='Replicate 1 posterior mean', y='Replicate 2 posterior mean',color = 'inoculum dilution')+
  theme(text=element_text(family='Helvetica',size=rel(4)),legend.text = element_text(size=rel(3)),legend.position = 'bottom') 







# Outliers from replicability plot (requires "preds" from above section) ---------------------
outliers_from_replicability_plot = preds%>% full_join(design) %>% 
  mutate(rep = ifelse(day %in% c('Jan11','Jan26'),'rep1','rep2')) %>% 
  subset(select = c(which,mean,upr,lwr,dilution,rep,gen1,gen2,col1,col2)) %>% 
  pivot_wider(names_from=rep,values_from=c(mean,upr,lwr)) %>% 
  filter(abs(mean_rep1 - mean_rep2) > 0.15) %>% mutate(which_outs = paste(gen1,gen2,col1,col2,dilution)) %>% pull(which_outs)



# these basically show that the competitions look like competitions (no weird aberrant fluorescence) which indicates
# these outcomes are really that noisy and the outliers should not be removed (can't choose which one to)
competition_data %>%   
  mutate(rep = ifelse(day %in% c('Jan11','Jan26'),'rep1','rep2'),
         ID_string = paste(gen1,gen2,col1,col2,dilution)) %>% filter(ID_string %in% outliers_from_replicability_plot) %>% 
  ggplot(aes(Hours)) + facet_grid(ID_string~rep) + scale_color_manual(values = c('darkturquoise','darkgreen','darkorange')) +
  geom_point(aes(y=ODeq1,color=col1)) + geom_point(aes(y=ODeq2,color=col2)) + theme_classic()

# Comparison to colony counts ---------------------------------------------

colony_count_data = read.csv('cleaned_datasets/colony_competition_raw_data.csv',row.names = 1)

colony_count_data %>% group_by(genA,genB,pair,Time,Date) %>% # grouping
  mutate(
    # convert cells to optical density using log-log relationship 
    # experimentally obtained
    A=exp(-14.5236887834642)*A^0.639273426672917,B=exp(-14.5236887834642)*B^0.639273426672917,
    A_over_B = log(A/B)) %>% subset(select=c(Date,Time,genA,genB,pair,A_over_B)) %>% 
  pivot_wider(names_from=Time,values_from=A_over_B,values_fn = list) %>% 
  mutate(s = map2_dbl(`0`,`24`, ~mean((.y-.x)/24))) %>% 
  subset(select = c(pair,s)) %>% group_by(pair) %>% nest %>% 
  mutate(colony_s = map_dbl(data,~mean(.x$s)),colony_se = map_dbl(data,~sd(.x$s)/length(.x$s))) %>% 
  ungroup %>% mutate(
    ## adding predicted effects
    preds = c(
      rvar(c(llg_fit_multK$draws(c('dc[3]')))),
      rvar(c(llg_fit_multK$draws(c('dc[2]')))),
      rvar(c(llg_fit_multK$draws(c('dc[2]'))))-rvar(c(llg_fit_multK$draws(c('dc[3]'))))),
    factor = preds/colony_s, meanpred = mean(preds), factorlabel = paste('X',round(mean(factor),2),'±',round(sd(factor),2),sep=" ")
  ) %>% 
  ggplot(aes(xdist=preds)) + stat_halfeye(normalize='xy') + facet_grid(~pair,scales="free_x") +
  geom_pointinterval(aes(x=colony_s,xmin=colony_s-colony_se,xmax=colony_s+colony_se,y=0.0),color='red4',size=5,linewidth=2) +
  geom_curve(aes(x=meanpred,xend=colony_s,y=0,yend=0.0),linetype='dashed',color = 'tomato3') + 
  geom_text(aes(x=colony_s*2,y=-0.2,label=factorlabel),color='red4',size=6,family='Helvetica')+
  theme_classic() + theme(text=element_text(family='Helvetica',size=rel(4)),strip.text.x = element_text(size=rel(5))) + 
  labs(x='selection coefficient',y='posterior density') + 
  geom_text(data=. %>% filter(pair == 'GFP vs OFP'), aes(x=-0.1,y=-0.4),label = 'scaling factor',color='red4',size=6,family='Helvetica') + 
  scale_x_continuous(labels=scales::label_number(accuracy = 0.01))


# how to calculate scaling factor
0.54/
  mean(diff(log(1.13/(1+((1-0.01)/0.01)*exp(-0.54*c(0:24))))))

# Epistatic terms  ---------------------------------------------------------

# mean non epistatic effect
#vectorized function to return absolute value vector/matrix
abs_sets <- function(vec){
  negs <- vec < 0
  vec[negs] <- vec[negs] * -1
  vec
}


meaneff = mean(abs_sets(llg_fit_multK$draws(paste('dc[',1:12,']',sep=""),format='matrix')))

as_tibble(llg_fit_multK$draws(paste('dc[',13:30,']',sep=""),format='draws_list')) %>% 
  mutate(ID = colnames(modmat1)[13:30]) %>% group_by(ID) %>% nest %>% ungroup %>% 
  mutate(alldat=map(data,unlist),
         rvar=map(alldat,rvar)) %>% subset(select=c(ID,rvar)) %>% unnest(rvar) %>% 
  separate(ID,into=c('col','gen'),sep=":") %>% 
  mutate(
    col=paste(gsub('col1','',col),'FP',sep=""),
    gen=gsub('gen1','',gen),
    gen=factor(gen,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                              "VIMR6.4","VIMR18.1","NDMR18")),
    col=factor(col, levels = c("CFP","GFP","OFP"))) %>% rowwise %>% 
  mutate(
    prop = ifelse(median(rvar) > 0,sum(rvar<0)/4000,sum(rvar>0)/4000),
    loc = quantile(rvar,probs=1)+0.02,
    sig = ifelse(prop<0.01,'wow','zow')
  ) %>% 
  ggplot(aes(xdist=rvar,y=(as.numeric(gen)))) + facet_grid(~col) + stat_halfeye(aes(fill=col),alpha=0.5,show.legend=F) +
  geom_vline(xintercept=0,show.legend=F) + scale_fill_manual(values = c('green4','orange3')) + theme_classic() +
  labs(x='epistatic fitness coefficient',y='Plasmid') + 
  scale_x_continuous(limits = c(-.13,0.12))+
  scale_color_manual(values = c('red4','grey20')) + scale_y_continuous(breaks = 2:10,limits = c(1,11))+
  geom_vline(xintercept=c(-meaneff,meaneff),linetype='11',color='grey40') +
  geom_text(aes(x=0.1,y=as.numeric(gen)+0.5,color=sig,
                label=paste(round(prop*100,2),'%',sep="")),size=6,show.legend=F) +
  theme(text = element_text(size=rel(5),family='Helvetica'),strip.text.x = element_text(size=rel(4.5)))





# Epistasis analysis for significant epistatic terms -------------------------------------------------

y_pred_full = as_tibble(llg_fit_multK$draws(paste('y_pred',sep=""),format='draws_list')) %>% 
  tibble(newdata %>% subset(select=-estimate_corr)) %>% group_by(gen1,col1) %>% nest %>% ungroup %>% 
  mutate(alldat=map(data,unlist),
         rvar=map(alldat,rvar)) %>% subset(select=c(gen1,col1,rvar)) %>% unnest(rvar) %>% mutate(pred = 'with epistasis')
y_pred_noepi = as_tibble(llg_fit_multK$draws(paste('y_pred_no_epi',sep=""),format='draws_list')) %>% 
  tibble(newdata %>% subset(select=-estimate_corr)) %>% group_by(gen1,col1) %>% nest %>% ungroup %>% 
  mutate(alldat=map(data,unlist),
         rvar=map(alldat,rvar)) %>% subset(select=c(gen1,col1,rvar)) %>% unnest(rvar) %>% mutate(pred = 'no epistatic coefficient')

significant_ones =  as_tibble(llg_fit_multK$draws(paste('dc[',13:30,']',sep=""),format='draws_list')) %>% 
  mutate(ID = colnames(modmat1)[13:30]) %>% group_by(ID) %>% nest %>% ungroup %>% 
  mutate(alldat=map(data,unlist),
         rvar=map(alldat,rvar)) %>% subset(select=c(ID,rvar)) %>% unnest(rvar) %>% 
  separate(ID,into=c('col','gen'),sep=":") %>% 
  mutate(
    col1=factor(gsub('col1','',col),levels=c('C','G','O')),
    gen1=factor(gsub('gen1','',gen),levels=levels(design$gen1))) %>% rowwise %>% 
  mutate(
    prop = ifelse(median(rvar) > 0,sum(rvar<0)/4000,sum(rvar>0)/4000),
    loc = quantile(rvar,probs=1)+0.01,
    sig = ifelse(prop<0.01,'wow','zow')) %>% subset(select=c(col1,gen1,sig))

mean_gr = mean(c(unlist(llg_fit_multK$draws('rbar',format='draws_list'))))
mean_cfp = mean(c(unlist(llg_fit_multK$draws('dc[1]',format='draws_list'))))
mean_gfp = mean_cfp + mean(llg_fit_multK$draws('dc[2]'))
mean_ofp = mean_cfp + mean(c(unlist(llg_fit_multK$draws('dc[3]',format='draws_list'))))


library(ggnewscale)
rbind(y_pred_full,y_pred_noepi) %>% #pivot_wider(names_from=pred,values_from=rvar) %>% 
  filter(col1 != 'C') %>%   full_join(significant_ones) %>% filter(sig == 'wow') %>% 
  mutate(col1 = paste(col1,'FP',sep=""),
         gen1=factor(gen1,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                     "VIMR6.4","VIMR18.1","NDMR18")),
         col1=factor(col1, levels = c("CFP","GFP","OFP"))) %>% 
  mutate(rvar=mean_gr + rvar) %>% 
  ggplot(aes(y=as.factor(as.numeric((gen1))))) + #facet_grid(~col1,scales="free_x") + 
  #stat_pointinterval(aes(group=pred),show.legend=F,color='grey50',size=2)+
  stat_slab(aes(xdist=rvar,color=pred),fill='transparent') + scale_color_manual(values = c('khaki3','tomato3'))+
  theme_classic()  +
  geom_vline(xintercept = mean_gr,linetype='12',color='grey20')+
  geom_vline(xintercept = c(
    quantile(c(unlist(llg_fit_multK$draws('rbar',format='draws_list'))),0.025),
    quantile(c(unlist(llg_fit_multK$draws('rbar',format='draws_list'))),0.975)
  ),linetype='14',color='grey20')+ theme(legend.position = 'bottom') + guides(color=guide_legend(nrow=2))+
  geom_vline(data=tibble(col1='GFP',intercept=(mean_gr+mean_gfp)),aes(xintercept=intercept),color='green4',linetype='12',linewidth=1)+
  geom_vline(data=tibble(col1='OFP',intercept=(mean_gr+mean_ofp)),aes(xintercept=intercept),color='darkorange',linetype='12',linewidth=1)+
  labs(x='maximal growth rate',y=NULL,color = NULL) +  new_scale_color() +
  geom_curve(data = . %>%   mutate(med = median(rvar)) %>% subset(select=-rvar) %>% pivot_wider(names_from=pred,values_from = med) %>% 
               #mutate(`epistasis type` = ifelse(abs(`with epistasis`-full_fit$summary('r')$mean) < abs(`no epistatic coefficient`-full_fit$summary('r')$mean),'type 1','type 2')),
               mutate(`epistasis type` = ifelse(`with epistasis` < `no epistatic coefficient`,'negative','positive')),
             aes(x=`no epistatic coefficient`, xend = `with epistasis`,yend = as.factor(as.numeric((gen1))),color=`epistasis type`),linewidth=1,
             arrow = arrow(type='closed',length = unit(0.03, "npc"))) + scale_color_manual(values = c('#39A293','#E1B265')) + labs(color=NULL)+
  guides(color=guide_legend(nrow=2))


rbind(y_pred_full_robust,y_pred_noepi_robust) %>% #pivot_wider(names_from=pred,values_from=rvar) %>% 
  filter(col1 != 'C') %>%   full_join(significant_ones) %>% filter(sig == 'wow') %>% 
  mutate(col1 = paste(col1,'FP',sep=""),
         gen1=factor(gen1,levels = c("Junk","CMP","VIMwt","VIMR1.1","VIMR2.1","VIMR2.2","VIMR6.1",
                                     "VIMR6.4","VIMR18.1","NDMR18")),
         col1=factor(col1, levels = c("CFP","GFP","OFP"))) %>% 
  mutate(rvar=(1+rvar)) %>% 
  ggplot(aes(y=as.factor(as.numeric((gen1))))) + facet_grid(~col1) + 
  #stat_pointinterval(aes(group=pred),show.legend=F,color='grey50',size=2)+
  stat_slab(aes(xdist=rvar,color=pred),fill='transparent') + scale_color_manual(values = c('khaki3','tomato3'))+
  theme_classic()  +
  geom_vline(xintercept = 1,linetype='12',color='grey20')+
  geom_vline(xintercept = c(
    1+ full_fit$summary('r')$mean - full_fit$summary('r')$q5,
    1+full_fit$summary('r')$q95 - full_fit$summary('r')$mean
  ),linetype='14',color='grey20')+ theme(legend.position = 'bottom') + guides(color=guide_legend(nrow=2))+
  labs(x='relative fitness',y=NULL,color = NULL) +  new_scale_color() +
  geom_curve(data = . %>%   mutate(med = median(rvar)) %>% subset(select=-rvar) %>% pivot_wider(names_from=pred,values_from = med) %>% 
               mutate(`epistasis type` = ifelse(abs(`with epistasis`-full_fit$summary('r')$mean) < abs(`no epistatic coefficient`-full_fit$summary('r')$mean),'type 1','type 2')),
             aes(x=`no epistatic coefficient`, xend = `with epistasis`,yend = as.factor(as.numeric((gen1))),color=`epistasis type`),linewidth=1,
             arrow = arrow(type='closed',length = unit(0.03, "npc"))) + scale_color_manual(values = c('#39A293','#E1B265')) + labs(color=NULL)+
  guides(color=guide_legend(nrow=2))






# Predict genotype growth rates  -------------------------------------------------------
#(useful for other applications with those genotypes e.g. cost correction)



ypred = as_tibble(llg_fit_multK$draws('y_pred',format='draws_list')) %>% tibble(newdata %>% subset(select=-estimate_corr)) %>% 
  group_by(gen1,col1) %>% nest %>% 
  mutate(alldat = map(data,~unlist(.x)),
         mean = map_dbl(alldat,mean),
         sd = map_dbl(alldat,sd),
         median = map_dbl(alldat,median),
         mad = map_dbl(alldat,mad),
         lwr = map_dbl(alldat,quantile,probs=0.025),
         upr = map_dbl(alldat,quantile,probs=0.975)) %>% subset(select=-data)




ypred %>% mutate(rvar = map(alldat,rvar)) %>% unnest(rvar) %>% 
  ggplot(aes(y=gen1,fill=col1,xdist=rvar)) +
  stat_halfeye(show.legend=F,point_size=.8,normalize='xy') + facet_grid(~col1) +
  geom_vline(xintercept=0) + scale_fill_manual(values = c('cyan3','green4','orange3')) + theme_classic() +
  labs(x='selection coefficient',y=NULL)





ypred %>% 
  mutate(meangr = map_dbl(alldat,~mean(.x+mean_gr)),
         sdgr =  map_dbl(alldat,~sd(.x+mean_gr)),
         lwrgr = map_dbl(alldat,~quantile(.x+mean_gr,probs=0.025)),
         uprgr = map_dbl(alldat,~quantile(.x+mean_gr,probs=0.975))
  ) %>% subset(select=-alldat) %>% 
  write.csv('cleaned_datasets/2307_sel_coefs_and_growth_rates.csv')


