setwd(dir = "/Users/dankehila/OneDrive - UBC/Dan/Thesis/Fluor fitness paper/2307 code/cleaned_datasets/")
library(tidyverse)


decay_dat <- read.csv("maturations_with_CMP_data.csv",row.names = 1) %>% 
  mutate(Hours=(Time)*24) %>% subset(select=-Time) %>% 
  pivot_longer(-c(Hours,Read),names_to = "Well") %>% 
  separate(Well,into = c("Row", "Column"), 
           sep = "(?<=[A-Za-z])(?=[0-9])",remove=F) %>% 
  filter(case_when(Row == "E" ~ Read == "CFP",
                   Row == "F" ~ Read == "GFP",
                   T ~ Read == "OFP")) %>% 
  group_by(Well) %>% arrange(Well,Hours) %>% 
  filter(as.integer(Column)<7) %>% 
  #GFP matures instantaneously so omitting it
  filter(Read != "GFP") %>% 
  mutate(fold_dilution = 2^as.integer(Column)) 


# using cubic splines as per Balleza et al. (see paper)
decay_dat %>% 
  group_by(Read,Column) %>% mutate(max=quantile(value,probs=0.999)) %>% 
  arrange(Read,Well,Hours) %>% 
  mutate(immature = 1-((value-first(value))/(max-first(value))),
         log_immature = log(immature)) %>% filter(log_immature >-Inf) %>% 
  nest %>% 
  mutate(css = map(data,~smooth.spline(.x$Hours,y=.x$log_immature,spar = 0.99)),
         csspred = map(css,function(x) exp(predict(x)$y))) %>% 
  unnest(c(data,csspred)) %>% 
  ggplot(aes(Hours,immature)) + facet_grid(~Read) + geom_point(aes(color=Read),alpha=0.1) +
  scale_y_continuous(trans=scales::log10_trans()) +
  geom_hline(yintercept = 0.5,linetype='dotted',col='grey60')+
  theme_classic() + labs(y="Fraction immature protein",color="fold dilution") + theme(
    text=element_text(family = 'Helvetica',size=rel(4)),strip.text=element_blank(),
    legend.position="none",
    legend.text=element_text(size=rel(4))) + scale_y_continuous(labels=scales::label_percent(),breaks = c(1,0.5,0)) +
  scale_color_manual(values=c("#3CA9D5","#CC5127")) +
  geom_line(aes(y=csspred,group=fold_dilution,color=Read),linewidth=1) 





# balleza method splines
decay_dat %>% 
  group_by(Read,fold_dilution,Column) %>% mutate(max=quantile(value,probs=0.999)) %>% 
  filter(value<max) %>% 
  arrange(Read,fold_dilution,Well,Hours) %>% 
  mutate(immature = log(1-((value-first(value))/(max-first(value))))) %>%
  filter(immature >-Inf) %>% 
  nest() %>% 
  mutate(css = map(data,~smooth.spline(.x$Hours,y=.x$immature,spar = 0.99)),
         hrs = list(seq(0,6,length.out=1e5)),
         csspred = map(css,function(x) exp(predict(object=x,x=seq(0,6,length.out=1e5))$y))) %>% 
  subset(select= c(fold_dilution,Read,Column,csspred,hrs)) %>%
  unnest(c(hrs,csspred)) %>% 
  # filter(fold_dilution < 16) %>% 
  # dplyr::filter(ifelse(Read =="CFP",Hours <  3, Hours> 1 & Hours < 3.329)) %>%filter(immature>0) %>%
  # dplyr::filter(Hours < 5 & ifelse(Read=="OFP",fold_dilution > 2,fold_dilution<16) )%>%
  group_by(Column,Read) %>%
  filter(abs(0.5-csspred)<0.0001 & abs(0.5-csspred) > 0) %>% 
  group_by(Read) %>% dplyr::summarize(n=n(),meanmat50=mean(hrs,na.rm=T),sd=sd(hrs,na.rm=T)/sqrt(n())) %>% 
  as.data.frame %>% mutate(across(where(is.double),round,digits=2))






