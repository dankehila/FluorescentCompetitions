setwd(dir = "make_datasets/")
library(tidyverse);library(broom)
`%notin%` <- Negate(`%in%`)
# Calibration files (background fluorescence/absorbance) for dynamic calibration dataset -------------------------------------------------------


bg_09_10 <- read.csv("220209_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% 
  filter(Colour == "blank") %>% group_by(Read) %>% dplyr::summarize(mean=min(value)) %>% as.data.frame()

bg_16 <- read.csv("220216_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% 
  filter(Colour == "blank") %>% group_by(Read) %>% dplyr::summarize(mean=min(value)) %>% as.data.frame()


bg_23 <- read.csv("220223_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(grepl("([1][7-9]|[2][0-4])",Well)) %>% subset(select=-Well) %>%  summarise_all(min) %>% t() %>% na.omit() %>% as.data.frame() %>% suppressWarnings()

bg_24 <- read.csv("220224_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(grepl("([1][7-9]|[2][0-4])",Well)) %>% subset(select=-Well) %>%summarise_all(min) %>% t() %>% na.omit() %>% as.data.frame() %>% suppressWarnings()

bg_0301 <- read.csv("220301_compete.csv",header=T) %>%   pivot_longer(-c(Read,Time),names_to = "Well") %>% filter(grepl("A8|A9|A10|A11|A12|A13",Well))  %>% 
  filter(Time ==min(Time)) %>% 
  group_by(Read) %>% dplyr::summarize(mean=min(value)) %>% as.data.frame()


calib_09_10 <- read.csv("220209_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% filter(Colour != "blank") %>% 
  pivot_wider(names_from=Read,values_from=value,values_fn=list) %>% unnest(cols=c(OD600,CFP,OFP)) %>% 
  mutate(CFP=CFP-bg_09_10[1,2],OFP=OFP-bg_09_10[3,2],OD600=OD600-bg_09_10[2,2]) %>%  filter(OD600<0.8) %>% 
  mutate(RFU=case_when(Colour=="T" ~ CFP, Colour == "O" ~ OFP)) %>% group_by(Colour) %>% do(fit=tidy(lm(OD600~RFU,data=.))) %>% unnest(fit) %>% 
  dplyr::select(c(Colour,term,estimate)) %>% pivot_wider(names_from=term,values_from=estimate) %>% rename("Intercept" = `(Intercept)`) %>% as.data.frame()

calib_16 <- read.csv("220216_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Colour,Genotype),names_to="column") %>% filter(Colour != "blank") %>% 
  pivot_wider(names_from=Read,values_from=value,values_fn=list) %>% unnest(cols=c(OD600,CFP,OFP,GFP)) %>% 
  mutate(CFP=CFP-bg_16[1,2],OFP=OFP-bg_16[4,2],OD600=OD600-bg_16[3,2],GFP=GFP-bg_16[2,2]) %>% filter(OD600<0.8) %>% 
  mutate(RFU=case_when(Colour=="T" ~ CFP, Colour == "O" ~ OFP, Colour == "G" ~ GFP)) %>% group_by(Colour) %>% do(fit=tidy(lm(OD600~RFU,data=.))) %>% unnest(fit) %>% 
  dplyr::select(c(Colour,term,estimate)) %>% pivot_wider(names_from=term,values_from=estimate) %>% rename("Intercept" = `(Intercept)`) %>% as.data.frame()



calib_23 <- read.csv("220223_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(!grepl("([1][7-9]|[2][0-4])",Well)) %>% 
  mutate(Colour=rep(c("T","G"),each=16),RFU=case_when(Colour == "T" ~ CFP-bg_23[2,1], T ~ GFP-bg_23[3,1]),OD600=OD600-bg_23[1,1]) %>% 
  group_by(Colour) %>% do(fit=tidy(lm(OD600~RFU,data=.))) %>% unnest(fit) %>% 
  dplyr::select(c(Colour,term,estimate)) %>% pivot_wider(names_from=term,values_from=estimate) %>% rename("Intercept" = `(Intercept)`) %>% as.data.frame()

calib_24 <- read.csv("220224_calib.csv",header=T) %>% filter(grepl("A|B",Well)) %>% 
  filter(!grepl("([1][7-9]|[2][0-4])",Well)) %>% 
  mutate(Colour=rep(c("T","G"),each=16),RFU=case_when(Colour == "T" ~ CFP-bg_24[2,1], T ~ GFP-bg_24[3,1]),OD600=OD600-bg_24[1,1]) %>% 
  group_by(Colour) %>% do(fit=tidy(lm(OD600~RFU,data=.))) %>% unnest(fit) %>% 
  dplyr::select(c(Colour,term,estimate)) %>% pivot_wider(names_from=term,values_from=estimate) %>% rename("Intercept" = `(Intercept)`) %>% as.data.frame()

# Load data for dynamic calibration ---------------------------------------------------------------




dat_09 <- read.csv("220209_glucose_curves.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>% 
  filter(!grepl("A|B|C",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>% 
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220210_des.csv",check.names = F)) %>%
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2)) %>% mutate(day="Feb9")


dat_10 <- read.csv("220210_glucose_curves.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>% 
  filter(!grepl("A|B|C",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>% 
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220210_des.csv",check.names = F)) %>%
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2))  %>% mutate(day="Feb10")


dat_16 <- read.csv("220216_glucose_curves.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>% 
  filter(grepl("L|M|N|O",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>% 
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220216_des.csv",check.names = F)) %>%
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2))  %>% mutate(day="Feb16")



dat_23 <- read.csv("220223_compete.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>%  filter(grepl("C|D|E|F|G|H|I|J",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>% 
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220223_des.csv",check.names = F)) %>% 
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2))  %>% mutate(day="Feb23")

dat_24 <- read.csv("220224_compete.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>%  filter(grepl("C|D|E|F|G|H|I|J",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>% 
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220224_des.csv",check.names = F)) %>% 
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2))  %>%
  mutate(glucose = glucose*2) %>% # doubled the serial dilution this day
  mutate(day="Feb24")

dat_0301 <- read.csv("220301_compete.csv") %>% 
  pivot_longer(-c(Read,Time),names_to = "Well") %>%  filter(!grepl("A|P",Well)) %>% 
  mutate(col=gsub("[^0-9]+", "",Well)) %>% filter(col %notin% c(1,24)) %>%
  mutate(Row=gsub("[0-9]+", "",Well)) %>% 
  full_join(.,read.csv("220301_des.csv",check.names = F)) %>% 
  mutate(glucose = case_when(
    col %in% c(2,3) ~ 0,col %in% c(4,5) ~ 2^-7,col %in% c(6,7) ~ 2^-6,
    col %in% c(8,9) ~ 2^-5,col %in% c(10,11) ~ 2^-4,col %in% c(12,13) ~ 2^-3,
    col %in% c(14,15) ~ 2^-2,col %in% c(16,17) ~ 2^-1,col %in% c(18,19) ~ 1, T ~ 2))  %>%
  mutate(glucose = glucose/2) %>% # halved the serial dilution this day
  mutate(day="Mar1")


monoculture_dat <- rbind(dat_09,dat_10,dat_16,dat_23,dat_24,dat_0301) %>% mutate(Hours = Time*24) %>% 
  filter(Type == "monoculture") %>% na.omit() %>% 
  pivot_wider(names_from=Read,values_from=value) %>% 
  mutate(RFU=case_when(Colour == "T" ~ CFP, Colour == "G" ~ GFP, T ~ OFP)) %>% 
  mutate(RFU = case_when(Colour == "G" & day == "Feb16" ~ RFU-bg_16[2,2],
                         Colour == "T" & day == "Feb16" ~ RFU-bg_16[1,2],
                         Colour == "T" & day %in% c("Feb9","Feb10") ~ RFU-bg_09_10[1,2],
                         Colour == "O" & day == "Feb16" ~ RFU-bg_16[4,2],
                         Colour == "O" & day %in% c("Feb9","Feb10") ~ RFU-bg_09_10[3,2],
                         Colour == "G" & day == "Feb23" ~ RFU-bg_23[3,1],
                         Colour == "G" & day == "Feb24" ~ RFU-bg_24[3,1],
                         Colour == "T" & day == "Feb23" ~ RFU-bg_23[2,1],
                         Colour == "T" & day == "Feb24" ~ RFU-bg_24[2,1],
                         Colour == "T" & day == "Mar1" ~ RFU-bg_0301[1,2],
                         Colour == "G" & day == "Mar1" ~ RFU-bg_0301[2,2],
                         Colour == "O" & day == "Mar1" ~ RFU-bg_0301[4,2])) %>% 
  mutate(OD600 = case_when(day %in% c("Feb9,Feb10") ~ (OD600-bg_09_10[2,2])/0.82,
                           day == "Feb23" ~ (OD600-bg_23[1,1])/0.82, day == "Feb24" ~ (OD600-bg_24[1,1])/0.82, day == "Mar1" ~ (OD600-0.045)/0.82,
                           T ~  (OD600-bg_16[3,2])/0.82)) %>% dplyr::select(-c(Time,CFP,GFP,OFP))


monoculture_dat %>% write.csv('../cleaned_datasets/monocultures_cleaned_all_glucose.csv')



# Make static calibration data --------------------------------------------



Jan27_dat <- read.csv("220127_calib.csv",check.names = F) %>% pivot_longer(-c(Read,Row),names_to="color") %>% 
  mutate(gen=recode(Row,"A"="CMP",
                    "B"="VIMwt",
                    "C"="Junk",
                    "D"="VIMR18.1",
                    "E"="VIMR1.1",
                    "F"= "VIMR6.4",
                    "G"= "VIMR6.1",
                    "H"= "VIMR2.2",
                    "I"="VIMR2.1",
                    "J"="NDMR18","K"="blank"))


Jan27_bg <- Jan27_dat %>% filter(gen=="blank") %>% group_by(Read) %>% dplyr::summarise(bg=mean(value)) %>% as.data.frame()

Jan27_dat <- Jan27_dat %>%  filter(gen!="blank") %>% dplyr::select(-Row) %>% pivot_wider(names_from=Read,values_fn=list) %>% unnest(cols = c(OD600, CFP, GFP, OFP) ) %>% 
  mutate(OD600=(OD600-Jan27_bg[3,2])/0.82,CFP=CFP-Jan27_bg[1,2],GFP=GFP-Jan27_bg[2,2],OFP=OFP-Jan27_bg[4,2]) %>% 
  mutate(RFU=case_when(color=="T" ~ CFP,color=="G"~GFP,T~OFP)) %>% dplyr::select(-c(CFP,GFP,OFP)) %>% mutate(day="Jan27")



feb1_dat <- read.csv("220201_overnight_calibration.csv") %>% 
  mutate(color=rep(rep(c("T","G","O"),each=40),4)) %>% 
  mutate(gen=rep(c(rep(c(1,4,6,8,10),8), # turquoise genotypes
                   rep(c(1,6,9,10,-1),8), # green genotypes
                   rep(c(1,2,8,9,-1),8)),4)) %>%   # orange genotypes
  mutate(gen=recode(gen,"1"="CMP",
                    "2"="VIMwt",
                    "3"="Junk",
                    "4"="VIMR18.1",
                    "5"="VIMR1.1",
                    "6"= "VIMR6.4",
                    "7"= "VIMR6.1",
                    "8"= "VIMR2.2",
                    "9"="VIMR2.1",
                    "10"="NDMR18","-1"="blank")) 


feb1_bg <- feb1_dat %>% filter(gen == "blank") %>% group_by(Read) %>% dplyr::summarise(bg=mean(value)) %>% as.data.frame()

feb1_dat <- feb1_dat %>% pivot_wider(names_from=Read) %>% filter(gen!="blank") %>% 
  mutate(OD600=(OD600-feb1_bg[3,2])/0.82,CFP=CFP-feb1_bg[1,2],GFP=GFP-feb1_bg[2,2],OFP=OFP-feb1_bg[4,2]) %>% 
  mutate(RFU=case_when(color=="T" ~ CFP,color=="G"~GFP,T~OFP)) %>% dplyr::select(-c(CFP,GFP,OFP,Well)) %>% mutate(day="Feb1")


## February 2 with old and new OFP wavelength
feb2_dat <- read.csv("220202_overnights.csv") %>% filter(!grepl("L|M|N|O|P",Well)) %>% 
  mutate(gen=rep(c(1,2,3,4,5,6,7,8,9,10,-1),24)) %>% 
  mutate(gen=recode(gen,"1"="CMP",
                    "2"="VIMwt",
                    "3"="Junk",
                    "4"="VIMR18.1",
                    "5"="VIMR1.1",
                    "6"= "VIMR6.4",
                    "7"= "VIMR6.1",
                    "8"= "VIMR2.2",
                    "9"="VIMR2.1",
                    "10"="NDMR18","-1"="blank")) %>% 
  mutate(color=rep(c("T","G","O"),each=88))

feb2_bg <- feb2_dat %>% filter(gen=="blank") %>% dplyr::summarise(across(OD600:OFP2, ~ mean(.x, na.rm = TRUE))) %>% as.vector()

feb2_dat <- feb2_dat %>% filter(gen!="blank") %>% 
  mutate(OD600=(OD600-as.numeric(feb2_bg[1]))/0.82,
         CFP=CFP-as.numeric(feb2_bg[2]),
         GFP=GFP-as.numeric(feb2_bg[3]),
         OFP=OFP-as.numeric(feb2_bg[4])) %>% 
  mutate(RFU=case_when(color=="T" ~ CFP,color=="G"~GFP,T~OFP)) %>% 
  dplyr::select(-c(CFP,GFP,OFP,OFP2,Well)) %>% mutate(day="Feb2") 



OD_RFU_dat_static <- full_join(Jan27_dat,feb1_dat,by = c("color", "gen", "OD600", "RFU", "day")) %>% full_join(.,feb2_dat) 


static_calibration = OD_RFU_dat_static %>%mutate(RFU=log(RFU),OD600=log(OD600)) %>% 
  group_by(color) %>% nest %>% 
  mutate(lm = map(data,~tidy(lm(OD600~RFU,data=.x)))) %>% unnest(lm) %>% 
  subset(select = c(color,term,estimate)) %>% 
  mutate(term = recode(term, '(Intercept)' = 'Intercept', 'RFU' = 'Slope')) %>% 
  pivot_wider(names_from=term,values_from=estimate) %>% as.data.frame 

# generating predicted values from log-log calibration to plot a line over the points
static_calib_preds = OD_RFU_dat_static %>% group_by(color) %>% dplyr::summarize(min=min(RFU),max=max(RFU)) %>% 
  full_join(static_calibration) %>% 
  mutate(RFU = map2(min,max,~seq(.x,.y,length.out=1000))) %>% unnest(RFU) %>% 
  mutate(OD600 = exp(Intercept)*RFU^Slope) %>% rename('Colour'='color') %>% 
  mutate(Colour = recode(Colour, "O" ='OFP',# "MKOkappa",
                         "G" = 'GFP',#"AausFP-1",
                         "T" = 'CFP',#"mTurquoise2"
  )) 



# # visualizing relationship
# OD_RFU_dat_static %>%# filter(OD600<1) %>% 
#   mutate(Colour = recode(color, "O" ='OFP',# "MKOkappa",
#                          "G" = 'GFP',#"AausFP-1",
#                          "T" = 'CFP',#"mTurquoise2"
#   )) %>% 
#   ggplot(aes(RFU,OD600,color=as.factor(as.numeric(as.factor(gen))))) + 
#   facet_wrap(~Colour,scales="free") + theme_classic()+
#   labs(y=expression("OD"[600]),color='Plasmid',shape="")+
#   geom_point(aes(shape = day),alpha=1) + scale_x_log10() + scale_y_log10()+
#   geom_line(data=static_calib_preds,color='grey15',linewidth=1) +
#   scale_colour_grey(start=0.3,end=0.8) +guides(shape = "none") 






write.csv(static_calibration,'../cleaned_datasets/static_calibration_curve.csv')
