setwd(dir = "make_datasets/")
library(tidyverse)
# Data loading ------------------------------------------------------------

## TURNING Data loading into a function that can be used for all the five plates

make_data = function (csv_file,design_file,background_file,which_order,which_orange) {
  csv_file %>% pivot_longer(-c(Time,Read),names_to="Well") %>% ## pivot_longer to long data format
    full_join(.,design_file,by="Well",relationship = "many-to-many") %>% filter(genotype !=-1) %>% ## Add color and genotype information from design (des) file and remove background (set to -1 in design file)
    group_by(Well) %>% mutate(Hours=Time*24) %>% dplyr::select(-Time) %>% ## SWitch to hours and remove time
    mutate(pair=case_when(any(color == "Orange" & genotype==0) ~ "CvG",
                          any(color == "Blue" & genotype==0) ~ "OvG",
                          any(color == "Green" & genotype==0) ~ "CvO"))%>% ## introduce column called "pair" which tells me which two colors are competing, e.g. Cyan vs orange will be "CvG"
    pivot_wider(names_from="Read",values_from="value") %>% ## creates three columns one for each color
    mutate(OD600=(OD600-background_file[1])/0.82,
           CFP=CFP-background_file[2],GFP=GFP-background_file[3],
           OFP=case_when(which_orange == "OFP" ~ OFP-background_file[4], # subtracting backgrounds from each column
                         T ~ OFP-background_file[5])) %>% 
    # the last case_when is for data number 5 because we found very similar but slightly better OFP wavelengths
    group_by(Well,color) %>% arrange(Well,Hours) %>% filter(genotype!=0) %>% ## Removing the genotype==0 which just tells me who isn't in the well (not interesting)
    # mutate(OD600=OD600/first(OD600)) %>% #optional normalization for OD600
    mutate(RFU=case_when(color=="Green" ~ GFP, 
                         color == "Blue" ~ CFP, T ~ OFP)) %>% ## Make a fluorescence column RFU and normalize OD to first timepoint
    mutate(genotype=case_when(which_order == "old" ~ recode(genotype,"1"="CMP", ## The order depends on whether the order is the new one or old one
                                                            "2"="VIMwt",
                                                            "3"="Junk",
                                                            "4"="VIMR18.1",
                                                            "5"="VIMR1.1",
                                                            "6"= "VIMR6.4",
                                                            "7"= "VIMR6.1",
                                                            "8"= "VIMR2.2",
                                                            "9"="VIMR2.1",
                                                            "10"="NDMR18"),
                              T ~ 
                                recode(genotype,"1"="CMP",
                                       "2"="Junk",
                                       "3"="VIMwt",
                                       "4"="VIMR1.1",
                                       "5"="VIMR2.1",
                                       "6"= "VIMR2.2",
                                       "7"= "VIMR6.1",
                                       "8"= "VIMR6.4",
                                       "9"="VIMR18.1",
                                       "10"="NDMR18"))) %>% ## Rename genotype column to actual names of genotypes
    pivot_wider(names_from=color,values_from=c(genotype,RFU)) %>% ### three colors splits genotype column into three and RFU column into three
    mutate(genA=case_when(pair == "CvG" ~ genotype_Blue,
                          pair == "CvO" ~ genotype_Blue, T ~ genotype_Orange),
           genB=case_when(pair == "CvG" ~ genotype_Green,
                          pair == "CvO" ~ genotype_Orange, T ~ genotype_Green),
           RFUA=case_when(pair == "CvG" ~ RFU_Blue,
                          pair == "CvO" ~ RFU_Blue, T ~ RFU_Orange),
           RFUB=case_when(pair == "CvG" ~ RFU_Green,
                          pair == "CvO" ~ RFU_Orange, T ~ RFU_Green)) %>% 
    dplyr::select(-c(genotype_Blue,genotype_Green,genotype_Orange,RFU_Blue,RFU_Orange,RFU_Green))
  
}

## The function takes in a dataset, a design file, a background file, asking whether the order is "old" or "new", and asking which orange wavelength was used (OFP or OFP2)

# Let's make the background file. 
bg =suppressWarnings(read.csv("plate_reader_bgs.csv")) # it gives a weird warning so I "suppress" it

# First four plates use a different OFP wavelength than the last one, this file will be just for them.
# We're first going to average all the different day's backgrounds assuming we don't need day-specific background
# Before running this code, make sure that bg columns are in the correct order that the function wants
bg = bg %>% dplyr::summarise(across(OD600:OFP2, ~ mean(.x,na.rm=T))) %>%  ## mean for each column from OD600 to OFP
  as.numeric(.)


#
des1<- read.csv("220111_des.csv",check.names = F) %>% pivot_longer(-c(Row,color),values_to="genotype") %>% 
  mutate(Well=paste(Row,name,sep="")) %>% dplyr::select(-c(Row,name))


dat1 <- make_data(read.csv("220111_dat.csv"),des1,bg,"old","OFP")

des2<- read.csv("220117_des.csv",check.names = F) %>% pivot_longer(-c(Row,color),values_to="genotype") %>% 
  mutate(Well=paste(Row,name,sep="")) %>% dplyr::select(-c(Row,name))

dat2 <- make_data(read.csv("220117_dat.csv"),des2,bg,"old","OFP")

des3<- read.csv("220126_des.csv",check.names = F) %>% pivot_longer(-c(Row,color),values_to="genotype") %>% 
  mutate(Well=paste(Row,name,sep="")) %>% dplyr::select(-c(Row,name))

dat3 <- make_data(read.csv("220126_dat.csv"),des3,bg,"new","OFP")

des4<- read.csv("220127_des.csv",check.names = F) %>% pivot_longer(-c(Row,color),values_to="genotype") %>% 
  mutate(Well=paste(Row,name,sep="")) %>% dplyr::select(-c(Row,name))

dat4 <- make_data(read.csv("220127_dat.csv"),des4,bg,"new","OFP")


raw_alldata =rbind(dat1 %>% mutate(day="Jan11",dilution=500),
                         dat2 %>% mutate(day="Jan19",dilution=500),
                         dat3 %>% mutate(day="Jan26",dilution=1250),
                         dat4 %>% mutate(day="Jan27",dilution=1250))

write.csv(raw_alldata,'../cleaned_datasets/raw_competition_dat_methodspaper.csv')
