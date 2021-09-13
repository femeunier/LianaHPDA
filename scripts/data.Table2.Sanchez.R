rm(list = ls())

library(dplyr)
library(tidyr)

data <- data.frame(Site = c(rep("PNM",4),rep("FTS",4)),
                   GF = rep(c("Liana","Liana","Tree","Tree"),2),
                   type = rep(c("m","sd"),4),
                   Cab = c(385.6,84.6,523.7,192.9,423.4,174.8,369.8,146.6),
                   Car = c(170.2,47.8,237.7,108.5,165.5,74.1,183.2,67.2),
                   WC = c(64.6,9.45,55.6,4.2,61.6,5.4,54.7,6.5),
                   ratio = c(0.341,0.0822,0.419,0.0604,0.471,0.049,0.477,0.0384),
                   LT = c(0.25,0.06,0.3,0.08,0.28,0.05,0.3,0.08),
                   SLA = c(14.6,6.6,10.5,4.2,9.4,2.6,7.3,1.7)) %>%
  group_by(Site,GF) %>%
  mutate(Cab = Cab*900/10000,
         Car = Car*(537 + 570)/2/10000,
         Cm = case_when(type == "m" ~ 1/(10*SLA),
                        type == "sd" ~ SLA[type == "sd"]/SLA[type == "m"]),
         # Cw = case_when(type == "m" ~ (1/ratio[type == "m"] - 1)*(1/(10*SLA[type == "m"])),
                        # type == "sd" ~ ratio[type == "sd"]/ratio[type == "m"])) %>%
         Cw = case_when(type == "m" ~ -(1/(10*SLA[type == "m"]))*(1 - 1/(1 - WC/100)),
                        type == "sd" ~ WC[type == "sd"]/WC[type == "m"])) %>%
         # Cw = case_when(type == "m" ~ (1/(10*SLA[type == "m"]))/((1 - WC/100)),
         #                type == "sd" ~ WC[type == "sd"]/WC[type == "m"])) %>%
  ungroup() %>%
  mutate(N = LT/LT[Site == "PNM" & GF == "Tree" & type == "m"]) %>%
  pivot_longer(cols = c("Cab","Car","WC","LT","SLA","Cw",'Cm',"N","ratio")) %>%
  pivot_wider(values_from = "value",
              names_from = "type") %>%
  ungroup() %>%
  mutate(sd = case_when(name %in% c("Cm","Cw") ~ sd*m,
                        TRUE ~ sd))

saveRDS(object = data,
        file = "./data/DATA_Sanchez_2009_Table2.RDS")

