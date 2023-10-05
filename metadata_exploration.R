library(readr)
library(dplyr)
library(ggplot2)

metadata <- read_delim("metadata_qiime.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

dx_ages <- data.frame(
  Subject_ID = c('T013815','E026079','E022137','E018113','E017751','E010629','E003989','T025418','E010937','E006574','E003251'),
  Age_at_Seroconversion = c(350,580,562,587,175,945,346,540,905,532,357),
  Age_at_T1D_Diagnosis = c(NA,NA,NA,NA,NA,NA,NA,879,959,1339,1168)
)

shannon <- read_delim("alpha-diversity-shannon.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE, col_names = c('SampleID', 'Shannon_Diversity'), skip = 1)

metadata <- metadata %>% 
  left_join(dx_ages, by = 'Subject_ID') %>%
  mutate(samp_before_SC = (Age_at_Collection < Age_at_Seroconversion)) %>%
  left_join(shannon, by = 'SampleID')

sum_table <- metadata %>%
  # Statistics at the individual level
  group_by(T1D_Status, Subject_ID, HLA_Risk_Class, Gender, Delivery_Route, Age_at_Seroconversion) %>%
  summarise(N_samples_per_individual = n(), N_samples_before_SC = sum(samp_before_SC, na.rm = T), .groups = 'drop') %>%
  group_by(T1D_Status) %>%
  summarise(N_samples_total = sum(N_samples_per_individual),
            N_individuals = n_distinct(Subject_ID),
            Avg_age_at_SC = mean(Age_at_Seroconversion) %>% round(2),
            SD_age_at_SC = sd(Age_at_Seroconversion) %>% round(2),
            Avg_samples_per_individual = mean(N_samples_per_individual) %>% round(2),
            SD_samples_per_individual = sd(N_samples_per_individual) %>% round(2),
            Avg_samples_pre_SC = mean(N_samples_before_SC) %>% round(2),
            SD_samples_pre_SC = sd(N_samples_before_SC) %>% round(2),
            Gender_perc_female = 100*sum(Gender == 'female')/n() %>% round(2),
            HLA_risk_3_perc = 100*sum(HLA_Risk_Class == 3)/n() %>% round(2),
            Veginal_delivery_perc = 100*sum(Delivery_Route == 'vaginal')/n() %>% round(2))
            
sum_table %>% tibble::column_to_rownames('T1D_Status') %>% t()

# Statistics at the sample level
metadata %>%
  group_by(T1D_Status) %>%
  summarise(N = n(), Avg_Total_Reads = mean(Total_Reads), Sd_Total_Reads = sd(Total_Reads))

bf_analysis <- metadata %>% 
  select(SampleID, Subject_ID, T1D_Status, Age_at_Collection, Exclusive_BF, BF, Infant_Formula, Solid_Food, Shannon_Diversity) %>%
  # mutate(BF_Status = ifelse((Exclusive_BF|BF) & (!Solid_Food) & (!Infant_Formula), 'BF-Exclusive', 'Other')) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & BF & (!Solid_Food) & Infant_Formula, 'BF+Formula', BF_Status)) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & (!BF) & Solid_Food & Infant_Formula, 'Formula+Solid', BF_Status)) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & BF & Solid_Food & (!Infant_Formula), 'BF+Solid', BF_Status)) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & BF & Solid_Food & Infant_Formula, 'BF+Formula+Solid', BF_Status)) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & (!BF) & Solid_Food & (!Infant_Formula), 'Solid-Only', BF_Status)) %>%
  # mutate(BF_Status = ifelse(BF_Status == 'Other' & (!BF) & (!Solid_Food) & Infant_Formula, 'Formula-Only', BF_Status)) %>%
  mutate(age_bin = cut(Age_at_Collection, breaks = c(0,90,180,270,360,450,540,630,720,810,900,Inf), ordered_result = TRUE))

tmp <- bf_analysis %>%
  group_by(age_bin, T1D_Status) %>%
  summarise(N_samples = n(),
            Sum_BF = sum(BF),
            Sum_BF_Exclusive = sum(Exclusive_BF),
            Sum_BF_Inclusive = sum(BF & !Exclusive_BF),
            Sum_No_BF = sum((!BF) & (!Exclusive_BF)), .groups = 'drop') %>%
  mutate(Perc_BF = Sum_BF/N_samples) %>%
  mutate(Perc_No_BF = Sum_No_BF/N_samples) %>%
  rowwise() %>% mutate(low_ci_BF =  prop.test(Sum_BF, N_samples, conf.level=0.95)$conf.int[1]) %>%
  rowwise() %>% mutate(high_ci_BF =  prop.test(Sum_BF, N_samples, conf.level=0.95)$conf.int[2]) %>%
  rowwise() %>% mutate(low_ci_No_BF = prop.test(Sum_No_BF, N_samples, conf.level=0.95)$conf.int[1]) %>%
  rowwise() %>% mutate(high_ci_No_BF = prop.test(Sum_No_BF, N_samples, conf.level=0.95)$conf.int[2]) 

ggplot(tmp, aes(x = age_bin, y = Perc_BF, fill = T1D_Status, group = T1D_Status)) +
  geom_point(color = 'black', position=position_dodge(width=0.5), shape = 21, size = 5) +
  geom_line(aes(color = T1D_Status), position=position_dodge(width=0.5), linewidth = 1.5) +
  geom_errorbar(aes(ymin = low_ci_BF, ymax = high_ci_BF), width = 0.4, position=position_dodge(width=0.5)) +
  geom_point(color = 'black', position=position_dodge(width=0.5), shape = 21, size = 5) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  xlab('Age bin [days]') +
  ylab('Percent of breastfed infants\n(exclusive + inclusive)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) +
  theme(legend.position = c(0.85,0.7), legend.box.background = element_rect(color = 'black', fill = 'white'))
ggsave("bf_over_time.png")  
  
tmp <- bf_analysis <- metadata %>% 
  mutate(age_bin = cut(Age_at_Collection, breaks = c(0,180,360,540,720,Inf), ordered_result = TRUE)) %>%
  group_by(age_bin, BF) %>%
  summarise(N_samples = n(),
            N_no_NA = sum(!is.na(Shannon_Diversity)),
            avg_shannon = mean(Shannon_Diversity, na.rm = TRUE),
            low_shannon = quantile(Shannon_Diversity, probs = 0.025, na.rm = TRUE),
            high_shannon = quantile(Shannon_Diversity, probs = 0.975, na.rm = TRUE))


ggplot(tmp, aes(x = age_bin, y = avg_shannon, fill = BF)) +
  geom_col(color = 'black', position=position_dodge(width=0.5), width = 0.7) +
  geom_errorbar(aes(ymin = low_shannon, ymax = high_shannon), width = 0.4, position=position_dodge(width=0.5)) +
  theme_classic() +
  xlab('Age bin [days]') +
  ylab('Shannon diversity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) +
  theme(legend.position = c(0.85,0.3), legend.box.background = element_rect(color = 'black', fill = 'white'))
ggsave("shannon_over_time.png")  
  
  
  
  
  
  