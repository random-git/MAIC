#Author: cz
library(dplyr)
library(boot)
library(survival)
library(ggplot2)
library(survminer)
library(flextable)
library(officer)
library(readxl)
library(maic)
library(sandwich)

set.seed(123)

#Pseudo AdaM data
adsl = read_xlsx("D:/Stat/Data/adsl.xlsx")
adrs = read_xlsx("D:/Stat/Data/adrs.xlsx")
adtte = read_xlsx("D:/Stat/Data/adtte.xlsx")

#The external study for comparison
#patient characteristics
df_compare = data.frame('ARM' =c(2),
                        'AGE' = c(50.2),
                        'SEX' = c(0.69))

#outcome treatment2 vs placebo
#effect size: 0.08 se:0.15 ci_lower: -0.51 ci_upper:0.45



df_all = adsl %>%
  full_join(adrs,by = "USUBJID")%>%
  full_join(adtte,by = "USUBJID")

df_index = df_all %>%
  select(USUBJID,AGE,SEX)%>%
  mutate(AGE_centered = AGE - df_compare$AGE,
         SEX_centered = SEX - df_compare$SEX)


#Overview of the data
df_all_summary = df_all %>%
  group_by(ARM)%>%
  summarize_at(c('AGE','SEX','AVAL','EVENT',"TIME"),mean, na.rm = T)

df_all_summary2 = df_all %>%
  summarize_at(c('AGE','SEX','AVAL','EVENT',"TIME"),mean, na.rm = T)%>%
  mutate(ARM = "ALL")%>%
  relocate(ARM) %>%
  rbind(df_all_summary)

#Primary endpoint
#Mean difference in ALC between treated (ARM=1) vs control (ARM=0)
t_test_result = t.test(AVAL~ARM,data=df_all)
index_outcome = data.frame(
  ARM = "Treatment 1 vs placebo (unadjusted by MAIC)",
  effect_size = as.numeric(t_test_result$estimate[2] - t_test_result$estimate[1]),
  se = (t_test_result$conf.int[2] - t_test_result$conf.int[1])/3.92,
  ci_lower = t_test_result$conf.int[1],
  ci_upper = t_test_result$conf.int[2]
)

# Set matching covariates
maic_dictionary = data.frame(
  'match.id' = c('AGE','SEX'),
  'target.variable' = c('AGE','SEX'),
  'index.variable' = c('AGE','SEX'),
  'match.type' = c('mean','mean')
)

#MAIC input matrix
maic_matrix = createMAICInput(
  index = df_index,
  target = df_compare,
  dictionary = maic_dictionary,
  matching.variables = c('AGE','SEX')
)

maic_weights = maicWeight(maic_matrix)

#Summarize the adjusted patients' characteristics using maic weights
char_summary = reportCovariates(
  index = df_index,
  target = df_compare,
  dictionary = maic_dictionary,
  matching.variables = c('AGE','SEX'),
  weights = maic_weights
)

#Recaled the weights
maic_weights_scaled = (maic_weights/sum(maic_weights))*dim(df_all)[1]

char_summary_scalled = reportCovariates(
  index = df_index,
  target = df_compare,
  dictionary = maic_dictionary,
  matching.variables = c('AGE','SEX'),
  weights = maic_weights_scaled
)

#Calcumate the effect sample size 
sum(maic_weights)^2/sum(maic_weights^2)

#linear model with MAIC weights
df_all = df_all %>%
  mutate(weights = maic_weights)

mod1 = lm(AVAL~ARM, data = df_all, weights = maic_weights)

#compare the weights results with unweights results - treatment vs placebo
index_outcome = index_outcome %>%
  add_row(
    ARM = "Treatment 1 vs placebo (adjusted by MAIC)",
    effect_size = coef(mod1)[[2]],
    se = sqrt(vcovHC(mod1)["ARM","ARM"]),
    ci_lower = effect_size-1.96*se,
    ci_upper = effect_size+1.96*se
  )

#Compared treatment
#outcome treatment2 vs placebo
#effect size: 0.08 se:0.15 ci_lower: -0.51 ci_upper:0.45

maic_comparison = rbind(
  data.frame(
    ARM = "Treatment 1 vs Treatment 2 (unadjusted by MAIC)",
    effect_size = index_outcome$effect_size[1] - 0.08,
    se = sqrt(index_outcome$se[1]^2+0.15^2)
  ) %>%
    mutate(ci_lower = effect_size - 1.96 *se,
           ci_upper = effect_size+1.96*se),
  
  data.frame(
    ARM = "Treatment 1 vs Treatment 2 (adjusted by MAIC)",
    effect_size = index_outcome$effect_size[2] - 0.08,
    se = sqrt(index_outcome$se[2]^2+0.15^2)
  ) %>%
    mutate(ci_lower = effect_size - 1.96 *se,
           ci_upper = effect_size+1.96*se)
)
