library(tidyverse)
library(bacondecomp)
library(haven)
library(lfe)
library(lmtest)
library(did)
library(didimputation)
library(DIDmultiplegt) 
library(broom)
library(fixest)

setwd("/Users/fernandasobrino/Documents/Huachicol")

base <- read.csv("base_all.csv")

#### Base Specification: 
base$time_to_treatment <- base$diffYears %>% replace(is.na(.), -1)

base <- base %>% mutate(across(c("lead3","lead2", "lead1", "lag0", 
                                 "lag1", "lag2", "lag3", "lag4",
                                 "lag5","lag6","lag7","lag8",
                                 "lag9"), function(x) as.integer(as.logical(x))))

did_simple <- felm(log_hom ~ post|CVEGEO + YEAR|0|CVEGEO, base, 
                   weights = base$pobtot,
                   keepCX = TRUE)
summary(did_simple)

event_study <- felm(log_hom~lead2 + lead3 +
                      lag0 + lag1 + lag2 + lag3 + lag4 + lag5 +
                      lag6 + lag7 + lag8 + lag9|CVEGEO+YEAR|0|CVEGEO,base,
                    weights = base$pobtot,
                    keepCX = TRUE)
summary(event_study)

plot_order <- c("lead3", "lead2", "lag0",
                "lag1", "lag2", "lag3",
                "lag4", "lag5", "lag6",
                "lag7", "lag8", "lag9")

leadslags_plot <- tibble(
  sd = c(event_study$cse[plot_order], 0),
  mean = c(coef(event_study)[plot_order], 0),
  label = c(-3, -2, 0, 1,2,3,4,5,6,7,8,9, -1)
)

leadslags_plot$label <- factor(leadslags_plot$label)
leadslags_plot %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = mean-1.96*sd,
                    ymax = mean+1.96*sd),
                colour = "#619CFF", width = 0.2) +
  scale_x_discrete(breaks = leadslags_plot$label) + 
  geom_hline(yintercept = 0,color="grey") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") + 
  theme_bw() +
  annotate("label", x = "2", y = 0.9,
           label ="Post(DD): 0.22**",size=3)+
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap') 



#### CS main specification: 
did_att <- att_gt(yname = "log_hom",
                  gname = "YearFirst",
                  idname = "CVEGEO",
                  tname = "YEAR",
                  clustervars   = "CVEGEO",
                  weightsname = 'pobtot',
                  panel = FALSE,
                  base_period = "varying",
                  data = base,
                  bstrap = TRUE,
                  cband = TRUE)
summary(did_att)


event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
summary(event_study_2)
ggdid(event_study_2)


leadslags_plot_2 <- tibble(
  sd = event_study_2$se.egt,
  mean = event_study_2$att.egt,
  ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
  ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
  label = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
)

leadslags_plot_2$label <- factor(leadslags_plot_2$label)

leadslags_plot_2 %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                colour = "#619CFF", width = 0.2) +
  geom_hline(yintercept = 0,color="gray") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
  theme_bw() +
  annotate("label", x = "2", y = 1.2,
           label ="Overall ATT's CS: 0.37*** ",size=3) +
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap')


## Los otros estimadores: 

static <- did_imputation(data = base, 
                         yname = "log_hom", 
                         gname = "YearFirst", tname = "YEAR", 
                         idname = "CVEGEO")

static

es <- did_imputation(data = base, 
                     yname = "log_hom", 
                     gname = "YearFirst",
                     tname = "YEAR", idname = "CVEGEO",
                     horizon=TRUE,
                     pretrends = TRUE,
                     wname = 'pobtot')

borusyak <- es %>%
  filter(as.numeric(term) >= -3 & as.numeric(term) <= 9) %>%
  select(-1,-5,-6) %>%
  mutate(term = as.numeric(term),
         Estimator = 'Borusyak, Jaravel, and Spiess')

event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)

leadslags_plot_2 <- tibble(
  std.error = event_study_2$se.egt,
  estimate = event_study_2$att.egt,
  ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
  ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
  term = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
) %>%
  mutate(Estimator = 'Callaway and  Sant’Anna')

borusyak <- borusyak %>%
  mutate(ci_lower = estimate - 1.96*std.error,
         ci_upper = estimate + 1.96*std.error)

all_estimates <- rbind(leadslags_plot_2, borusyak)

all_estimates[is.na(all_estimates)] = 0
all_estimates$term <- factor(all_estimates$term)

all_estimates %>%
  ggplot(aes(x=term, y=estimate, color=Estimator)) +
  geom_point(aes(shape=Estimator), position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=ci_lower, 
                    ymax=ci_upper),position=position_dodge(0.5),
                width=0.25) + 
  geom_hline(yintercept = 0,color="grey") + 
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
  labs(y = "log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap') + 
  theme_bw() +
  theme(legend.position="bottom") + 
  annotate("label", x = "2", y = 1.35,
           label ="Overall ATT's CS:  0.37* \n BJS mean estimate: 0.44*" ,size=3) +
  
  ylim(c(-1.5, 1.5)) 

## Covariates: 
did_att <- att_gt(yname = "log_hom",
                  gname = "YearFirst",
                  idname = "CVEGEO",
                  tname = "YEAR",
                  clustervars   = "CVEGEO",
                  weightsname = 'pobtot',
                  panel = FALSE,
                  xformla = ~pobhlindigena + pobanalf + gradoescolar,
                  base_period="varying",
                  data = base,
                  bstrap = TRUE,
                  cband = TRUE)
                  #panel = FALSE)
summary(did_att)

event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
summary(event_study_2)
ggdid(event_study_2)


leadslags_plot_2 <- tibble(
  sd = event_study_2$se.egt,
  mean = event_study_2$att.egt,
  ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
  ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
  label = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
)

leadslags_plot_2$label <- factor(leadslags_plot_2$label)

leadslags_plot_2 %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                colour = "#619CFF", width = 0.2) +
  geom_hline(yintercept = 0,color="gray") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
  theme_bw() +
  annotate("label", x = "2", y = 1.2,
           label ="Overall ATT's CS: 0.35*** ",size=3) +
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap')




### Mecanismos 
### 1. Number of cartesls involved in hauchicol
base <- base %>% mutate(otros_carteles = num_cartels - carteles_huachicol)

did_att <- att_gt(yname = "carteles_huachicol",
                  gname = "YearFirst",
                  idname = "CVEGEO",
                  tname = "YEAR",
                  clustervars   = "CVEGEO",
                  panel = FALSE,
                  base_period = "varying",
                  data = base,
                  bstrap = TRUE,
                  cband = TRUE)
summary(did_att)


event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
summary(event_study_2)
ggdid(event_study_2)


leadslags_plot_2 <- tibble(
  sd = event_study_2$se.egt,
  mean = event_study_2$att.egt,
  ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
  ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
  label = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
)

leadslags_plot_2$label <- factor(leadslags_plot_2$label)

leadslags_plot_2 %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                colour = "#619CFF", width = 0.2) +
  geom_hline(yintercept = 0,color="gray") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
  theme_bw() +
  annotate("label", x = "2", y = 2,
           label ="Overall ATT's CS: 0.49*** ",size=3) +
  labs(y="Cartels involved in oil theft",
       x = 'Year since first illegal tap')





## 2. States with and without military action 

ents_calderon <- c(2, 8, 10, 12, 16, 19, 25, 28)

## Con operativo Militar: 
base_operativo <- base %>% filter(ENT %in% ents_calderon)

did_simple <- felm(log_hom ~ post|CVEGEO + YEAR|0|CVEGEO, base_operativo, 
                   weights = base_operativo$pobtot,
                   keepCX = TRUE)
summary(did_simple)

mod_twfe = feols(log_hom ~ post|                    ## Other controls
                   YEAR + CVEGEO,                             
                 cluster = ~CVEGEO, 
                 weights = base_sin_operativo$pobtot,
                 data = base_sin_operativo)
summary(mod_twfe)

mod_twfe = feols(log_hom ~ i(time_to_treatment, treated, ref = -1) |                    ## Other controls
                   YEAR + CVEGEO,                             
                 cluster = ~CVEGEO,
                 weights = base_sin_operativo$pobtot,
                 data = base_sin_operativo)

iplot(mod_twfe, 
      xlab = 'Months since earthquake')


event_study <- felm(log_hom~lead2 + lead3 +
                      lag0 + lag1 + lag2 + lag3 + lag4 + lag5 +
                      lag6 + lag7 + lag8 + lag9|CVEGEO+YEAR|0|CVEGEO,base_operativo,
                    weights = base_operativo$pobtot,
                    keepCX = TRUE)
summary(event_study)

plot_order <- c("lead3", "lead2", "lag0",
                "lag1", "lag2", "lag3",
                "lag4", "lag5", "lag6",
                "lag7", "lag8", "lag9")

leadslags_plot <- tibble(
  sd = c(event_study$cse[plot_order], 0),
  mean = c(coef(event_study)[plot_order], 0),
  label = c(-3, -2, 0, 1,2,3,4,5,6,7,8,9, -1)
)

leadslags_plot$label <- factor(leadslags_plot$label)
leadslags_plot %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = mean-1.96*sd,
                    ymax = mean+1.96*sd),
                colour = "#619CFF", width = 0.2) +
  scale_x_discrete(breaks = leadslags_plot$label) + 
  geom_hline(yintercept = 0,color="grey") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") + 
  theme_bw() +
  annotate("label", x = "2", y = 1,
           label ="Post(Average): 0.25",size=3)+
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap') 





did_att <- att_gt(yname = "log_hom",
                  gname = "YearFirst",
                  idname = "CVEGEO",
                  tname = "YEAR",
                  clustervars   = "CVEGEO",
                  weightsname = 'pobtot',
                  panel = FALSE,
                  base_period = "varying",
                  data = base_operativo,
                  bstrap = TRUE,
                  cband = TRUE)
summary(did_att)


event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
summary(event_study_2)
ggdid(event_study_2)

leadslags_plot_2 <- tibble(
  sd = event_study_2$se.egt,
  mean = event_study_2$att.egt,
  ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
  ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
  label = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
)

leadslags_plot_2$label <- factor(leadslags_plot_2$label)

leadslags_plot_2 %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                colour = "#619CFF", width = 0.2) +
  geom_hline(yintercept = 0,color="gray") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
  theme_bw() +
  annotate("label", x = "2", y = 1.4,
           label ="Overall ATT's CS: 0.31 ",size=3) +
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap')

### Sin Operativo militar 

base_sin_operativo <- base %>% filter(!(ENT %in% ents_calderon))

did_simple <- felm(log_hom ~ post|CVEGEO + YEAR|0|CVEGEO, base_sin_operativo, 
                   weights = base_sin_operativo$pobtot,
                   keepCX = TRUE)
summary(did_simple)

event_study <- felm(log_hom~lead2 + lead3 +
                      lag0 + lag1 + lag2 + lag3 + lag4 + lag5 +
                      lag6 + lag7 + lag8 + lag9|CVEGEO+YEAR|0|CVEGEO,base_sin_operativo,
                    weights = base_sin_operativo$pobtot,
                    keepCX = TRUE)
summary(event_study)

plot_order <- c("lead3", "lead2", "lag0",
                "lag1", "lag2", "lag3",
                "lag4", "lag5", "lag6",
                "lag7", "lag8", "lag9")

leadslags_plot <- tibble(
  sd = c(event_study$cse[plot_order], 0),
  mean = c(coef(event_study)[plot_order], 0),
  label = c(-3, -2, 0, 1,2,3,4,5,6,7,8,9, -1)
)

leadslags_plot$label <- factor(leadslags_plot$label)
leadslags_plot %>%
  ggplot(aes(x = label, y = mean)) +
  geom_point(color = '#619CFF') +
  geom_errorbar(aes(ymin = mean-1.96*sd,
                    ymax = mean+1.96*sd),
                colour = "#619CFF", width = 0.2) +
  scale_x_discrete(breaks = leadslags_plot$label) + 
  geom_hline(yintercept = 0,color="grey") +
  geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") + 
  theme_bw() +
  #annotate("label", x = "2", y = 1,
           #label ="Post(Average): 0.25",size=3)+
  labs(y="log Homicides x 100,000 inhabitants",
       x = 'Year since first illegal tap') 

did_att <- att_gt(yname = "log_hom",
                  gname = "YearFirst",
                  idname = "CVEGEO",
                  tname = "YEAR",
                  clustervars   = "CVEGEO",
                  weightsname = 'pobtot',
                  panel = FALSE,
                  base_period = "varying",
                  data = base_sin_operativo,
                  bstrap = TRUE,
                  cband = TRUE)
summary(did_att)


event_study_2 <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
summary(event_study_2)
ggdid(event_study_2)

# leadslags_plot_2 <- tibble(
#   sd = event_study_2$se.egt,
#   mean = event_study_2$att.egt,
#   ci_lower = event_study_2$att.egt - event_study_2$crit.val.egt * event_study_2$se.egt,
#   ci_upper = event_study_2$att.egt + event_study_2$crit.val.egt * event_study_2$se.egt,
#   label = c(-3, -2,-1, 0, 1,2,3,4,5, 6, 7, 8 , 9)
# )
# 
# leadslags_plot_2$label <- factor(leadslags_plot_2$label)
# 
# leadslags_plot_2 %>%
#   ggplot(aes(x = label, y = mean)) +
#   geom_point(color = '#619CFF') +
#   geom_errorbar(aes(ymin = ci_lower,
#                     ymax = ci_upper),
#                 colour = "#619CFF", width = 0.2) +
#   geom_hline(yintercept = 0,color="gray") +
#   geom_vline(xintercept = which(levels(leadslags_plot$label) == "0"), color = "grey") +
#   theme_bw() +
#   annotate("label", x = "2", y = 1.3,
#            label ="Overall ATT's CS: 0.25*** ",size=3) +
#   labs(y="log Homicides x 100,000 inhabitants",
#        x = 'Year since first illegal tap')


### Robustness Checks: ir quitanado estados 

perform_analysis <- function(data, state_to_remove) {
  data_filtered <- data %>% filter(ESTADO != state_to_remove)
  did_att <- att_gt(yname = "log_hom",
                    gname = "YearFirst",
                    idname = "CVEGEO",
                    tname = "YEAR",
                    clustervars = "CVEGEO",
                    weightsname = 'pobtot',
                    panel = FALSE,
                    base_period = "varying",
                    data = data_filtered,
                    bstrap = TRUE,
                    cband = TRUE)
  event_study <- aggte(did_att, type = "dynamic", na.rm = TRUE,
                       min_e = -3, max_e = 9)
  leadslags_plot <- tibble(
    sd = event_study$se.egt,
    mean = event_study$att.egt,
    label = seq(min(event_study$egt), max(event_study$egt))
  )
  
  return(leadslags_plot)
}

unique_states <- unique(base$ESTADO)

results_list <- list()

for (state in unique_states) {
  result <- perform_analysis(base_operativo, state)
  result <- result %>% mutate(state = state)  
  results_list[[state]] <- result
}

combined_results <- bind_rows(results_list)

combined_results <- combined_results %>% filter(label == 1)


ggplot(combined_results, aes(x = factor(state), y = mean)) +
  geom_point(color = "#619CFF") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "#619CFF") +
  labs(title = "Robustness: leaving one state out at a time",
       x = "State taken out",
       y = "CS Event-study (average post period)") +
  scale_x_discrete(drop = FALSE) +  # Ensure all x-axis values are shown
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # Rotate x-axis labels for better readability
  )
