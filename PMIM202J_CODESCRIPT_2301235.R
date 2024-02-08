#Connecting to Db2 
source("P:/okekec/PMIM202J/login.R")

# Loading Libraries for analysis 
library(tidyverse)
library(lubridate)
library(survival)
library(scales)
library(caret)
library(rpart)
library(survminer)

############  NOTICE #################
# LINE 202 CONTAINS THE CLEANED DATATSET USED FOR ANALYSIS

# WLGP DATA 
Gp_data <- sqlQuery(channel, "
     SELECT ALF_PE , GNDR_CD , WOB, EVENT_CD , EVENT_DT
     FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
      WHERE EVENT_CD LIKE 'H33%'AND EVENT_DT >= '2020-01-01' 
      AND EVENT_DT <= '2022-12-31'
              ")

#Defining severity  --- severe == 0, else 1 and removing duplicates 
Gp_data <- Gp_data %>%
  distinct(ALF_PE, .keep_all = TRUE) %>%
  mutate(severity_status = 
           if_else(EVENT_CD == 'H333.' 
                   | EVENT_CD == 'H334.' | EVENT_CD == 'H335.'| 
                EVENT_CD == 'H333' | EVENT_CD == 'H334' 
                | EVENT_CD == 'H335' , 0,1))

#adding the age using the wob 
w_o_b <- as.numeric(as.Date(Gp_data$WOB))
END_DATE <- as.numeric(as.Date('2022-12-31'))
Age <- round(END_DATE - w_o_b, 0)
Age <- round(Age/365.25, 0)
Gp_data$Age <- Age


# sex 0 == male, 1 == female 
Gp_data <- Gp_data %>%
  mutate(Sex = ifelse(GNDR_CD == 1, 0, 1)) %>%
  select(-c(GNDR_CD, WOB)) %>%
  select(c(ALF_PE, EVENT_CD, EVENT_DT, Age, Sex, severity_status))

#####################COVID #################

gp_covid <- sqlQuery(channel, "
     SELECT ALF_PE , EVENT_CD AS EVENT_C 
     FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
      WHERE EVENT_CD = 'A795.'AND EVENT_DT >= '2020-01-01'
              ")

gp_covid <- gp_covid %>%
  distinct(ALF_PE, .keep_all = TRUE)

# Covid status 0 == no covid , 1 == covid
general_tab <- Gp_data %>%
  left_join(gp_covid, by = 'ALF_PE', relationship = 'many-to-many') %>%
  mutate(covid_status = ifelse(is.na(EVENT_C), 0 , 1)) %>%
  select(-EVENT_C)

table(general_tab$covid_status)

# Death for Censoring 

WDSD <- sqlQuery(channel, "
     SELECT ALF_PE , DOD 
     FROM SAIL1281V.WDSD_AR_PERS_20220307 wap 
     WHERE DOD >= '2020-01-01'
              ")

WDSD <- WDSD %>%
  distinct(ALF_PE, .keep_all = TRUE)

general_tab <- general_tab %>%
  left_join(WDSD, by = 'ALF_PE')

################# SMOKING ############################
#Read code for current smokers 137R. 

smoking <- sqlQuery(channel, "
     SELECT ALF_PE , EVENT_CD AS EVENT_C2
     FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
     WHERE EVENT_CD = '137R.' AND EVENT_DT >= '2020-01-01'
              ")

smoking <- smoking %>%
  distinct(ALF_PE, .keep_all = TRUE)
 

general_tab <- general_tab %>%
  left_join(smoking, by = 'ALF_PE') %>%
  mutate(smoking = ifelse(is.na(EVENT_C2), 0, 1)) %>%
  select(-EVENT_C2)

table(general_tab$smoking)

###################### Obesity ##########################

obesity <- sqlQuery(channel, "
     SELECT ALF_PE , EVENT_CD AS EVENT_C3
     FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
     WHERE EVENT_CD = 'C380.' AND EVENT_DT >= '2020-01-01'
              ")

obesity <- obesity %>%
  distinct(ALF_PE, .keep_all = TRUE)


general_tab <- general_tab %>%
  left_join(obesity, by = 'ALF_PE') %>%
  mutate(obesity = ifelse(is.na(EVENT_C3), 0, 1)) %>%
  select(-EVENT_C3)

table(general_tab$obesity)

################ SINUSITIS ################################

sinusitis <- sqlQuery(channel, "
    SELECT ALF_PE , EVENT_CD AS EVENT_C4
    FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
    WHERE EVENT_CD LIKE 'H13%' AND EVENT_DT >= '2020-01-01'
              ")

sinusitis <- sinusitis %>%
  distinct(ALF_PE, .keep_all = TRUE)

general_tab <- general_tab %>%
  left_join(sinusitis, by = 'ALF_PE') %>%
  mutate(sinusitis = ifelse(is.na(EVENT_C4), 0, 1)) %>%
  select(-c(EVENT_C4))

table(general_tab$sinusitis)

############ ALLERGIC RHINITIS ###############################
a_rhinitis <- sqlQuery(channel, "
    SELECT ALF_PE , EVENT_CD AS EVENT_C5
    FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
    WHERE EVENT_CD LIKE 'H17%' AND EVENT_DT >= '2020-01-01'
              ")

a_rhinitis <- a_rhinitis %>%
  distinct(ALF_PE, .keep_all = TRUE)

general_tab <- general_tab %>%
  left_join(a_rhinitis, by = 'ALF_PE') %>%
  mutate(a_rhinitis = ifelse(is.na(EVENT_C5), 0, 1)) %>%
  select(-EVENT_C5)

table(general_tab$a_rhinitis)

####### Exploratory Analysis ################################

str(general_tab)


summary(general_tab)


####################### Cleaning ###########################
# To calculate time-to-event (asthma attack), we have to restrict 
# the dataset to the study  period where study start date  
# is 2020-01-01 and study end date is 2022-12-31

general_tab <- general_tab %>%
  mutate(START_DATE = as.Date("2020-01-01")) %>%
  mutate(END_DATE = as.Date("2022-12-31"))

# consideration for those with error (OUTLIERS) in dates above the 
# end of study period "2022-12-31"
temp <- general_tab %>%
  mutate(EVENT_DT = if_else(EVENT_DT > as.Date('2021-12-31'), 
                            as.Date('2022-12-31'), 
                           EVENT_DT))

gen_tab <- temp %>%
  filter(Age < 100)

#Creating the time variable for the study 
gen_tab <-  gen_tab %>%
  mutate(study_period = ifelse(is.na(EVENT_DT),
                               as.integer(difftime(
                                 END_DATE , START_DATE, units = "days")),
                               as.integer(difftime(
                                 EVENT_DT, START_DATE, units = "days")))) %>%
  select(-c(EVENT_DT, START_DATE, END_DATE))

max(gen_tab$study_period)

# checking once more for NA
gen_tab<- gen_tab %>%
   filter(!is.na(ALF_PE))
 
# saving as a csv file 

#write.csv(gen_tab, file = 'asthma_data.csv', row.names = FALSE)

#IMporting the cleaned dataset for analysis 

asthma_data <- read_csv("asthma_data.csv")


#################### decision tree ##################################
asthma_data1 <- asthma_data

par(mfrow = c(1,1), xpd = NA)

asthma_data1$Age <- scale(asthma_data1$Age)


decision_tree <- rpart(severity_status ~ Age + Sex + covid_status 
                       + smoking + obesity + sinusitis + a_rhinitis, 
                       data = asthma_data1, 
                       control = rpart.control(minsplit = 10, cp = 0.0007))
plot(decision_tree)
text(decision_tree, digits = 2, use.n = TRUE)


printcp(decision_tree)

summary(decision_tree)

############# Logistics Regression ################################

model <- glm(severity_status ~ Age + Sex + covid_status + smoking + obesity + 
               sinusitis + a_rhinitis, 
             data = asthma_data, family = binomial(link = "logit"))

summary(model)


model2 <- glm(severity_status ~ Age + Sex + smoking + obesity + 
                sinusitis + a_rhinitis, 
              data = asthma_data, family = binomial(link = "logit"))

summary(model2)

# is there an interaction between these factors ?

model3 <- glm(severity_status ~ Age * Sex * covid_status , 
              data = asthma_data, family = binomial)

summary(model3)

############## Survival Analysis #############################

#creating the survival object of severity and time 
surv_obj <- with(asthma_data, Surv(study_period, severity_status))
surv_obj

surv_model <- survfit(surv_obj~ 1)

print(surv_model)

# plotting the surv model 
cox_plot <- ggsurvplot(surv_model, data = asthma_data)
cox_plot

mytable <- table(asthma_data$severity_status, asthma_data$Sex)
row.names(mytable) <- c("Severe", "Not Severe")
colnames(mytable) <- c("Male", "Female")
mytable

# Survival based on sex
model71 <- survfit(surv_obj~ Sex, data = asthma_data)

print(model71)

g_plot <- ggsurvplot(model71, data = asthma_data)
g_plot

#saving the plot 

png("Survival_plot.png", width = 6, height = 4, units = 'in', res = 300)

print(g_plot$plot)

dev.off()


# setting up a cox regression model

cox_model <- coxph(surv_obj ~ Age + Sex + covid_status + smoking + obesity + 
                  sinusitis + a_rhinitis , data = asthma_data)

summary(cox_model)

################# TIME-T0-EVENT BETWEEN AGE GROUPS ################
# GROUP A == 0 - 20 
# GROUP B == 21 - 40
# GROUP C == 41 - 60
# GROUP D == > 60

asthma_data <- asthma_data %>%
  mutate(Age_Group = 'NA') %>%
  mutate(Age_Group = if_else(Age >= 0 & Age <= 20 , 'Group A', Age_Group)) %>%
  mutate(Age_Group = if_else(Age >= 21 & Age <= 40 , 'Group B', Age_Group)) %>%
  mutate(Age_Group = if_else(Age >= 41 & Age <= 60 , 'Group C', Age_Group)) %>%
  mutate(Age_Group = if_else(Age >= 61 , 'Group D', Age_Group))

# plotting the surv model 

mytable <- table(asthma_data$severity_status, asthma_data$Age_Group)
row.names(mytable) <- c("Severe Asthma", "Asthma")
colnames(mytable) <- c("Group A", "Group B", "Group C", "Group D")
mytable

model72 <- survfit(surv_obj~ Age_Group , data = asthma_data)

print(model72)

ggsurvplot(model72, data = asthma_data)


############### END OF PROJECT EXPLORATORY ANALYSIS #############

cat("The total number of people with asthma in wales is ", nrow(asthma_data), '\n')

cat("The min and max age of people in this study is ", min(asthma_data$Age),
    "and", max(asthma_data$Age), '\n')

cat("The total number of people with severe asthma in wales is ",
    sum(asthma_data$severity_status == 0), '\n')

cat("The total number of astha patients with Covid is ",
    sum(asthma_data$covid_status == 1), '\n')

cat("The total number of astha patients that smokes is ",
    sum(asthma_data$smoking == 1), '\n')


cat("The total number of astha patients with obesity is ",
    sum(asthma_data$obesity == 1), '\n')
  
cat("The total number of astha patients with sinusitis is ",
    sum(asthma_data$sinusitis == 1), '\n')

cat("The total number of astha patients with a_rhinitis is ",
    sum(asthma_data$a_rhinitis == 1), '\n')



