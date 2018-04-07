#######################################
####### Survival Analysis ##
#######################################
# TODO: Clean data
# TODO: Figure out what time to event data to use. Both the patient and 
#       the follow up clinical info have slightly different days to death data
#       Also I think there's some progression free time data hiding somewhere
# TODO: Determine clinical variables for inclusion in the Cox model
# TODO: Get code ready

#######################################
######### Libraries ############
#######################################
library("TCGAbiolinks")
library(survival)

#######################################
######### Data Preparation ######
#######################################
#Load the clinical data
clinical.query <- GDCquery(project = "TCGA-OV", 
                           data.category = "Clinical")

# TODO: Clean, figure out what's useful
clinical <- GDCprepare_clinic(clinical.query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(clinical.query, clinical.info = "drug")
clinical.fup <- GDCprepare_clinic(clinical.query, clinical.info = "follow_up")
clinical.admin <- GDCprepare_clinic(clinical.query, clinical.info = "admin")

#######################################
######### Kaplan Meier ######
#######################################

surv.vars <- c("bcr_patient_barcode", "days_to_death", "days_to_last_followup",
               "age_at_initial_pathologic_diagnosis", "ethnicity")
surv.data <- clinical[, surv.vars]
# Define the time to event. If we have the days to death, use that. 
# Otherwise use their last follow up
surv.data$EventTime <- ifelse(!is.na(surv.data$days_to_death), 
                              surv.data$days_to_death, 
                              surv.data$days_to_last_followup)
# Create an indicator variable for whether the event is death or censoring
surv.data$censored <- ifelse(!is.na(surv.data$days_to_death), 
                              FALSE, 
                              TRUE)
# Fit the KM curves
km.surv.fit <- survfit(Surv(surv.data$EventTime, surv.data$censored) ~ surv.data$ethnicity, conf.type = "plain")
# Example KM Plot
plot(km.surv.fit, main=expression(paste("Kaplan-Meier-estimate ", hat(S)[g](t), " for groups g")),
     xlab="t", ylab="Survival", lwd=2, col=1:3)
legend(x="topright", col=1:3, lwd=2, legend=LETTERS[1:3])


#######################################
######### Cox Model ######
#######################################

# Cox model based on age. Just for example
cox.age <- coxph(data = surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis + ethnicity)
# Reference group for ethnicity is either missing or nonwhite, not sure
summary(cox.age)
