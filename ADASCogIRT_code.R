# Install required packages - mirt and lme4
install.packages('mirt')
install.packages('lme4')

# load requried packages
library(mirt)
library(lme4)

##############################################################################
# SECTION 1: TRAINING THE IRT MODEL USED FOR CALCULATING PATIENT LATENT TRAITS 
##############################################################################

# load example cross-sectional data to build IRT model (included on Git)
load(file = 'Cross_section_ADAS_Data.RData')

# Training IRT model using mirt package
# define factor loadings (F1-memory, F2-language, F3-praxis)
# All 3 traits allowed to be correlated with each other 
# (defined by COV = F1*F2, F2*F3, F1*F3)
MIRT_F3 <- mirt.model('F1 = 1, 23-29, 35, 30
                      F2 = 2-8, 10-13, 31-34
                      F3 = 14-22
                      COV = F1*F2, F2*F3, F1*F3')

# Define type of individual ADAS-Cog items 
Item_type <- c('graded',                                     # Q1 
               '2PL','2PL','2PL','2PL','2PL','2PL','3PL',    # Q2 Objects
               '2PL','2PL','2PL','2PL','2PL',                # Q2 fingers
               '2PL','2PL','2PL',                            # Q3
               '2PL','2PL','3PL',                            # Q4
               '2PL','2PL','2PL',                            # Q5
               '2PL','2PL','2PL','2PL','3PL','2PL','2PL',    # Q6
               'graded','graded','graded',                   
               'graded','graded','graded')


# Define interaction terms between patient variables and ADAS-Cog items
# Defined based on measurement invariance analysis (refer to published paper)
Q25_Gender_Recog  <- mat.or.vec(35,1);  Q25_Gender_Recog[3]<-1
Q210_Gender_Recog <- mat.or.vec(35,1); Q210_Gender_Recog[6]<-1
Q46_Gender_Recog  <- mat.or.vec(35,1);  Q46_Gender_Recog[19]<-1
Q11_Gender_Recog  <- mat.or.vec(35,1);  Q11_Gender_Recog[34]<-1

Itemdesign_F3 <- data.frame(Q25_Gender_Recog, Q210_Gender_Recog, 
                            Q46_Gender_Recog, Q11_Gender_Recog)

# redefine gender as 0,1 
Cross_Covariates$Gender <- Cross_Covariates$Gender-1

# Train Cross-sectional IRT model 
# ADAS_Cross is cross-sectional ADAS-Cog data
# Cross_Covariates are patient covariate terms
Cross_Model_F3 <- 
  mixedmirt(data = ADAS_Cross, 
            covdata = Cross_Covariates,
            model = MIRT_F3,  
            fixed =~ Gender:Q25_Gender_Recog +  Gender:Q210_Gender_Recog +
              Gender:Q46_Gender_Recog +  Gender:Q11_Gender_Recog, 
            itemtype = Item_type, 
            itemdesign = Itemdesign_F3,  
            SE = TRUE, GenRandomPars = TRUE, rotate = "oblimin",
            technical = list(NCYCLES = 4000, MAXQUAD= 20000))


################################################################
# SECTION 2: CALCULATING PATIENT LATENT TRAITS USING IRT MODEL 
################################################################

# GLMER data preparation 
# load the cognitive dysfunction scale parameters
load(file = 'Cross_IRT_Model_Results.RData')
load(file = 'Longitudinal_ADAS_Data.RData')
dat <- HU_ADAS_Long
# compile patient covariates
dat_Covariates <- HU_Covariates_Long
dat_Covariates$ARM <- as.character(dat_Covariates$ARM)
dat_Covariates$VISCODE <- as.character(dat_Covariates$VISCODE)
dat_Covariates$Gender <- as.character(dat_Covariates$Gender)
dat_Covariates$RID <- as.character(dat_Covariates$RID) # unique patient ID

# restructure data
Theta1_Slopes <- c(); Theta2_Slopes <- c(); Theta3_Slopes <- c();
Theta_Intercept <- c(); Item_Responses <- c(); 
Patient_EDU <- c(); Patient_Gender <- c(); Patient_VISCODE <- c();
Patient_RID <- c(); Patient_Visits <- c(); Patient_Arm <- c(); 
Patient_Age <- c();
Item_Model <- c(1:35)
for(it in 1:length(Item_Model)){
  item <- Item_Model[it]
  uniq_item <- sort(unique(dat[,item]))
  for(unq in 2:(length(uniq_item))){
    # item responses as separate rows
    Item_Responses <- c(Item_Responses, (dat[,item]>=uniq_item[unq])*1)
    # slopes associated with memory (theta1), language (theta2), and 
    # praxis (theta3) latent traits
    Theta1_Slopes <- c(Theta1_Slopes, 
                       (F3_Slopes_Not_CM[item, 1]*(1-dat_Covariates$CM_Med) +
                          F3_Slopes_CM[item, 1]*dat_Covariates$CM_Med))
    Theta2_Slopes <- c(Theta2_Slopes, 
                       (F3_Slopes_Not_CM[item, 2]*(1-dat_Covariates$CM_Med) +
                          F3_Slopes_CM[item, 2]*dat_Covariates$CM_Med))
    Theta3_Slopes <- c(Theta3_Slopes, 
                       (F3_Slopes_Not_CM[item, 3]*(1-dat_Covariates$CM_Med) +
                          F3_Slopes_CM[item, 3]*dat_Covariates$CM_Med))
    
    Patient_RID <- c(Patient_RID, dat_Covariates$RID)
    Patient_Arm <- c(Patient_Arm, dat_Covariates$ARM)
    Patient_Visits <- c(Patient_Visits, dat_Covariates$VIS)
    Patient_Age <- c(Patient_Age, dat_Covariates$AGE)
    Patient_VISCODE <- c(Patient_VISCODE, dat_Covariates$VISCODE)
    Patient_Gender <- c(Patient_Gender, dat_Covariates$Gender)
    Patient_EDU <- c(Patient_EDU, dat_Covariates$EDU)
    # Check the pattern of encoding of gender (COMMON MISTAKE)
    M_e3 <- which(dat_Covariates$Gender==0) 
    F_e3 <- which(dat_Covariates$Gender==1)
    temp_intercept <- mat.or.vec(length(dat[,item]),1) + NA
    temp_intercept[M_e3] <- F3_Intercepts_Male[item,(unq-1)]
    temp_intercept[F_e3] <- F3_Intercepts_Female[item,(unq-1)]
    # item intercepts
    Theta_Intercept <- c(Theta_Intercept, temp_intercept) 
  }
}
# define progression rate terms
Theta1_Rate <- Theta1_Slopes*Patient_Visits
Theta2_Rate <- Theta2_Slopes*Patient_Visits
Theta3_Rate <- Theta3_Slopes*Patient_Visits

Patient_Arm[Patient_Arm=='Placebo'] <- 'APlacebo'

# Create a data frame for usage with glmer
GLMER_F3_Data <- data.frame(Theta1_Slopes, Theta2_Slopes, Theta3_Slopes,
                            Theta1_Rate, Theta2_Rate, Theta3_Rate,
                            Patient_Arm, Patient_RID, Patient_VISCODE, 
                            Patient_EDU, Patient_Age, Patient_Gender,
                            Item_Responses, Theta_Intercept, Patient_Visits, 
                            HU200 = (Patient_Arm=='HU200ug')*1, 
                            HU400 = (Patient_Arm=='HU400ug')*1)

# standardize patient education and age
GLMER_F3_Data$Patient_EDU <- (GLMER_F3_Data$Patient_EDU-
                                mean(GLMER_F3_Data$Patient_EDU))/
  sd(GLMER_F3_Data$Patient_EDU)
GLMER_F3_Data$Patient_Age <- (GLMER_F3_Data$Patient_Age-
                                mean(GLMER_F3_Data$Patient_Age))/
  sd(GLMER_F3_Data$Patient_Age)

# estimate latent traits of patients using glmer 
HU_GLMER_F3 <- glmer(Item_Responses ~ 0 + 
                       Theta1_Slopes + Theta2_Slopes + Theta3_Slopes + 
                       Theta1_Rate + Theta2_Rate + Theta3_Rate + 
                       I(Theta1_Slopes*HU200) + 
                       I(Theta2_Slopes*HU200) + 
                       I(Theta3_Slopes*HU200) + 
                       I(Theta1_Rate*HU200) + 
                       I(Theta2_Rate*HU200) + 
                       I(Theta3_Rate*HU200) + 
                       I(Theta1_Slopes*HU400) + 
                       I(Theta2_Slopes*HU400) + 
                       I(Theta3_Slopes*HU400) + 
                       I(Theta1_Rate*HU400) + 
                       I(Theta2_Rate*HU400) + 
                       I(Theta3_Rate*HU400) + 
                       I(Theta1_Rate*Patient_Age) + 
                       I(Theta2_Rate*Patient_Age) + 
                       I(Theta3_Rate*Patient_Age) +                              
                       I(Theta1_Slopes*Patient_Age) + 
                       I(Theta2_Slopes*Patient_Age) + 
                       I(Theta3_Slopes*Patient_Age) +
                       (0 + Theta1_Slopes + Theta2_Slopes + Theta3_Slopes + 
                          Theta1_Rate + Theta2_Rate + Theta3_Rate|Patient_RID), 
                     data = GLMER_F3_Data, 
                     offset = Theta_Intercept, family = binomial, verbose = 2,
                     control = glmerControl(optimizer="bobyqa"))





