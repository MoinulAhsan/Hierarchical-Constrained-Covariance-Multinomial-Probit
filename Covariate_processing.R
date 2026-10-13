

covariate_processing <- function(covariate) {


x1<-data.frame(covariate)
n<-nrow(x1)


#Creating dummy variables


# Create dummy variables for HEALTH (baseline: Fair/poor dropped)
HEALTH_mat<-matrix(0,nrow(as.matrix(x1$HEALTH)),ncol=(length(unique(x1$HEALTH))-3))
for(i in 1:(length(unique(x1$HEALTH))-3))
{
  HEALTH_mat[,i]<-as.numeric(x1$HEALTH==i)
}

colnames(HEALTH_mat)<-paste("HL",1:(length(unique(x1$HEALTH))-3),sep="")



# Create dummy variables for DOCTOR (baseline:MORE THAN 2 YEARS AGO dropped)
DOCTOR_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$DOCTOR))-2))
for(i in 1:(length(unique(x1$DOCTOR))-2))
{
  DOCTOR_mat[,i]<-as.numeric(x1$DOCTOR==i)
}

colnames(DOCTOR_mat)<-paste("DV",1:(length(unique(x1$DOCTOR))-2),sep="")


# BMI as a continuous covariate
BMI_mat<-as.matrix(x1$BMI)
colnames(BMI_mat)<-"BMI"

####standardize BMI column##
BMI_mat <- (BMI_mat - mean(BMI_mat, na.rm = TRUE)) / sd(BMI_mat, na.rm = TRUE)



# Create dummy variables for employment status (baseline: NOT EMPLOYED dropped)
EMPSTATUS_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$EMPSTATUS))-2))
for(i in 1:(length(unique(x1$EMPSTATUS))-2))
{
  EMPSTATUS_mat[,i]<-as.numeric(x1$EMPSTATUS==i)
}

colnames(EMPSTATUS_mat)<-paste("EM",1:(length(unique(x1$EMPSTATUS))-2),sep="")





# Create dummy variables for marital status (baseline: NEVER MARRIED dropped)
MARITAL_mat <- matrix(0, nrow = nrow(x1), ncol = 2)

# assign dummy values
MARITAL_mat[, 1] <- as.numeric(x1$MARITAL == 1)
MARITAL_mat[, 2] <- as.numeric(x1$MARITAL %in% c(2, 3, 4)) # SEPARATED, DIVORCED, and WIDOWED combined as previously married

# name columns
colnames(MARITAL_mat) <- c("CM", "PM")




# Create dummy variables for last grade completed (baseline: (graduate degree)=SOME GRADUATE SCHOOL/GRADUATE DEGREE (PHD, MD, JD, MASTER'S DEGREE) dropped)
EDUCATION_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$EDUCATION))-3))
for(i in 1:(length(unique(x1$EDUCATION))-3))
{
  EDUCATION_mat[,i]<-as.numeric(x1$EDUCATION==i)
}

# combined LESS THAN HIGH SCHOOL GRADUATE, HIGH SCHOOL GRADUATE, and TECHNICAL/TRADE SCHOOL as <college degree
# combined SOME COLLEGE, COLLEGE GRADUATE (ASSOCIATES/COMMUNITYCOLLEGE, BA, OR BS) as college degree
EDUCATION_mat<-cbind(rowSums(EDUCATION_mat[,c(1,2,3)]), rowSums(EDUCATION_mat[,c(4,5)]) )

colnames(EDUCATION_mat)<-paste("ED",1:(length(unique(x1$EDUCATION))-6),sep="")




# Create dummy variables for age group (baseline:61 - 70,>=70 as >=61 dropped)
AGEGRP_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$AGEGRP))-3))
for(i in 1:(length(unique(x1$AGEGRP))-3))
{
  AGEGRP_mat[,i]<-as.numeric(x1$AGEGRP==i)
}

# combined 18 - 30 , 31 - 40 =18-40 (AG1), and  41 - 50 , 51-60 = 41-60 (AG2)
AGEGRP_mat<-cbind(rowSums(AGEGRP_mat[,c(1,2)]), rowSums(AGEGRP_mat[,c(3,4)]) )

colnames(AGEGRP_mat)<-paste("AG",1:(length(unique(x1$AGEGRP))-5),sep="")




# Create dummy variables for living place (baseline: RURAL AREA dropped)
LIVEPL_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$LIVEPL))-2))
for(i in 1:(length(unique(x1$LIVEPL))-2))
{
  LIVEPL_mat[,i]<-as.numeric(x1$LIVEPL==i)
}

colnames(LIVEPL_mat)<-paste("LP",1:(length(unique(x1$LIVEPL))-2),sep="")




# Create dummy variables for Race (baseline: OTHER dropped)
RACETH_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$RACETH))-2))
for(i in 1:(length(unique(x1$RACETH))-2))
{
  RACETH_mat[,i]<-as.numeric(x1$RACETH==i)
}

colnames(RACETH_mat)<-paste("RE",1:(length(unique(x1$RACETH))-2),sep="")




# Create dummy variables for Income (baseline: $75,000 TO UNDER $100,000 / $100,000 TO UNDER $150,000 / $150,000 OR MORE as  >75000 dropped)
INCOME_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$INCOME))-4))
for(i in 1:(length(unique(x1$INCOME))-4))
{
  INCOME_mat[,i]<-as.numeric(x1$INCOME==i)
}

#combined UNDER $10,00/ $10,000 TO UNDER $20,000/ $20,000 TO UNDER $30,000 as <30000 (HI1), and
# combined $30,000 TO UNDER $40,000 / $40,000 TO UNDER $50,000 / $50,000 TO UNDER $75,000 as 30000-75000 (HI2)
INCOME_mat<-cbind(rowSums(INCOME_mat[,c(1,2,3)]), rowSums(INCOME_mat[,c(4,5,6)]) )

colnames(INCOME_mat)<-paste("HI",1:(length(unique(x1$INCOME))-8),sep="")



# Create dummy variables for gender (baseline: FEMALE dropped)
GENDER_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$GENDER))-2))
for(i in 1:(length(unique(x1$GENDER))-2))
{
  GENDER_mat[,i]<-as.numeric(x1$GENDER==i)
}

colnames(GENDER_mat)<-"ML"



# Create dummy variables for census region (baseline: WEST dropped)
CENSUS_REGION_mat<-matrix(0,nrow=n,ncol=(length(unique(x1$CENSUS_REGION))-1))
for(i in 1:(length(unique(x1$CENSUS_REGION))-1))
{
  CENSUS_REGION_mat[,i]<-as.numeric(x1$CENSUS_REGION==i)
}

colnames(CENSUS_REGION_mat)<-paste("CR",1:(length(unique(x1$CENSUS_REGION))-1),sep="")





dummy_list<-list(HEALTH_mat , DOCTOR_mat ,BMI_mat ,EMPSTATUS_mat , MARITAL_mat , EDUCATION_mat,
                 AGEGRP_mat,LIVEPL_mat, RACETH_mat, INCOME_mat ,GENDER_mat, CENSUS_REGION_mat )


final_mat <- do.call(cbind, dummy_list)





X2<-as.matrix(final_mat)
p<-ncol(X2)



#######################
#Check the multicollinearity  based on VIF
MAX_VIF<-c()

for(j in 1:p)
{
  
  VIF<-c()
  for(i in 1:(p))
  {
    aa<-lm(X2[,i]~X2[,-c(i)])
    R2<-summary.lm(aa)$r.squared
    VIF[i]<- 1/(1-R2)
    
  }
  
  if(sum(as.numeric(VIF>5))< 1) {break}
  
  MAX_VIF[j]<-which.max(VIF)
  
  X2<-X2[,-MAX_VIF[j]]
  
  p<-ncol(X2)
  
}


# removed RE1 because of the multicollinearity 
if(is.null(MAX_VIF)==F)
{
  final_mat<-final_mat[,-(MAX_VIF)]
}


return(final_mat)
}
