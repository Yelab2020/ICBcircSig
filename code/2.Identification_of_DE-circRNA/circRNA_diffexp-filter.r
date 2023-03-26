library(nlme)
library(tidyverse)
library(data.table)
library(ggplot2)
library(magrittr)
library(gridExtra)
library(grid)
library(gridBase)
library(DESeq2)
library(dplyr)
library(miscTools)


#lme function
runlme.My = function(circID,datRPM,datTraits,label){
  oneRPM <- datRPM %>% 
    dplyr::select(ID,circID) %>% 
    dplyr::rename(circID=circID)
	mode(oneRPM$circID)='numeric' 

  one_dat <- oneRPM %>% 
    left_join(datTraits, by="ID")
  one_dat$label=label;
  one_dat$label <- relevel(factor(one_dat$label), ref ="T1") 

  flme <- lme(asin(circID) ~ label , data=one_dat, random = ~1|Annotate, na.action = na.exclude)
  slme <- summary(flme)
  effect_lme <- rownames(slme$tTable)
  
  output_lme <- as_tibble(slme$tTable) %>% 
    add_column(effect_lme) %>% 
    dplyr::select(effect_lme, everything()) %>% 
    dplyr::filter(effect_lme=="labelT2") 
  
  output_lme
}

#===================input 1.clinicl data
PRJEB_Ann <- read.delim("/work/dy/project/2020CircRNA/temp_github/data/SampleClinicInformation.txt")
preBenefitID <- PRJEB_Ann[PRJEB_Ann$Benefit=="Benefit" & PRJEB_Ann$PrePost=="PRE",]$SRR
preNonBenefitID <- PRJEB_Ann[PRJEB_Ann$Benefit=="NonBenefit" & PRJEB_Ann$PrePost=="PRE",]$SRR
postBenefitID <- PRJEB_Ann[PRJEB_Ann$Benefit=="Benefit" & PRJEB_Ann$PrePost=="EDT",]$SRR
postNonBenefitID <- PRJEB_Ann[PRJEB_Ann$Benefit=="NonBenefit" & PRJEB_Ann$PrePost=="EDT",]$SRR

PRJEB_Ann_LME=PRJEB_Ann
PRJEB_Ann_LME=dplyr::rename(PRJEB_Ann_LME,'ID'=SRR)
PRJEB_Ann_LME$SRR=PRJEB_Ann_LME$ID

#===================input 2.circRNA data
all_CircDataMF_LME=readr::read_rds("/work/dy/project/2020CircRNA/temp_github/data/circRNA_exp.rds");

all_CircDataMF_LME_exp=t(all_CircDataMF_LME[,as.character(PRJEB_Ann$SRR)]/max(all_CircDataMF_LME[,as.character(PRJEB_Ann$SRR)]))
all_CircDataMF_LME_exp=cbind(ID=PRJEB_Ann$SRR, as.data.frame(all_CircDataMF_LME_exp))

SampleList = list(list(preNonBenefitID,preBenefitID),list(postNonBenefitID,postBenefitID))

#==LME 
LME_result <-  tibble::tibble(Class=c("Pre_BenNonBen","Post_BenNonBen"),
                                  SampleList = SampleList)
times=1
LME_result <- LME_result %>% dplyr::mutate(LME_diff_NoAdjust = purrr::map(.x=SampleList,function(.x){
		# define data format of input files		
		#.x=list(preNonBenefitID,preBenefitID);
		print(times)		
		times=times+1
		
		.x1=as.character(unlist(.x))
		temp_datRPM <- as_tibble(all_CircDataMF_LME_exp[.x1,])
		temp_datTraits <- right_join(PRJEB_Ann_LME,data.frame(.x1), by=c("ID" = ".x1"))
		label=c(rep('T2',length(.x[[1]])),rep('T1',length(.x[[2]])))
	
		DE_PD1s <- data.frame(circID=NA, Beta=NA, StdErr=NA, Pval=NA)
		numCirc <- ncol(temp_datRPM)-1
		circV <- colnames(temp_datRPM)

		for (i in 1:numCirc){ 
			circID <- circV[i+1]
			result <- try(runlme.My(circID,temp_datRPM,temp_datTraits,label), silent=T)
			if(length(result) == 6){
				Beta_Diagnosis <- result$Value
				StdErr_Diagnosis <- result$Std.Error
				Pval_Diagnosis <- result$`p-value`
			} else{
				#cat('Error in LME for circRNA', circID, '\n')
				#cat('Setting Beta value=0, StdErr=Inf, and Pval=1\n')
				Beta_Diagnosis <- 0
				StdErr_Diagnosis <- 0
				Pval_Diagnosis <- 1
				}
			DE_PD1s[i, "circID"] <- circID
			DE_PD1s[i, "Beta"] <- Beta_Diagnosis
			DE_PD1s[i, "StdErr"] <- StdErr_Diagnosis
			DE_PD1s[i, "Pval"] <- Pval_Diagnosis
    
		}
	DE_PD1s$BH_adjPval = p.adjust(DE_PD1s$Pval, method="BH", n=nrow(DE_PD1s))
    return(DE_PD1s)
}))

		
		
		