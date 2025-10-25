# ==== ABOUT ==== 
# R-code for reproducing the result of the manuscript titled
#
#  Gendered educational disparities and in-home smoking and use of electronic tobacco/nicotine 
#  devices among cohabiting couples: findings from a Hungarian cross-sectional survey
#
# Authors: <anonymized>
#
# License: CC-BY-NC-ND 4.0


# ==== ENVIRONMENT, PACKAGES, COMMON FOLDERS/FILES AND FUNCTIONS ====
rm(list = ls())   #clear environment
if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!") #check R environment

fn_input <- "./2023_dataset.RDS" #downloadable by researchers from KDKD

# packages
library("dplyr")
library("gtsummary")
library("parallel")
library("precrec")
library("pwr")        
library("DescTools")  
library("effectsize") 
library("coin")      
library("rstatix")    
library("Matrix")    
library("glmnet")    
library("carData")    
library("car")        
library("ggplot2")
library("ggpubr")
library("fastshap")   
library("margins")

#special auxiliary function
  #df = input data frame of eval tables
  p025_p975_df <- function(res_df)
  {
    
    #checking input
    if (!is(res_df, "data.frame"))  {stop("res_df is not a data.frame")}

    #set up result data frame
    sum_res <- data.frame(variable_name = character(),
                          p05 = numeric(),
                          p95 = numeric()
                          )
    
    for (col_name in colnames(res_df)) {
      if (is.numeric(res_df[[col_name]])) {
        p025 <- as.numeric(quantile(res_df[[col_name]], 0.025, na.rm = TRUE))
        p500 <- as.numeric(quantile(res_df[[col_name]], 0.50, na.rm = TRUE))
        p975 <- as.numeric(quantile(res_df[[col_name]], 0.975, na.rm = TRUE))
        
        sum_res <- rbind(sum_res, data.frame(variable_name = col_name, 
                                             p025 = p025,
                                             p500 = p500,
                                             p975 = p975))
      }
      
    }
   return(sum_res)
  } 
  
# pred wrapper for Shapley
  pred_wrapper <- function(object, newdata) {
    predict(object, newx=newdata, s = "lambda.min") }

# multiple prop test  
prop_test <- function(x, y) {
  table_data <- table(x, y)
  results <- lapply(1:ncol(table_data), function(col_idx) {
    group <- colnames(table_data)[col_idx]
    successes <- table_data[2, col_idx]  
    total <- colSums(table_data)[col_idx]
    test_result <- prop.test(successes, total)
    data.frame(
      Group = group,
      Proportion = test_result$estimate,
      Lower_CI = test_result$conf.int[1],
      Upper_CI = test_result$conf.int[2]
    )
  })
  do.call(rbind, results)
}

# jacc for allowing in-home smoking/vaping w/o both banned group
calculate_jaccard_ci <- function(vec1, vec2, n_bootstrap = 1000, conf_level = 0.95) {
  jaccard_binary <- function(v1, v2) {
    M11 <- sum(v1 & v2)
    M01 <- sum(!v1 & v2)
    M10 <- sum(v1 & !v2)
    return(M11 / (M11 + M01 + M10))
  }
  
  original_jaccard <- jaccard_binary(vec1, vec2)
  n <- length(vec1)
  
  bootstrap_jaccards <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(n, n, replace = TRUE)
    boot_vec1 <- vec1[idx]
    boot_vec2 <- vec2[idx]
    bootstrap_jaccards[i] <- jaccard_binary(boot_vec1, boot_vec2)
  }
  
  ci <- quantile(bootstrap_jaccards, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))
  
  return(list(
    original_jaccard = original_jaccard,
    ci_lower = ci[1],
    ci_upper = ci[2],
    bootstrap_samples = bootstrap_jaccards
  ))
}


# ==== LOAD INPUT  ====
# read RDS
  survey_raw <- readRDS(fn_input)

# convert to data frame
  df_survey <- as.data.frame(survey_raw)

# check the rows/columns
  if (dim(df_survey)[1] != 1500) stop("Input data is incorrect!")
  if (dim(df_survey)[2] != 195) stop("Input data is incorrect!")

# close connection and delete temporary vars
  rm(survey_raw, fn_input)
  
  
# ==== VALIDATE DATA  ====

# a relevant non-quota variable checked
# in 2021 a representative face-to-face study estimated proportion of Tobacco users 95% MT 33-35,2%
  at_prop_tob_users = as.numeric(df_survey$H3[as.numeric(df_survey$H3)<6]) < 4
  prop.test( x = sum(at_prop_tob_users), 
             n = length(at_prop_tob_users))
#35-40% is acceptable as more then a year passed
  prop.test( x = c(sum(at_prop_tob_users), 2387) ,  
             n = c(length(at_prop_tob_users), 7000 ),
             conf.level = 0.95)
  
  
# ==== FILTER PARTICIPANTS ====
#has a romantic partner and they are cohabiting with their romantic partner 
#and there are no more adults living in the household. 
  f1 <- (df_survey$par=="Igen, együtt élünk.") #respondent have a partner and they live together
  sum(f1, na.rm=TRUE)
  f2 <- (df_survey$SZD2==2) #only two adult persons live in the household
  sum(f2, na.rm=TRUE)

  filt<- f1&f2
  filt[is.na(filt)] <- FALSE 

  sum(filt) #number of selected couples

#var E19 is the cohabiting couples dataset
  E19 <- df_survey[filt,]

#handling critical non-responses
  E19$cri_fil <- (E19$HE1=="Nem tudom") | E19$HE1==("Nem akarok válaszolni") |
    (E19$HE2=="Nem tudom") | E19$HE2==("Nem akarok válaszolni") |
    (E19$H3=="Nem tudom") | E19$H3==("Nem akarok válaszolni") |
    (E19$HP3=="Nem tudom") | E19$HP3==("Nem akarok válaszolni")

  sum(E19$cri_fil) #31 person do not reply to critical question

#var F19 is the cohahibting couples who provided answers to all critical questions
  F19 <- subset(E19, cri_fil==FALSE)

#cleaning unneccassary variables and data frames
  rm(df_survey, E19)


# ==== GENERATE VARIABLES ====
#OUTPUT
#y1 is ALLOWING inhome smoking
  F19$y1 <- (F19$HE1 %in% c("Igen, bárhol.", "Igen, de nem bárhol.") )
#y2 is ALLOWING inhome vaping
  F19$y2 <- (F19$HE2 %in% c("Igen, bárhol.", "Igen, de nem bárhol.") )

#INPUT
#x1 = male age,
  
  #if the partner age is unkown, imputation of respondent age
  #assumed typo correction
  F19$P2[F19$P2 == 2013 & !is.na(F19$P2)] <- 2003
  
  sum(is.na(F19$P2))  #16 cases
  F19$P2[is.na(F19$P2)] = F19$szulev[is.na(F19$P2)] 
  
  F19$x1 <- ifelse(F19$SZD1 == "Férfi.", 
                   (2022 - F19$szulev), 
                   (2022 - F19$P2) )
#x2 = female age,
  F19$x2 <- ifelse(F19$SZD1 == "Nő.", 
                     (2022 - F19$szulev), 
                     (2022 - F19$P2) )

#x3 = difference between male minus female age in years
  F19$x3 <- F19$x1 - F19$x2

#x4 = presence of minors living the household
  sum(is.na(F19$SZD2g)) #missing
  F19$x4 <- (F19$SZD2g>0)
  F19$x4[is.na(F19$x4)] <- FALSE #imputation of median

#x5 = living in village
  summary(F19$SZD6a) # no missings
  F19$x5 <- (F19$SZD6a == "falu ")
  
#x6 = living in counties with less than 31K GDP per capita in USD in 2022
  summary(F19$SZD6b) # no missings
  counties_less_31kGDPpercap = c("Nógrád megye ", 
                                 "Szabolcs-Szatmár-Bereg megye ", 
                                 "Békés megye ",
                                 "Somogy megye ",
                                 "Jász-Nagykun-Szolnok megye ",
                                 "Borsod-Abaúj-Zemplén megye ",
                                 "Baranya megye ")
  F19$x6 <- (F19$SZD6b %in% counties_less_31kGDPpercap)
  
#FINANCIAL STATUS
#x7 = financial situation is OK or „just about OK”
  summary(F19$SZD9) # 3 missings, handled at x8
  
#x7 = financial situation is deprived
  F19$x7 <- ifelse(F19$SZD9 %in% c("hónapról-hónapra anyagi gondjai(k) van/vannak, vagy "
                                   ,"nélkülözések között él/élnek? "),1,0)

#EDUCATION
  summary(F19$SZD4) #no missing
  summary(F19$P3)   #no missing
  
#EDU
  tmp_male_col <- ifelse(F19$SZD1 == "Férfi.", 
                   (F19$SZD4 %in% c("doktori (PhD/DLA) vagy magasabb", "egyetemi diploma (MA/MSc)",
                                    "főiskolai diploma (BA/BSc)")), 
                   (F19$P3 %in% c("doktori (PhD/DLA) vagy magasabb", "egyetemi diploma (MA/MSc)",
                                  "főiskolai diploma (BA/BSc)")))
  

  tmp_female_col <- ifelse(F19$SZD1 == "Nő.", 
                   (F19$SZD4 %in% c("doktori (PhD/DLA) vagy magasabb", "egyetemi diploma (MA/MSc)",
                                    "főiskolai diploma (BA/BSc)")), 
                   (F19$P3 %in% c("doktori (PhD/DLA) vagy magasabb", "egyetemi diploma (MA/MSc)",
                                  "főiskolai diploma (BA/BSc)")))
  
  tmp_male_high <- ifelse(F19$SZD1 == "Férfi.", 
                    (F19$SZD4 %in% c("érettségi")), 
                    (F19$P3 %in% c("érettségi")))
  
  tmp_female_high <- ifelse(F19$SZD1 == "Nő.", 
                    (F19$SZD4 %in% c("érettségi")), 
                    (F19$P3 %in% c("érettségi")))
  
#x8 = dual high school or higher 
  F19$x8 <- (tmp_male_col==TRUE  |  tmp_male_high==TRUE) & (tmp_female_col==TRUE  | tmp_female_high==TRUE)

#x9 = just female is a high/shool college graduate, male is lower  
  F19$x9 <- (tmp_female_col==TRUE  |  tmp_female_high==TRUE) & (tmp_male_col==FALSE & tmp_male_high==FALSE)
  
#x10 = just male high school or higher       
  F19$x10 <- (tmp_male_col==TRUE  |  tmp_male_high==TRUE) & (tmp_female_col==FALSE & tmp_female_high==FALSE)

#x11 = male is lower high scool, female is lower high-school (REF)
  F19$x11 <- (tmp_female_col==FALSE  &  tmp_female_high==FALSE) & (tmp_male_col==FALSE & tmp_male_high==FALSE)

#SMOKING
  tmp_male_smk <- ifelse(F19$SZD1 == "Férfi.", 
                       (as.numeric(as.numeric(F19$H3)<3)), 
                       (as.numeric(as.numeric(F19$HP3)<3)) )
  
  tmp_female_smk <- ifelse(F19$SZD1 == "Nő.", 
         (as.numeric(as.numeric(F19$H3)<3)), 
         (as.numeric(as.numeric(F19$HP3)<3))  )
  
#x13 = dual smoker couple: male smoking and female smoking
  F19$x12 = (tmp_male_smk==TRUE) & (tmp_female_smk==TRUE)
  
#x13 = just female smoking
  F19$x13 <- (tmp_male_smk==FALSE) & (tmp_female_smk==TRUE)

#x14 = just male smoking
  F19$x14 <- (tmp_male_smk==TRUE) & (tmp_female_smk==FALSE)

#x15 = no-smoker couple (REF)
  F19$x15 <- (tmp_male_smk==FALSE) & (tmp_female_smk==FALSE)
  

#VAPING
  tmp_male_vap <- ifelse(F19$SZD1 == "Férfi.", 
                    (as.numeric(as.numeric(F19$H3)==2)) + (as.numeric(as.numeric(F19$H3)==3)), 
                    (as.numeric(as.numeric(F19$HP3)==2)) + (as.numeric(as.numeric(F19$HP3)==3))  )
  
  tmp_female_vap <- ifelse(F19$SZD1 == "Nő.", 
                    (as.numeric(as.numeric(F19$H3)==2)) + (as.numeric(as.numeric(F19$H3)==3)), 
                    (as.numeric(as.numeric(F19$HP3)==2)) + (as.numeric(as.numeric(F19$HP3)==3))  )
  
#x17 = dual, male vaping and female vaping
  F19$x16 = (tmp_male_vap==TRUE) & (tmp_female_vap==TRUE)
  
#x17 = just female vaping
  F19$x17 <- (tmp_male_vap==FALSE) & (tmp_female_vap==TRUE)
  
#x18 = just male vaping
  F19$x18 <- (tmp_male_vap==TRUE) & (tmp_female_vap==FALSE)

#x19 = non-vaper couple
  F19$x19 <- (tmp_male_vap==FALSE) & (tmp_female_vap==FALSE)       


#G19 contains only the variables selected to analyze
  at_var_names = c("y1","y2",
                   "x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                   "x11","x12","x13","x14","x15","x16","x17","x18","x19")
  G19 <- subset(F19, select = at_var_names)
  
#remove unusued vars
rm(F19)


# ==== DESCRIPTIVE - TABLE1 ====
# variable Labels
  at_var_labels <- list(
    y1  ~ "Allowing in-home smoking",
    y2  ~ "Allowing in-home vaping",
    x1  ~ "Male age (in years)",
    x2  ~ "Female age (in years)",
    x3  ~ "Age difference (male minus female) (in years)",
    x4  ~ "Presence of minors",
    x5  ~ "Living in village",
    x6  ~ "Living in counties with less than 31K GDP per capita",
    x7  ~ "Financial situation is deprived",
    x8 ~ "Dual high/college graduate couple",
    x9 ~ "Just female high/college graduate",
    x10 ~ "Just male high/college graduate",
    x11 ~ "None high/college graduates",
    x12 ~ "Dual smoker couple",
    x13 ~ "Just female smoking",
    x14 ~ "Just male smoking",
    x15 ~ "No-smoker couple",
    x16 ~ "Dual vaper couple",
    x17 ~ "Just female vaping",
    x18 ~ "Just male vaping",
    x19 ~ "Non-vaper couple"
  )
  
  at_var_types <- list(
    y1  ~ "dichotomous",
    y2  ~ "dichotomous",
    x1  ~ "continuous",
    x2  ~ "continuous",
    x3  ~ "continuous",
    x4  ~ "dichotomous",
    x5  ~ "dichotomous",
    x6  ~ "dichotomous",
    x7  ~ "dichotomous",
    x8  ~ "dichotomous",
    x9  ~ "dichotomous",
    x10 ~ "dichotomous",
    x11 ~ "dichotomous",
    x12 ~ "dichotomous",
    x13 ~ "dichotomous",
    x14 ~ "dichotomous",
    x15 ~ "dichotomous",
    x16 ~ "dichotomous",
    x17 ~ "dichotomous",
    x18 ~ "dichotomous",
    x19 ~ "dichotomous"
    )

#how to report
  at_report <- list(all_continuous() ~ "{median} ({p25} - {p75})", 
                    all_dichotomous() ~ "{n} ({p}%)")


#full table
  tbl_full <- gtsummary::tbl_summary(G19, 
                                     type  = at_var_types,
                                     label = at_var_labels,
                                     statistic = at_report,
                                     digits = list(all_categorical() ~ c(0, 1)))
  
  
  table1 = tbl_full
  table1


# ==== BIVARIATE ANALYSES - TABLE 2, FIGURE 1 and FIGURE 2 ====
#prevalence of smoking ban and vaping ban and difference:
  p1 <- table(G19$y1)[2]/sum(table(G19$y1))
  p2 <- table(G19$y2)[2]/sum(table(G19$y2))
  
  prop.test(x = c(table(G19$y1)[2]), n = sum(table(G19$y1)))
  prop.test(x = c(table(G19$y2)[2]), n = sum(table(G19$y2)))

  prop.test(x = c(table(G19$y1)[2], table(G19$y2)[2]), n = c(sum(table(G19$y1)), sum(table(G19$y2))))
  mcnemar.test(table(G19$y1, G19$y2)) #since paired data

  h1<-pwr::ES.h(p1,p2)
  h1

  pwr.p.test(h=h1,n=sum(table(G19$y1)),sig.level=0.05,alternative="two.sided")

#prevalence of smoking ban and vaping ban and difference -- Table 2
  tbl2_n = table(G19$y1, G19$y2)
  tbl2_p = prop.table(tbl2_n) 
  table2 = tbl2_p
  chisq.test(tbl2_n)
  effectsize::oddsratio(tbl2_n, ci = 0.95)
  SomersDelta(tbl2_n, direction = "column", conf.level = 0.95)
  
  #jaccard
  tbl2_n <- matrix(c(444, 57, 23, 99), nrow = 2, byrow = TRUE)
  calculate_jaccard_ci(G19$y1, G19$y2)
  
#Figure 1
  Fig12DF = subset(G19, select= c("y1", "y2", "x8", "x9", "x10", "x11"))
  Fig12DF$label <- case_when(
      Fig12DF$x8 == 1 ~ "Dual secondary or higher",
      Fig12DF$x9 == 1 ~ "Female-only secondary or higher",
      Fig12DF$x10 == 1 ~ "Male-only secondary or higher",
      Fig12DF$x11 == 1 ~ "Dual lower than secondary"
    )
  
  Fig12DF$label <- factor(Fig12DF$label, levels = c("Dual secondary or higher", 
                                                "Female-only secondary or higher", 
                                                "Male-only secondary or higher", 
                                                "Dual lower than secondary"))
  
  Fig1_comparisons <- list( c("Dual secondary or higher", "Male-only secondary or higher"), 
                          c("Dual secondary or higher", "Dual lower than secondary"), 
                          c("Female-only secondary or higher", "Male-only secondary or higher"), 
                          c("Female-only secondary or higher", "Dual lower than secondary") )
  
  Fig2_comparisons <- list( c("Dual secondary or higher", "Male-only secondary or higher"), 
                            c("Dual secondary or higher", "Dual lower than secondary") )
  
  Fig12DF$y1n <- as.numeric(Fig12DF$y1)
  Fig12DF$y2n <- as.numeric(Fig12DF$y2)
  
  #raw tables for Fig 1, 2
  table(Fig12DF$label , Fig12DF$y1n)
  table(Fig12DF$label , Fig12DF$y2n)
  
  Fig1_mcp <- pairwise.prop.test(x = as.matrix( table(Fig12DF$label , Fig12DF$y1n)) ,
                     p.adjust.method = "BH")
  Fig2_mcp <- pairwise.prop.test(x = as.matrix( table(Fig12DF$label , Fig12DF$y2n)) ,
                                 p.adjust.method = "BH")
  
  Fig1_padj_vals <- c(Fig1_mcp$p.value[2,1], Fig1_mcp$p.value[3,1], Fig1_mcp$p.value[2,2], Fig1_mcp$p.value[3,2])
  Fig2_padj_vals <- c(Fig2_mcp$p.value[2,1], Fig2_mcp$p.value[3,1], Fig2_mcp$p.value[2,2], Fig2_mcp$p.value[3,2])
  
  Fig1_padj_vals
  Fig2_padj_vals


jpeg("Figure1.jpeg",
       width = 12,      # Width in inches
       height = 12,     # Height in inches
       units = "in",    # Units (inches)
       res = 300,       # Resolution in DPI ( 
       quality = 100)   # Highest quality (0-100)
  
  p <- ggerrorplot(Fig12DF, x = "label", y = "y1n", add = "mean", error.plot = "errorbar",
                   order = c("Dual secondary or higher", 
                             "Female-only secondary or higher", 
                             "Male-only secondary or higher", 
                             "Dual lower than secondary"),
                   palette = "jco") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    labs(x = "Educational composition", y = "Ratio of in-home smoking") +
    stat_compare_means(comparisons = Fig1_comparisons,
                       method = "t.test" , #z-test and t-test does not differ after n>500
                       label.y = c(0.7, 0.6, 0.5, 0.4),
                       label = "p.signif",
                       aes(label = paste0("p=", Fig1_padj_vals)) ) #it does not work :-()
  #manually write padj_vals to the figure
  print(p)
dev.off()

jpeg("Figure2.jpeg",
     width = 12,      # Width in inches
     height = 12,     # Height in inches
     units = "in",    # Units (inches)
     res = 300,       # Resolution in DPI ( 
     quality = 100)   # Highest quality (0-100)

p <- ggerrorplot(Fig12DF, x = "label", y = "y2n", add = "mean", error.plot = "errorbar",
                 order = c("Dual secondary or higher", 
                           "Female-only secondary or higher", 
                           "Male-only secondary or higher", 
                           "Dual lower than secondary"),
                 palette = "jco") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Educational composition", y = "Ratio of in-home EDTNAD-use") +
  stat_compare_means(comparisons = Fig2_comparisons,
                     method = "t.test" , #z-test and t-test does not differ after n>500
                     label.y = c(0.7, 0.6, 0.5, 0.4),
                     label = "p.signif",
                     aes(label = paste0("p=", Fig2_padj_vals)) ) #it does not work :-()
#manually write padj_vals to the figure
print(p)
dev.off()

  Fig1_padj_vals
  Fig2_padj_vals
#note: have to adjust asterisks to the adjusted p-values at Fig 2 manually


# ==== MULTIVARIATE ANALYSES ====
#removing reference categories
at_xvar_labels <- list(
#  y1  ~ "Allowing in-home smoking",
#  y2  ~ "Allowing in-home vaping",
  x1  ~ "Male age (in years)",
  x2  ~ "Female age (in years)",
  x3  ~ "Age difference (male minus female) (in years)",
  x4  ~ "Presence of minors",
  x5  ~ "Living in village",
  x6  ~ "Living in counties with less than 31K GDP per capita",
  x7  ~ "Financial situation is deprived",
  x8 ~ "Dual high/college graduate couple",
  x9 ~ "Just female high/college graduate",
  x10 ~ "Just male high/college graduate",
#  x11 ~ "None high/college graduates",
  x12 ~ "Dual smoker couple",
  x13 ~ "Just female smoking",
  x14 ~ "Just male smoking",
#  x15 ~ "No-smoker couple",
  x16 ~ "Dual vaper couple",
  x17 ~ "Just female vaping",
  x18 ~ "Just male vaping"
#  x19 ~ "Non-vaper couple"
)

  at_xvar_names <- sapply(at_xvar_labels, function(x) as.character(x)[2])
  
#Variable names
#output variables
  y1var <- c("y1") #allowing smoking at home
  y2var <- c("y2") #allowing vaping at home

#input variables
  xvarsname <- gsub('[ -]', "_", gsub('["\'~,]', "", at_xvar_labels))
  xvars <- at_xvar_names
  n_vars = length(xvars)

#data preparation for analysis, in matrix format
  x1s <-  scale(as.matrix(G19[,xvars]))
  y1  <-  as.matrix(G19[,y1var ])
  y2  <-  as.matrix(G19[,y2var]) 

#LASSO 1st run, saved for plotting later
  sm <- glmnet(x1s, y1, 
               alpha = 1, 
               family = "binomial") #smoking lasso modell
  vm <- glmnet(x1s, y2, 
               alpha = 1, 
               family = "binomial") #vaping lasso modell

#LASSO 100 times
  k = 100
  shap_k = 10 #defalult value for fastshap::explain()
  
#set up results data frame
  #smoking ban
  smk_col_names = c("smk_ln_lambda", 
                    "smk_dev",
                    "smk_ROC_AUCf",
                    "smk_PR_AUCf",
                    "smk_intercept",
                    xvarsname)
  res_smk <- data.frame(matrix(ncol = length(smk_col_names), nrow = k))
  colnames(res_smk) = smk_col_names
  
  #vaping ban
  vap_col_names = c("vap_ln_lambda", 
                    "vap_dev",
                    "vap_ROC_AUCf",
                    "vap_PR_AUCf",
                    "vap_intercept",
                    xvarsname)
  res_vap <- data.frame(matrix(ncol = length(vap_col_names), nrow = k))
  colnames(res_vap) = vap_col_names
  
  #smoking shapley results (effect sizes) medians
  smk_shap_meds <- NULL
  vap_shap_meds <- NULL
  
#lasso cycle 
  for (c1 in 1:k) {
    
    #parallelize cv.glmnet for smoking and vaping models
    results <- mclapply(1:2, function(i) {
      if (i == 1) {
        cv_sm <- cv.glmnet( x1s, y1, 
                           alpha = 1, 
                           family = "binomial", 
                           nfolds = 4,    
                           type.measure = "deviance",
                           relax = FALSE)
        sm_lam <- cv_sm$lambda.min  
    
        sm1 <- glmnet(x1s, y1, 
                      alpha  = 1, 
                      family = "binomial", 
                      lambda = sm_lam)
        
        smk_shap <- fastshap::explain(
          object = sm1,
          X = x1s,          
          nsim = shap_k,       
          pred_wrapper = pred_wrapper)
        
        return(list(cv = cv_sm, model = sm1, shap = smk_shap, lam=sm_lam))
      } else {
        cv_vm <- cv.glmnet(x1s, y2, 
                           alpha  = 1, 
                           family = "binomial", 
                           nfolds = 4,    
                           type.measure = "deviance",
                           relax = FALSE)
        vm_lam <- cv_vm$lambda.min  
        
        vm1 <- glmnet(x1s, y2, 
                      alpha  = 1, 
                      family = "binomial", 
                      lambda = vm_lam)
        
        vap_shap <- fastshap::explain(
          object = vm1,
          X = x1s,          
          nsim = shap_k,       
          pred_wrapper = pred_wrapper)
        
        return(list(cv = cv_vm, model = vm1, shap=vap_shap, lam=vm_lam))
      }
    }, mc.cores = 2) # Use 2 cores
    
    # extract results
    cv_sm <- results[[1]]$cv
    sm1   <- results[[1]]$model
    smk_shap <- results[[1]]$shap
    sm_lam <- results[[1]]$lam
    
    cv_vm <- results[[2]]$cv
    vm1   <- results[[2]]$model
    vap_shap <- results[[2]]$shap
    vm_lam <- results[[2]]$lam
    
    # store results for smoking
    res_smk[c1, "smk_ln_lambda"] <- log(sm_lam)

    #eval done by precrec, but the focus is those, who allow!
    tmp_AUCS <- evalmod(
      scores = (predict(sm1, newx = x1s, s = sm_lam, type = "response")), 
      labels = y1 
      )
    
    res_smk[c1, "smk_dev"]      <- cv_sm$cvm[cv_sm$lambda == sm_lam]
    res_smk[c1, "smk_ROC_AUCf"] <- precrec::auc(tmp_AUCS)$aucs[1]
    res_smk[c1, "smk_PR_AUCf"]  <- precrec::auc(tmp_AUCS)$aucs[2]
    
    tmp_coefs_sm <- as.data.frame(t(matrix(coef(sm1))))
    colnames(tmp_coefs_sm) <- c("smk_intercept", xvarsname)
    res_smk[c1, 5:(5+n_vars)] <- tmp_coefs_sm
    
    # store results for vaping
    res_vap[c1, "vap_ln_lambda"] <- log(vm_lam)

    tmp_AUCS <- evalmod(
      scores = (predict(vm1, newx = x1s, s = vm_lam, type = "response")), 
      labels = y2 
    )
    res_vap[c1, "vap_dev"]  <- cv_vm$cvm[cv_vm$lambda == vm_lam]
    res_vap[c1, "vap_ROC_AUCf"] <- precrec::auc(tmp_AUCS)$aucs[1]
    res_vap[c1, "vap_PR_AUCf"]  <- precrec::auc(tmp_AUCS)$aucs[2]
    
    tmp_coefs_vm <- as.data.frame(t(matrix(coef(vm1))))
    colnames(tmp_coefs_vm) <- c("vap_intercept", xvarsname)
    res_vap[c1, 5:(5+n_vars)] <- tmp_coefs_vm
    
    #shapley
    smk_shap_median = apply(abs(smk_shap), 2, median)
    vap_shap_median = apply(abs(vap_shap), 2, median)
    
    smk_shap_meds = rbind(smk_shap_meds, smk_shap_median)
    vap_shap_meds = rbind(vap_shap_meds, vap_shap_median)  
    
    cat("process\t = " , round( ((c1/k)*100) ) , "%\n")
         
    Sys.sleep(0.01)
} #end of the lasso cycle, results might slightly differ due to stohasticity
  
#table 3
  table3a = p025_p975_df(res_smk)
  table3a$smk = paste0(round(table3a$p500,2), " [", round(table3a$p025,2), "–", round(table3a$p975,2), "]")
  table3a$smk_ast = ifelse((table3a$p025<0.0 & table3a$p975<0.0) | (table3a$p025>0.0 & table3a$p975>0.0), 
                           paste0(table3a$smk,"*"), table3a$smk)
  table3a$smk_shap_mm = round(c(NA,NA,NA,NA,NA, apply(abs(smk_shap_meds), 2, median) ),3)

  table3b = p025_p975_df(res_vap)
  table3b$vap = paste0(round(table3b$p500,2), " [", round(table3b$p025,2), "–", round(table3b$p975,2), "]")
  table3b$vap_ast = ifelse((table3b$p025<0.0 & table3b$p975<0.0) | (table3b$p025>0.0 & table3b$p975>0.0), 
                           paste0(table3b$vap,"*"), table3b$vap)
  table3b$vap_shap_mm = round(c(NA,NA,NA,NA,NA,apply(abs(vap_shap_meds), 2, median) ),3)
    
  table3    = cbind(table3a[,c(1,6,7)], table3b[,6:7])
