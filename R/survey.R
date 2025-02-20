#pv_survey_oo <- read_csv("~/Policy/CAMG/SolarPVReport/PVBESS_calibrater/pv_survey_oo.csv")
#pv_survey <- read_csv("~/Policy/CAMG/SolarPVReport/PVBESS_calibrater/pv_survey.csv")
#bill_values <- read_csv("~/Policy/CAMG/SolarPVReport/PVBESS_calibrater/bills.csv")
#pv_qanda <- read_csv("~/Policy/CAMG/SolarPVReport/PVBESS_calibrater/pv_qanda.csv")

#remove household profile
#pv_survey1 <- pv_survey1 %>% dplyr::select(-qi)

#dim(pv_survey_oo)


#' recode_bill
#'
#' converts to annual bill feature (q_ab), dropping min and max bills (q14 & q15)
#'
#' @param pv_data_in input survey data
#'
#' @returns a new survey dataframe with annual bills range 1~ < 900, 2~900-1350, 3~1350 to 1650, 4~1650-2100,5~> 2100
#' @export
#'
#' @examples
recode_bill <- function(pv_data_in){

  pv_data_in <- pv_data_in %>% dplyr::inner_join(bill_values %>% dplyr::rename("q14"=response_code,"highest_bill"=bill))
  pv_data_in <- pv_data_in %>% dplyr::inner_join(bill_values %>% dplyr::rename("q15"=response_code,"lowest_bill"=bill))
  pv_data_in <- pv_data_in %>% dplyr::mutate(lowest_bill= suppressWarnings(as.numeric(lowest_bill)),highest_bill=suppressWarnings(as.numeric(highest_bill)), annual_bill=3*(lowest_bill+highest_bill))
  pv_data_in <- pv_data_in %>% dplyr::mutate(q_ab = dplyr::case_when(annual_bill < 900~1,annual_bill>=900 & annual_bill < 1350~2,
                                              annual_bill>=1350 & annual_bill < 1650~3,
                                             annual_bill>=1650 & annual_bill < 2100~4,
                                             annual_bill > 2100~5))
  #pv_data_in <- pv_data_in %>% dplyr::mutate(q_ab = dplyr::case_when(annual_bill < 1000~1,annual_bill>=1000 & annual_bill < 2000~2,
   #                                          annual_bill > 2000~3))
  pv_data_in <- pv_data_in %>% dplyr::mutate(q_ab = tidyr::replace_na(q_ab,6))
  pv_data_in <- pv_data_in %>% dplyr::select(-q14,-q15,-lowest_bill,-highest_bill,-annual_bill)
  return(pv_data_in)
}

#' feature_select
#'
#' removes selected features from pv survey dataset before calibration
#' if recode_bills=T an annual bill (q_ab) feature is introduced and min-max bills are dropped
#'
#' @param pv_data_in pv_survey_oo or pv_survey
#' @param recode_bills if true an annual bill feature q_ab is introduced with values
#' @param drop_lowestbills if true q15 is dropped
#'
#' @returns reduced survey dataframe
#' @export
#'
#' @examples
feature_select <- function(pv_data_in, recode_bills=T, drop_lowestbills=F){

  pv_data_out <- pv_data_in %>% dplyr::filter(q262 != 1) #remove current adopters
  pv_data_out <- pv_data_out %>% dplyr::select(-serial,-q261,-q262,-q263,-q264,-q265) #remove non-PV ZET ownership
  #remove social grade and tariff questions
  pv_data_out <- pv_data_out %>% dplyr::select(-qd,-q40,-q41,-q43,-q44) #remove social grade, tariff questions
  #remove additional enviro questions
  pv_data_out <- pv_data_out %>% dplyr::select(-q30_3,-q30_4,-q30_5,-q30_6)
  #recode bills
  if(recode_bills) pv_data_out <- recode_bill(pv_data_out)
  if(drop_lowestbills & !recode_bills)  pv_data_out <- pv_data_out %>% dplyr::select(-q15)
  pv_data_out %>% return()

}


#' transform_to_utils
#'
#' Transforms survey Likert scores for likelihood of pv adoption (1..5) to utilities in range -1,1.
#' The transformation depends on s (utility uncertainty) and epsilon (an initial hypothetical bias estimate).
#' The error in the hypithetical bias estimate is corrected at the micro-calibration stage
#'
#' @param pv_data_in survey data including qsp22_7 (Likert scores for adopting pv)
#' @param s utility uncertainty (default 0.15)
#' @param epsilon degree of hypothetical bias (default 0.7, 1 = no bias). Further hypothetical bias correction is applied at ABM tuning stage.
#' @param s utility uncertainty (standard deviation) default 0.15
#' @param epsilon survey hypothetical bias (zero bias corresponds to 1) default 0.7
#'
#' @return modified pv survey dataframe
#' @export
#'
#' @examples
transform_to_utils <-function(pv_data_in,s=0.15,epsilon=0.7){

  #pv_data1 <- pv_data %>% dplyr::select(-ID)
  util_mapping <- map_likertscores_to_utilities(s,epsilon)
  pv_data1 <- pv_data_in %>% dplyr::rowwise() %>% dplyr::mutate(u = util_mapping[[q46_5,"dU"]]) %>% dplyr::select(-q46_5)
  return(pv_data1)
}


#' find_optimum_rounds_from_crossvalidation
#'
#' helper function to find optimum learning complexity for for given learning rate and tree depth
#'
#' @param pv_data_in solar PV survey dataset e.g. pv_data_oo
#' @param learning_rate eta parameter (default 0.02)
#' @param tree_depth maximum tree depth (default 5)
#' @param k_crossvalidation k-fold cross validation

#'
#' @return n_opt
#' @export
#'
#' @examples
find_optimum_rounds_from_crossvalidation <- function(pv_data_in, learning_rate=0.02, tree_depth=5, k_crossvalidation=5){

  #pv_data1 <- pv_data %>% dplyr::select(-ID)
  #pv.train <- xgboost::xgb.DMatrix(as.matrix(pv_util[,-dim(pv_util)[2]]),label=as.vector(pv_util$u), missing=NA)
  pv_train <- xgboost::xgb.DMatrix(as.matrix(pv_data_in %>% dplyr::select(-u)),label=as.vector(pv_data_in$u), missing=NA)

  #if(!train_on_utilities) #train on Likert scores
  # pv.train <- xgboost::xgb.DMatrix(as.matrix(pv_data1[,-dim(pv_data1)[2]]),label=as.vector(pv_data1$qsp22_7-1), missing=NA)
  #
  paramlist <- list(booster="gbtree",
                    tree_method = "exact",
                    eta=learning_rate,
                    max_depth=tree_depth,
                    gamma=0,
                    subsample=0.9,
                    colsample_bytree = 0.9,
                    objective="reg:squarederror",
                    eval_metric="rmse"
                    #objective="multi:softprob",
                    #eval_metric = "mlogloss"
  )

  bst <- xgboost::xgb.cv(params=paramlist,pv_train,nrounds=500,nfold=k_crossvalidation)

  cv_data <- bst["evaluation_log"] %>% as.data.frame() %>% tibble::as_tibble()
  nopt <- cv_data[,"evaluation_log.test_rmse_mean"][[1]] %>% as.numeric() %>% which.min()
  print(paste("optimal nrounds",nopt))
  return(nopt)

}


#' get_boosted_tree_model
#'
#' creates a cross-validated boosted tree regression model from pv survey data
#'
#' @param pv_data_in input survey data
#' @param learning_rate eta. typical value 0.02
#' @param tree_depth tree depth typical value 5
#' @param k_crossvalidation k-fold cross validation typical value 5
#' @param complexity_factor "over-fitting" enhancement relative to optimal model complexity from cross-validation. Values in range 1-1.5.
#'
#'
#' @return xgboost model
#' @export
#'
#' @examples
#'
get_boosted_tree_model <- function(pv_data_in, learning_rate=0.02, tree_depth=5, k_crossvalidation=5,complexity_factor = 1){

  #if("ID" %in% names(pv_data_in)) pv_data_in <- pv_data_in %>% dplyr::select(-ID)
  #pv_util <- transform_to_utils(pv_data_in,s,epsilon)
  pv_train <- suppressWarnings(xgboost::xgb.DMatrix(as.matrix(pv_data_in %>% dplyr::select(-u)),label=as.vector(pv_data_in$u), missing=NA))

  #if(!train_on_utilities) #train on Likert scores
  # pv.train <- xgboost::xgb.DMatrix(as.matrix(pv_data1[,-dim(pv_data1)[2]]),label=as.vector(pv_data1$qsp22_7-1), missing=NA)

  paramlist <- list(booster="gbtree",
                    tree_method = "exact",
                    eta=learning_rate,
                    max_depth=tree_depth,
                    gamma=0,
                    subsample=0.9,
                    colsample_bytree = 0.9,
                    objective="reg:squarederror",
                    eval_metric="rmse"
                    #objective="multi:softprob",
                    #eval_metric = "mlogloss"
  )
  n_opt <- find_optimum_rounds_from_crossvalidation(pv_data_in,learning_rate,tree_depth,k_crossvalidation)
  bst <- xgboost::xgboost(data=pv_train,params=paramlist,pv_train,nrounds=complexity_factor*n_opt)
  return(bst)
}


#bst <- get_boosted_tree_model(transform_to_utils(feature_select(pv_survey_oo,recode_bills=F,drop_lowestbills = T),epsilon=0.7))

#' get_shap_scores
#'
#' @param pv_data_in the pv_survey dataset
#' @param bst xgboost model
#'
#' @returns dataframe long format
#' @export
#'
#' @examples
get_shap_scores <- function(pv_data_in,bst){

  pv_data_long <- pv_data_in
  pv_data_long$ID <- 1:dim(pv_data_in)[1]
  pv_data_long <- pv_data_long %>% tidyr::pivot_longer(-ID,names_to="question_code",values_to="response_code")
  pv_data_long <- pv_data_long %>% dplyr::left_join(pv_qanda,by=c("question_code","response_code"))

  shap_scores <- predict(bst, as.matrix(pv_data_in %>% dplyr::select(-u)), predcontrib = TRUE, approxcontrib = F) %>% tibble::as_tibble()
  shap_scores$ID <- 1:dim(shap_scores)[1]
  shap_scores_long <- tidyr::pivot_longer(shap_scores,-ID,values_to="shap","names_to"="question_code")
  #add predictions
  preds <- shap_scores_long %>% dplyr::group_by(ID) %>% dplyr::summarise(u_predicted=sum(shap)) #includes BIAS
  #preds$actual <- pv_data$qsp22_7
  shap_scores_long1 <- shap_scores_long  %>% dplyr::inner_join(preds,by="ID")
  shap_scores_long1$u_actual <- sapply(pv_data_in$u, rep, dim(pv_data_in)[2]) %>% as.vector()
  #shap_scores_long1$pred <- shap_scores_long1$pred + 1 #+ dplyr::filter(shap_scores,name=="BIAS")$value
  #pv_data1$ID <- 1:dim(pv_data1)[1]
  #shap_scores_long1 <- shap_scores_long1 %>% dplyr::inner_join(pv_data1)
  shap_scores_long1 <-  shap_scores_long1 %>% dplyr::left_join(pv_data_long)
  return(shap_scores_long1)
}

#shap_scores_long <- get_shap_scores(transform_to_utils(feature_select(pv_survey_oo,recode_bills=F,drop_lowestbills = T),epsilon=0.7),bst)


#' get_abm_calibration
#'
#' returns abm partial utilities for the selected features (q9_1 and qsp21) and barrier (theta) terms in ABM model for each agent. Results are expressed as mean
#' partial utilities corresponding to each survey response and individual weights for each agent. A regularisation parameter can be
#' used to ensure that model weights for financial and social variables are > 1, or some negative weights can be tolerated.
#'
#' @param shap_scores_long individual shap scores by feature (output from get_shap_scores)
#' @param stat statistic - median (default) or mean
#' @param regularisation regularisation 0=none, 1= full, > 1 over, < 1 under
#'
#'
#' @return data frame giving partial utilities for abstracted model features and residual (theta) terms and individual weights
#' @export
#'
#' @examples
get_abm_calibration <- function(shap_scores_long, stat="mean",regularisation=1){
  #
  shap_scores <- shap_scores_long %>% dplyr::select(-question,-response)
  #q45 social
  #q1,q2,q3,q5,q10,q14,q15 financial and roof constraint
  social_code <- "q45"
  #finance_codes <- c("q_ab")
  finance_codes <- "q14"
  u_theta <- shap_scores %>% dplyr::filter(!(question_code %in% c(finance_codes,social_code))) %>% dplyr::group_by(ID) %>% dplyr::summarise(shap=sum(shap))
  u_theta$question_code <- "theta"
  u_theta$response_code <- NA

  shap_scores_abm <- shap_scores %>% dplyr::filter(question_code %in% c(finance_codes,social_code)) %>% dplyr::select(-u_predicted,-u_actual)
  #shap_scores_abm <- shap_scores %>% dplyr::filter(question_code %in% c(finance_codes,social_code)) %>% dplyr::select(-u_actual,-shap)

  shap_scores_abm <- shap_scores_abm %>% dplyr::bind_rows(u_theta) %>% dplyr::arrange(ID)
  shap_scores_abm <- shap_scores_abm %>% dplyr::rename("du"=shap)
  #regularise
  #regularise
  min_shap <- shap_scores_abm %>% filter(question_code != "theta") %>% group_by(question_code) %>% slice_min(du) %>% select(question_code,du) %>% rename("du_min"=du)
  theta_shift <- min_shap %>% pull(du_min) %>% sum()
  min_shap <- min_shap %>% bind_rows(tibble(question_code="theta",du_min = -theta_shift))

  shap_scores_abm <- shap_scores_abm %>% inner_join(min_shap) %>% mutate(du=du-regularisation*du_min) %>% select(-du_min)



  if(stat=="mean") {shap_scores_mean <- shap_scores_abm %>% dplyr::group_by(question_code,response_code) %>% dplyr::summarise(du_average=mean(du))
  shap_scores_abm <- shap_scores_abm %>% dplyr::inner_join(shap_scores_mean)}
  if(stat=="median") {shap_scores_median <- shap_scores_abm %>% dplyr::group_by(question_code,response_code) %>% dplyr::summarise(du_average=median(du))
  shap_scores_abm <- shap_scores_abm %>% dplyr::inner_join(shap_scores_median)}

  shap_scores_abm <- shap_scores_abm %>% dplyr::mutate(weight=du/du_average)
  return(shap_scores_abm)
}

#get_abm_calibration(shap_scores_long)


#' get_model_weights
#'
#' The agent weights for financial, social and barrier terms
#'
#' @param shap_scores_long shaps scores (partial utilities)
#' @param stat median (default) or mean
#' @param regularisation 0 none, 1 full
#'
#' @return a dataframe with colums ID w_q9_1  w_qsp21 W_theta
#' @export
#'
#' @examples
get_model_weights <- function(shap_scores_long, stat="mean",regularisation = 1){

  shap_scores_abm <- get_abm_calibration(shap_scores_long,stat, regularisation)

  weights_abm <- shap_scores_abm %>% tidyr::pivot_wider(id_cols=c(-du,-response_code,-du_average),values_from=weight,names_from=question_code)
  names(weights_abm)[2:ncol(weights_abm)] <- paste("w_",names(weights_abm)[2:ncol(weights_abm)],sep="")
  return(weights_abm)
}


#' get_empirical_partial_utilities
#'
#' partial (dis) utilities for pv adoption derived from survey
#'
#'
#' @param shap_scores_long shap scores (partial utilities)
#' @param stat median (default) or mean
#' @param regularisation none 0, full 1, over > 1
#'
#' @return dataframe
#' @export
#'
#' @examples
get_empirical_partial_utilities <- function(shap_scores_long,stat="median", regularisation=1){

  shap_scores_abm <- get_abm_calibration(shap_scores_long,stat,regularisation=1)
  if(stat=="mean") partial_utils <- shap_scores_abm %>% dplyr::group_by(question_code,response_code) %>% dplyr::summarise(du_average=mean(du))
  if(stat=="median") partial_utils <- shap_scores_abm %>% dplyr::group_by(question_code,response_code) %>% dplyr::summarise(du_average=median(du))
  return(partial_utils)
}



#' pbeta_util
#'
#' Beta function probability generalised to range -1,1
#'
#' @param x real in range -1 to 1
#' @param shape1 "a" shape parameter
#' @param shape2 "b" shape parameter
#'
#' @return generalised beta function value
#' @export
#'
#' @examples
pbeta_util <- function(x,shape1,shape2){
  #beta function generalised to -2,1
  return(pbeta((x+1)/2,shape1,shape2))
}

#' dbeta_util
#'
#' Beta function distribution generalised to interval -1,1
#'
#' @param x in -1,1
#' @param shape1 "a" shape parameter
#' @param shape2 "b" shape parameter
#'
#' @return function value
#' @export
#'
#' @examples
dbeta_util <- function(x,shape1,shape2){
  #beta function generalised to -2,1
  return(stats::dbeta((x+1)/2,shape1,shape2))
}

#' probs_from_shape_params
#'
#' mapping between shape parameters of generalised Beta function and total probability that value > 0
#'
#' @param s standard deviation of Beta distribution
#'
#' @return dataframe
#' @export
#'
#' @examples
probs_from_shape_params <- function(s){

  df <- tidyr::tibble()
  for( a in seq(0.15,300, by=0.1)){
    f <- function(b)  {4*a*b-((a+b+1)*(a+b)^2)*s^2} #its a polynomial
    f1 <- function(b) {4*a*b-b^3*s^2-3*a*b^2*s^2-b^2*s^2-3*a^2*b*s^2-2*a*b*s^2-a^3*s^2-a^2*s^2}
    #  #find roots of cubic polynomial in b given s (sd) and a
    f.roots <- polyroot(c(-a^3*s^2-a^2*s^2,-3*a^2*s^2-2*a*s^2+4*a,-3*a*s^2-s^2,-s^2)) %>% Re()
    #  #f.roots
    i.roots <- polyroot(c(-a^3*s^2-a^2*s^2,-3*a^2*s^2-2*a*s^2+4*a,-3*a*s^2-s^2,-s^2)) %>% Im() %>% round(5)
    #
    #  #real roots only
    b <- f.roots[which(i.roots==0)]
    #positive real roots only
    b <- b[b>0]
    for(b1 in b)
      #  #print(paste("b=",b1,"mean=", 2*a/(a+b1)-1, "  sd=",2*sqrt(a*b1/((a+b1+1)*(a+b1)^2)),"prob=",round(1-pbeta_util(0,a,b1),2) ))
      df <- dplyr::bind_rows(df,tidyr::tibble(s=s,a=a,b=b1,mean= 2*a/(a+b1)-1,prob=round(1-pbeta_util(0,a,b1),3) ))
  }
  return(df)
}


#' shape_params_from_prob
#'
#' returns shape parameters corresponding to an adoption probability (i.e. area of generalised beta distribution > 0)
#'
#' @param df.in dataframe produced by probs_from_shape_params()
#' @param prob value in range 0,1
#'
#' @return shape paramaters (a and b)
#' @export
#'
#' @examples
shape_params_from_prob <- function(df.in,prob){
  #
  return(c(approx(df.in$prob,df.in$a,prob,ties="mean")$y,approx(df.in$prob,df.in$b,prob,ties="mean")$y))
}

#' map_likertscores_to_utilities
#'
#' maps likert scores 1..5 to adoption probabilities and then to expected utility.
#' Epsilon controls the degree and sign of the hypothetical bias correction needed when tuning the ABM
#'
#' @param s agent utility uncertainty (standard deviation of utility at time of survey)
#' @param epsilon probability scale factor (to allow for hypothetical bias) 1= no hypothetical bias. default 0.75.
#'
#' @return a dataframe of shape parameters and expected utility values corresponging to likert probabilities
#' @export
#'
#' @examples
map_likertscores_to_utilities <- function(s=0.15,epsilon=0.75){

  probs <- epsilon*c(0.1,0.3,0.5,0.7,0.9) #hard-wired
  df <- probs_from_shape_params(s)
  distrib_params <- tidyr::tibble()
  for(prob in probs)
    distrib_params <- distrib_params %>% dplyr::bind_rows(tidyr::tibble(s=s,epsilon=epsilon,prob=prob,a=shape_params_from_prob(df,prob)[1],b=shape_params_from_prob(df,prob)[2]))
  distrib_params <- distrib_params %>% dplyr::mutate(dU = 2*a/(a+b)-1)
  distrib_params$likert <- 1:5
  return(distrib_params)
}


