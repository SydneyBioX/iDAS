
#' Interpretable differential abundance analysis (three-way analysis)
#'
#' @param Z A matrix/dataframe of omics or gene expression data, row as sample.
#' @param f1 A vector of factor 1 variables.
#' @param f2 A vector of factor 2 variables.
#' @param f3 A vector of factor 3 variables.
#' @param random A vector of random effect term of ANOVA analysis,
#' by default is NULL, which means the model doesn't include random effect term.
#' @param test_func Testing function used, either stats::lm or lme4::lmer. By default is "lm".
#' @param Sig_cutoff  No effect test significance level is defined by a fraction value
#' to indicate when ordering the p-values
#' and defining the top X% as the significance level. If both Sig_cutoff and Sig
#' are set, the algorithm will, by default, use Sig_cutoff to determine
#' the significance level.
#' @param Sig No effect test significance level is defined directly
#' by a fraction value, by default is 0.05
#' @param Int Interaction effect test significance level,by default is 0.01
#' @param F1 Factor 1 effect test significance level, by default is 0.01
#' @param F2 Factor 2 effect test significance level, by default is 0.01
#' @param F3 Factor 3 effect test significance level, by default is 0.01
#' @param F1F2 F1 and F2 interaction effect test significance level, by default is 0.01
#' @param F1F3 F1 and F3 interaction effect test significance level, by default is 0.01
#' @param F2F3 F2 and F3 interaction effect test significance level, by default is 0.01
#' @param threeways three-way effect test significance level, by default is 0.01
#' @param adj_method Pvalue adjust method. See p.adjust. By default is "BH".
#' @param f1name The column name of factor 1, by default is F1.
#' @param f2name The column name of factor 2, by default is F2.
#' @param f3name The column name of factor 3, by default is F3.
#' @param randomname The column name of random effect term, by default is Random.
#'
#' @return A list of hypothesis test outcome, P_mat is the pvalue matrix of all tests,
#' S_mat is the statistics matrix of all test, cls_df is the Classification data frame of all tests.
#' @export
#' @importFrom stats anova formula p.adjust
#' @importFrom lme4 lmer
#' @examples #res=iDAS_3F(Z = X,
#' #f1 = all.timepoint, f2 = all.pcellstats, f3 = all.pcelltype, random = all.pid, test_func = "lmer",
#' #Sig_cutoff = 0.02,Int = 0.01,F1 = 0.01,F2=0.02,F3=0.01,adj_method = "BH")

iDAS_3F=function(Z,f1,f2,f3,random=NULL,test_func="lm",
               Sig_cutoff=0.02, Sig = 0.05,Int= 0.01, F1 = 0.01,F2 = 0.01,F3=0.01,
               F1F2=0.01,F1F3=0.01,F2F3=0.01,
               threeways=0.01,
               adj_method="BH",f1name=NULL,f2name=NULL,f3name=NULL,randomname=NULL){

  P_mat =S_mat <- matrix(NA, nrow = ncol(Z), ncol = 9)
  if(test_func=="lm"&is.null(random)){

    factor_tmp=check_factor_name3(f1name,f2name,f2name,randomname,f1,f2,f3,random)

    lm_full <- paste("Y~(", factor_tmp$f1_tmp, ")*(", factor_tmp$f2_tmp,
                     ")*(", factor_tmp$f3_tmp, ")",
                     sep = "")


    lm_int_alt2way <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp,  ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp,")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")",
                            sep = "")

    lm_int_altf1f2 <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")",
                            sep = "")

    lm_int_altf1f3 <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")",
                            sep = "")


    lm_int_altf2f3 <- paste("Y~(", factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")",
                            sep = "")

    lm_int_null <- paste("Y~(", factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp,
                         ")+(", factor_tmp$f3_tmp,
                         ")", sep = "")

    lm_f1 <- paste("Y~(", factor_tmp$f1_tmp, ")", sep = "")
    lm_f2 <- paste("Y~(", factor_tmp$f2_tmp, ")", sep = "")
    lm_f3 <- paste("Y~(", factor_tmp$f3_tmp, ")", sep = "")

  }else if(test_func=="lmer"&(!is.null(random))){

    factor_tmp=check_factor_name3(f1name,f2name,f2name,randomname,f1,f2,f3,random)

    lm_full <- paste("Y~(", factor_tmp$f1_tmp, ")*(", factor_tmp$f2_tmp,
                     ")*(", factor_tmp$f3_tmp,
                     ")+(1|",factor_tmp$random_tmp,")", sep = "")

    lm_int_alt2way <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp,  ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp,")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")+(1|",factor_tmp$random_tmp,")",
                            sep = "")

    lm_int_altf1f2 <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")+(1|",factor_tmp$random_tmp,")",sep = "")

    lm_int_altf1f3 <- paste("Y~(", factor_tmp$f2_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")+(1|",factor_tmp$random_tmp,")",sep = "")

    lm_int_altf2f3 <- paste("Y~(", factor_tmp$f1_tmp, "):(", factor_tmp$f3_tmp, ")+(",
                            factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp, ")+(",
                            factor_tmp$f3_tmp,")+(1|",factor_tmp$random_tmp,")",sep = "")

    lm_int_null <- paste("Y~(", factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp,
                         ")+(", factor_tmp$f3_tmp,
                         ")+(1|",factor_tmp$random_tmp,")",sep = "")

    lm_f1 <- paste("Y~(",  factor_tmp$f1_tmp, ")+(1|",
                   factor_tmp$random_tmp,")",sep = "")
    lm_f2 <- paste("Y~(", factor_tmp$f2_tmp, ")+(1|",
                   factor_tmp$random_tmp,")", sep = "")

    lm_f3 <- paste("Y~(", factor_tmp$f3_tmp, ")+(1|",
                   factor_tmp$random_tmp,")", sep = "")

  }else{
    print("wrong!")
  }

  # Sig vs Not sig
  p_sig =s_sig= c()
  sig_level=list()
  if(test_func=="lm"&is.null(random)){
    calc0 <- "Y~1"
    func=get("lm")
    calc1 <- lm_full
  }else{
    calc0 <- "Y~1+(1|random)"
    func=get("lmer")
    calc1 <- lm_full
  }
  for (i in 1:ncol(Z)) {
    if(test_func=="lm"&is.null(random)){
      dat <- data.frame(Y = Z[, i], f1, f2,f3)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      dat <- data.frame(Y = Z[, i], f1, f2,f3,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_sig[i] <-p
    s_sig[i]=S
  }
  if (is.null(adj_method)) {
    p_sig_adj <- p_sig
  }else {
    p_sig_adj <- p.adjust(p_sig, method = adj_method)
  }

  P_mat[, 1] <- p_sig_adj
  S_mat[,1]=s_sig
  if(!is.null(Sig_cutoff)){
    cutoff=sort(p_sig_adj)[length(p_sig_adj)*Sig_cutoff] # adjusted p value 1 problem
    cls_df <- data.frame(Sig0 = ifelse(p_sig_adj > cutoff,"non-sig","sig"))
    idx_sig <- which(cls_df$Sig0=="sig")
    if(length(idx_sig)>1){
      Z_sig <- Z[, idx_sig]
    }else{
      Z_sig <- data.frame(Z[, idx_sig] )
      colnames(Z_sig)=colnames(Z)[idx_sig]
      rownames(Z_sig)=rownames(Z)
    }
  }else if(Sig){
    cls_df <- data.frame(Sig0 = ifelse(p_sig_adj > Sig,"notSig","Sig"))
    idx_sig <- which(cls_df$Sig0=="Sig")
    if(length(idx_sig)>1){
      Z_sig <- Z[, idx_sig]
    }else{
      Z_sig <- data.frame(Z[, idx_sig] )
      colnames(Z_sig)=colnames(Z)[idx_sig]
      rownames(Z_sig)=rownames(Z)
    }

  }

  # Int vs not Int
  p_int_tmp=s_int_tmp   <- c()
  for (i in 1:ncol(Z_sig)) {

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_int_null
      calc1 <- lm_full
      dat <- data.frame(Y = Z_sig[, i], f1, f2,f3)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_int_null
      calc1 <- lm_full
      dat <- data.frame(Y = Z_sig[, i], f1, f2,f3,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }

    p_int_tmp[i] <- p
    s_int_tmp[i] <- S

  }

  if (is.null(adj_method)) {
    P_mat[idx_sig,2]=p_int_tmp
  }else{
    P_mat[idx_sig,2]=p.adjust(p_int_tmp,method =adj_method )
  }
  S_mat[idx_sig,2]=s_int_tmp

  cls_df[idx_sig,"Sig1"]= ifelse(P_mat[idx_sig,2] > Int,"notInt","Int")
  idx_int <- which(cls_df$Sig1=="Int")
  idx_notint = which(cls_df$Sig1=="notInt")
  if(length(idx_int)>1){
    Z_int <- Z[, idx_int]
  }else{
    Z_int <- data.frame(Z[, idx_int] )
    colnames(Z_int)=colnames(Z)[idx_int]
    rownames(Z_int)=rownames(Z)
  }
  if(length(idx_notint)>1){
    Z_notint = Z[, idx_notint]
  }else{
    Z_notint <- data.frame(Z[, idx_notint] )
    colnames(Z_notint)=colnames(Z)[idx_notint]
    rownames(Z_notint)=rownames(Z)
  }

  # f1,f2,f3

  p_f1_tmp=p_f2_tmp=p_f3_tmp=c()
  s_f1_tmp=s_f2_tmp=s_f3_tmp=c()
  for (i in 1:ncol(Z_notint)) {
    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_f1
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_f1
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }

    p_f1_tmp[i] <- p
    s_f1_tmp[i] <- S


    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_f2
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_f2
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f2_tmp[i] <- p
    s_f2_tmp[i] <- S

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_f3
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_f3
      calc1 <- lm_int_null
      dat <- data.frame(Y = Z_notint[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f3_tmp[i] <- p
    s_f3_tmp[i] <- S

  }


  if (is.null(adj_method)) {
    P_mat[idx_notint, 3] = p_f1_tmp
    P_mat[idx_notint, 4] = p_f2_tmp
    P_mat[idx_notint, 5] = p_f3_tmp
  }else {
    P_mat[idx_notint, 3] = p.adjust(p_f1_tmp,method = adj_method)
    P_mat[idx_notint, 4] = p.adjust(p_f2_tmp,method = adj_method)
    P_mat[idx_notint, 5] = p.adjust(p_f3_tmp,method = adj_method)
  }


  S_mat[idx_notint, 3] = s_f1_tmp
  S_mat[idx_notint, 4] = s_f2_tmp
  S_mat[idx_notint, 5] = s_f3_tmp

  cls_notint_tmp <- rep("Add", ncol(Z_notint))
  cls_notint_tmp[P_mat[idx_notint, 3] >F1 & P_mat[idx_notint, 4] <= F2 & P_mat[idx_notint, 5] <= F3]="F1"
  cls_notint_tmp[P_mat[idx_notint, 3] <=F1 & P_mat[idx_notint, 4] > F2 & P_mat[idx_notint, 5] <= F3]="F2"
  cls_notint_tmp[P_mat[idx_notint, 3] <=F1 & P_mat[idx_notint, 4] <= F2 & P_mat[idx_notint, 5] > F3]="F3"


  cls_df[idx_notint,"notInt"]= cls_notint_tmp


  p_2ways_tmp=s_2ways_tmp=c()
  for (i in 1:ncol(Z_int)) {
    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_int_alt2way
      calc1 <- lm_full
      dat <- data.frame(Y = Z_int[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_int_alt2way
      calc1 <- lm_full
      dat <- data.frame(Y = Z_int[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_2ways_tmp[i] <- p
    s_2ways_tmp[i] <- S
  }

  if (is.null(adj_method)) {
    P_mat[idx_int, 6] = p_2ways_tmp
  }else {
    P_mat[idx_int, 6] = p.adjust(p_2ways_tmp,method = adj_method)
  }

  S_mat[idx_int, 6] = s_2ways_tmp

  cls_df[idx_int,"Int"]=ifelse(P_mat[idx_int, 6] < threeways,"threeways","twoways")

  idx_threeways <- which(cls_df$Int=="threeways")
  idx_twoways = which(cls_df$Int=="twoways")
  if(length(idx_twoways)>1){
    Z_twoways = Z[, idx_twoways]
  }else{
    Z_twoways = data.frame(Z[, idx_twoways] )
    colnames(Z_twoways)=colnames(Z)[idx_twoways]
    rownames(Z_twoways)=rownames(Z)
  }

  if(length(idx_threeways)>1){
    Z_threeways = Z[, idx_threeways]
  }else{
    Z_threeways = data.frame(Z[, idx_threeways] )
    colnames(Z_threeways)=colnames(Z)[idx_threeways]
    rownames(Z_threeways)=rownames(Z)
  }


  p_f1f2_tmp=p_f1f3_tmp=p_f2f3_tmp=c()
  s_f1f2_tmp=s_f1f3_tmp=s_f2f3_tmp=c()
  for (i in 1:ncol(Z_twoways)) {
    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_int_altf1f2
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_int_altf1f2
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f1f2_tmp[i] <- p
    s_f1f2_tmp=S

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_int_altf2f3
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_int_altf2f3
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f2f3_tmp[i] <- p
    s_f2f3_tmp[i]=S

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_int_altf1f3
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_int_altf1f3
      calc1 <- lm_int_alt2way
      dat <- data.frame(Y = Z_twoways[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f1f3_tmp[i] <- p
    s_f1f3_tmp[i] <- S

  }

  if (is.null(adj_method)) {
    P_mat[idx_twoways, 7] = p_f1f2_tmp
    P_mat[idx_twoways, 8] = p_f2f3_tmp
    P_mat[idx_twoways, 9] = p_f1f3_tmp
  }else {
    P_mat[idx_twoways, 7] = p.adjust(p_f1f2_tmp,method = adj_method)
    P_mat[idx_twoways, 8] = p.adjust(p_f2f3_tmp,method = adj_method)
    P_mat[idx_twoways, 9] = p.adjust(p_f1f3_tmp,method = adj_method)
  }

  S_mat[idx_twoways, 7] = s_f1f2_tmp
  S_mat[idx_twoways, 8] = s_f2f3_tmp
  S_mat[idx_twoways, 9] = s_f1f3_tmp


  cls_twoways_tmp <- rep("twowaycombinations", ncol(Z_twoways))
  cls_twoways_tmp[P_mat[idx_twoways, 7] <=F1F2 & P_mat[idx_twoways, 8] > F2F3 & P_mat[idx_twoways, 9] > F1F3]="F1F2"
  cls_twoways_tmp[P_mat[idx_twoways, 7] > F1F2 & P_mat[idx_twoways, 8] <= F2F3 & P_mat[idx_twoways, 9] >  F1F3]="F2F3"
  cls_twoways_tmp[P_mat[idx_twoways, 7] > F1F2 & P_mat[idx_twoways, 8] > F2F3 & P_mat[idx_twoways, 9] <= F1F3]="F1F3"


  cls_df[idx_twoways,"twoways"]= cls_twoways_tmp

  cls_df$varname <- colnames(Z)
  colnames(P_mat) =colnames(S_mat) <- c("Sig0", "Intornotint", "F1",
                       "F2","F3","twowaysorthreeways","F1F2","F2F3","F1F3")
  rownames(P_mat)=  rownames(S_mat) =rownames(cls_df)  <- colnames(Z)
  return(list(P_mat=P_mat,S_mat=S_mat,cls_df=cls_df))
}

#' Title
#'
#' @param f1name a string of factor 1's name
#' @param f2name a string of factor 2's name
#' @param f3name a string of factor 3's name
#' @param randomname a string of random effect term's name
#' @param f1 A vector of factor 1 variables.
#' @param f2 A vector of factor 2 variables.
#' @param f3 A vector of factor 3 variables.
#' @param random A vector of random effect term variables.
#'
#' @return A list of each factor's name and values
#'
#' @examples  #factor_tmp=check_factor_name3(f1name,f2name,f2name,randomname,f1,f2,f3,random)

check_factor_name3=function(f1name,f2name,f3name,randomname,f1,f2,f3,random){
  if (is.null(f1name)) {
    f1name = "f1"
  }
  if (is.null(f2name)) {
    f2name = "f2"
  }
  if (is.null(f3name)) {
    f3name = "f3"
  }
  if (is.null(randomname)) {
    randomname = "random"
  }

  if (is.factor(f1)) {
    f1 <- data.frame(f1)
    colnames(f1) <- f1name
  }else if (is.vector(f1)) {
    f1 <- data.frame(f1)
    colnames(f1) <- f1name
  }
  f1name <- colnames(f1)
  if (is.factor(f2)) {
    f2 <- data.frame(f2)
    colnames(f2) <- f2name
  }else if (is.vector(f2)) {
    f2 <- data.frame(f2)
    colnames(f2) <- f2name
  }
  f2name <- colnames(f2)

  if (is.factor(f3)) {
    f3 <- data.frame(f3)
    colnames(f3) <- f3name
  }else if (is.vector(f3)) {
    f3 <- data.frame(f3)
    colnames(f3) <- f3name
  }
  f3name <- colnames(f3)


  if(is.null(random)){
    random=NULL
  }else{
    if (is.factor(random)) {
      random <- data.frame(random)
      colnames(random) <- randomname
    }else if (is.vector(random)) {
      random <- data.frame(random)
      colnames(random) <- randomname
    }
  }
  randomname=colnames(random)
  if (length(f1name) > 1) {
    f1_tmp <- paste(f1name, collapse = "+")
  }else {
    f1_tmp <- f1name
  }
  if (length(f2name) > 1) {
    f2_tmp <- paste(f2name, collapse = "+")
  }else {
    f2_tmp <- f2name
  }

  if (length(f3name) > 1) {
    f3_tmp <- paste(f3name, collapse = "+")
  }else {
    f3_tmp <- f3name
  }

  if (length(randomname) > 1) {
    random_tmp <- paste(randomname, collapse = "+")
  }else {
    random_tmp <- randomname
  }
  return(list(f1_tmp=f1_tmp,f2_tmp=f2_tmp,f3_tmp=f3_tmp,
              random_tmp=random_tmp,f1=f1,f2=f2,f3=f3,random=random))
}

