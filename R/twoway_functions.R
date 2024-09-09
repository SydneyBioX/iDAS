
#' Interpretable differential abundance analysis (two-way analysis)
#'
#' @param Z A matrix/dataframe of omics or gene expression data, row as sample.
#' @param f1 A vector of factor 1 variables.
#' @param f2 A vector of factor 2 variables.
#' @param random A vector of random effect term of ANOVA analysis,
#' by default is NULL, which means the model doesn't include random effect term.
#' @param test_func Testing function used, either stats::lm or lme4::lmer. By default is "lm".
#' @param Sig_cutoff No effect test significance level is defined by a fraction value
#' to indicate when ordering the p-values
#' and defining the top X% as the significance level. If both Sig_cutoff and Sig
#' are set, the algorithm will, by default, use Sig_cutoff to determine
#' the significance level.
#' @param Sig No effect test significance level is defined directly
#' by a fraction value, by default is 0.05
#' @param Int Interaction effect test significance level,by default is 0.01
#' @param F1 F1 effect test significance level, by default is 0.01
#' @param F2 F2 effect test significance level, by default is 0.01
#' @param adj_method Pvalue adjust method. See p.adjust. By default is "BH".
#' @param f1name The column name of factor 1, by default is F1.
#' @param f2name The column name of factor 2, by default is F2.
#' @param randomname The column name of random effect term, by default is Random.
#'
#' @return A list of hypothesis test outcome, P_mat is the pvalue matrix of all tests,
#' S_mat is the statistics matrix of all test, cls_df is the Classification data frame of all tests.
#' @export
#' @importFrom stats anova formula p.adjust
#' @importFrom lme4 lmer
#' @examples # res=iDAS_2F(Z= X,
#' #f1=pcelltype,f2=pcell_stats,random=NULL,test_func="lm",
#' #Sig_cutoff=0.02, Sig = 0.1,Int= 0.01, F1 = 0.01,F2 = 0.01,
#' #adj_method="BH",f1name=NULL,f2name=NULL,randomname=NULL)

iDAS_2F=function(Z,f1,f2,random=NULL,test_func="lm",
              Sig_cutoff=0.02, Sig = 0.05,Int= 0.01, F1 = 0.01,F2 = 0.01,
              adj_method="BH",f1name=NULL,f2name=NULL,randomname=NULL){

  P_mat =S_mat<- matrix(NA, nrow = ncol(Z), ncol = 4)
  if(test_func=="lm"&is.null(random)){

    factor_tmp=check_factor_name(f1name,f2name,randomname,f1,f2,random)
    lm_full <- paste("Y~(", factor_tmp$f1_tmp, ")*(", factor_tmp$f2_tmp, ")",
                     sep = "")
    lm_int_alt <- paste("Y~(", factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp, ")",
                        sep = "")
    lm_int_null <- paste("Y~(", factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp,
                         ")", sep = "")
    lm_f1 <- paste("Y~(", factor_tmp$f1_tmp, ")", sep = "")
    lm_f2 <- paste("Y~(", factor_tmp$f2_tmp, ")", sep = "")

  }else if(test_func=="lmer"&(!is.null(random))){

    factor_tmp=check_factor_name(f1name,f2name,randomname,f1,f2,random)
    lm_full <- paste("Y~(",  factor_tmp$f1_tmp, ")*(", factor_tmp$f2_tmp,
                     ")+(1|",factor_tmp$random_tmp,")", sep = "")
    lm_int_alt <- paste("Y~(",  factor_tmp$f1_tmp, "):(", factor_tmp$f2_tmp,
                        ")+(1|",factor_tmp$random_tmp,")", sep = "")
    lm_int_null <- paste("Y~(",  factor_tmp$f1_tmp, ")+(", factor_tmp$f2_tmp,
                         ")+(1|",factor_tmp$random_tmp,")",sep = "")
    lm_f1 <- paste("Y~(",  factor_tmp$f1_tmp, ")+(1|",
                   factor_tmp$random_tmp,")",sep = "")
    lm_f2 <- paste("Y~(", factor_tmp$f2_tmp, ")+(1|",
                   factor_tmp$random_tmp,")", sep = "")

  }else{
    print("wrong!")
  }

  lm_formula=c("full" = lm_full, "int_null" = lm_int_null,
               "int_alt" = lm_int_alt,
               "f1" = lm_f1, "f2" = lm_f2)

  p_sig =s_sig= c()
  sig_level=list()

  if(test_func=="lm"&is.null(random)){
    calc0 <- "Y~1"
    func=get("lm")
    calc1 <- lm_formula[1]
  }else{
    calc0 <- "Y~1+(1|random)"
    func=get("lmer")
    calc1 <- lm_formula[1]
  }

  for (i in 1:ncol(Z)) {

    if(test_func=="lm"&is.null(random)){
      dat <- data.frame(Y = Z[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{

      dat <- data.frame(Y = Z[, i], f1, f2,random)
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
  cls <- rep("sig", ncol(Z))
  if(!is.null(Sig_cutoff)){
    cutoff=sort(p_sig_adj)[length(p_sig_adj)*Sig_cutoff] # adjusted p value 1 problem
    cls[p_sig_adj > cutoff] <- "non-sig"
    cls_df <- data.frame(Sig0 = cls)
    idx_sig <- which(cls=="sig")
    if(length(idx_sig)>1){
      Z_tmp <- Z[, idx_sig]
    }else{
      Z_tmp <- data.frame(Z[, idx_sig] )
      colnames(Z_tmp)=colnames(Z)[idx_sig]
      rownames(Z_tmp)=rownames(Z)
    }
  }else if(Sig){
    cls[p_sig_adj > Sig] <- "non-sig"
    cls_df <- data.frame(Sig0 = cls)
    idx_sig <- which(cls=="sig")
    if(length(idx_sig)>1){
      Z_tmp <- Z[, idx_sig]
    }else{
      Z_tmp <- data.frame(Z[, idx_sig] )
      colnames(Z_tmp)=colnames(Z)[idx_sig]
      rownames(Z_tmp)=rownames(Z)
    }

  }


  p_int_tmp=s_int_tmp  <- rep(NA, ncol(Z_tmp))
  p_f1_tmp=s_f1_tmp <- rep(NA, ncol(Z_tmp))
  p_f2_tmp=s_f2_tmp <- rep(NA, ncol(Z_tmp))

  for (i in 1:ncol(Z_tmp)) {

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_formula["int_null"]
      calc1 <- lm_formula["full"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_formula["int_null"]
      calc1 <- lm_formula["full"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }

    p_int_tmp[i] <- p
    s_int_tmp[i] <- S

    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_formula["f2"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_formula["f2"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }

    p_f1_tmp[i] <- p
    s_f1_tmp[i] <- S


    if(test_func=="lm"&is.null(random)){
      calc0 <- lm_formula["f1"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F")
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-1]
    }else{
      calc0 <- lm_formula["f1"]
      calc1 <- lm_formula["int_null"]
      dat <- data.frame(Y = Z_tmp[, i], f1, f2,random)
      M0 <- func(formula(calc0), data = dat)
      M1 <- func(formula(calc1), data = dat)
      anova_res <- anova(M0, M1, test = "F",refit=FALSE)
      p <- anova_res[2,ncol(anova_res)]
      S <-  anova_res[2,ncol(anova_res)-2]
    }
    p_f2_tmp[i] <- p
    s_f2_tmp[i] <- S
  }


  if (is.null(adj_method)) {
    P_mat[idx_sig, 2] = p_int_tmp
    P_mat[idx_sig, 3] = p_f1_tmp
    P_mat[idx_sig, 4] = p_f2_tmp
    S_mat[idx_sig, 2] = s_int_tmp
    S_mat[idx_sig, 3] = s_f1_tmp
    S_mat[idx_sig, 4] = s_f2_tmp
  }else {
    P_mat[idx_sig, 2] = p.adjust(p_int_tmp,method = adj_method)
    P_mat[idx_sig, 3] = p.adjust(p_f1_tmp,method = adj_method)
    P_mat[idx_sig, 4] = p.adjust(p_f2_tmp,method = adj_method)
    S_mat[idx_sig, 2] = s_int_tmp
    S_mat[idx_sig, 3] = s_f1_tmp
    S_mat[idx_sig, 4] = s_f2_tmp
  }

  cls_int_tmp <- rep("Add", ncol(Z_tmp))
  cls_int_tmp[p_int_tmp < Int] = "Int"
  names(cls_int_tmp) <- colnames(Z_tmp)
  cls_tmp <- rep(NA, ncol(Z))
  idx_add <- which(cls_int_tmp == "Add")
  cls_tmp1 <- cls_int_tmp[idx_add]
  p_f2_tmp1 <- p_f2_tmp[idx_add]
  p_f1_tmp1 <- p_f1_tmp[idx_add]
  cls_tmp1[(p_f2_tmp1 <= F2) & (p_f1_tmp1 > F1)] = "F2"
  cls_tmp1[(p_f2_tmp1 >F2) & (p_f1_tmp1 <= F1)] = "F1"
  cls_int_tmp[idx_add] <- cls_tmp1
  cls_tmp[idx_sig] <- cls_int_tmp
  cls_df$Sig1 <- cls_tmp
  cls_df$varname <- colnames(Z)
  colnames(P_mat) = colnames(S_mat) <- c("Sig", "Interaction", "F1",
                       "F2")

  rownames(P_mat) =rownames(S_mat) <- colnames(Z)
  rownames(cls_df) <- colnames(Z)
  return(list(P_mat=P_mat,S_mat=S_mat,cls_df=cls_df))
}






#' Check the iDAS_2F input factors' name
#'
#' @param f1name a string of factor 1's name
#' @param f2name a string of factor 2's name
#' @param randomname a string of random effect term's name
#' @param f1  A vector of factor 1 variables.
#' @param f2  A vector of factor 2 variables.
#' @param random  A vector of random effect term variables.
#'
#' @return A list of each factor's name and values
#'
#' @examples  #factor_tmp=check_factor_name(f1name,f2name,randomname,f1,f2,random)
check_factor_name=function(f1name,f2name,randomname,f1,f2,random){
  if (is.null(f1name)) {
    f1name = "f1"
  }
  if (is.null(f2name)) {
    f2name = "f2"
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

  if (length(randomname) > 1) {
    random_tmp <- paste(randomname, collapse = "+")
  }else {
    random_tmp <- randomname
  }
  return(list(f1_tmp=f1_tmp,f2_tmp=f2_tmp,random_tmp=random_tmp,f1=f1,f2=f2,random=random))
}
