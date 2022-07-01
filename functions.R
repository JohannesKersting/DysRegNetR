library(DESeq2)
library(edgeR)
library(psych)
library(MASS)
library(ggplot2)
library(ggpubr)
library(methods)
library(data.table)




################################################################################
# Util
################################################################################

getCounts <- function(regulations, matrix){
  lapply(regulations, function(regulation){
    reg <- strsplit(regulation, ":")
    source <- reg[[1]][1]
    target <- reg[[1]][2]
    return(list(source=matrix[source,], target = matrix[target,]))
  })
}

get_model_params <- function(design, meta_cols, network_cols, con_col, shared_covariates){
  
  covariates <- meta_cols[meta_cols %in% all.vars(design)]
  tg_col <- as.character(design[[2]])
  tf_col <- network_cols[1:2][network_cols[1:2] != tg_col]
  
  tf_index <- which(tf_col==all.vars(design))
  tf_coef <- attr(terms(design),"term.labels")[tf_index - 1]
  tf_call <- attr(terms(design),"variables")[[tf_index+1]]
  
  cov_coefs <- sapply(covariates, function(x){
    index <- which(x==all.vars(design))
    attr(terms(design),"term.labels")[index - 1]
  })
  
  list("tg_col"=tg_col, "tf_col"=tf_col, "covariates"=covariates, "shared_covariates" = shared_covariates, "cov_coefs" = cov_coefs, "tf_coef"=tf_coef, "tf_call"=tf_call, "con_col"=con_col)
  
}

get_model_data <- function(tg_counts, tf_counts, meta, size_factors, model_params){
  
  model_data <- data.frame(offset = size_factors)
  
  model_data[[model_params$tg_col]] <- tg_counts
  model_data[[model_params$tf_col]] <- tf_counts / size_factors
  
  model_data <- cbind(model_data, meta[model_params$covariates])
  model_data <- cbind(model_data, meta[model_params$con_col])
  
  model_data
}

get_sd <- function(model, predicted){
  if(!is.null(model$theta)){
    sd <- sqrt(predicted + (predicted ^ 2) / model$theta)
    return(sd)
  }
  else{
    return(sqrt(predicted))
  }
}

get_transformed_tf_counts <- function(regulations, model_params, counts, size_factors){
  tfs <- unique(sapply(regulations, function(regulation){
  strsplit(regulation, ":")[[1]][1]
  }))
  
  result <- t(sapply(tfs, function(tf){
    data <- list()
    data[[model_params$tf_col]] <- (counts[tf,])/size_factors
    eval(model_params$tf_call, envir = data)
  }))
  
  return(result)
}

get_pca <- function(counts, size_factors, normalized = T){
  
  if(normalized){
    counts <- sweep(counts, 2, size_factors, FUN = '/')
  }
  
  prcomp(t(counts), scale. = T)
  
}

get_qq_data <- function(residuals, regulation, residual_types, condition_type){
  
  if(length(residual_types)!=1){
    
    if(!is.list(residuals)){
      stop("multiple residuals types require a named list, containing the individual residual matrices")
    }
    qq_data <- rbindlist(lapply(residual_types, function(residual_type){
      data.table("residual"=residuals[[residual_type]][regulation,], "condition" = condition_type, "residual_type" = residual_type)
    }))
  }
  
  else{
    if(is.list(residuals)){
      qq_data <- data.table("residual"=residuals[[residual_types]][regulation,], "condition" = condition_type, "residual_type" = residual_types)
    }else{
      qq_data <- data.table("residual"=residuals[regulation,], "condition" = condition_type, "residual_type" = residual_types)
    }
  }
  
  return(qq_data)
}

################################################################################
# Normalization
################################################################################

get_size_factors <- function(counts, norm_method){
  
  if(norm_method == "mrn") {
    
    return (mrn_size_factors(counts))
    
  } else if(norm_method == "tmm"){
    
    return (tmm_size_factors(counts))
    
  } else if(norm_method == "uq"){
    
    return (uq_size_factors(counts))
    
  } else if(norm_method == "cpm"){
    
    return (cpm_size_factors(counts))
    
  } else if(norm_method == "raw"){
    
    return (raw_size_factors(counts))
    
  } else {
    
    stop(paste("unknown norm_method ",norm_method))
    
  }
}

mrn_size_factors <- function(counts){
  col_data <- data.frame(sample=c(rep("sample",times = ncol(counts))))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts,col_data, ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  DESeq2::sizeFactors(dds)
}

tmm_size_factors <- function(counts){
  norm_factors <- edgeR::calcNormFactors(counts, method="TMM")
  cpm_size_factors(counts, norm_factors)
}

uq_size_factors <- function(counts){
  norm_factors <- edgeR::calcNormFactors(counts, method="upperquartile")
  cpm_size_factors(counts, norm_factors)
}

cpm_size_factors <- function(counts, norm_factors = 1){
  (colSums(counts) * norm_factors) / 1000000
}

raw_size_factors <- function(counts){
  size_factors <- rep(1, ncol(counts))
  names(size_factors) <- colnames(counts)
  size_factors
}



################################################################################
# Network filtering
################################################################################

pre_process_network <- function(network, gene_ids, silent = F){
  
  n_regs <- nrow(network)
  network <- network[network[[1]]!=network[[2]],]
  
  if(nrow(network) < n_regs && !silent){
    print(paste(n_regs-nrow(network), "self-regulating pairs were removed from the network"))
  }
  
  n_regs <- nrow(network)
  network <- unique(network)
  
  if(nrow(network) < n_regs && !silent){
    print(paste(n_regs-nrow(network), "duplicate regulations were removed from the network"))
  }
  
  n_regs <- nrow(network)
  network <- network[network[[1]] %in% gene_ids & network[[2]] %in% gene_ids,]
  
  if(nrow(network) < n_regs && !silent){
    print(paste(n_regs-nrow(network), "regulations with unmatchable gene ids were removed from the network"))
  }
  
  return(network)
}


network_cor_test <- function(x, y, cor_method){
  if(cor_method == "spearman"){
    cor_test <- psych::corr.test(x, y, adjust="none", method = "spearman")
  } else if(cor_method == "log-pearson"){
    cor_test <- psych::corr.test(log(x), log(y), adjust="none", method = "pearson")
  } else{
    stop(paste("unknown cor_method", cor_method))
  }
  c("cor"=cor_test$r,"pval"=cor_test$p)
}


filter_network <- function(network, counts, size_factors, cor_method = c("log-pearson","spearman"), min_cor, max_cor_pval){
  
  
  if(length(size_factors == 1) && is.na(size_factors)) 
    stop("size_factors need to be estimated before correlation filtering")
  
  cor_method <- match.arg(cor_method)
  
  cor_tests <- t(apply(network, 1, function(regulation){
    network_cor_test(counts[regulation[1],] / size_factors, counts[regulation[2],] / size_factors, cor_method)
    
  }))
  
  network$cor <- cor_tests[,"cor"]
  network$pval <- cor_tests[,"pval"]
  network$use <- abs(network$cor) >= min_cor & network$pval <= max_cor_pval
  
  return(network)
}


################################################################################
# Model building
################################################################################

get_models <- function(network, counts, size_factors, design, meta, con_col, model_params, 
                          cor_filter = T, 
                          model_type = c("negative-binomial","poisson"), 
                          cooks_filter = T,
                          case_for_disp = T,
                          silent = T
                          ){
  
  model_type <- match.arg(model_type)

  if(length(size_factors == 1) && is.na(size_factors)) 
    stop("size_factors need to be estimated before model building")
  
  if(cor_filter && is.null(network$use))
    stop("network wasn't filtered by correlation, turn off with cor_filter = FALSE")
  
  
  model_use <- T
  if(!is.null(network$use) && cor_filter){
    
    model_use <- network$use
    
    if(sum(model_use) == 0){
      print("no regulation passed the correlation filtering")
      return(list())
    }
  }
  
  if(model_type == "negative-binomial"){
    model_function <- nb_model
  } else{
    model_function <- p_model
  }
  
  covariates <- model_params$covariates
  tg_col <- model_params$tg_col
  tf_col <- model_params$tf_col
  
  if(model_type == "negative-binomial" && case_for_disp){
    
    train_offset <- size_factors
    train_counts <- counts
    train_meta <- meta
    
    design <- update(design,  paste("~ . +",con_col,"-",model_params$tf_coef, "+", model_params$tf_coef, ":", con_col))
    
    for(covariate in covariates){
      if(!(covariate %in% model_params$shared_covariates)){
        design <- update(design,  paste("~ . -",model_params$cov_coefs[covariate], "+", model_params$cov_coefs[covariate], ":", con_col))
      }
    }
    
  } else {
    
    is_control <- !meta[[con_col]]
    
    train_offset <- size_factors[is_control]
    train_counts <- counts[,is_control]
    train_meta <- meta[is_control,]
    
  }
  

  
  
  result <- apply(network[model_use,], 1, function(regulation){
    
    
    tg_train_counts <- train_counts[regulation[[tg_col]],]
    tf_train_counts <- train_counts[regulation[[tf_col]],]
    
    train_data <- get_model_data(tg_counts = tg_train_counts, 
                                 tf_counts = tf_train_counts, 
                                 meta = train_meta, 
                                 size_factors = train_offset, 
                                 model_params = model_params)
    
    model <- get_model(data = train_data, design = design, model_function = model_function,
                       model_params  = model_params, cooks_filter = cooks_filter, silent = silent)
    
    
    if(cooks_filter && !is.null(model)){
      model$n_outlier <- nrow(train_data) - length(model$y)
    }
    
    slim_model(model)
    
  
  })
  
  
  
  names(result) <- paste(network[model_use,][[tf_col]], network[model_use,][[tg_col]], sep = ":")
  result
}


get_model <- function(data, design, model_function, model_params, cooks_filter = T, no_poisson = T, silent = F) {
  
  use_sample <- T
  build <- T
  outliers <- c()
  force_outlier_poisson <- F
  
  if(no_poisson && identical(model_function, nb_model)){
    force_outlier_poisson <- T
  }
  
  while(build){
    
    data <- data[use_sample,]
    
    if(silent){
      model <- suppressWarnings(model_function(data = data, design = design, silent = silent))
    } else{
      model <- model_function(data = data, design = design, silent = silent)
    }
    
    if(cooks_filter && !is.null(model)){
      
      #clipped_abs_resids <- pmin(sqrt(nrow(data)), abs(resid(model, type="pearson")))
      disp <- sum(resid(model, type="pearson")^2) / df.residual(model)
      dist <- cooks.distance(model, dispersion = disp)
      cutoff <- qf(0.99, nrow(data), df.residual(model), lower.tail=T)
      
      if(any(dist>cutoff) || (is.null(model$theta) && force_outlier_poisson)){
        
        build <- T
        max_index <- which.max(dist)
        
        use_sample <- rep(T, nrow(data))
        use_sample[max_index] <- F
        outliers <- c(outliers,rownames(data)[max_index])
        
      }else{
        build <- F
      }
    } else{
      build <- F
    }
    
  }
  
  if(! is.null(model)){
    model$outliers <- outliers
    model$tf_mean <- mean(eval(model_params$tf_call, envir = data[data[[model_params$con_col]]==F,])) 
  }
  
  return(model)
  
}


nb_model <- function(data, design, silent = F){
  
  nb_design <- update(design, ~ . + offset(log(offset)))
  control <- glm.control(maxit=50)
  
  model <- tryCatch(glm.nb(nb_design, data, control = control), error=function(e) {
    if(!silent){
      print(e) 
      print("")
    }
    NULL})
  if(is.null(model)){
    return(p_model(data, design, silent))
  }
  
  return(model)
}


p_model <- function(data, design, silent = F){
  
  model <- tryCatch(glm(design, family ="poisson", data, offset=log(offset)), error=function(e) {
    if(!silent){
      print(e) 
      print("")
    }
    NULL})
  
  return(model)
}


slim_model <- function(model) {
  model$y = c()
  model$model = c()
  
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()  
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  
  
  model$family$variance = c()
  model$family$dev.resids = c()
  model$family$aic = c()
  model$family$validmu = c()
  model$family$simulate = c()
  #attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model$offset = c()
  
  model
}

gof_filter <- function(models, max_statistic = 1.5, min_pval = 0.001){
tests <- gof_test(models)
return(models[ !is.nan(tests[,"statistic"]) & tests[,"statistic"] <= max_statistic & tests[,"pval"] >= min_pval])
}

gof_test <- function(models){
  
  t(sapply(models, function(model){
    
    statistic <- NaN
    pval <- NaN
    
    if(!is.null(model)){
      dev <- deviance(model)
      df <- model$df.residual
      statistic <- dev/df
      pval <- 1 - pchisq(dev, df)
    }
    
    return(c("statistic"=statistic, "pval" = pval))
  }))
  

  
}

################################################################################
# Residuals
################################################################################

get_residuals <- function(models, counts, size_factors, samples, model_params, meta,
                          residual_type = c("deviance", "pearson", "anscombe", "raw"),
                          condition = c("all", "case", "control"),
                          regulations = NULL){
  
  residual_type <- match.arg(residual_type)
  condition <- match.arg(condition)
  
  if(is.null(regulations)){
    regulations <- names(models)
  }
  
  tg_col <- model_params$tg_col
  tf_col <- model_params$tf_col
  
  if(condition == "case"){
    counts <- counts[,samples]
    size_factors <- size_factors[samples]
    meta <- meta[samples,]
  } else if(condition == "control"){
    counts <-  counts[,!samples]
    size_factors <- size_factors[!samples]
    meta <- meta[!samples,]
  } 
  
  result <- t(sapply(regulations, function(regulation){
    
    model <- models[[regulation]]
    
    if(!is.null(model)){
      
      reg <- strsplit(regulation, ":")
      source <- reg[[1]][1]
      target <- reg[[1]][2]
      
      tg_counts <- counts[target,]
      tf_counts <- counts[source,]
      
      data <- get_model_data(tg_counts = tg_counts, 
                             tf_counts = tf_counts, 
                             meta = meta, 
                             size_factors = size_factors, 
                             model_params = model_params)
      
      data[[model_params$con_col]] <- F
      
      y_hat <- predict(model, newdata=data, type = "response")

      y <- tg_counts
      
      if(is.null(model$theta)){
        residuals <- p.resid(y, y_hat, residual_type)
      } else {
        residuals <- nb.resid(y, y_hat, model$theta, residual_type)
      }
      
    } else {
      residuals <- rep(NaN,length(size_factors))
      names(residuals) <- names(size_factors)
    }
    
    
    return(residuals)
    
  }))
  
  rownames(result) <- regulations
  return(result)
   
}


nb.resid <- function(y, y_hat, theta, type){
  
  if(type=="pearson"){
    
    return(nb.resid_pearson(y, y_hat, theta))
    
  } else if(type=="deviance"){
    
    return(nb.resid_deviance(y, y_hat, theta))
    
  } else if(type=="anscombe"){
    
    return(nb.resid_anscombe(y, y_hat, theta))
    
  } else if(type=="raw"){
    
    return(resid_raw(y, y_hat))
    
  } else {
    
    return(NULL)
    
  }
}

nb.resid_pearson <- function(y, y_hat, theta){
  (y - y_hat) / sqrt(y_hat + y_hat * y_hat * (1/theta))
}

nb.resid_deviance <- function(y, y_hat, theta){
  residuals <- 
    sign( y - y_hat ) *
    sqrt(
      2 * (
        y * log( y / y_hat ) - 
          ( y + theta ) * log((y + theta) / (y_hat + theta))
      )
    )
  return(residuals)
}

nb.resid_anscombe <- function(y, y_hat, theta){
  a <- 1/theta
  
  top1 <-
    (3/a) * (
      (1 + a * y)^(2/3) -
        (1 + a * y_hat)^(2/3)
    )
  top2 <- 3 * ( y^(2/3) - y_hat^(2/3) )
  bot <- 2 * (y_hat + a * y_hat^2)^(1/6)
  
  (top1 +  top2) / bot 
}


p.resid <- function(y, y_hat, type){
  
  if(type=="pearson"){
    
    return(p.resid_pearson(y, y_hat))
    
  } else if(type=="deviance"){
    
    return(p.resid_deviance(y, y_hat))
    
  } else if(type=="anscombe"){
    
    return(p.resid_anscombe(y, y_hat))
    
  } else if(type=="raw"){
    
    return(resid_raw(y, y_hat))
    
  } else {
    
    return(NULL)
    
  }
}

p.resid_pearson <- function(y, y_hat){
  (y - y_hat) / sqrt(y_hat)
}

p.resid_deviance <- function(y, y_hat){
  residuals <- 
    sign(y-y_hat) * 
    sqrt(
      2 * (
        y * log(y / y_hat) - (y - y_hat)
      )
    )
  return(residuals)
}

p.resid_anscombe <- function(y, y_hat){
  (3/2) * ( y^(2/3) * y_hat^(-1/6) - y_hat^(1/2) ) 
}

resid_raw <- function(y, y_hat){
  y - y_hat
}



################################################################################
# Results
################################################################################

normality_test <- function(residuals){
  
  if(!is.matrix(residuals)){
    residuals <- matrix(residuals, nrow = 1)
  }
  
  t(apply(residuals, 1, function(x){

    test <- c(NaN, NaN)
    if(sum(!is.na(x)) > 3){
      test <- shapiro.test(x)
      test <- c(test["statistic"], test["p.value"])
    } 
    test <- as.numeric(test)
    names(test) <- c("statistic", "pval")
    test
  }))
  
}


resid_z_scores <- function(resid_control, resid_case){
  (resid_case - mean(resid_control)) / sd(resid_control)
}

get_z_scores <- function(residuals, samples, condition = c("all", "case", "control")){
  condition <- match.arg(condition)
  
  if(condition == "control"){
    means <- apply(residuals, 1, mean, na.rm = T)
    sds <- apply(residuals, 1, sd, na.rm = T)
    
    return((residuals - means) / sds)
    
  } else {
    
    means <- apply(residuals[,!samples, drop = F], 1, mean)
    sds <- apply(residuals[,!samples, drop = F], 1, sd)
    
    if(condition == "case"){
      
      return((residuals[,samples, drop = F] - means) / sds)
      
    }else{
      
      return((residuals - means) / sds)
      
    }
  }
}


filter_z_scores <- function(z_scores,  models, model_params, 
                           bonferroni_alpha = 0.01, transformed_tf_counts = NULL,
                           criterion = c("direction", "leverage", "none")){
  
  criterion <- match.arg(criterion)
  
  result <- t(sapply(rownames(z_scores), function(regulation){
    
   row <- z_scores[regulation,]  
   p_vals <- pnorm(abs(row), lower.tail = F)*2
   p_vals <- p.adjust(p_vals, method = "bonferroni")
   
   if(criterion != "none"){
     
     model <- models[[regulation]]
     
     if(is.null(model)){
       return(row)
     }
     
     coef <- coefficients(model)[model_params$tf_coef]
     if(is.na(coef)){
       coef <- coefficients(model)[paste0(model_params$con_col,"FALSE:",model_params$tf_coef)]
     }
     if(is.na(coef)){
       stop("can not find tf coefficient")
     }
     
     direction <- sign(coef)
     row[direction * row > 0] <- 0

     
   }
   
   if(criterion == "leverage"){
     source <-  strsplit(regulation, ":")[[1]][1]
     row[transformed_tf_counts[source,] < model$tf_mean] <- 0
   }
   
   row[p_vals >= bonferroni_alpha] <- 0
   row
   
  }))
  
  rownames(result) <- rownames(z_scores)
  result
  
}


filter_residuals <- function(residuals,  models, counts, size_factors, meta, model_params, 
                            bonferroni_alpha = 0.01, transformed_tf_counts = NULL,
                            criterion = c("direction", "leverage", "none")){
  
  
  criterion <- match.arg(criterion)
  tg_col <- model_params$tg_col
  tf_col <- model_params$tf_col
  
  result <- t(sapply(rownames(residuals), function(regulation){

    model <- models[[regulation]]
    row <- residuals[regulation,]  
    
    if(is.null(model)){
      return(row)
    }
    

    reg <- strsplit(regulation, ":")
    source <- reg[[1]][1]
    target <- reg[[1]][2]
    
    tg_counts <- counts[target,]
    tf_counts <- counts[source,]
    
    data <- get_model_data(tg_counts = tg_counts, 
                           tf_counts = tf_counts, 
                           meta = meta, 
                           size_factors = size_factors, 
                           model_params = model_params)
    
    data[[model_params$con_col]] <- F
    
    y_hat <- predict(model, newdata=data, type = "response")
    
    if(!is.null(model$theta)){
      cdf_lower <- pnbinom(tg_counts, size = model$theta , mu = y_hat)
      cdf_upper <- 1 - pnbinom(tg_counts-1, size = model$theta , mu = y_hat)
    }
    else{
      cdf_lower <- ppois(tg_counts, lambda = y_hat)
      cdf_upper <- 1 - ppois(tg_counts-1, lambda = y_hat)
    }

    
    p_vals <- 2*pmin(0.5, cdf_lower, cdf_upper)
    
    
    p_vals <- p.adjust(p_vals, method = "bonferroni")
    
    if(criterion != "none"){
      

      coef <- coefficients(model)[model_params$tf_coef]
      if(is.na(coef)){
        coef <- coefficients(model)[paste0(model_params$con_col,"FALSE:",model_params$tf_coef)]
      }
      if(is.na(coef)){
        stop("can not find tf coefficient")
      }
      
      direction <- sign(coef)
      row[direction * row > 0] <- 0
      
      
    }
    
    if(criterion == "leverage"){
      source <-  strsplit(regulation, ":")[[1]][1]
      row[transformed_tf_counts[source,] < model$tf_mean] <- 0
    }
    
    row[p_vals >= bonferroni_alpha] <- 0
    row
    
  }))
  
  rownames(result) <- rownames(residuals)
  result
  
}

################################################################################
# Plots
################################################################################


count_scatter <- function(regulations, count_matrix, size_factors, samples, model_params, meta , models, 
                          condition = c("all", "case", "control"),
                          norm_tg = T,
                          norm_tf = T,
                          titles = NULL,
                          show_fit = T,
                          show_sd = F,
                          log_tg = F,
                          log_tf = F,
                          alpha = 1,
                          sd_alpha = 0.2,
                          case_outlier = NULL,
                          control_outlier = T,
                          interaction_model = F,
                          ...
                          ){
  
  condition = match.arg(condition)
  
  
  color_values <- c() 
  shape_values <- c()
  
  
  if(condition == "case"){
    
    counts <- getCounts(regulations, count_matrix[,samples])
    size_factors <- size_factors[samples]
    meta <- meta[samples,]
    type <- "case"
    color_values <- c(color_values, "case" = "#F0E442")
    shape_values <- c(shape_values, "case" = 16)
    
  } else if(condition == "control"){
    
    counts <- getCounts(regulations, count_matrix[,!samples])
    size_factors <- size_factors[!samples]
    meta <- meta[!samples,]
    type <- "control"
    color_values <- c(color_values, "control" = "#56B4E9")
    shape_values <- c(shape_values, "control" = 16)
    
  } else{
    
    counts <- getCounts(regulations, count_matrix)
    type <- rep("control", length(samples))
    type[samples] <- "case"
    color_values <- c(color_values, "control" = "#56B4E9",
      "case" = "#F0E442")
    shape_values <- c(shape_values, "control" = 16, "case" = 16)
    
  }
  
  if(condition != "case" && control_outlier){
    color_values <- c(color_values, "control outlier" = "#0072B2")
    shape_values <- c(shape_values, "control outlier" = 4)

  } 
  if(condition != "control" && !is.null(case_outlier)){
    color_values <- c(color_values ,"case outlier" = "#D55E00")
    shape_values <- c(shape_values, "case outlier" = 4)
  } 
  if(show_fit){
    color_values <- c(color_values, "fitted" = "#009E73") 
    shape_values <- c(shape_values, "fitted" = 3)
  }
  
  
  if(is.null(titles)){
    titles <- regulations
  }

  names(counts) <- regulations
  
  plot_list <- lapply( c(1:length(regulations)), function(i){
    
    regulation <- regulations[i]
    
    reg <- strsplit(regulation, ":")
    labs <- paste(reg[[1]], "counts")

    source_counts <- counts[[regulations[i]]]$source
    target_counts <- counts[[regulations[i]]]$target
    
    patient_ids <- names(source_counts)
    
    data <- data.table(
      source_counts=source_counts,
      target_counts=target_counts,
      type = type
    )
    
    model <- models[[regulation]]
    sd_data <- NULL
    
    if(!is.null(model) && control_outlier && condition != "case"){
      data[patient_ids %in% model$outliers, type:="control outlier"]
    }
    
    if(!is.null(model) && !is.null(case_outlier) && condition != "control"){
      reg_case_outlier <- colnames(case_outlier)[case_outlier[regulation,] != 0]
      data[patient_ids %in% reg_case_outlier, type:="case outlier"]
    }
    
      
    if(show_fit && !is.null(model)){
      
      fit_data <- get_model_data(tg_counts = target_counts, 
                                 tf_counts = source_counts,
                                 meta = meta,
                                 size_factors = size_factors,
                                 model_params = model_params
                                 )
      
      if(!interaction_model){
        fit_data[[model_params$con_col]] <- F
      }
      
      predicted <- predict(model, newdata=fit_data, type="response")
      predicted_data <- data.table(target_counts = predicted, source_counts = source_counts , type = "fitted")
      data <- rbind(data, predicted_data)
      
      if(show_sd){
        sd <- get_sd(model, predicted)
        sd_data <- data.table(source_counts=source_counts, type = "fitted")
        sd_data[,ymin:= predicted - sd]
        sd_data[,ymax:= predicted + sd]
      }
      
    }
    
    
    
    if(norm_tg){
      data[ ,target_counts :=(target_counts)/size_factors]
      if(!is.null(sd_data)){
        sd_data$ymin<- sd_data$ymin/size_factors
        sd_data$ymax <- sd_data$ymax/size_factors
      }
    }
    
    if(norm_tf){
      data[ ,source_counts :=(source_counts)/size_factors]
      if(!is.null(sd_data)){
        sd_data[ ,source_counts :=(source_counts )/size_factors]
      }
    }

    
    g <- ggplot(data, aes(x=source_counts, y=target_counts, color = type, shape = type)) +
      geom_point(alpha=alpha) +
      labs(title = titles[i], x = labs[1], y= labs[2], color = "Type", shape = "Type" ) +
      theme_bw() +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values)
    
    if(!is.null(sd_data)){
      g <- g +
        geom_errorbar(data=sd_data, aes(ymin=ymin, ymax=ymax, x=source_counts), alpha = sd_alpha)
    }
    
    if(log_tg){
      g <- g + scale_y_log10()
    }
    
    if(log_tf){
      g <- g + scale_x_log10()
    }
    
    return(g)
  })
  return(plot_list)
}

plot_pca <- function(pca_data, samples, condition = c("all", "case", "control"), alpha = 0.5, scree = F, ...){
  
  condition <- match.arg(condition)
  
  color_values <- c() 
  
  if(condition == "case"){
    type = "case"
    color_values <- c(color_values, "case" = "#D55E00")
    
    if(length(samples) == nrow(pca_data$x)){
      pca_data$x <- pca_data$x[samples,]
    }
  }
  else if( condition == "control"){
    color_values <- c(color_values, "control" = "#0072B2")
    type = "control"
    
    if(length(samples) == nrow(pca_data$x)){
      pca_data$x <- pca_data$x[!samples,]
    }
  }
  else{
    
    if(length(samples)!= nrow(pca_data$x)){
      stop("Sample number doesn't match the number of samples in the pca data")
    }
    
    type <- rep("control", length(samples))
    type[samples] <- "case"
    color_values <- c(color_values, "control" = "#0072B2",
                      "case" = "#D55E00")
  }
  
  var <- pca_data$sdev^2
  var_per <- round((var/sum(var))*100, 1)
  
  plot_data <- data.table("PC1" = pca_data$x[,1], "PC2" = pca_data$x[,2], "type" = type)
  
  g <- ggplot(plot_data, aes(x = PC1, y = PC2, color = type)) +
    geom_point(alpha = alpha) +
    xlab(paste("PC1 -", var_per[1], "%"))+
    ylab(paste("PC2 -", var_per[2], "%"))+
    scale_color_manual(values = color_values)+
    labs(color = "Condition")+
    theme_bw()
    
  if(scree){
    
    var_per <- var_per[var_per > 0.001]
    s <- ggplot(data.frame(var = var_per, pc = 1:length(var_per)), aes(y=var, x = pc)) + 
      geom_bar(stat="identity") +
      theme_bw()
    
    g <- ggarrange(g, s,
               labels = c("A", "B"),
               ncol = 2)
  }
  
  return(g)
}


plot_qq <- function(regulations, residuals, samples, 
                    condition = c("all", "case", "control"), 
                    residual_types = c("deviance", "pearson", "anscombe", "raw"), 
                    titles = NULL,
                    sn_line = T,
                    ...
                    ){
  
  condition <- match.arg(condition)
  residual_types <- match.arg(residual_types, several.ok = T)
  
  
  if(is.null(regulations)){
    if(is.matrix(residuals)){
      regulations <- rownames(residuals)
    } else{
      regulations <- rownames(residuals[[1]])
    }
  }
  
  
  if(is.null(titles)){
    titles <- regulations
  }
  
  
  color_values <- c() 
  
  if(condition == "case"){
    condition_type <- "case"
    color_values <- c(color_values, "case" = "#D55E00")
    
  }
  else if( condition == "control"){
    color_values <- c(color_values, "control" = "#0072B2")
    condition_type <- "control"
    
  }
  else{
    
    condition_type <- rep("control", length(samples))
    condition_type[samples] <- "case"
    color_values <- c(color_values, "control" = "#0072B2",
                      "case" = "#D55E00")
  }
  
  lapply(1:length(regulations), function(i){
    regulation <- regulations[i]
    
    plot_data <- get_qq_data(residuals, regulation, residual_types, condition_type)
    
    g <- ggplot(plot_data, aes(sample=residual, color = condition)) +
      stat_qq()+
      stat_qq_line()+  
      geom_vline(xintercept=0)+
      labs(title = titles[i], x = "Normal distribution", y="Residuals", color = "Condition")+
      scale_color_manual(values = color_values)+
      facet_wrap(vars(residual_type), scales = "free")+
      theme_bw()
    
    if(sn_line){
      g <- g + geom_abline(slope=1,intercept=0)
    }
    g
  })
  
}

norm_test_plot <- function(residuals, samples,                     
                           condition = c("all", "case", "control"), 
                           residual_types = c("deviance", "pearson", "anscombe", "raw"),  annotate = T, ...){
  
  condition <- match.arg(condition)
  residual_types <- match.arg(residual_types, several.ok = T)
  
  color_values <- c() 
  
  if(condition == "case"){
    color_values <- c(color_values, "case" = "#D55E00")
    
  }
  else if( condition == "control"){
    color_values <- c(color_values, "control" = "#0072B2")
    
  }
  else{
    color_values <- c(color_values, "control" = "#0072B2",
                      "case" = "#D55E00")
  }
  
  plot_data <- rbindlist(lapply(residual_types, function(residual_type){
    
    if(condition == "all"){
      control_tests <- as.data.table(normality_test(residuals[[residual_type]][,!samples]))
      case_tests <- as.data.table(normality_test(residuals[[residual_type]][,samples]))
      
      control_tests$condition_type <- "control"
      case_tests$condition_type <- "case"
      
      tests <- rbind(control_tests, case_tests)
    } else {
      tests <- as.data.table(normality_test(residuals[[residual_type]]))
      tests$condition_type <- condition
    }

    tests$residual_type  <- residual_type
    tests
  }))
  
  stat_plot <- ggplot(plot_data,aes(y=statistic, x = residual_type, color = condition_type)) +
    geom_boxplot() +
    scale_color_manual(values = color_values)+
    labs(x = "Residual type", y = "Statistics", color = "Condition")+
    theme_bw()
  
  pval_plot <- ggplot(plot_data,aes(y=pval, x = residual_type, color = condition_type)) +
    geom_boxplot() +
    scale_color_manual(values = color_values)+
    labs(x = "Residual type", y = "P-values", color = "Condition")+
    theme_bw()
  
  plot <- ggarrange(stat_plot, pval_plot, labels = c("A","B"), nrow = 1, common.legend = T)
  if(annotate){
    plot <- annotate_figure(plot, top = text_grob("Shapiro-Wilk test of normality", face = "bold", size = 14))
  }
 plot
}