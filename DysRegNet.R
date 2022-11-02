

source("functions.R")



setClass("DysRegNetR",
         slots = c(counts = "matrix",
                   design = "formula",
                   norm_method = "character",
                   size_factors = "numeric",
                   network = "data.frame",
                   meta = "data.frame",
                   counts_add_one = "logical",
                   con_col = "character",
                   
                   cor_filter = "logical",
                   cor_method = "character",
                   min_cor = "numeric",
                   max_cor_pval = "numeric",
                   
                   model_type = "character",
                   case_for_disp = "logical",
                   shared_covariates = "character",
                   cooks_filter = "logical",
                   
                   gof_filter = "logical",
                   max_gof_stat = "numeric",
                   min_gof_pval = "numeric",
                   
                   norm_filter = "logical",
                   min_norm_stat = "numeric",
                   min_norm_pval = "numeric",
                   
                   rating_method = "character",
                   residual_type = "character",
                   outlier_criterion = "character",
                   bonferroni_alpha = "numeric",
                   
                   model_params = "list",
                   models = "list",
                   n_removed_models = "numeric"
                   )
)



setValidity("DysRegNetR", function(object){
  
  errors = character()
  
  if(any(is.na(object@counts)) ||
     !is.numeric(object@counts) || 
     any(object@counts %% 1 != 0) || 
     any(object@counts < 0)){
    msg = "counts must contain non-negative integers"
    errors = c(errors, msg)
  }
  
  if(is.null(rownames(object@counts))){
    msg = "counts must have gene-ids as row names"
    errors = c(errors, msg)
  }
  
  if(is.null(colnames(object@counts))){
    msg = "counts must have patient-ids as column names"
    errors = c(errors, msg)
  }
  
  if(ncol(object@network) != 2){
    msg = "network must have exactly two columns"
    errors = c(errors, msg)
  } else if (!any( network[[1]] %in% rownames(object@counts) & 
                   network[[2]] %in% rownames(object@counts) )
             ){
    msg = "none of the pairs in the network can be matched to gene-ids in the row names of counts"
    errors = c(errors, msg)
  }
  
  if(!all(colnames(object@network) %in% all.vars(object@design))){
    msg = "network column names must be part of the design formula"
    errors = c(errors, msg)
  }
  
  if(!any( colnames(object@network) == object@design[[2]] )){
    msg = "one of the network column names must be the outcome of the design formular"
    errors = c(errors, msg)
  }
  
  if(!object@norm_method %in% c("mrn", "tmm", "uq", "cpm", "raw")){
    msg = 'norm_method must be one of "mrn", "tmm", "uq", "cpm", "raw"'
    errors = c(errors, msg)
  }
  
  if(length(object@con_col) != 1){
    msg = 'con_col must be a character vector with only one element'
    errors = c(errors, msg)
  } else if(object@con_col %in% all.vars(object@design)){
    msg = paste0("the con_col '", object@con_col ,"' must not be part of the design formula")
    errors = c(errors, msg)
  } else if (is.null(object@meta[[object@con_col]])){
    msg = paste0("con_col '", object@con_col ,"' doesn't match to a column name in meta")
    errors = c(errors, msg)
  } else if(!is.numeric(meta[[object@con_col]]) && !is.logical(meta[[object@con_col]])){
    msg = paste0("the column '", object@con_col ,"' in meta must be numeric or logical")
    errors = c(errors, msg)
  } else if(is.numeric(meta[[object@con_col]]) && !all(meta[[object@con_col]] %in% c(0,1))){
    msg = paste0("the numeric column '", object@con_col ,"' in meta must only contain 1's (case samples) and 0's (control samples)")
  }
  
  if(!all(all.vars(object@design) %in% c(colnames(object@network), colnames(object@meta)))){
    msg = "all variables in the design formula must be a column name in network or meta "
    errors = c(errors, msg)
  }
  
  if(length(object@shared_covariates)!=0){
    if( !all(object@shared_covariates %in% all.vars(object@design))){
      msg = "all covariates in shared_covariates must appear in the design formula "
      errors = c(errors, msg)
    }
    if(any(object@shared_covariates %in% colnames(object@network))){
      msg = "the outcome and the transcription factor covariate can't be in shared_covariates"
      errors = c(errors, msg)
    }
  }
  
  if(is.null(rownames(object@meta))){
    msg = "meta must have patient-ids as rownames"
    errors = c(errors, msg)
  } else if(length(rownames(object@meta)) != length(colnames(object@counts)) ||
            any(rownames(object@meta) != colnames(object@counts))){
    msg = "row names of meta must match the column names of counts regarding content and sequence"
    errors = c(errors, msg)
  }

  
  if(length(object@cor_filter)!=1){
    msg = 'cor_filter must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@cooks_filter)!=1){
    msg = 'cooks_filter must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@counts_add_one)!=1){
    msg = 'counts_add_one must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@case_for_disp)!=1){
    msg = 'case_for_disp must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@min_cor)!=1 || object@min_cor > 1 || object@min_cor < 0){
    msg = 'min_cor must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  if(length(object@max_cor_pval)!=1 || object@max_cor_pval > 1 || object@max_cor_pval < 0){
    msg = 'max_cor_pval must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  if(length(object@gof_filter)!=1){
    msg = 'gof_filter must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@max_gof_stat)!=1 || object@max_gof_stat < 0){
    msg = 'max_gof_stat must be a single positive numeric'
    errors = c(errors, msg)
  }
  
  if(length(object@min_gof_pval)!=1 || object@min_gof_pval > 1 || object@min_gof_pval < 0){
    msg = 'min_gof_pval must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  
  if(length(object@norm_filter)!=1){
    msg = 'norm_filter must be TRUE or FALSE'
    errors = c(errors, msg)
  }
  
  if(length(object@min_norm_stat)!=1 || object@min_norm_stat < 0 || object@min_norm_stat > 1){
    msg = 'min_norm_stat must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  if(length(object@min_norm_pval)!=1 || object@min_norm_pval > 1 || object@min_norm_pval < 0){
    msg = 'min_norm_pval must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  
  if(length(object@bonferroni_alpha)!=1 || object@bonferroni_alpha > 1 || object@bonferroni_alpha < 0){
    msg = 'bonferroni_alpha must be a single numeric between 0 and 1'
    errors = c(errors, msg)
  }
  
  if(length(object@size_factors) !=1 ){
    if(length(object@size_factors) != nrow(meta)){
      msg = 'the size_factors must have the same length as the number of samples'
      errors = c(errors, msg)
      
    } else if(is.null(names(size_factors))){
      msg = 'size_factors must be a named vector'
      errors = c(errors, msg)
    
    } else if(any(names(object@size_factors) != rownames(meta))){
      msg = 'size_factors must be a named vector, the names must match the rownames of meta'
      errors = c(errors, msg)
      
    } else if(any(object@size_factors == 0)){
      msg = 'size_factors must not include zeros'
      errors = c(errors, msg)
    }

  }
  
  if(length(errors) == 0) TRUE else errors
  
})



DysRegNetR <- function(counts,
                       network,
                       meta,
                       counts_add_one = TRUE,
                       con_col = "condition",
                       design = tg ~ log(tf), 
                       norm_method = c("mrn", "tmm", "uq", "cpm", "raw"),
                       size_factors = NaN,
                       cor_filter = TRUE,
                       cor_method = c("log-pearson", "spearman"),
                       min_cor = 0.5,
                       max_cor_pval = 0.05,
                       model_type = c("negative-binomial","poisson"),
                       cooks_filter = TRUE,
                       case_for_disp = TRUE,
                       shared_covariates = character(),
                       gof_filter = TRUE,
                       max_gof_stat = 1.5,
                       min_gof_pval = 0.001,
                       norm_filter = FALSE,
                       min_norm_stat = 0.9,
                       min_norm_pval = 0.001,
                       rating_method = c("z-scores","model"),
                       residual_type = c("deviance","pearson","anscombe"),
                       outlier_criterion = c("direction","leverage","none"),
                       bonferroni_alpha = 0.01
                      ){
  
  norm_method <- match.arg(norm_method)
  cor_method <- match.arg(cor_method)
  model_type <- match.arg(model_type)
  residual_type <- match.arg(residual_type)
  rating_method <- match.arg(rating_method)
  outlier_criterion <- match.arg(outlier_criterion)
  
  drnr <- methods::new("DysRegNetR",
              counts = counts,
              design = design,
              network = network,
              meta = meta,
              counts_add_one = counts_add_one,
              con_col = con_col,
              norm_method = norm_method,
              size_factors = size_factors,
              cor_filter = cor_filter,
              cor_method = cor_method,
              min_cor = min_cor,
              max_cor_pval = max_cor_pval,
              model_type = model_type,
              cooks_filter = cooks_filter,
              case_for_disp = case_for_disp,
              shared_covariates = shared_covariates,
              gof_filter = gof_filter,
              max_gof_stat = max_gof_stat,
              min_gof_pval = min_gof_pval,
              norm_filter = norm_filter,
              min_norm_stat = min_norm_stat,
              min_norm_pval = min_norm_pval,
              rating_method = rating_method,
              residual_type = residual_type,
              outlier_criterion = outlier_criterion,
              bonferroni_alpha = bonferroni_alpha,
              n_removed_models = 0
             )
  
  if(drnr@counts_add_one){
    drnr@counts <- drnr@counts + 1
  }

  drnr@model_params <- get_model_params(design=drnr@design,
                                        meta_cols = colnames(drnr@meta), 
                                        network_cols = colnames(drnr@network),
                                        con_col = drnr@con_col,
                                        shared_covariates = drnr@shared_covariates)
  
  drnr@meta[[drnr@con_col]] <- as.logical(drnr@meta[[drnr@con_col]])
  drnr@network <- pre_process_network(drnr@network, rownames(drnr@counts))
  
  return(drnr)
}



setGeneric("run", function(object) standardGeneric("run"))
setMethod("run", signature("DysRegNetR"), function(object){
  
  if(length(object@size_factors) == 1 && is.na(object@size_factors)){
    print("estimating size factors")
    object <- normalize(object)
  } else {
    print("size factors already exist, estimation will be skipped")
  }
  
  if(object@cor_filter){
    print("filtering the network based on correlation")
    object <- corFilter(object)
  }
  
  print("building models")
  object <- buildModels(object)
  
  if(object@gof_filter){
    print("filtering models by goodness of fit")
    object <- gofFilter(object)
  }
  
  
  return(object)
})


setGeneric("normalize", 
           function(object, norm_method = object@norm_method) standardGeneric("normalize"),
           signature = "object")
              
setMethod("normalize", signature("DysRegNetR"), function(object, norm_method){
  
  object@size_factors <- get_size_factors(object@counts, norm_method)
  return(object)
})


setGeneric("corFilter", function(object) standardGeneric("corFilter"))
setMethod("corFilter", signature("DysRegNetR"), function(object){
  
  is_control <- !object@meta[[object@con_col]]
  counts <- object@counts[,is_control]
  size_factors <- object@size_factors[is_control]
  
  object@network <- filter_network(network = object@network, 
                               counts = counts, 
                               size_factors = size_factors, 
                               cor_method = object@cor_method, 
                               min_cor = object@min_cor,
                               max_cor_pval = object@max_cor_pval)
  return(object)
})


setGeneric("buildModels", function(object,  silent = T) standardGeneric("buildModels"), signature = "object")
setMethod("buildModels", signature("DysRegNetR"), function(object, silent){
  
  object@models <- get_models(
             network = object@network, 
             counts = object@counts, 
             size_factors = object@size_factors, 
             design = object@design,
             meta = object@meta, 
             con_col = object@con_col,
             model_params = object@model_params,
             cor_filter = object@cor_filter,
             model_type = object@model_type,
             cooks_filter = object@cooks_filter,
             case_for_disp = object@case_for_disp,
             silent = silent
             )
  
  return(object)
  
})


setGeneric("gofFilter", function(object) standardGeneric("gofFilter"), signature = "object")
setMethod("gofFilter", signature("DysRegNetR"), function(object){
  
  n_models <- length(object@models)
  if(n_models==0){
    stop("no models built yet")
  }
  
  
  object@models <- gof_filter(object@models, max_statistic = object@max_gof_stat, min_pval = object@min_gof_pval)
  object@n_removed_models <- n_models - length(object@models)
  
  return(object)
  
})


setGeneric("getResiduals", 
           function(object, 
                    residual_type = object@residual_type,
                    condition = c("all", "case", "control"),
                    regulations = NULL,
                    norm_filter = object@norm_filter,
                    min_norm_stat = object@min_norm_stat,
                    min_norm_pval = object@min_norm_pval,
                    ...
                    ) standardGeneric("getResiduals"),
           signature = "object")
setMethod("getResiduals",signature("DysRegNetR"), function(object, residual_type, condition, regulations, norm_filter, min_norm_stat, min_norm_pval, ...){
  

  condition <- match.arg(condition)
  
  if(norm_filter){
    if(condition=="case"){
      stop("Shapiro-Wilk test of normality works only, if condition is 'all' or 'control'")
    }
  }
  
  if(length(object@models)==0){
    stop("no models found, call run() first")
  }

  residuals <- get_residuals(models = object@models,
                counts = object@counts,
                size_factors = object@size_factors,
                samples = object@meta[, object@con_col],
                model_params = object@model_params,
                meta = object@meta,
                residual_type = residual_type,
                condition = condition,
                regulations = regulations
                )
  
  if(norm_filter){
    if(condition=="control"){
      tests <- normality_test(residuals)
    } else{
      tests <- normality_test(residuals[,!object@meta[[object@con_col]]])
    }
    residuals <- residuals[ !is.na(tests[,"statistic"]) & tests[,"statistic"] >= min_norm_stat & tests[,"pval"] >= min_norm_pval,,drop=F]
    
  }
  
  return(residuals)
  
})
      

setGeneric("getZScores", 
           function(object, 
                    residual_type = object@residual_type,
                    condition = c("all", "case", "control"),
                    regulations = names(object@models), 
                    ...
           ) standardGeneric("getZScores"),
           signature = "object")
setMethod("getZScores",signature("DysRegNetR"), function(object, residual_type, condition, regulations, ...){
  
  condition = match.arg(condition)
  resid_condition <- "all"
  
  if(condition == "control"){
    resid_condition <- "control"
  }
  
  residuals <- getResiduals(object, residual_type, resid_condition, regulations, ...)
  
  get_z_scores(residuals = residuals, samples = object@meta[, object@con_col], condition)
  
})

setGeneric("getPCA", 
           function(object, condition = c("all", "case", "control"), normalized = T, ...) standardGeneric("getPCA"), 
           signature = "object")
setMethod("getPCA", signature("DysRegNetR"), function(object, condition, normalized, ...){
  
  condition <- match.arg(condition)
  
  if(normalized && length(object@size_factors) == 1 && is.na(object@size_factors)){
    stop("calculate size_factors first (normalize())")
  }
  
  samples <- object@meta[[object@con_col]]
  
  if(condition == "case"){
    counts <- object@counts[,samples]
    size_factors <- object@size_factors[samples]
  } else if(condition == "control"){
    counts <- object@counts[,!samples]
    size_factors <- object@size_factors[!samples]
  } else{
    counts <- object@counts
    size_factors <- object@size_factors
  }
  

  get_pca(counts = counts, 
          size_factors = size_factors, 
          normalized = normalized)
})



setGeneric("getResults", 
           function(object, 
                    rating_method = object@rating_method,
                    residual_type = object@residual_type,
                    bonferroni_alpha = object@bonferroni_alpha,
                    criterion = object@outlier_criterion,
                    regulations = names(object@models),
                    ...
           ) standardGeneric("getResults"),
           signature = "object")
setMethod("getResults",signature("DysRegNetR"), function(object, rating_method, residual_type, bonferroni_alpha, criterion, regulations, ...){
  
  rating_method <- match.arg(rating_method, choices = c("z-scores","model"))
  
  if(rating_method == "z-scores"){
    result <- getZScores(object = object, residual_type = residual_type, condition = "case", regulations = regulations, ...)
  }
  
  else{
    result <- getResiduals(object = object, residual_type = residual_type, condition = "case", regulations = regulations, ...)
  }

  
  transformed_tf_counts <- NULL
  
  if(criterion=="leverage"){
    
    transformed_tf_counts <- get_transformed_tf_counts(regulations = regulations, 
                                                       model_params = object@model_params, 
                                                       counts = object@counts[,colnames(result)],
                                                       size_factors = object@size_factors[colnames(result)])
  }
  
  if(rating_method == "z-scores"){
      result <- filter_z_scores(z_scores = result,
                              models = object@models,
                              model_params  = object@model_params,
                              bonferroni_alpha = bonferroni_alpha,
                              transformed_tf_counts = transformed_tf_counts,
                              criterion = criterion)
  }
  else{
    result <- filter_residuals(residuals = result,
                               models = object@models,
                               counts = object@counts[,colnames(result)],
                               size_factors = object@size_factors[colnames(result)],
                               meta = object@meta[colnames(result),],
                               model_params = object@model_params,
                               bonferroni_alpha = bonferroni_alpha,
                               transformed_tf_counts = transformed_tf_counts,
                               criterion = criterion)
  }
  

  
  return(result)
})



setMethod("plot", signature(x="DysRegNetR"), 
          function(x, y, type = c("count_scatter", "pca", "qq", "norm_test"), 
                   pca_data = NULL, 
                   residuals = NULL, 
                   residual_types = c("deviance", "pearson", "anscombe", "raw"),
                   ...){
            
  type = match.arg(type)
  
  if(type == "count_scatter"){
    return(count_scatter(count_matrix=x@counts, 
                         size_factors=x@size_factors,
                         samples = x@meta[,x@con_col],
                         model_params = x@model_params,
                         meta = x@meta,
                         models = x@models,
                          ...)
           )
    

  } else if (type == "pca"){
    
    if(is.null(pca_data)){
      pca_data <- getPCA(x, ...)
    } else{
      print("using existing pca_data")
    }
    
    return(plot_pca(pca_data = pca_data, samples = x@meta[,x@con_col], ...))
    
    
  } else if (type == "qq"){
    
    residual_types <- match.arg(residual_types, several.ok = T)
    
    if(is.null(residuals)){
      residuals <- lapply(residual_types, function(residual_type){
        getResiduals(x, residual_type = residual_type, norm_filter = F, ...)
      })
      names(residuals) <- residual_types
      
    } else{
      print("using existing residuals")
    }
    
    return(plot_qq(residuals = residuals, samples = x@meta[,x@con_col], residual_types = residual_types, ...))
    
  } else if (type == "norm_test"){
    
    residual_types <- match.arg(residual_types, several.ok = T)
    
    residuals <- lapply(residual_types, function(residual_type){
      getResiduals(x, residual_type = residual_type, norm_filter = F, ...)
    })
    names(residuals) <- residual_types
    
    return(norm_test_plot(residuals = residuals, samples = x@meta[,x@con_col], residual_types = residual_types, ...))
  }
})



setMethod("show", signature("DysRegNetR"), function(object){
            
            # samples
            total_samples <- nrow(object@meta)
            case_samples <- sum(object@meta[[object@con_col]])
            control_samples <- total_samples - case_samples
            
            # regulations
            model_params <- object@model_params
          
            regulations <- nrow(object@network)
            tfs <- length(unique(object@network[,model_params$tf_col]))
            tgs <- length(unique(object@network[,model_params$tg_col]))
            
            if(!is.null(object@network$use)){
              filtered_network <- object@network[object@network$use,]
              if(length(filtered_network)==0){
                filtered_regulations <- 0
                filtered_tfs <- 0
                filtered_tgs <- 0
              }
              else{
                filtered_regulations <- nrow(filtered_network)
                filtered_tfs <- length(unique(filtered_network[,model_params$tf_col]))
                filtered_tgs <- length(unique(filtered_network[,model_params$tg_col]))
              }
              regulation_report <- c("Regulations:", regulations,"(",filtered_regulations, "after filtering )\n",
                                     "Unique transcription factors:", tfs, "(",filtered_tfs, "after filtering )\n",
                                     "Unique target genes:", tgs, "(",filtered_tgs, "after filtering )\n\n")
              
            } 
            
            else{
              regulation_report <- c("Regulations:", regulations, "\n",
                                     "Unique transcription factors:", tfs, "\n",
                                     "Unique target genes:", tgs, "\n\n")
            }
            
            # models 
            if(length(object@models)==0){
              model_report <- "No models built yet\n\n"
            } else {
              total_models <- length(object@models)
              null_models <- sum(sapply(object@models,is.null))
              warn_models <- sum(sapply(object@models, function(model){!is.null(model) && !is.null(model$th.warn)}))
              nb_models <- sum(sapply(object@models, function(model){!is.null(model) && !is.null(model$theta)}))
              p_models <- total_models - (null_models + nb_models)
              n_removed_models <- object@n_removed_models
              
              directions <- sapply(object@models, function(model){
                if(is.null(model)){
                  return(0)
                }
                
                coef <- coefficients(model)[model_params$tf_coef]
                if(is.na(coef)){
                  coef <- coefficients(model)[paste0(model_params$con_col,"FALSE:",model_params$tf_coef)]
                }
                if(is.na(coef)){
                  stop("can not find tf coefficient")
                }
                return(sign(coef))
              })
              
              model_report <- c(
                                "Total models:", total_models,"\n",
                                "Models failed building:", null_models,"\n",
                                "Models with warnings:", warn_models,"\n",
                                "Models removed by GoF filtering:", n_removed_models ,"\n",
                                "Negative-binomial models:", nb_models,"\n",
                                "Poisson models:", p_models,"\n",
                                "Activating regulations:", sum(directions>0),"\n",
                                "Repressing regulations:", sum(directions<0),"\n\n"
                                )

              
              if(object@cooks_filter){
                
                outlier_table <- table(sapply(object@models, function(model){
                  if(!is.null(model$n_outlier)){
                    return(model$n_outlier)
                  }
                  else{
                    return(-1)
                  }
                }))
                
                outlier_table <- paste(capture.output(show(outlier_table)), collapse = "\n")
                model_report <- c(model_report, "Table of removed outliers during training:\n")
                model_report <- c(model_report, outlier_table)
                model_report <- c(model_report, "\n")

              }
            }
            
            
            cat(is(object)[[1]], "object\n\n",
                "Total samples:", total_samples, "\n",
                "Case samples:", case_samples, "\n",
                "Control samples:", control_samples, "\n\n",
                regulation_report,
                model_report)
          })
