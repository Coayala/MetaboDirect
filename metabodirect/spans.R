# Calculate SPANS Score for a Number of Normalization Methods
# Based on the spans functions from the pmartR package (https://github.com/pmartR/pmartR)
# 

modified_spans_procedure <- function(omicsData, norm_fn = c('max', 'mean', 'median', 'minmax', 'sum', 'zscore'), 
                                     subset_fn = c("all", "los", "ppp"), params = NULL, group = NULL, n_iter = 1000, 
                                     sig_thresh = 0.0001, nonsig_thresh = 0.5, 
                                     min_nonsig = 20, min_sig = 20, ...){
  .spans_procedure(omicsData, norm_fn = norm_fn, subset_fn = subset_fn, params = params, group = group, 
                   n_iter = n_iter, sig_thresh = sig_thresh, nonsig_thresh = nonsig_thresh, min_nonsig = min_nonsig, 
                   min_sig = min_sig, ...)
}

.spans_procedure <- function(omicsData, norm_fn = c('max', 'mean', 'median', 'minmax', 'sum', 'zscore'), 
                             subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"), 
                             params = NULL, group = NULL, n_iter = 1000, sig_thresh = 0.0001, nonsig_thresh = 0.5, 
                             min_nonsig = 20, min_sig = 20, max_nonsig = NULL, max_sig = NULL, location_thresh = 0.05, 
                             scale_thresh = 0.05, verbose = TRUE, parallel = TRUE){ 
  
  edata_cname = get_edata_cname(omicsData)
  fdata_cname = get_fdata_cname(omicsData)
  nsamps = attributes(omicsData)$data_info$num_samps
  
  # error checks
  if(!inherits(omicsData, c("pepData", "proData"))) stop("omicsData must be of class 'pepData', or 'proData'")
  
  # params defaults if none specified
  if(is.null(params)){
    params <- list("los" = list(0.05,0.1,0.2,0.3), 
                   "ppp" = list(0.1,0.25,0.50,0.75))
    
    for(name in names(params)){
      if(!(name %in% subset_fn)){
        params[[name]] <- NULL
      }
    }
  }
  
  # simple function to check if an element of params contains all values between 0 and 1
  checkvals <- function(list){
    sapply(list, function(x){
      isTRUE(all(x >= 0) & all(x<=1))
    })
  }
  
  if(!isTRUE(all(names(params) %in% c("los", "ppp")))) stop("params must be a named list with names of the normalization functions to be tested, one or more of: 'los', 'ppp'")
  
  # subset function parameter value checks
  if("los" %in% subset_fn){
    if(is.null(params$los)) stop("params must contain a list element named 'los'")
    if(!all(checkvals(params$los))) stop("each element of 'los' in params must be a numeric value between 0 and 1")
  }
  if("ppp" %in% subset_fn){
    if(is.null(params$ppp)) stop("params must contain a list element named 'ppp'")
    if(!all(checkvals(params$ppp))) stop("each element of 'ppp' in params must be a numeric value between 0 and 1")
  }

  # assign group variable or internally call group_designation() if user specified a grouping column in f_data
  omicsData <- group_designation(omicsData, main_effects = group)
  group_df = attr(omicsData, "group_DF")
  reorder = match(colnames(omicsData$e_data)[-which(colnames(omicsData$e_data) == edata_cname)], as.character(group_df[,fdata_cname]))
  group = group_df[reorder,]$Group
  
  # misc input checks
  if(!inherits(attr(omicsData, "group_DF"), "data.frame")) stop("the omicsData object must have a grouping structure, usually set by calling group_designation() on the object")
  if(any(c(min_nonsig, min_sig) < 1) | any(c(min_nonsig, min_sig) > nrow(omicsData$e_data))) stop("min_nonsig and min_sig must be an integer value no greater than the number of observed biomolecules")
  if(!is.null(max_sig)){
    if(!is.numeric(max_sig) | max_sig <= min_sig) stop('max_sig must be an integer value greater than min_sig')
  }
  if(!is.null(max_nonsig)){
    if(!is.numeric(max_nonsig) | max_nonsig <= min_nonsig) stop('max_nonsig must be an integer value greater than min_nonsig')
  }
  if(!all(subset_fn %in% c("all", "los", "ppp", "complete", "rip", "ppp_rip"))) stop("subset_fn must be a character vector containing more than one of the elements 'all', 'los', 'ppp', 'complete', 'rip', 'ppp_rip'")
  if(is.null(attr(omicsData, "group_DF"))) stop("omicsData object must have a grouping structure set by calling group_designation()")
  
  ### end main error checking ###
  
  # get indices of significant and nonsignificant p-values
  kw_pvals <- kw_rcpp(omicsData$e_data %>% dplyr::select(-edata_cname) %>% as.matrix(), as.character(group))
  
  # initial storage of both vectors
  sig_inds <- (kw_pvals <= sig_thresh & !is.na(kw_pvals))
  nonsig_inds <- (kw_pvals >= nonsig_thresh & !is.na(kw_pvals))
  
  ## these two while loops change the p-value threshold until at least min_sig and min_nonsig indices are selected
  
  # low to high p-values from KW test
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>% unique() %>% sort()
  iter <- 1
  # while we dont have enough indices, add the indices corresponding to the 1:iter lowest p-values
  while(sum(sig_inds) < min_sig){
    sig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }
  
  # same thing but for high p-values
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>% unique() %>% sort(decreasing = TRUE)
  iter <- 1
  while(sum(nonsig_inds) < min_nonsig){
    nonsig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }
  
  # if user set a maximum for the number of significant and nonsignificant molecules and we are over the maximum
  # randomly sample max_sig indices from the indices of significant molecules ...
  if(!is.null(max_sig)){
    if(sum(sig_inds) > max_sig){
      make_false <- sample(which(sig_inds), sum(sig_inds) - max_sig)
      sig_inds[make_false] <- FALSE
    }
  }
  # ... and max_nonsig indices from the indices of non-significant molecules.
  if(!is.null(max_nonsig)){
    if(sum(nonsig_inds) > max_nonsig){
      make_false <- sample(which(nonsig_inds), sum(nonsig_inds) - max_nonsig)
      nonsig_inds[make_false] <- FALSE
    }
  }
  
  # get a vector of n_iter sample sizes for randomly selecting peptides to determine normalization factors
  scaling_factor <- sum(!is.na(omicsData$e_data %>% dplyr::select(-edata_cname)))/100
  select_n <- ceiling(runif(n_iter, nsamps/scaling_factor, 100)*scaling_factor) - nsamps
  
  ### produce a list with all combinations of subset functions, normalization functions, and parameters ###
  all_calls <- list()
  
  # for each normalization/subset method combination...
  for(nf in norm_fn){
    for(sf in subset_fn){
      #...if it is a normalization function that doesn't have parameters specified, just append a list of the normalization and subset function names...
      if(is.null(params[[sf]])){
        all_calls[[length(all_calls)+1]] <- list(norm_fn = nf, subset_fn = sf)
      }
      #...otherwise append the same list with a parameters variable as well
      else{
        for(par in params[[sf]]){
          if(sf == "los") temp_par <- list("los" = par)
          if(sf == "ppp") temp_par <- list("ppp" = par)
          all_calls[[length(all_calls)+1]] <- list(norm_fn = nf, subset_fn = sf, params = temp_par)
        }
      }
    }
  }
  
  if(length(all_calls) < 2) stop("Your input parameters did not result in more than 1 normalization method.")
  
  n_methods <- length(all_calls)
  
  ### STEP 0:  create random distribution ####
  
  # set up parallel backend
  if(parallel){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else foreach::registerDoSEQ()
  
  # get a median significant and non-significant p-value for n_iter iterations
  background_distribution <- foreach::foreach(i = 1:n_iter, .packages = "pmartR",
                                              .export = c('spans_make_distribution', 'data_normalization', 'kw_rcpp')) %dopar% {
    # function defined at bottom of this script
    spans_make_distribution(omicsData, norm_fn, sig_inds, nonsig_inds, select_n[i])
  }

  # # make empirical cdfs based on vectors of n_iter median p-values
  sig_cdf <- sapply(background_distribution, function(el){el[[1]]}) %>% ecdf()
  nonsig_cdf <- sapply(background_distribution, function(el){el[[2]]}) %>% ecdf()
  
  if(verbose) print("Finished creating background distribution, moving to method candidate selection")
  
  #### STEP 1:  
  # determine which methods (subset function + normalization function + parameters combination) will be assessed in step 2.
  # returned list contains information on the method applied and a T/F value for whether it passed to step 2.

  which_spans <- foreach::foreach(i = 1:length(all_calls), .packages = c("pmartR", 'tidyverse'), 
                                  .export = c('data_normalization', 'kw_rcpp')) %dopar% {
    
    el <- all_calls[[i]]
    norm_object <- data_normalization(omicsData$e_data, norm_method = el$norm_fn, norm_subset = el$subset_fn, subset_parameter = el$params)
    
    p_scale <- kw_rcpp(t(as.matrix(norm_object$scale_param)), group = as.character(group))
    
    if(!is.null(norm_object$location_param)){
      p_location <- kw_rcpp(t(as.matrix(norm_object$location_param)), group = as.character(group))
      if(any(c(p_location, p_scale) < c(location_thresh, scale_thresh))){
        res <- list(passfail = FALSE, step1_pvals = c(p_location, p_scale))
      } else {
        res <- list(passfail = TRUE, step1_pvals = c(p_location, p_scale))
      }
    } else if(p_scale < scale_thresh) {
      res <- list(passfail = FALSE, step1_pvals = c(NA, p_scale))
    } else {
      res <- list(passfail = TRUE, step1_pvals = c(NA, p_scale))
    }
    
    res <- c(el, res)
    
    return(res)
    
  }
  
  if(verbose) print("Finished method candidate selection, proceeding to score selected methods.")

  # STEP 2: Score each method that passed step 1 by normalizing the full data and getting median Kruskal-Wallis p-values for significant and nonsignificant peptides

  scores <- foreach::foreach(el = which_spans, .packages = "pmartR", .export = c('data_normalization', 'kw_rcpp')) %dopar% {
    if(el$passfail){
      norm_data <- data_normalization(omicsData$e_data, norm_method = el$norm_fn, norm_subset = el$subset_fn, subset_parameter = el$params)
      abundance_matrix <- norm_data$norm_matrix

      sig_score <- -log10(median(kw_rcpp(abundance_matrix[sig_inds,], group = as.character(group)), na.rm = TRUE))
      non_sig_score <- log10(median(kw_rcpp(abundance_matrix[nonsig_inds,], group = as.character(group)), na.rm = TRUE))

      score <- (sig_cdf(sig_score) + nonsig_cdf(non_sig_score))/2

      return(list(score, sig_cdf(sig_score), nonsig_cdf(non_sig_score)))
    }
    else return(list(NA, NA, NA))
  }

  if(verbose) print("Finished scoring selected methods")

  # create dataframe with selected methods
  spansres_obj <- data.frame("subset_method" = character(n_methods), "normalization_method" = character(n_methods), "SPANS_score" = numeric(n_methods),
                             "parameters" = character(n_methods), "passed_selection" = logical(n_methods),
                             stringsAsFactors = FALSE, check.names = FALSE)

  extra_info <- data.frame("subset_method" = character(n_methods), "normalization_method" = character(n_methods), "parameters" = character(n_methods),
                           "location_p_value" = numeric(n_methods), "scale_p_value" = numeric(n_methods),
                           "F_log_HSmPV" = numeric(n_methods), "F_log_NSmPV" = numeric(n_methods),
                           "SPANS_score" = numeric(n_methods), stringsAsFactors = FALSE, check.names = FALSE)

  # populate the dataframe from which_spans
  for(i in 1:n_methods){
    ss <- which_spans[[i]]$subset_fn
    norm <- which_spans[[i]]$norm_fn
    score <- scores[[i]][[1]]
    params <- which_spans[[i]]$params %>% unlist() %>% as.character() %>% paste(collapse = ";")
    p_loc <- which_spans[[i]]$step1_pvals[1]
    p_scale <- which_spans[[i]]$step1_pvals[2]
    F_HSmPV <- scores[[i]][[2]]
    F_NSmPV <- scores[[i]][[3]]
    pass_fail <- which_spans[[i]]$passfail

    # store into row of df
    spansres_obj[i,] <- list(ss, norm, score, params, pass_fail)
    extra_info[i, ] <- list(ss, norm, params, p_loc, p_scale, F_HSmPV, F_NSmPV, score)
  }

  spansres_obj <- dplyr::arrange(spansres_obj, dplyr::desc(SPANS_score))
  extra_info <- dplyr::arrange(extra_info, dplyr::desc(SPANS_score)) %>% dplyr::select(-SPANS_score)

  attr(spansres_obj, "method_selection_pvals") <- extra_info
  attr(spansres_obj, "group_vector") = group
  attr(spansres_obj, "significant_thresh") = sig_thresh
  attr(spansres_obj, "nonsignificant_thresh") = nonsig_thresh
  attr(spansres_obj, "n_not_significant") = sum(nonsig_inds)
  attr(spansres_obj, "n_significant") = sum(sig_inds)
  attr(spansres_obj, "location_threshold") = location_thresh
  attr(spansres_obj, "scale_thresh") = scale_thresh
  class(spansres_obj) <- c("SPANSRes", "data.frame")

  return(spansres_obj)
  
}

#' Creates the list of median p-values used to make the background distribution used to compute the SPANS score in step 2.
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'lipidData', or 'metabData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param norm_fn a character vector of normalization methods to choose from.  Current options are 'mean', 'median', 'zscore', and 'mad'.
#' @param sig_inds significant peptide indices (row indices) based on a Kruskal-Wallis test on the un-normalized data
#' @param nonsig_inds non-significant peptide indices (row indices) based on a Kruskal-Wallis test on the un-normalized data
#' @param select_n number of peptide by sample indices in the data to randomly select to determine normalization parameters
#' 
#' @return a list with 2 elements. The median of highly significant p-values, and the median of nonsignificant p-values.
#'  These are obtained from a SINGLE Kruskal-Wallis test on data normalized by scale/location factors determined from a randomly selected subset of peptides and normalization method
spans_make_distribution <- function(omicsData, norm_fn, sig_inds, nonsig_inds, select_n){
  edata_cname = get_edata_cname(omicsData)
  group_df = attr(omicsData, "group_DF")
  
  # put the samples in the group vector in the same order as in the columns of e_data
  reorder = match(colnames(omicsData$e_data)[-which(colnames(omicsData$e_data) == edata_cname)], as.character(group_df[,get_fdata_cname(omicsData)]))
  group_vector = as.character(group_df[reorder,]$Group)
  nsamps = attributes(omicsData)$data_info$num_samps
  
  # need a matrix to pass to kw_rcpp
  abundance_matrix <- omicsData$e_data %>% dplyr::select(-edata_cname) %>% as.matrix()
  
  # indices vector that will be used for subsetting
  inds <- NULL
  
  # for each sample, randomly select an index to include for that sample, this ensures each sample gets at least 1 observation
  # remember matrices are stored as vectors, we select elements using a single number
  for(j in 1:nsamps){
    forced_ind <- sample(which(!is.na(abundance_matrix[,j])), 1)*j
    inds <- c(inds, forced_ind)
  }
  
  # randomly assign the rest of the indices
  inds <- c(inds, sample(setdiff(which(!is.na(abundance_matrix)), inds), select_n))
  
  # randomly select a normalization method
  rand_norm <- sample(norm_fn ,1)
  
  # Normalize random distributions using own normalization methods
  
  abundance_matrix <- data_normalization(cbind(omicsData$e_data %>% dplyr::select(edata_cname), replace(abundance_matrix, -inds, NA)), rand_norm)
  # # get normalization parameters from ubsetted matrix.  
  # # normalize_global_matrix is not intended to be used outside this function, it returns the location and (if applicable) scale parameters for normalization.
  # norm_params <- normalize_global_matrix(cbind(omicsData$e_data %>% dplyr::select(edata_cname), replace(abundance_matrix, -inds, NA)),
  #                                        rand_norm, apply_norm = FALSE)
  # 
  # # apply the normalization
  # lapply(1:ncol(abundance_matrix), function(i){
  #   abundance_matrix[,i] <<- abundance_matrix[,i] - norm_params$location[i]
  #   if(!is.null(norm_params$scale)) abundance_matrix[,i] <<- abundance_matrix[,i]/norm_params$scale[i]
  # })
  
  # run Kruskal-Wallis on the normalized dataset
  hs_mpv <- kw_rcpp(abundance_matrix$norm_matrix[sig_inds,], group_vector)
  ns_mpv <- kw_rcpp(abundance_matrix$norm_matrix[nonsig_inds,], group_vector)
  
  # NA values are from peptides with no observations in 1 group
  hs_mpv <- hs_mpv[!is.na(hs_mpv)]
  ns_mpv <- ns_mpv[!is.na(ns_mpv)]
  
  # log transformed median p-values
  hs_mpv <- -log10(median(hs_mpv))
  ns_mpv <- log10(median(ns_mpv))
  
  return(list(hs_mpv, ns_mpv))
}

#' Gets the parameters for the highest ranked methods from spans.
#' 
#' @param SPANSRes_obj an object of the class SPANSRes obtained by calling \code{spans_procedure()}
#' 
#' @return A list of lists, where there are multiple sublists only if there were ties for the top SPANS score.  Each sublist contains named elements for the subset and normalization methods, and the parameters used for the subset method. \cr
#' 
#' @examples 
#' 
#' library(pmartR)
#' library(pmartRdata)
#' 
#' data(pep_object)
#' 
#' # data must be log transformed and grouped
#' myobject <- edata_transform(pep_object, data_scale = "log2")
#' myobject <- group_designation(myobject, main_effects = "Condition")
#' 
#' spans_result <- spans_procedure(myobject)
#' 
#' # a list of the parameters for any normalization procedure with the best SPANS score
#' best_params <- get_spans_params(spans_result)
#' 
#' # extract the arguments from the first list element
#' subset_fn = best_params[[1]]$subset_fn
#' norm_fn = best_params[[1]]$norm_fn
#' params = best_params[[1]]$params
#' 
#' # pass arguments to normalize global
#' norm_object <- normalize_global(omicsData = myobject, subset_fn = subset_fn, norm_fn = norm_fn, params = params)
#' 
#' @export 
get_spans_params <- function(SPANSRes_obj, sort_by_nmols = FALSE){
  
  if(all(is.na(SPANSRes_obj$SPANS_score))) stop("No methods were selected for scoring, there is no 'best' set of parameters to return.")
  
  # get rows that are tied for top score
  best_df <- SPANSRes_obj %>% 
    dplyr::top_n(1, wt = SPANS_score)
  
  if(sort_by_nmols){
    best_df <- best_df %>%
      dplyr::top_n(1, mols_used_in_norm) 
  }
  
  ## populate a list with the subset method, normalization method, and subset parameters.
  params <- vector("list", nrow(best_df))
  
  for(i in 1:nrow(best_df)){
    params[[i]][["subset_fn"]] = best_df[i, "subset_method"]
    params[[i]][["norm_fn"]] = best_df[i, "normalization_method"]
    
    pars_from_df = as.numeric(strsplit(best_df[i,"parameters"], ";")[[1]])
    
    if(params[[i]][["subset_fn"]] == "ppp_rip"){
      params[[i]][["params"]] = list(ppp_rip = list(ppp = pars_from_df[1], rip = pars_from_df[2]))
    }
    else if(params[[i]][["subset_fn"]] == "all"){
      params[[i]][["params"]] = list(NULL)
    }
    else if(params[[i]][["subset_fn"]] %in% c("los", "rip", "ppp")){
      sublist = list()
      sublist[[params[[i]]["subset_fn"]]] <- pars_from_df
      params[[i]][["params"]] <- sublist
    }
  }
  
  return(params)
}

kw_rcpp <- function(mtr, group) {
  .Call('_pmartR_kw_rcpp', PACKAGE = 'pmartR', mtr, group)
}

data_normalization <- function(edata, norm_method, norm_subset = 'all', subset_parameter = NULL){
  
  if(is.null(subset_parameter)){
    subset_parameter <- 1
  } else {
    subset_parameter <- as.numeric(subset_parameter)
  }
  id_list <- edata[, 1]
  
  # Data subsetting
  
  if(norm_subset == 'ppp'){
    min_observations <- as.integer(subset_parameter * ncol(edata))
    observations <- rowSums(!is.na(edata[-1]))
    observations <- observations > min_observations
    selected_peaks = edata$Mass[observations]
  } else if (norm_subset == 'los'){
    ntop_peaks <- as.integer(subset_parameter * nrow(edata))
    selected_peaks <- list()
    for(sample in 1:ncol(edata)){
      peaks <- edata[, sample]
      names(peaks) <- id_list
      ord_peaks <- sort(abs(peaks), decreasing = TRUE)
      selected_peaks <- names(ord_peaks[1:ntop_peaks])
    }
    selected_peaks <- unique(selected_peaks)
  } else {
    selected_peaks <- id_list
  }
  
  subset_edata <- edata[edata$Mass %in% selected_peaks,]
  subset_edata <- subset_edata[-1]
  edata <- edata[-1]
  
  # Calculating normalization parameters and scales
  sample_mean <- apply(subset_edata, 2, mean, na.rm = TRUE)
  sample_min <- apply(subset_edata, 2, min, na.rm = TRUE)
  sample_max <- apply(subset_edata, 2, max, na.rm = TRUE)
  sample_median <- apply(subset_edata, 2, median, na.rm = TRUE)
  sample_std <- apply(subset_edata, 2, sd, na.rm = TRUE)
  sample_sum <- apply(subset_edata, 2, sum, na.rm = TRUE)
  
  edata <- t(edata)
  
  norm_data = list(location_param = NA, scale_param = NA, norm_matrix = NA)
  
  if(norm_method == 'mean'){
    norm_data$location_param <- sample_mean
    norm_data$scale_param <- (sample_max - sample_min)
  } else if(norm_method == 'median'){
    norm_data$location_param <- sample_median
    norm_data$scale_param <- (sample_max - sample_min)
  } else if(norm_method == 'zscore'){
    norm_data$location_param <- sample_mean
    norm_data$scale_param <- sample_std
  } else if(norm_method == 'sum'){
    norm_data$location_param <- NULL
    norm_data$scale_param <- sample_sum
  } else if(norm_method == 'max'){
    norm_data$location_param <- NULL
    norm_data$scale_param <- sample_max
  } else if(norm_method == 'minmax'){
    norm_data$location_param <- sample_min
    norm_data$scale_param <- (sample_max - sample_min)
  }
  
  if(!is.null(norm_data$location_param)){
    norm_data$norm_matrix <- t((edata - norm_data$location_param) / norm_data$scale_param)
  } else{
    norm_data$norm_matrix <- t(edata / norm_data$scale_param)
  }
  
  
  return(norm_data)
}
