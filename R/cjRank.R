############################
# cjRank class
############################

#' Example conjoint data in cjRank
#'
#' @name ukraine
#' @docType data
#' @description A subset of the data used in Dill, Janina, Howlett, Marnie, & 
#' Muller-Crepon, Carl (2022). At Any Cost: How Ukrainians Think about 
#' Self-Defense Against Russia. _American Journal of Political Science_, forthcoming. 
#' @keywords data
NULL


#' @title Attribute feature ranking and nested marginal means for conjoint experiments
#' 
#' @name cjRank
#' @importFrom R6 R6Class
#' @importFrom plyr join
#' @importFrom fixest feols
#' @importFrom stats as.formula
#' @importFrom stats quantile
#' @export
cjRank <- R6::R6Class(classname = "cjRank",
                      public = 
                        list(
                          ## PUBLIC DATA
                          #' @field  outcome
                          #' Vector with binary outcome data
                          outcome = NULL, 
                          
                          #' @field  task_id
                          #' Unique identifier by task
                          task_id = NULL, 
                          
                          #' @field  resp_id
                          #' Unique respondent identifier
                          resp_id = NULL, 
                          
                          #' @field  attr_cat
                          #' Attributes as categorical variables
                          attr_cat = NULL, 
                          
                          #' @field  cluster_var
                          #' Variable to cluster standard errors on
                          cluster_var = NULL, 
                          
                          #' @field  features
                          #' Binary feature matrix
                          features = NULL, 
                          
                          #' @field  feature_key
                          #' Feature labels from categorical attributes
                          feature_key = NULL, 
                          
                          #' @field  variance_mat
                          #' Matrix with task-level attribute variance
                          variance_mat = NULL, 
                          
                          #' @field  order
                          #' Order of features
                          order = NULL, 
                          
                          #' @field  rank
                          #' Rank at which a task is decided
                          rank = NULL, 
                          
                          #' @field  min_obs
                          #' Minimum number of observations to continue ranking
                          min_obs = NULL, 
                          
                          
                          ## PUBLIC METHODS
                          
                          ### Initialize object
                          #' @param data \code{data.frame} with long-format conjoint data, one row per profile in task. Only two profiles per task. 
                          #'
                          #' @param outcome.var String indicating the binary outcome variable. No missings allowed.
                          #' @param task_id.var String indicating an ID-variable that identifies unique tasks, 
                          #' each with exactly two profiles, nested within respondents. No missings allowed.
                          #' @param attr_cat.var Character vector indicating the names of categorical attribute variables, 
                          #' each with two or more values. No missings allowed.
                          #' @param resp_id.var String indicating an ID-variable that identifies unique respondents, 
                          #' each with one or multiple tasks. No missings allowed.
                          #' @param cluster_var String indicating a variable that identifies clusters for the clustering of standard errors. 
                          #' No missings allowed.
                          #' @param min_obs Minimum number of observations to continue ranking and computation of nested marginal means.
                          #'
                          #' @description Initialize a new cjRank object
                          initialize = function(data, 
                                                outcome.var,
                                                task_id.var, 
                                                attr_cat.var, 
                                                resp_id.var = NULL,
                                                cluster_var = NULL,
                                                min_obs = 100){
                            
                            # Checks
                            task.tab <- table(data[,c(task_id.var)])
                            stopifnot(
                              "Missing data" = 
                                nrow(na.omit(data[, c(outcome.var, task_id.var, attr_cat.var, resp_id.var, cluster_var)])) ==
                                nrow(data),
                              "Specified variable missing in data" = 
                                all(c(outcome.var, task_id.var, attr_cat.var, resp_id.var, cluster_var) %in%
                                      colnames(data)),
                              "Less observations in data than min_obs" = min_obs <= nrow(data),
                              "Only binary outcome data" = all(data[,outcome.var] %in% c(0,1)),
                              "No duplicated outcome within task" = nrow(unique(data[,c(task_id.var,outcome.var)])) == nrow(data),
                              "Tasks not nested within respondents" = length(unique(data[,c(task_id.var)])) == 
                                nrow(unique(data[,c(task_id.var, resp_id.var), drop = FALSE])),
                              "Data has tasks with more than 2 profiles" = !any(task.tab > 2 & task.tab != 0),
                              "Data has tasks with less than 2 profiles" = !any(task.tab == 1),
                              TRUE
                            )
                            
                            # Transfer
                            self$outcome <- data[,outcome.var]
                            self$task_id <- data[,task_id.var]
                            self$attr_cat <- data[,attr_cat.var]
                            if(!is.null(resp_id.var)){
                              self$resp_id <- data[,resp_id.var]
                            }
                            self$cluster_var <- data[,cluster_var]
                            self$min_obs <- min_obs
                            
                            # Make features
                            self$features <- do.call(cbind, lapply(attr_cat.var, 
                                                                   function(v){
                                                                     catvar2binvars(self$attr_cat[,v],
                                                                                    prefix = v) 
                                                                   }))
                            self$feature_key <- unlist(lapply(attr_cat.var, function(v){
                              sort(unique(self$attr_cat[,v]))
                            }))
                            names(self$feature_key) <- colnames(self$features)
                            
                            # Make variance matrix
                            self$variance_mat <- do.call(cbind, lapply(attr_cat.var, 
                                                                       function(v){
                                                                         taskvariance(self$attr_cat[,v],
                                                                                      prefix = v,
                                                                                      self$task_id) 
                                                                       }))
                            
                          },
                          
                          #' Get order of attribute features
                          #' 
                          #' @description Applies the hierarchical ranking algorithm described in Dill, Howlett, and Muller-Crepon (2023) using the 
                          #' \code{min_obs} cutoff to stop the derivation of ranks for small samples and lowly ranked features.
                          #' 
                          #' @return Returns the order of features f_a of attributes a in A. 
                          #' The last two remaining features per attribute are equally ranked and only one is included in the order. 
                          #' See \code{cjRank$get_rank()} to get the full ranking of all features.
                          get_order = function(){
                            
                            # Run function
                            order <- get_order(data.bin = self$features,
                                               data.var = self$variance_mat,
                                               outcome = self$outcome,
                                               taskid = self$task_id, 
                                               attr.vars = colnames(self$attr_cat),
                                               min_obs = self$min_obs)
                            
                            # Return
                            self$order <- order
                            return(order)
                          } ,
                          
                          
                          #' Get ranking of attribute features
                          #' 
                          #' @param order Defaults to the order computed by cjRank. But users can provide alternative orders. 
                          #'
                          #' @description Derives the ranking of all features from the ordering returned by the
                          #'  algorithm described in Dill, Howlett, and Muller-Crepon (2023). Equally ranked
                          #'   features receive the same rank. 
                          #' Features that have not been differentially ranked due to the \code{min_obs} threshold 
                          #' receive the largest (i.e., lowest) rank of the ranked features \code{+1}.
                          #' 
                          #' @return Returns the rank of features f_a of attributes a in A. 
                          get_ranks = function(order = NULL){
                            vars <- colnames(self$features)
                            
                            ## Get order
                            if(is.null(order) & is.null(self$order)){
                              order <- self$get_order()
                            } else if(is.null(order)){
                              order <- self$order 
                            }
                            
                            ## Get full feature ranking
                            ranks <- sapply(vars, function(v){
                              x <- which(order == v)
                              if(length(x) == 0){
                                a <- gsub("\\.[0-9]+$", "", v)
                                o.a <- gsub("\\.[0-9]+$", "", order)
                                f.a <- gsub("\\.[0-9]+$", "", vars)
                                if(sum(o.a == a) == sum(f.a == a) - 1){
                                  x <- which(gsub("\\.[0-9]+$", "", order) ==
                                               a)
                                  x <- x[length(x)]
                                }
                              }
                              if(length(x) == 0){
                                x <- length(order) + 1
                              }
                              x
                            })
                            
                            ## Return
                            return(ranks)
                          },
                          
                          
                          #' Bootstrap ordering
                          #'
                          #' @param iter Number of iterations. Defaults to 100
                          #' @param respondents Resample respondents (with replacement)
                          #'
                          #' @description Applies a non-parametric bootstrap by resampling (with replacement) either respondents or tasks. 
                          #' Can be used to derive uncertainty estimates for the ordering in the observed data. 
                          #' 
                          #' @return A matrix with \code{iter} number of rows, each containing one set of ordering derived from the bootstrapped data.
                          #' if the orderings have different lengths due to \code{min_obs}, remaining cells in the matrix are filled with NAs.
                          #' 
                          get_order_bs = function(iter = 100, respondents = TRUE){
                            ## Variable over which to bootstrap
                            if(is.null(self$resp_id)){
                              warning("Respondent ID is missing. Bootstrapping over tasks instead.")
                              cat()
                              sample.var <- unique(self$task_id)
                              full.var <- self$task_id
                            } else {
                              warning("Bootstrapping over respondents.")
                              sample.var <- unique(self$resp_id)
                              full.var <- self$resp_id
                            }
                            
                            ## Bootstrap
                            order_bs_ls <- lapply(seq_len(iter), function(i){
                              ## Sample
                              new.ids <- sample(sample.var, length(sample.var), replace = T)
                              
                              ## Expand to full sample
                              new.df <- plyr::join(data.frame(id = new.ids),
                                                   data.frame(id = full.var,
                                                              obs = seq_along(full.var)), 
                                                   by = "id", match = "all", type = "left")
                              
                              ## Get order
                              get_order(data.bin = self$features[new.df$obs,],
                                        data.var = self$variance_mat[new.df$obs,],
                                        outcome = self$outcome[new.df$obs],
                                        taskid = self$task_id[new.df$obs],
                                        attr.vars = colnames(self$attr_cat),
                                        min_obs = self$min_obs)
                            })
                            
                            ## Fill up
                            num <- max(sapply(order_bs_ls, length))
                            order_bs_ls <- lapply(order_bs_ls, function(x){
                              c(x, rep(NA, num - length(x)))
                            })
                            order_bs <- do.call(rbind, order_bs_ls)
                            
                            ## Return
                            return(order_bs)
                          },
                          
                          #' Compute bootstrapped confidence intervals for feature ranking
                          #'
                          #' @param iter How many iterations. Must be integer.
                          #' @param respondents Logical. Bootstrap over respondents? Defaults to TRUE. 
                          #' @param lower_bound Percentile for lower bound of confidence interval. Defaults to .025.
                          #' @param upper_bound Percentile of upper bound of confidence interval. Defaults to .975. 
                          #' @param quantile_type Passed to \code{stats::quantile(type = quantile_type)}. 
                          #' Type of quantile algorithm used. 
                          #' Can make a difference with few iterations. 
                          #'
                          #' @return A data.frame with confidence intervals
                          #'  of the ranking of each attribute feature. 
                          get_ranks_bs = function(iter = 1000, respondents = TRUE,
                                                 lower_bound = .025, 
                                                 upper_bound = .975,
                                                 quantile_type = 1){
                            
                            ## Bootstrap ordering
                            order.boot <- self$get_order_bs(iter = iter, 
                                                            respondents = respondents)
                            
                            ## Get ranks from orders
                            ranks.boot <- apply(order.boot, 1, function(r){
                              self$get_ranks(order = na.omit(r))
                            })
                            
                            ## Summarize
                            sum.boot <- t(apply(ranks.boot, 1, function(res){
                              c(mean = mean(res),median =  median(res), 
                                lower_bound = stats::quantile(res, 
                                                              probs = c(lower_bound), 
                                                              type = quantile_type),
                                upper_bound = stats::quantile(res, 
                                                              probs = c(upper_bound), 
                                                              type = quantile_type))
                            }))
                            colnames(sum.boot) <- c("mean", "median", "lower_bound", "upper_bound")
                            
                            ## Put into table
                            sum.boot <- data.frame(attribute = gsub("\\.[0-9]+$", "", 
                                                                    rownames(sum.boot)),
                                                   feature = rownames(sum.boot),
                                                   feature_key = self$feature_key[rownames(sum.boot)],
                                                   rank = self$get_ranks(),
                                                   sum.boot)
                            rownames(sum.boot) <- NULL
                            
                            ## Return
                            return(sum.boot)
                          },
                          
                          #' Get rank at which a decision for profiles in a task is made
                          #'
                          #' @param order  Defaults to the order computed by cjRank. But users can provide alternative orders. 
                          #'
                          #' @return Returns a vector with the rank at which decisions for a profile has been 
                          #' made if the ordering of features is strictly applied.
                          decision_rank = function(order = NULL){
                          
                            ## Get order
                            if(is.null(order) & is.null(self$order)){
                              order <- self$get_order()
                            } else if(is.null(order)){
                              order <- self$order 
                            }
                            
                            
                            ## Get data
                            data <- self$features
                            
                            
                            ## Evaluate
                            eval.mat <- NULL
                            drop = c()
                            rank <- rep(NA, nrow(self$features))
                            for(v in order){
                              # Which is the next one to drop
                              drop <- c(drop,
                                        v)
                              var <- paste0(gsub("\\.[0-9]+$", "", drop), ".var")
                              
                              # Drop data
                              loss.obs <- apply(self$features[, drop, drop = F] == 1 &
                                                  self$variance_mat[, var, drop = F] == 1, 1,
                                                any)
                              drop.pairs <- self$task_id[loss.obs]
                              
                              # Save rank of dropped ones
                              these.obs <- self$task_id %in% drop.pairs & is.na(rank)
                              rank[these.obs] <- which(order == v)
                            }
                            
                            # NA as last rank
                            rank[is.na(rank)] <- max(rank, na.rm = T) + 1
                            
                            # Return
                            self$rank <- rank
                            return(rank)
                          },
                          
                          
                          
                          #' F-Statistics of contrasts between subsequent ranks
                          #'
                          #' @param order Defaults to the order computed by cjRank. But users can provide alternative orders. 
                          #'
                          #' @return Returns a data.frame() with the results of Wald-Test for differences in the 
                          #' estimates of marginal means for subsequent ranks.
                          get_rank_fstats = function(order = NULL){
                            
                            ## Get order
                            if(is.null(order) & is.null(self$order)){
                              order <- self$get_order()
                            } else if(is.null(order)){
                              order <- self$order 
                            }
                            
                            
                            # Get data
                            data <- data.frame(outcome = self$outcome,
                                               self$features,
                                               clustervar = self$cluster_var)
                            
                            # Get rank at which obs is decided
                            data$rank <- self$decision_rank()
                            
                            # Make contrasts
                            eval_ranks <- seq_along(order)
                            test.df <- do.call(rbind, lapply(eval_ranks, function(i){
                              # Model 
                              
                              ## All features
                              attr.vars <- colnames(self$features)
                              attr.vars <- attr.vars[!attr.vars %in% 
                                                       order[seq_len(min(c(i, length(order))))]]
                              
                              ## Restricted variables
                              restr.vars <- paste0("I(",
                                                   attr.vars, "*",
                                                   rep(paste0("(rank >", i, ")"), 1), ")")
                              
                              ## Formula
                              f <-  as.formula(paste("outcome ~ ", paste0("I(rank >", i, ") +"), 
                                                     paste(c(order[i], 
                                                             attr.vars, 
                                                             restr.vars),
                                                           collapse = " + ")))
                              
                              # Estimate and suppress collinearity warnings
                              mf <- suppressMessages(fixest::feols(f,
                                                  data = data[data$rank >= i,],
                                                  cluster = "clustervar"))
                              
                              # Test
                              test <- fixest::wald(mf, 
                                                   keep = paste0("rank >"), 
                                                   cluster = "clustervar", 
                                                   print = FALSE)
                              
                              # Return
                              data.frame(rank_i = i, 
                                         N = mf$nobs,
                                         var = order[i],
                                         fstat = test$stat,
                                         pval = test$p
                              )
                            }))
                            return(test.df)
                          } ,
                          
                          #' Compute nested marginal means
                          #' 
                          #' @description Computed nested marginal means as described in Dill, Howlett, and Muller-Crepon (2023). 
                          #' For a given rank, nested marginal means are computed by subsetting the data to those observations for
                          #' which decisions are made at numerically higher (less important) ranks, hence the name "nested". 
                          #' The first set of estimates at "rank 0" therefore consists of baseline marginal mean estimates computed for the full sample. 
                          #'
                          #' @param order Defaults to the order computed by cjRank. But users can provide alternative orders. 
                          #' @param adj_coocurrence Adjust for co-occurence by dropping tasks without variance on an
                          #'  attribute from the sample when computing the respective marginal means. 
                          #' Defaults to TRUE to avoid biases. 
                          #'
                          #' @return A data.frame() with nested marginal means (coefficients and clustered standard errors) 
                          #' organized by rank, attribute and attribute feature. 
                          #' \code{N_full} denotes the number of observations/profiles remaining at a a given rank, 
                          #' \code{N_est} is the number of observation used to compute a given marginal mean after co-occurence adjustment.
                          #' \code{coef} and \code{se} contain the marginal mean and its (clustered) standard error. 
                          get_nested_margmean = function(order = NULL, adj_coocurrence = TRUE){
                            ## Get order
                            if(is.null(order) & is.null(self$order)){
                              order <- self$get_order()
                            } else if(is.null(order)){
                              order <- self$order 
                            }
                            
                            
                            # Get data
                            data <- data.frame(outcome = self$outcome,
                                               self$features,
                                               clustervar = self$cluster_var)
                            
                            # Get rank at which obs is decided
                            data$rank <- self$decision_rank()
                            
                            ## Make contrasts
                            result <- do.call(rbind, lapply(seq_len(length(order)-1), function(i){
                              do.call(rbind, lapply(colnames(self$attr_cat), function(a){
                                # Formula
                                attr <- colnames(self$features)[
                                  gsub("\\.[0-9]+$", "", colnames(self$features)) == a
                                ]
                                attr <- attr[!attr %in% order[seq_len(i-1)]]
                                if(length(attr) <= 1){
                                  return(NULL)
                                }
                                f <- as.formula(paste("outcome", " ~ 0 + ", 
                                                      paste(attr,
                                                            collapse = " + ")))
                                
                                # Co-occurrence adjustment 
                                if(adj_coocurrence){
                                  keep = (self$variance_mat[,paste0(a, ".var")] == 1)
                                } else {
                                  keep = rep(TRUE, nrow(data))
                                }
                                
                                # Estimate
                                m <- fixest::feols(f,
                                                   data = data[keep &
                                                                 data$rank >= i, ],
                                                   cluster = "clustervar")
                                
                                # Extract results
                                data.frame(N_full = sum(data$rank >= i),
                                           N_est = m$nobs,
                                           rank = i - 1, 
                                           rankvar = c("Baseline", order)[i],
                                           attribute = a,
                                           feature = names(m$coefficients),
                                           feature_key = self$feature_key[names(m$coefficients)],
                                           coef = m$coefficients,
                                           se = m$se)
                              }))
                            }))
                            
                            # Delete rownames
                            row.names(result) <- NULL
                            
                            # Return
                            return(result)
                          },
                          
                          
                          
                          
                          #' Compute marginal means within ranks
                          #' 
                          #' @description Computed marginal means within ranks as described in Dill, Howlett, and Muller-Crepon (2023). 
                          #' For a given rank, nested marginal means are computed by subsetting the data to those observations for
                          #' which decisions are made at that specific rank. 
                          #' The estimates for rank 1 therefore show marginal means for all tasks in which the most important feature varies. 
                          #' Rank 2 shows results for tasks in which the most important feature does not vary, but the second-most important does, and so on. 
                          #'
                          #' @param order Defaults to the order computed by cjRank. But users can provide alternative orders. 
                          #' @param adj_coocurrence Adjust for co-occurence by dropping tasks without variance on an
                          #'  attribute from the sample when computing the respective marginal means. 
                          #' Defaults to TRUE to avoid biases. 
                          #'
                          #' @return A data.frame() with nested marginal means (coefficients and clustered standard errors) 
                          #' organized by rank, attribute and attribute feature. 
                          #' \code{N_full} denotes the number of observations/profiles decided at a a given rank, 
                          #' \code{N_est} is the number of observation used to compute a given marginal mean after co-occurence adjustment.
                          #' \code{coef} and \code{se} contain the marginal mean and its (clustered) standard error. 
                          get_within_rank_margmean = function(order = NULL, adj_coocurrence = TRUE){
                            ## Get order
                            if(is.null(order) & is.null(self$order)){
                              order <- self$get_order()
                            } else if(is.null(order)){
                              order <- self$order 
                            }
                            
                            
                            # Get data
                            data <- data.frame(outcome = self$outcome,
                                               self$features,
                                               clustervar = self$cluster_var)
                            
                            # Get rank at which obs is decided
                            data$rank <- self$decision_rank()
                            
                            ## Make contrasts
                            result <- do.call(rbind, lapply(seq_len(length(order)-1), function(i){
                              print(i)
                              do.call(rbind, lapply(colnames(self$attr_cat), function(a){
                                print(a)
                                # Formula
                                attr <- colnames(self$features)[
                                  gsub("\\.[0-9]+$", "", colnames(self$features)) == a
                                ]
                                attr <- attr[!attr %in% order[seq_len(i-1)]]
                                if(length(attr) <= 1){
                                  return(NULL)
                                }
                                f <- as.formula(paste("outcome", " ~ 0 + ", 
                                                      paste(attr,
                                                            collapse = " + ")))
                                
                                # Co-occurrence adjustment 
                                if(adj_coocurrence){
                                  keep = (self$variance_mat[,paste0(a, ".var")] == 1)
                                } else {
                                  keep = rep(TRUE, nrow(data))
                                }
                                
                                # Estimate
                                m <- fixest::feols(f,
                                                   data = data[keep &
                                                                 data$rank == i, ],
                                                   cluster = "clustervar")
                                
                                # Extract results
                                data.frame(N_full = sum(data$rank == i),
                                           N_est = m$nobs,
                                           rank = i, 
                                           rankvar = order[i],
                                           attribute = a,
                                           feature = names(m$coefficients),
                                           feature_key = self$feature_key[names(m$coefficients)],
                                           coef = m$coefficients,
                                           se = m$se)
                              }))
                            }))
                            
                            # Delete rownames
                            row.names(result) <- NULL
                            
                            # Return
                            return(result)
                          }
                          
                        ))





#' Encode categorical variable to series of binary indicators
#'
#' @param catvar Categorical variable vector
#' @param prefix name of variable / prefix to be used for naming of output
#'
#' @return Matrix of binary indicators, one column per category
#' @keywords internal
catvar2binvars <- function(catvar, prefix){
  # Unique values
  uni.vals <- sort(unique(catvar))
  
  # Binary matrix
  binmat <- do.call(cbind, lapply(uni.vals, function(v){
    as.integer(catvar == v)
  }))
  
  # Name it
  colnames(binmat) <- paste0(prefix, ".", seq_along(uni.vals))
  
  # Return
  return(binmat)
}

#' Compute binary variance indicator
#'
#' @param catvar Categorical variable vector
#' @param prefix name of variable / prefix to be used for naming of output
#' @param taskid Categorical task ID vector
#'
#' @return A vector indicating variance of attribute in a task (0,1)
#' @keywords internal
taskvariance <- function(catvar, prefix, taskid){
  # Unique task x value combinations
  uni.comb <- unique(data.frame(var = catvar,
                                taskid = taskid))
  
  
  # Duplicated taskids
  dupl.tasks <- unique(uni.comb$taskid[duplicated(uni.comb$taskid)])
  
  # Result
  res.df <- data.frame(res = as.integer(taskid %in% dupl.tasks))
  colnames(res.df) <- paste0(prefix, ".var")
  
  # Return
  return(res.df)
}


#' Function or compute feature ordering
#'
#' @param data.bin Matrix of binary features
#' @param data.var Matrix showing variance by attribute within tasks
#' @param outcome Outcome vector
#' @param taskid Task ID vector
#' @param attr.vars Character vector of categorical attribute names 
#' @param min_obs Minimum number of observations needed for ranking.
#'
#' @return Ordering
#' @importFrom fixest feols
#' @importFrom stats as.formula
#' @keywords internal
get_order <- function(data.bin,data.var,
                      outcome,
                      taskid, 
                      attr.vars,
                      min_obs){
  
  # Locals to keep track of things
  
  ## Dropped features (those higher ranked)
  drop <- c()
  
  ## Order of features
  order <- c()
  
  ## predictions based on ranked features
  pred <- rep(NA, nrow(data.bin))
  
  ## Tasks dropped
  drop.tasks <- c()
  
  ## Matrix to store evaluation results in
  eval.mat <- NULL
  
  # Transform data to binary feature indicators
  features <- colnames(data.bin)
  
  # Data for analysis
  data.rank <- data.frame(outcome = outcome, 
                          taskid = taskid, 
                          data.bin,
                          data.var)
  
  # Loop over all features (vars) until none left
  while(length(features) > 0 &
        min_obs <= sum(!data.rank$taskid %in% drop.tasks)){
    
    # Estimate marginal means per feature
    model.ls <- lapply(attr.vars, function(a){
      
      # Get attribute features
      these.vars <- features[
        gsub("\\.[0-9]+$", "", features) == a
      ]
      these.vars <- these.vars[!these.vars %in% order]
      if(length(these.vars) <= 1){
        return(NULL)
      }
      
      # Formula
      f <- as.formula(paste("outcome", " ~ 0 + ", 
                            paste(these.vars,
                                  collapse = " + "),
                            ""))
      
      # Subset data to variant pairs only
      this.df <-  data.rank[!data.rank$taskid %in% drop.tasks &
                              data.rank[, paste0(a, ".var")] == 1,]
      
      
      # Estimate
      m <- fixest::feols(f,
                         data = this.df)
    })
    names(model.ls) <- attr.vars
    
    # Which is the next ranked feature to be dropped in next round?
    coefs <- lapply(model.ls, function(m){
      if(is.null(m)){
        c()
      } else {
        m$coefficients
      }
    })
    names(coefs) <- NULL
    coefs <- unlist(coefs)
    
    ## Name of highest ranked
    new.drop <- names(coefs)[which.max(abs(coefs-.5))]
    
    ## Model of highest ranked
    new.drop.m <- which(sapply(attr.vars, function(a){
      gsub("\\.[0-9]+$", "", new.drop) == a
    }))
    
    ## Expand rank to contrasting feature if only 2 are left in that attribute
    if(length(model.ls[[new.drop.m]]$coefficients) == 2){
      new.drop <- c(new.drop, 
                    names(model.ls[[new.drop.m]]$coefficients)[
                      names(model.ls[[new.drop.m]]$coefficients) != new.drop
                    ])
    }
    
    ## Update order and drop vectors
    order <- c(order, new.drop[1])
    drop <- c(drop, 
              new.drop) 
    var <- sapply(drop, function(d){
      colnames(data.var)[which(sapply(colnames(data.var), function(v){
        gsub("\\.[0-9]+$", "", d) == 
          gsub("\\.var$", "", v)
      }))]
    })
    
    # Drop features to analyze in next rounds
    features <- features[!features %in% drop]
    
    # Drop data
    loss.obs <- apply(data.rank[, drop, drop = F] == 1 & 
                        data.rank[, var, drop = F] == 1, 1, 
                      any)
    drop.tasks <- data.rank$taskid[loss.obs]
    
  }
  
  # Return
  return(order)
}




