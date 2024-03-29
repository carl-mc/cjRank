% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cjRank.R
\name{cjRank}
\alias{cjRank}
\title{Attribute feature ranking and nested marginal means for conjoint experiments}
\description{
Attribute feature ranking and nested marginal means for conjoint experiments

Attribute feature ranking and nested marginal means for conjoint experiments
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{outcome}}{Vector with binary outcome data}

\item{\code{task_id}}{Unique identifier by task}

\item{\code{resp_id}}{Unique respondent identifier}

\item{\code{attr_cat}}{Attributes as categorical variables}

\item{\code{cluster_var}}{Variable to cluster standard errors on}

\item{\code{features}}{Binary feature matrix}

\item{\code{feature_key}}{Feature labels from categorical attributes}

\item{\code{variance_mat}}{Matrix with task-level attribute variance}

\item{\code{order}}{Order of features}

\item{\code{rank}}{Rank at which a task is decided}

\item{\code{min_obs}}{Minimum number of observations to continue ranking}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-cjRank-new}{\code{cjRank$new()}}
\item \href{#method-cjRank-get_order}{\code{cjRank$get_order()}}
\item \href{#method-cjRank-get_ranks}{\code{cjRank$get_ranks()}}
\item \href{#method-cjRank-get_order_bs}{\code{cjRank$get_order_bs()}}
\item \href{#method-cjRank-get_ranks_bs}{\code{cjRank$get_ranks_bs()}}
\item \href{#method-cjRank-decision_rank}{\code{cjRank$decision_rank()}}
\item \href{#method-cjRank-get_rank_fstats}{\code{cjRank$get_rank_fstats()}}
\item \href{#method-cjRank-get_nested_margmean}{\code{cjRank$get_nested_margmean()}}
\item \href{#method-cjRank-get_within_rank_margmean}{\code{cjRank$get_within_rank_margmean()}}
\item \href{#method-cjRank-clone}{\code{cjRank$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-new"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-new}{}}}
\subsection{Method \code{new()}}{
Initialize a new cjRank object
Get order of attribute features
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$new(
  data,
  outcome.var,
  task_id.var,
  attr_cat.var,
  resp_id.var = NULL,
  cluster_var = NULL,
  min_obs = 100
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{data.frame} with long-format conjoint data, one row per profile in task. Only two profiles per task.}

\item{\code{outcome.var}}{String indicating the binary outcome variable. No missings allowed.}

\item{\code{task_id.var}}{String indicating an ID-variable that identifies unique tasks, 
each with exactly two profiles, nested within respondents. No missings allowed.}

\item{\code{attr_cat.var}}{Character vector indicating the names of categorical attribute variables, 
each with two or more values. No missings allowed.}

\item{\code{resp_id.var}}{String indicating an ID-variable that identifies unique respondents, 
each with one or multiple tasks. No missings allowed.}

\item{\code{cluster_var}}{String indicating a variable that identifies clusters for the clustering of standard errors. 
No missings allowed.}

\item{\code{min_obs}}{Minimum number of observations to continue ranking and computation of nested marginal means.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_order"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_order}{}}}
\subsection{Method \code{get_order()}}{
Applies the hierarchical ranking algorithm described in Dill, Howlett, and Muller-Crepon (2023) using the 
\code{min_obs} cutoff to stop the derivation of ranks for small samples and lowly ranked features.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_order()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Returns the order of features f_a of attributes a in A. 
The last two remaining features per attribute are equally ranked and only one is included in the order. 
See \code{cjRank$get_rank()} to get the full ranking of all features.
Get ranking of attribute features
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_ranks"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_ranks}{}}}
\subsection{Method \code{get_ranks()}}{
Derives the ranking of all features from the ordering returned by the
 algorithm described in Dill, Howlett, and Muller-Crepon (2023). Equally ranked
  features receive the same rank. 
Features that have not been differentially ranked due to the \code{min_obs} threshold 
receive the largest (i.e., lowest) rank of the ranked features \code{+1}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_ranks(order = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{Defaults to the order computed by cjRank. But users can provide alternative orders.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns the rank of features f_a of attributes a in A. 
Bootstrap ordering
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_order_bs"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_order_bs}{}}}
\subsection{Method \code{get_order_bs()}}{
Applies a non-parametric bootstrap by resampling (with replacement) either respondents or tasks. 
Can be used to derive uncertainty estimates for the ordering in the observed data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_order_bs(iter = 100, respondents = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{iter}}{Number of iterations. Defaults to 100}

\item{\code{respondents}}{Resample respondents (with replacement)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A matrix with \code{iter} number of rows, each containing one set of ordering derived from the bootstrapped data.
if the orderings have different lengths due to \code{min_obs}, remaining cells in the matrix are filled with NAs.

Compute bootstrapped confidence intervals for feature ranking
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_ranks_bs"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_ranks_bs}{}}}
\subsection{Method \code{get_ranks_bs()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_ranks_bs(
  iter = 1000,
  respondents = TRUE,
  lower_bound = 0.025,
  upper_bound = 0.975,
  quantile_type = 1
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{iter}}{How many iterations. Must be integer.}

\item{\code{respondents}}{Logical. Bootstrap over respondents? Defaults to TRUE.}

\item{\code{lower_bound}}{Percentile for lower bound of confidence interval. Defaults to .025.}

\item{\code{upper_bound}}{Percentile of upper bound of confidence interval. Defaults to .975.}

\item{\code{quantile_type}}{Passed to \code{stats::quantile(type = quantile_type)}. 
Type of quantile algorithm used. 
Can make a difference with few iterations.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame with confidence intervals
 of the ranking of each attribute feature. 
Get rank at which a decision for profiles in a task is made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-decision_rank"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-decision_rank}{}}}
\subsection{Method \code{decision_rank()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$decision_rank(order = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{Defaults to the order computed by cjRank. But users can provide alternative orders.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns a vector with the rank at which decisions for a profile has been 
made if the ordering of features is strictly applied.
F-Statistics of contrasts between subsequent ranks
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_rank_fstats"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_rank_fstats}{}}}
\subsection{Method \code{get_rank_fstats()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_rank_fstats(order = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{Defaults to the order computed by cjRank. But users can provide alternative orders.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns a data.frame() with the results of Wald-Test for differences in the 
estimates of marginal means for subsequent ranks.
Compute nested marginal means
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_nested_margmean"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_nested_margmean}{}}}
\subsection{Method \code{get_nested_margmean()}}{
Computed nested marginal means as described in Dill, Howlett, and Muller-Crepon (2023). 
For a given rank, nested marginal means are computed by subsetting the data to those observations for
which decisions are made at numerically higher (less important) ranks, hence the name "nested". 
The first set of estimates at "rank 0" therefore consists of baseline marginal mean estimates computed for the full sample.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_nested_margmean(order = NULL, adj_coocurrence = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{Defaults to the order computed by cjRank. But users can provide alternative orders.}

\item{\code{adj_coocurrence}}{Adjust for co-occurence by dropping tasks without variance on an
 attribute from the sample when computing the respective marginal means. 
Defaults to TRUE to avoid biases.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame() with nested marginal means (coefficients and clustered standard errors) 
organized by rank, attribute and attribute feature. 
\code{N_full} denotes the number of observations/profiles remaining at a a given rank, 
\code{N_est} is the number of observation used to compute a given marginal mean after co-occurence adjustment.
\code{coef} and \code{se} contain the marginal mean and its (clustered) standard error. 
Compute marginal means within ranks
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-get_within_rank_margmean"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-get_within_rank_margmean}{}}}
\subsection{Method \code{get_within_rank_margmean()}}{
Computed marginal means within ranks as described in Dill, Howlett, and Muller-Crepon (2023). 
For a given rank, nested marginal means are computed by subsetting the data to those observations for
which decisions are made at that specific rank. 
The estimates for rank 1 therefore show marginal means for all tasks in which the most important feature varies. 
Rank 2 shows results for tasks in which the most important feature does not vary, but the second-most important does, and so on.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$get_within_rank_margmean(order = NULL, adj_coocurrence = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{Defaults to the order computed by cjRank. But users can provide alternative orders.}

\item{\code{adj_coocurrence}}{Adjust for co-occurence by dropping tasks without variance on an
 attribute from the sample when computing the respective marginal means. 
Defaults to TRUE to avoid biases.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame() with nested marginal means (coefficients and clustered standard errors) 
organized by rank, attribute and attribute feature. 
\code{N_full} denotes the number of observations/profiles decided at a a given rank, 
\code{N_est} is the number of observation used to compute a given marginal mean after co-occurence adjustment.
\code{coef} and \code{se} contain the marginal mean and its (clustered) standard error.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cjRank-clone"></a>}}
\if{latex}{\out{\hypertarget{method-cjRank-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cjRank$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
