# cjRank: Conjoint Attribute Feature Ranking Methods

## Introduction

This package helps researchers to compute attribute feature rankings,
nested marginal means for conjoint experiments, and within-rank marginal
means as described by Dill, Howlett and Müller-Crepon (AJPS, 2023) and
briefly outlined below.

When using the method and/or package, please cite:

Dill, Janina, Howlett, Marnie, & Müller-Crepon, Carl (2022). At Any
Cost: How Ukrainians Think about Self-Defense Against Russia. *American
Journal of Political Science*, forthcoming.

## The method

In conjoint experiments, respondents may at times overwhelmingly respond
in categorical terms to some attribute features. This occurs when they
(almost) always accept or reject a profile with given attribute level,
irrespective of the values of other attributes. Such a categorical logic
implies a ranking of (un)desirable features so that a strategy
characterized by the most resisted (desired) feature $f_1$ across all
attributes and levels is rejected (accepted) . If $f_1$ characterizes
either none or both profiles in a task, choices are guided by
categorical reactions to the second-ranked feature $f_2$, etc. Such
decision-making does not contradict the assumptions underlying conjoint
experiments and typical Average Marginal Component Effect (AMCE)
estimates ([Hainmüller et al.,
2014](https://doi.org/10.1093/pan/mpt024)). Yet, AMCEs and other average
quantitites such as Marginal Means ([Leeper et al.,
2020](https://doi.org/10.1017/pan.2019.30)) depict such decision-making
inadequately as they average over preference directions and intensities
([Abramson et al., 2022](https://doi.org/10.1111/ajps.12714)) across all
tasks. The latter can be very high for highly-ranked features in
categorical decision-making.

But what is the ranking of attribute features? As we explain in the
article, we answer use a heuristic approach: We start choosing the
first-ranked feature $f_1$ as that with a co-occurrence adjusted
marginal mean closest to either 0 or 1, being the feature with the
greatest predictive power over respondents’ choices. A feature with a
marginal mean predicts respondents’ choices no better than a coin toss
while marginal means of 0 or 1 signal perfect predictive power since
respondents always decide in favor or against profiles with that
attribute feature if given a choice. We then identify the second-ranked
feature $f_2$, but using only strategy pairs in which $f_1$ is either
absent or invariant. For this sub-sample we proceed as before,
estimating *nested marginal means* to delineate $f_2$. Again only
keeping pairs without variation in $f_2$, we proceed in the same manner
until all features are ranked or we run out of statistical power.
Alongside the nested marginal means we can also compute *within-rank*
marginal means which tell us how important features are at a certain
rank. I.e., for the second ranked feature $f_2$, we can compute marginal
means for all attributes across all tasks in which 1) $f_2$ varies but
2) the first ranked feature $f_1$ does not vary.

Standard errors of the feature ranks and ordering can be derived by
bootstrapping clustered at the respondent- or task-level. Naturally, the
subset of the data to compute each consecutive rank shrinks, at times
rapidly. Users should therefore pay attention not to overinterpret the
relative ranking among lowly-ranked attributes.

## Installation

You can directly download and install the cjRank package from GitHub.
Upon installation, the package should automatically install all
necessary R dependencies. If this is not the case, the user may see an
error and have install the following modules manually.

``` r
# Download pspm package
library(devtools)
install_github(repo = "carl-mc/cjRank")
```

## Getting started

Users get started by preparing a “long” dataset of their conjoint
experiment in which each row corresponds to one profile seen by one
respondent. A set of columns should correspond to one categorical
attribute with two or more levels from which one was randomly chosen, as
indicated by the value the variable takes for a given profile.
Co-occurrence of the same attribute level (a.k.a. feature) within a
subset of tasks facilitates the ranking as it leads to larger samples
for distinguishing between lower-ranked features. The dataset must have
one column with a binary outcome variable, denoting the outcome of a
forced choice between exactly two profiles in a task. The ranking method
does not work for non-binary outcomes and is not prepared for tasks with
more than two profiles (though this extension is, in theory, possible).
Each dataset should also have one column with a unique respondent ID and
one column with a unique task ID, nested within respondents. The package
includes a mock dataset (taken as a subset of the data in the original
article) for illustrative purposes.

``` r
# Packages
library(cjRank)
library(ggplot2)

# Load data
data("ukraine")

# Define objects with varible names

## Categorical attribute variables
attr.vars <- paste0("attr", 1:5, "_level")

## Binary outcome variable
outcome <- "choice"

## Unique respondent ID
respid <- "id"

## Unique task x respondent ID
taskid <- "resp_pair"
```

## Setting up a cjRank object

To initialize the computation of a feature ranking, we first set up a
`cjRank` object (working through the `R6` package). All operations will
be made through this object which acts as a data container that ensures
the consistency of our operations. Note that upon initialization, we
must set a `min_obs` parameter, which indicates the minimum number of
observations needed to compute the next-ranked feature. In addition tho
the bootstrapping of standard errors, setting this parameter to a
reasonable high value avoids inferences from small subsets of the data.

``` r
ukr_rank <- cjRank$new(data = ukraine, 
                   outcome.var = outcome,
                   task_id.var = taskid,
                   attr_cat.var = attr.vars,
                   resp_id.var = respid,
                   min_obs = 100,
                   cluster_var = respid)
```

## Derive order and ranking of attribute features

We can first query the ordering of features as long as the data subsets
the computation relies on remains larger than `min_obs`. the
`feature_key` entry of our `cjRank`object retains the original entries
of the categorical variable. Since the ordering vector of attribute
features does not include all features (the two last features of an
attribute receive the same rank and one is exclude from the ordering
vector), we can also query the rank each attribute feature receives.

``` r
order <- ukr_rank$get_order()
ukr_rank$feature_key[order]
```

    ##                   attr5_level.3                   attr1_level.1 
    ## "Russian-controlled government"                "Full integrity" 
    ##                   attr1_level.2                   attr2_level.2 
    ##                  "minus Crimea"                        "24'000"

``` r
ukr_rank$get_ranks()
```

    ## attr1_level.1 attr1_level.2 attr1_level.3 attr2_level.1 attr2_level.2 
    ##             2             3             3             5             4 
    ## attr2_level.3 attr3_level.1 attr3_level.2 attr3_level.3 attr4_level.1 
    ##             5             5             5             5             5 
    ## attr4_level.2 attr4_level.3 attr5_level.1 attr5_level.2 attr5_level.3 
    ##             5             5             5             5             1

## Bootstrapp the ordering

To assess the uncertainty of the ranking, we can bootstrap, i.e.,
resampling respondents or tasks in our data with replacement.

``` r
set.seed(1)
ukr_rank$get_ranks_bs(iter = 100)
```

    ## Warning in self$get_order_bs(iter = iter, respondents = respondents):
    ## Bootstrapping over respondents.

    ##      attribute       feature                   feature_key rank mean median
    ## 1  attr1_level attr1_level.1                Full integrity    2 2.40      2
    ## 2  attr1_level attr1_level.2                  minus Crimea    3 3.94      4
    ## 3  attr1_level attr1_level.3         minus Donbas & Crimea    3 3.37      3
    ## 4  attr2_level attr2_level.1                        12'000    5 5.37      5
    ## 5  attr2_level attr2_level.2                        24'000    4 4.82      5
    ## 6  attr2_level attr2_level.3                         6'000    5 4.52      5
    ## 7  attr3_level attr3_level.1                        12'000    5 5.40      5
    ## 8  attr3_level attr3_level.2                        24'000    5 5.13      5
    ## 9  attr3_level attr3_level.3                         6'000    5 5.49      5
    ## 10 attr4_level attr4_level.1                      Low (5%)    5 5.45      5
    ## 11 attr4_level attr4_level.2                Moderate (10%)    5 5.42      5
    ## 12 attr4_level attr4_level.3                     None (0%)    5 5.03      5
    ## 13 attr5_level attr5_level.1                 Full autonomy    5 5.20      5
    ## 14 attr5_level attr5_level.2         Negotiated neutrality    5 5.28      5
    ## 15 attr5_level attr5_level.3 Russian-controlled government    1 1.33      1
    ##    lower_bound upper_bound
    ## 1            1           4
    ## 2            3           6
    ## 3            1           6
    ## 4            4           6
    ## 5            3           6
    ## 6            2           6
    ## 7            4           6
    ## 8            4           6
    ## 9            5           6
    ## 10           5           6
    ## 11           4           6
    ## 12           3           6
    ## 13           2           6
    ## 14           3           6
    ## 15           1           2

## Nested marginal means

With our ranking of features, we can proceed to producing marginal means
estimates for every rank.

``` r
nested_mm <- ukr_rank$get_nested_margmean(adj_coocurrence = TRUE)
head(nested_mm)
```

    ##   N_full N_est rank  rankvar   attribute       feature           feature_key
    ## 1   1000   650    0 Baseline attr1_level attr1_level.1        Full integrity
    ## 2   1000   650    0 Baseline attr1_level attr1_level.2          minus Crimea
    ## 3   1000   650    0 Baseline attr1_level attr1_level.3 minus Donbas & Crimea
    ## 4   1000   650    0 Baseline attr2_level attr2_level.1                12'000
    ## 5   1000   650    0 Baseline attr2_level attr2_level.2                24'000
    ## 6   1000   650    0 Baseline attr2_level attr2_level.3                 6'000
    ##        coef         se
    ## 1 0.7186147 0.03039793
    ## 2 0.4529148 0.03025753
    ## 3 0.2959184 0.03213631
    ## 4 0.5504587 0.03547077
    ## 5 0.3750000 0.03026231
    ## 6 0.5817308 0.03580560

## F-statistics of differences between nested marginal means

Next, we can estimate whether the marginal means computed for the
various ranks are indeed different from each other using an omnibus test
and it’s F-Statistic (see [Leeper et al.,
2020](https://doi.org/10.1017/pan.2019.30))

``` r
ukr_rank$get_rank_fstats()
```

    ##   rank_i    N           var    fstat         pval
    ## 1      1 1000 attr5_level.3 6.901297 2.690179e-11
    ## 2      2  548 attr1_level.1 4.278365 4.032525e-06
    ## 3      3  292 attr1_level.2 3.412311 3.169684e-04
    ## 4      4  184 attr2_level.2 2.099840 2.712386e-02

## Within-rank marginal means

Lastly, we can estimate within-rank marginal means for every rank. These
answer the question of the effect of a feature – say $f_2$ – on choices
for tasks in which higher ranked feature (i.e., $f_1$) do not appear or
are invariant.

``` r
within_mm <- ukr_rank$get_within_rank_margmean(adj_coocurrence = TRUE)
head(within_mm)
```

    ##   N_full N_est rank       rankvar   attribute       feature
    ## 1    452   286    1 attr5_level.3 attr1_level attr1_level.1
    ## 2    452   286    1 attr5_level.3 attr1_level attr1_level.2
    ## 3    452   286    1 attr5_level.3 attr1_level attr1_level.3
    ## 4    452   276    1 attr5_level.3 attr2_level attr2_level.1
    ## 5    452   276    1 attr5_level.3 attr2_level attr2_level.2
    ## 6    452   276    1 attr5_level.3 attr2_level attr2_level.3
    ##             feature_key      coef         se
    ## 1        Full integrity 0.6893204 0.04663568
    ## 2          minus Crimea 0.4500000 0.05452622
    ## 3 minus Donbas & Crimea 0.3253012 0.05278393
    ## 4                12'000 0.5600000 0.05135576
    ## 5                24'000 0.4090909 0.05098955
    ## 6                 6'000 0.5227273 0.05498856

## Feedback, comments, questions

We are very grateful for any bug reports, feedback, questions, or
contributions to this package. Please report any issues here or write to
c.a.muller-crepon \[at\] lse.ac.uk .
