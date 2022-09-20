--- 
title: "The Principle of R Package GSClassifier"
author: "Weibin Huang"
date: "2022-09-20"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
csl: briefings-in-bioinformatics.csl
link-citations: yes
geometry:
  - top=1in
  - left=1in
  - right=1in
  - bottom=1in
fontsize: 12pt
linestretch: 1.5
linkcolor: "blue"
description: "The Principle of GSClassifier"
---


<!-- 

# Welcome {-}

<!--

[GSClassifier](https://github.com/huangwb8/GSClassifier) is an R package for modeling and identification of Gene Expression Profiles (GEPs) subtypes. The detail of **GSClassifier** package usage had been demonstrated in [Github WiKi](https://github.com/huangwb8/GSClassifier/wiki). Here, we propose to introduce the principle of GSClassifier, including flowchart, **top scoring pairs (TSP)** algorithm, and batch effect control. 

emoji: https://github.com/rstudio/blogdown/issues/171

-->


## Basic information {-}

**The Priciple of GSClassifier** is a book for users of R package **GSClassifier** who want to know the most details. If you're looking for the PDF edition, you can find it at <a rel="pdf" href="https://github.com/huangwb8/GSClassifier.principle/blob/master/docs/GSClassifier.principle.pdf">here</a>.

+ [**GSClassifier**](https://github.com/huangwb8/GSClassifier) is an R-based comprehensive classification tool for subtypes modeling and personalized calling based on pure transcriptomics. It could be used for precision medicine, such as cancer diagnosis.
+ The inspiration of **GSClassifier** come from [ImmuneSubtypeClassifier](https://github.com/CRI-iAtlas/ImmuneSubtypeClassifier), an R package for classification of PanCancer immune subtypes based on the work of Gibbs et al [@RN160; @RN315].
+ Lots of surprising features in **GSClassifier** as follows: 
  + Optimized for just `one sample`
  + Available for modeling and calling of brand-new `GEPs-based subtypes` in any diseases (cancers)
  + No limitation of the amount of `gene signatures`(≥1) or `subtypes`(≥2)
  + `Normalization insensitive` due to the use of  the individual `gene rank matrix`
  + More ensemble and repeatable modeling process
  + More optimizations in the parallel computing
  + New useful functions as supplements
+ **ATTENTION!** In the future, there might be third-party contributors in `GSClassifier` platform, with some useful models for specific usages. If you use models provided by these people, **you had better know more details as possible**, including **designs, data sources, destinations, training scripts and limitations** of models, expecially those from studies under peer-review.
+ **MORE PROJECTS**:
  + [**The Principle of GSClassifier**](https://huangwb8.github.io/GSClassifier.principle/): A eBook with more details about GSClassifier package
  + [**luckyModel**](https://github.com/huangwb8/luckyModel): Model ensemble for third-party lucky series, such GSClassifier
  

## License {-}

**GSClassifier** is released under the Apache-2.0 license. See [LICENSE](https://github.com/huangwb8/GSClassifier/blob/master/license.txt) for details.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>

The technical documentation, as a whole, is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

The code contained in this book is simultaneously available under the [MIT license](https://opensource.org/licenses/MIT); this means that you are free to use it in your own packages, as long as you cite the source.

## Installation {-}

`RStudio` is one of the best Integrated Development Environments (IDE) in R programming. If you're struggling in R-GUI, it is recommanded to turn to [RStudio](https://www.rstudio.com/).

For installation of **GSClassifier**, please run these commands in an R environment: 

```R
# Install "devtools" package
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install dependencies
if (!requireNamespace("luckyBase", quietly = TRUE))
    devtools::install_github("huangwb8/luckyBase")

# Install the "GSClassifier" package
if (!requireNamespace("GSClassifier", quietly = TRUE))
    devtools::install_github("huangwb8/GSClassifier")
```

In the future, a stable `GSClassifier` version might be sent to [`CRAN`](https://cran.r-project.org/). Still beta.


## [Mirror](https://gitee.com/huangwb8/GSClassifier) {-}

For some special countries or regions, users could also try:

```R
# Install dependencies
install.packages("https://gitee.com/huangwb8/luckyBase/repository/archive/Primary?format=tar.gz", repos=NULL, method="libcurl")

# Install the "GSClassifier" package
install.packages("https://gitee.com/huangwb8/GSClassifier/repository/archive/Primary?format=tar.gz", repos=NULL, method="libcurl")
```

## Change log {-}

### Version 0.1.9 {-}

+ Optimize function verbose
+ Optimize for a routine scenario: one gene set and two subtypes
+ Optimize the strategy of automatic parameters selection for modeling training with R package `caret`
+ Interact with external models from the [luckyModel](https://github.com/huangwb8/luckyModel) package

### Version 0.1.8 {-}

+ Primary public version of `GSClassifier`
+ Apache License, Version 2.0
+ Friendly wiki-based tutorial
+ Platform for developers

## TODO {-}

+ More medical fields included, such as in Pan-cancer uitility
+ Advanced methods (such as artificial intelligence) for enhanced robustness
+ Unsupervised learning for de-novo classification based on intrinsic frames of omics instead of human knowledges
+ Multi-omics exploration and support
+ More friendly characteristics for developers and contributors
+ Web application for newbies of R programing

## Other Projects {-}

You may also be interested in:

* __"[GSClassifier](https://github.com/huangwb8/GSClassifier)"__ A comprehensive classification tool based on pure transcriptomics for precision medicine.

* __"[luckyModel](https://github.com/huangwb8/luckyModel)"__ Model ensemble for third-party lucky series, such GSClassifier.

-->
