



# Suggestions for GSClassifier model developers


## About

+ The book `R packages` is a straightaway and useful reference book for R developers. The free-access website of `R packages` is [https://r-pkgs.org/](https://r-pkgs.org/). As a developer of R, if you haven't hear about it, it's strongly recommanded to just read it. Hadley Wickham, the main author of the book, is an active R developer and have led some master works like `ggplot2` and `plyr`.

+ With `GSClassifier` package, it could be easy for users to build a model only with certain gene sets and transcriptomics data. If you are interesting in sharing your model, `GSClassifier` also provides a simple methodology for this vision. In this section, let's see how to achieve it!

First, load the package


```r
library(GSClassifier)
# 载入需要的程辑包：luckyBase
```

## Available models

With `GSClassifier_Data()`, all models supported in the current `GSClassifier` package would showed.


```r
GSClassifier_Data()
# Available data:
# Usage example:
#   ImmuneSubtype.rds 
#   PAD.train_20200110.rds 
#   PAD.train_20220916.rds 
#   PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier")) 
#   ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
```

For more details of `GSClassifier_Data()`, just:

```
?GSClassifier_Data()
```

Set `model=F`, all `.rds` data would be showed:


```r
GSClassifier_Data(model = F)
# Available data:
# Usage example:
#   general-gene-annotation.rds 
#   ImmuneSubtype.rds 
#   PAD.train_20200110.rds 
#   PAD.train_20220916.rds 
#   testData.rds 
#   PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier")) 
#   ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
```



## Components of a GSClassifier model

Currently, a GSClassifier model and related product environments is designed as a `list` object. Let's take `PAD.train_20210110`(also called `PADi`) as an example.


```r
PADi <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier")) 
```


This picture shows the components of `PADi`:

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{./fig/GSClassifier-model} 

}

\caption{Details of a GSClassifier model}(\#fig:GSClassifierModel)
\end{figure}

As shown, a typical `GSClassifier` model is consist of four parts (with different colors in the picture):

+ `1. ens`: 
  + `Repeat`: productive parameters of `GSClassifier` models
  + `Model`: `GSClassifier` models. Here, `PADi` had 20 models from different subs of the training cohorts
+ `2. scaller`: 
  + `Repeat`: productive parameters of the `scaller` model, which was used for `BestCall` calling
  + `Model`: the `scaller` model
+ `3. geneAnnotation`: a data frame containing gene annotation information
+ `4. geneSet`: a list contains several gene sets

Thus, you can assemble your model like:

```
model <- list()

# bootstrap models based on the training cohort
model[['ens']] <- <Your model for subtypes calling>

# Scaller model
model[['scaller']] <- <Your scaller for BestCall calling>

# a data frame contarining gene annotation for IDs convertion
model[['geneAnnotation']] <- <Your gene annotation>

# Your gene sets
model[['geneSet']] <- <Your gene sets>

saveRDS(model, 'your-model.rds')

```

More tutorials for model establishment, please go to [markdown tutorial](https://github.com/huangwb8/GSClassifier/wiki/Model-establishment) or [html tutorial](http://htmlpreview.github.io/?https://raw.githubusercontent.com/wiki/huangwb8/GSClassifier/Model-establishment.html).

## Submit models to luckyModel package 

Considering most users of `GSClassifier` might have no need for lots of models, We divided the model storage feature into an new ensembl package called [**luckyModel**](https://github.com/huangwb8/luckyModel). Don't worry, the usage is very easy!

If you want to summit your model, you should apply for a contributor of `luckyModel` first. Then, just send the model (`.rds`) into the `inst/extdata/<project>` path of `luckyModel`. After audit, your branch would be accepted and available for the users. 

The name of your model must be the format as following:

```
# <project>
GSClassifier

# <creator>_<model>_v<yyyymmdd>:
HWB_PAD_v20211201.rds
```

## Repeatablility of models

For repeatablility, you had better submit a `.zip` or `.tar.gz` file that containing the information of your model. Here are some suggestions:  

+ `<creator>_<model>_v<yyyymmdd>.md`

  + **Destinations**: Why you develop the model

  + **Design**: The evidence for gene sigatures, et al

  + **Data sources**: The data for model training and validating, et al

  + **Applications**: Where to use your model

  + **Limintations**: Limitation or improvement direction of your model

+ `<creator>_<model>_v<yyyymmdd>.R`: The code you used for model training and validating. 

+ `Data-of-<creator>_<model>_v<yyyymmdd>.rds` (Optional): Due to huge size of omics data, it's OK for you not to submit the raw data.

:cupid: Welcome your contributions!

## Gene Annotation

For convenience, we provided a general gene annotation dataset for different genomics:


```r
gga <- readRDS(system.file("extdata", "general-gene-annotation.rds", package = "GSClassifier"))
names(gga)
# [1] "hg38" "hg19" "mm10"
```

I believe they're enough for routine medicine studies.

Here, take a look at `hg38`:


```r
hg38 <- gga$hg38
head(hg38)
#           ENSEMBL       SYMBOL  ENTREZID
# 1 ENSG00000223972      DDX11L1 100287102
# 3 ENSG00000227232       WASH7P      <NA>
# 4 ENSG00000278267    MIR6859-1 102466751
# 5 ENSG00000243485 RP11-34P13.3      <NA>
# 6 ENSG00000284332    MIR1302-2 100302278
# 7 ENSG00000237613      FAM138A    645520
```

With this kind of data, it's simple to customize your own gene annotation (take `PADi` as examples):


```r

tGene <- as.character(unlist(PADi$geneSet))
geneAnnotation <- hg38[hg38$ENSEMBL %in% tGene, ]
dim(geneAnnotation)
# [1] 32  3
```

Have a check:


```r
head(geneAnnotation)
#              ENSEMBL SYMBOL ENTREZID
# 353  ENSG00000171608 PIK3CD     5293
# 1169 ENSG00000134686   PHC2     1912
# 2892 ENSG00000134247 PTGFRN     5738
# 3855 ENSG00000117090 SLAMF1     6504
# 3858 ENSG00000117091   CD48      962
# 4043 ENSG00000198821  CD247      919
```

This `geneAnnotation` could be the `model[['geneAnnotation']]`.

Also, we use a function called `convert` to do gene ID convertion.


```r
luckyBase::convert(c('GAPDH','TP53'), 'SYMBOL', 'ENSEMBL', hg38)
# [1] "ENSG00000111640" "ENSG00000141510"
```

Note: the `luckyBase` package integrates lots of useful tiny functions, you could explore it sometimes.



