

# Discussion

In this section, we would discuss some key topics about **GSClassifier**, including **Missing value imputation (MVI)**, **Batch effect**, **hyperparameters**, and so on.

## Packages


```r

# Install "devtools" package
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install dependencies
if (!requireNamespace("luckyBase", quietly = TRUE))
  devtools::install_github("huangwb8/luckyBase")

# Install the "**GSClassifier**" package
if (!requireNamespace("GSClassifier", quietly = TRUE))
  devtools::install_github("huangwb8/GSClassifier")

# Install the "pacman" package
if (!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
  library(pacman)
} else {
  library(pacman)
}

# Load needed packages
packages_needed <- c(
  "readxl",
  "ComplexHeatmap",
  "GSClassifier",
  "rpart",
  "tidyr",
  "reshape2",
  "ggplot2")
for(i in packages_needed){p_load(char=i)}
```

Here is the environment of R programming:


```
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.936 
# [2] LC_CTYPE=Chinese (Simplified)_China.936   
# [3] LC_MONETARY=Chinese (Simplified)_China.936
# [4] LC_NUMERIC=C                              
# [5] LC_TIME=Chinese (Simplified)_China.936    
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
# [1] ggplot2_3.3.6        reshape2_1.4.4       tidyr_1.2.0         
# [4] rpart_4.1.16         GSClassifier_0.1.22  luckyBase_0.1.0     
# [7] ComplexHeatmap_2.4.3 readxl_1.4.0         pacman_0.5.1        
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3     ggsignif_0.6.3       rjson_0.2.21        
#   [4] ellipsis_0.3.2       class_7.3-20         rprojroot_2.0.3     
#   [7] circlize_0.4.15      GlobalOptions_0.1.2  fs_1.5.2            
#  [10] clue_0.3-57          rstudioapi_0.13      ggpubr_0.4.0        
#  [13] listenv_0.8.0        remotes_2.4.2        prodlim_2019.11.13  
#  [16] fansi_1.0.3          lubridate_1.8.0      codetools_0.2-18    
#  [19] splines_4.0.3        doParallel_1.0.17    cachem_1.0.6        
#  [22] knitr_1.30           pkgload_1.2.4        jsonlite_1.8.0      
#  [25] pROC_1.18.0          caret_6.0-92         broom_1.0.0         
#  [28] cluster_2.1.3        png_0.1-7            compiler_4.0.3      
#  [31] backports_1.4.1      assertthat_0.2.1     Matrix_1.2-18       
#  [34] fastmap_1.1.0        cli_3.3.0            htmltools_0.5.2     
#  [37] prettyunits_1.1.1    tools_4.0.3          gtable_0.3.0        
#  [40] glue_1.6.2           dplyr_1.0.9          Rcpp_1.0.8.3        
#  [43] carData_3.0-5        cellranger_1.1.0     vctrs_0.4.1         
#  [46] nlme_3.1-149         iterators_1.0.14     timeDate_3043.102   
#  [49] xfun_0.33            gower_1.0.0          stringr_1.4.0       
#  [52] globals_0.15.1       ps_1.4.0             testthat_3.1.0      
#  [55] lifecycle_1.0.1      devtools_2.4.3       rstatix_0.7.0       
#  [58] future_1.26.1        MASS_7.3-53          scales_1.2.0        
#  [61] ipred_0.9-12         parallel_4.0.3       RColorBrewer_1.1-3  
#  [64] yaml_2.3.5           memoise_2.0.1        stringi_1.7.6       
#  [67] desc_1.4.1           randomForest_4.6-14  foreach_1.5.2       
#  [70] hardhat_1.1.0        pkgbuild_1.3.1       lava_1.6.10         
#  [73] shape_1.4.6          tuneR_1.4.0          rlang_1.0.2         
#  [76] pkgconfig_2.0.3      evaluate_0.15        lattice_0.20-41     
#  [79] purrr_0.3.4          recipes_0.2.0        processx_3.7.0      
#  [82] tidyselect_1.1.2     parallelly_1.32.0    plyr_1.8.7          
#  [85] magrittr_2.0.3       bookdown_0.21        R6_2.5.1            
#  [88] generics_0.1.2       DBI_1.1.3            pillar_1.7.0        
#  [91] withr_2.5.0          survival_3.3-1       abind_1.4-5         
#  [94] nnet_7.3-17          tibble_3.1.7         future.apply_1.9.0  
#  [97] crayon_1.5.1         car_3.1-0            xgboost_1.6.0.1     
# [100] utf8_1.2.2           rmarkdown_2.14       GetoptLong_1.0.5    
# [103] usethis_2.1.3        data.table_1.14.2    callr_3.7.0         
# [106] ModelMetrics_1.2.2.2 digest_0.6.29        stats4_4.0.3        
# [109] signal_0.7-7         munsell_0.5.0        sessioninfo_1.2.2
```

## Missing value imputation (MVI)

<!--
Missing Value Imputation (MVI); 
+ Missing Value type: Missing Completely at Random (MCAR), Missing at Random (MAR), and Not Missing at Random (NMAR).
+ The main reason of missing value in the scenario of **GSClassifier**?
+ Why is missing value imputation needed in **GSClassifier** modeling?
+ Routine algorithms of missing value imputation for DNA array (mRNA expression)?
+ Missing value imputation in subtype calling?
+ The flowchart of quantile algorithm?
+ zero == Cutout; https://arxiv.org/abs/1708.04552; In this paper, we show that the simple regularization technique of randomly masking out square regions of input during training, which we call cutout, can be used to improve the robustness and overall performance of convolutional neural networks. Not only is this method extremely easy to implement, but we also demonstrate that it can be used in conjunction with existing forms of data augmentation and other regularizers to further improve model performance. We evaluate this method by applying it to current state-of-the-art architectures on the CIFAR-10, CIFAR-100, and SVHN datasets, yielding new state-of-the-art results of 2.56%, 15.20%, and 1.30% test error respectively.
+ Experiments: leave genes out of validation cohorts and then test the performance of **GSClassifier**-PADi.
-->

Due to reasons like low expression/weak signal, contamination of microarray surfaces, inappropriate manual operations, insufficient resolution or systematic errors during the laboratory process [@RN387; @RN389; @RN382], **missing value** in high-imput genetic data is common. Generally, tiny missing value could be just dealed with case deletion, while the biological discovery might be damaged when the missing rate tops 15% [@RN392; @RN386]. Currently, lots of methods, including statistic-based or machine learning-based methods (Figure \@ref(fig:mvi01)), had been developed for **missing value imputation (MVI)** [@RN386]. Wang et al [@RN384] categorized MVI methods into simple (zeros or average),biology knowledge-, global learning-, local learning-, hybrid-based methods. In order to satisfy the working conditions of xgboost [@xgboost] functions (`xgb.train`, `xgboost`, and `xgb.cv`) in GSClassifer, the missing value in expression matrix must be deleted or imputation.

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{./fig/mvi-01} 

}

\caption{Missing value imputation methods reviewed by Hasan et al.}(\#fig:mvi01)
\end{figure}

In **PAD** project, several strategies were applied to reduce the impact of missing values as possible. First, both **PIAM** and **PIDG** in **PAD** project were curated GEPs that were not be missing in over 80% gastric cancer datasets. Here we showed the actual distribution of missing value across samples in gastric cancer datasets we used.


```r
# Data
testData <- readRDS(
  system.file("extdata", 
              "testData.rds", 
              package = "GSClassifier")
  )
expr_pad <- testData$PanSTAD_expr_part

# Missing value
expr_pad_na <- apply(expr_pad, 2, 
                     function(x) sum(is.na(x))/length(x))
expr_pad_na_df <- data.frame(
  sample = names(expr_pad_na),
  prob = as.numeric(expr_pad_na),
  stringsAsFactors = F
)
```

As shown in Figure \@ref(fig:mvi02), the percentage of all samples in gastric cancer datasets we used is lower than 8%.


```r
# ggplot
ggplot(data = expr_pad_na_df, 
       aes(x = sample, y = prob)) + 
  geom_bar(stat = 'identity', color = mycolor[3]) + 
  scale_y_continuous(labels=scales::percent) + 
  labs(x = 'Samples in gastric cancer cohorts', 
       y = 'Percentage of missing value') + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
```

\begin{figure}

{\centering \includegraphics[width=0.8\linewidth]{Discussion_files/figure-latex/mvi02-1} 

}

\caption{The distribution of missing value across gastric cancer samples.}(\#fig:mvi02)
\end{figure}

Second, we did conduct some MVI strategy to deal with data before model training in **GSClassifier**. Due to the low missing rate of our experimental data, we just **set missing value as zero** during model training and subtype identification in the early version of PADi (**PAD.train.v20200110**). The model seemed to be robust in both the internal cohort and external cohorts, and greatly predicted the response to immune checkpoint inhibitors (ICIs) in advanced gastric cancer.

In the new version of PADi (**PAD.train.v20220917**), we designed the so-called **quantile** algorithm for random MVI during **PADi** model training, which also seemed to work well for PADi model training.

<!--

We supposed that ...
Add a flowchart for quantile algorithm

-->

Here, we demonstrated the principle of **quantile** algorithm in the simulated dataset:


```r
# Simulated data
x <- read_xlsx('./data/simulated-data.xlsx', sheet = 'RNA')
expr0 <- as.matrix(x[,-1])
rownames(expr0) <- as.character(as.matrix(x[,1])); rm(x)

# MVI with Quantile algorithm
expr <- expr0
na.pos <- apply(expr,2,is.one.na)
set.seed(478); seeds <- sample(1:ncol(expr)*10, sum(na.pos), replace = F)
tSample <- names(na.pos)[na.pos]
quantile_vector <- (1:1000)/1000
for(i in 1:length(tSample)){ # i=1
  
  sample.i <- tSample[i]
  expr.i <- expr[, sample.i]
  expr.i.max <- max(expr.i, na.rm = T)
  expr.i.min <- min(expr.i, na.rm = T)
  set.seed(seeds[i]);
  
  # Details of quantile algorithm
  expr.i[is.na(expr.i)] <-
    expr.i.min +
    (expr.i.max-expr.i.min) * sample(quantile_vector,
                                     sum(is.na(expr.i)),
                                     replace = T)
  expr[, sample.i] <- expr.i
}
  

# Report
cat('RNA expression:', '\n')
print(expr0)
cat('\n')
cat('RNA expression without NA value:', '\n')
print(expr)
# RNA expression: 
#       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
# Gene1    0.51    0.52    0.60    0.21    0.30    0.40
# Gene2    0.52    0.54    0.58    0.22    0.31    0.35
# Gene3    0.53    0.60    0.61      NA    0.29    0.30
# Gene4    0.21    0.30    0.40    0.51    0.52    0.60
# Gene5    0.22    0.31    0.35    0.52    0.54    0.58
# Gene6    0.23    0.29    0.30    0.53      NA    0.61
# Gene7    0.10    0.12    0.09    0.11    0.12    0.14
# 
# RNA expression without NA value: 
#       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
# Gene1    0.51    0.52    0.60 0.21000 0.30000    0.40
# Gene2    0.52    0.54    0.58 0.22000 0.31000    0.35
# Gene3    0.53    0.60    0.61 0.43256 0.29000    0.30
# Gene4    0.21    0.30    0.40 0.51000 0.52000    0.60
# Gene5    0.22    0.31    0.35 0.52000 0.54000    0.58
# Gene6    0.23    0.29    0.30 0.53000 0.32622    0.61
# Gene7    0.10    0.12    0.09 0.11000 0.12000    0.14
```

Look at the new matrix via heatmap, where the clustering result is not obviously disturbed after MVI:


```r
Heatmap(t(scale(t(expr))), name = "Z-score", column_title = "After MVI")
```



\begin{center}\includegraphics[width=0.6\linewidth]{Discussion_files/figure-latex/unnamed-chunk-5-1} \end{center}

Due to missing value might damage the integrity of biological information, we explored **how much the number of missing value in one sample impacts subtype identification via PADi**. The steps are as following: (i) we used quantile algorithm to do MVI in the internal validation cohort of gastric cancer; (ii) we randomly masked different proportion of genes as zero expression; (iii) we calculated the relative multi-ROC [@pROC] (masked data vs. MVI data). In **GSClassifier**, we developed a function called **mv_tolerance** to complete the task.

(i) Load the internal validation cohort:


```r
# Internal validation cohort
testData <- readRDS(
  system.file("extdata", "testData.rds", package = "GSClassifier")
  )
expr_pad <- testData$PanSTAD_expr_part
modelInfo <- modelData(
  design = testData$PanSTAD_phenotype_part,
  id.col = "ID",
  variable = c("platform", "PAD_subtype"),
  Prop = 0.7,
  seed = 19871
)
validInform <- modelInfo$Data$Valid
expr_pad_innervalid <- expr_pad[,validInform$ID]
```

(ii) Missing value tolerance analysis:


```r
# Time-consuming
mvt <- mv_tolerance(
  X = expr_pad_innervalid,
  gene.loss = c(2, 4, 6, 8, 10, 12),
  levels = c(1, 2, 3, 4),
  model = "PAD.train_20220916",
  seed = 487,
  verbose = T
)
```



(iii) multi-ROC analysis:


```r
# Data
mvt_auc <- mvt$multiAUC
mvt_auc_df <- data.frame()
for(i in 1:length(mvt_auc)){ # i=1
  df.i <- data.frame(
    x = as.integer(Fastextra(names(mvt_auc)[i], '=', 2)),
    y = as.numeric(mvt_auc[[i]]$auc),
    stringsAsFactors = F
  )
  mvt_auc_df <- rbind(mvt_auc_df, df.i)
}
  
```


```r
ggplot(mvt_auc_df, aes(x,y)) +
  geom_point() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) + 
  stat_smooth(formula = y ~ x,method = 'glm') +
  labs(x = 'No. of missing value', 
       y = 'Relative AUC in multi-ROC analysis') +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
```

\begin{figure}

{\centering \includegraphics[width=0.6\linewidth]{Discussion_files/figure-latex/mvi03-1} 

}

\caption{The association between the number of missing value and subtype identification performance.}(\#fig:mvi03)
\end{figure}

As showed in Figure \@ref(fig:mvi03), there is linear negative correlation between the number of missing value (mising rate ranges from 6.25% to 37.5%) and subtype identification performance of **PADi** model. One of the reasons might be that PIAM/PIDG were small GEPs, so little gene loss might significantly impact the performance of **PADi**. By the way, there is no mising value in PIAM/PIDG of the 'Kim2018' cohort, an external validation cohort for ICIs therapy response prediction via **PADi**. Nonetheless, we still used **zero strategy** during subtype identification of **PADi** if any missing value exist, because randomization might make the result unstable, which is not suitable for clinical decision.

<!--
Nonetheless, we still used **zero strategy** during subtype identification of **PADi** if any missing value exist, because randomization might make the result unstable, which is not suitable for clinical decision. 
-->

In conclusion, zero or quantile strategy could be applied for MVI before **GSClassifier** model training. However, missing value should be avoided as possible in subtype identification for missing value really damage the performance of **GSClassifier**. Nonetheless, due to low-input GEPs used in **PADi** model (No. of Gene=32), it's easy to avoid missing value in clinical practice.

## Batch effect {#batch-effect}

<!--

gene-pair batch effect in Google Scholar

relative expression orderings (REO)

single-sample enrichment analysis 

Batch effect control & gene number

-->


**TSP** was widely applied to control batch effects in transciptomic data [@RN369; @RN367; @RN368; @RN364; @RN363; @RN362; @RN366; @RN365]. Still, we tested whether **TSP** is a robust method for batch effect control in real-world data. As demonstrated in Figure \@ref(fig:be01), the obvious batch effects across gastric cancer datsets were significantly reduced after **TSP** normalization.

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{./fig/bactch-effect-01} 

}

\caption{Batch effects across gastric cancer cohorts. All gene pairs were used because subtype vectors were not specified. Top: Raw expression of all genes across samples. Middle: Raw expression of PIAM and PIDG across samples. Bottom: TSP of PIAM and PIDG across samples.}(\#fig:be01)
\end{figure}

In order to confirmed the association between **gene counts** in modeling and batch effect control via **TSP** normalization, we selected random genes with counts ranging 4, 8, 20, 40, and 80 for TSP matrix establishment. As shown in Figure \@ref(fig:be02), **TSP** normalization works greatly in different gene counts for batch effect control compared with raw expression matrix.

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{./fig/bactch-effect-02} 

}

\caption{Batch effects of random genes across gastric cancer cohorts.  All gene pairs were used because subtype vectors were not specified. Gene counts 4, 8, 20, 40, and 80 were detected. Data of set difference were not available because only one gene set were applied.}(\#fig:be02)
\end{figure}


<!--
Batch effect reduction of TSP could be explained in three parts:

+ **Binned expression**. Regardless of what raw distribution of a numeric vector is, it would be converted into a binned vector with exactly the same distribution under a fixed quantile vector.

+ **Pair/set difference**.
-->

## Hyperparameters

<!--
+ modelData(): Prop, seed

+ fitEnsembleModel():

  - n=100
  - sampSize=0.7
  - sampSeed = 2020
  - na.fill.seed = 443
  - breakVec=c(0, 0.25, 0.5, 0.75, 1.0)
  - ptail=0.8/2
  
  - https://www.youtube.com/watch?v=f3ryHJ05h5k
  
  - https://xgboost.readthedocs.io/en/latest/parameter.html
  - Parameters for Tree Booster
  - params = list(max_depth = 10,
                  eta = 0.5,
                  nrounds = 100,
                  nthread = 10,
                  nfold=5)
                  
  所以，在时间非常宝贵的情况下，你必须应用自 己的判断力，并使用小数据样本来找出合适的调优参数，否则硬盘空间可能不足
  nrounds ：最大迭代次数（最终模型中树的数量）。
  colsample_bytree ：建立树时随机抽取的特征数量，用一个比率表示，默认值为1（使
用100%的特征）。
  min_child_weight ：对树进行提升时使用的最小权重，默认为1。
  eta ：学习率，每棵树在最终解中的贡献，默认为0.3。
  gamma ：在树中新增一个叶子分区时所需的最小减损。
  subsample ：子样本数据占整个观测的比例，默认值为1（100%）。
  max_depth ：单个树的最大深度。
-->


### Number of SubModel {#topicSubmodel}


Test

















