# BayesImpute: a Bayesian imputation method for single-cell RNA-seq data

Siqi Chen, Ruiqing Zheng, Luyi Tian, Fang-Xiang Wu, Min Li * 2021-04-18

## Introduction

Single-cell RNA-sequencing (scRNA-seq) data suffer from a large number of zeros. Such dropout events hinder the downstream data analyses. We propose BayesImpute, a statistical algorithm to impute dropouts in scRNA-seq data. BayesImpute first identifies likely dropouts based on expression rate and coefficient of variation of genes within cell subpopulation, and then constructs the posterior distribution and utilizes the posterior mean to impute dropout values. 

## Workflow

![](E:\chensiqi\ＣＳＱ\Bayesmodel\FIG\Workflow.png)

## Overview

- The first step is preprocessing, which includes gene filtering, normalization, and logarithmization;

- The second step is dimensionality reduction and preliminary clustering;

- The third step is the identification of dropout events;

- The last step is the imputation of dropout events.

  

  

## Usage

- import  imputation.py into  bayesmodel.py
- run the  bayesmodel.py

