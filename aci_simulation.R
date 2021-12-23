######install libraries#######
# if (!require(devtools)){
#   install.packages("devtools")
# }
#if (!require(oro.nifti)){install.packages("oro.nifti")}
#if (!require(oro.dicom)){install.packages("oro.dicom")}
#if (!require(fslr)){devtools::install_github("muschellij2/fslr")}
#if (!require(cmaker)){devtools::install_github("stnava/cmaker")}
#if (!require(ITKR)){devtools::install_github("stnava/ITKR")}
#if (!require(ANTsR)){devtools::install_github("stnava/ANTsR")}
#if (!require(extrantsr)){devtools::install_github("muschellij2/extrantsr")}
#if (!require(AnalyzeFMRI)){install.packages("AnalyzeFMRI")}
#if (!require(WaveletComp)) {install.packages("WaveletComp")}
# if (!require(dplyr)) {install.packages("dplyr", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
# if (!require(tidyr)) {install.packages("tidyr", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
# if (!require(foreach)) {install.packages("foreach", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
# if (!require(doParallel)) {install.packages("doParallel", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
# if (!require(RColorBrewer)) {install.packages("RColorBrewer", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
# #if (!require(ggplot2)) {install.packages("ggplot2")}
# #if (!require(wavelets)) {install.packages("wavelets")}
# #if (!require(wavethresh)){install.packages("wavethresh")}
# if (!require(waveslim)){install.packages("waveslim", lib = "/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5")}
#library(devtools);
#library(oro.nifti);library(oro.dicom);library(fslr); library(ITKR); library(ANTsR);library(extrantsr);library(AnalyzeFMRI);library(cmaker);library(scales); 
#library(WaveletComp)
#library(wavelets); 
#library(wavethresh); 
.libPaths(c(.libPaths(),"/storage/home/h/hqd1/R/x86_64-redhat-linux-gnu-library/3.5"))
library(waveslim)
library(dplyr)
#library(ggplot2); 
library(tidyr)
library(RColorBrewer)
jBrewColors = brewer.pal(n = 9, name = "BuPu"); color_contrast = c(jBrewColors[2],jBrewColors[9])
shadesOfGrey <- colorRampPalette(c("grey0", "grey100"))
library(foreach)
library(doParallel) #In addition, foreach supports a parallelizable operator %dopar% from the doParallel package. This allows each iteration through the loop to use different cores or different machines in a cluster. Here, we demonstrate with using all the cores on the current machine


setwd("/storage/home/h/hqd1/NeuroProject")
source("Functions.R")

registerDoParallel(cores=20)
######
arg_pbs = commandArgs(trailingOnly=TRUE)
it = as.numeric(arg_pbs[1]);id = as.numeric(arg_pbs[2]); ft = toString(arg_pbs[3]); cat(id, "\t",ft, "\n")

#3d array simulation with no head movement#####
set.seed(it)
# load("/storage/home/h/hqd1/NeuroProject/slice_30_hm_t25_array.RData")
# x4 = slice_30_hm_t25_array
x4 = array(rnorm(64*64*64, 0, 0), dim = c(64,64,64))
x4[, , 50] = (as.matrix(expand.grid(1:64, 1:64, 47), ncol = 3) %>%
                mvtnorm::dmvnorm(mean = c(30,30,47), sigma = diag(c(40,40,1))) %>%
                array(dim = c(64,64,1)))*10000
# x4[, , 10] = sample(c(rep(0, 2096), rep(seq(10,10.4, by = 0.1), 500)), 64*64, replace = FALSE) %>%
# array(dim = c(64,64,1))

x4[11:50, 20, 20] = c(seq(0,10,length.out = 20),seq(10,0,length.out = 20))
x4[10, 10, 30] = 10
x4[43, 43, 21:40] = c(seq(0,10, length.out = 12),seq(9.8,0, length.out = 8))
x4 = x4 + array(rnorm(64*64*64, 0, 0.1), dim = c(64,64,64))
# load(paste0("/storage/home/h/hqd1/NeuroProject/x4_despike_thrice.rdata"))
# x4=x4_despike_thrice
# x4 = array(rnorm(64*64*64, 0, 0), dim = c(64,64,64))
# x4[, , 50] = (as.matrix(expand.grid(1:64, 1:64, 47), ncol = 3) %>%
#                 mvtnorm::dmvnorm(mean = c(30,30,47), sigma = diag(c(40,40,1))) %>%
#                 array(dim = c(64,64,1)))*10000
# # x4[, , 10] = sample(c(rep(0, 2096), rep(seq(10,10.4, by = 0.1), 500)), 64*64, replace = FALSE) %>%
# #   array(dim = c(64,64,1))
# 
# x4[11:50, 20, 20] = c(seq(0,10,length.out = 20),seq(10,0,length.out = 20))
# x4[10, 10, 30] = 10
# x4[43, 43, 21:40] = c(seq(0,10, length.out = 10),seq(9.8,0, length.out = 10))
# x4 = x4 + array(rnorm(64*64*64, 0, 1), dim = c(64,64,64))
x4_modwt = waveslim::modwt.3d(x4, "la8", J = 3)
#####denoise by hard thresholding ######
# x4_dwt = waveslim::dwt.3d(x4,"la8", J = 4)
# sigma_hat = median(abs(x4_dwt$HHH1))/0.6745
# #estimate the decay of noise
# nchars = nchar(names(x4_modwt))
# unique_ft = unique(substring(names(x4_modwt), 1, nchars-1))
# nscale = max(substring(names(x4_modwt), nchars, nchars))
# T = 3*sigma_hat #recommended by Mallat 2008 11.3.2
# c = T/2^(-0.65) #-0.65 is found by looking at the decay of mean magnitute of noise across scales
# T_adjusted = c()
# for (i in 1:nscale){
#   T_adjusted[i] = c*(2^i)^(-0.65) 
# }
x4_threshold_modwt = x4_modwt
nchars = nchar(names(x4_modwt))
unique_ft = unique(substring(names(x4_modwt), 1, nchars-1))
names(unique_ft) = unique_ft #so that below, maxima will retain the filter name
nscale = max(substring(names(x4_modwt), nchars, nchars))

for (filter in unique_ft[-length(unique_ft)]) {
  for(i in 1:nscale){
    ft_modwt = x4_modwt[paste0(filter,i)][[1]]
    Th = sqrt(median((ft_modwt - median(ft_modwt))^2))
    x4_threshold_modwt[paste0(filter,i)][[1]][abs(x4_modwt[paste0(filter,i)][[1]]) < 3*Th] = 0
  }
  #x4_threshold_modwt[paste0(ft,1)][[1]][abs(x4_modwt[paste0(ft,1)][[1]]) < T] = 0
}
x4_denoise = waveslim::imodwt.3d(x4_threshold_modwt)
x4_denoise_modwt = waveslim::modwt.3d(x4_denoise, "la8")


shift = sapply(1:3, function(i){
  la8_phase(substring(ft,i,i), as.numeric(substring(ft,4,4)))
})
x4_denoise_modwt_ft_aligned = align_objects(x4_denoise_modwt[ft][[1]], x4, max_shift = shift, min_shift = shift)
#Lipschitz estimates
assign(paste0("maxima_", ft,"_",it),local_max_3d(x4_denoise_modwt_ft_aligned, w = 1, jump = 2, nclus = 20))
# #real data with head motion#####
# load("/storage/home/h/hqd1/NeuroProject/slice_30_hm_t25_array.RData")
# slice_30_modwt = waveslim::modwt.3d(slice_30_hm_t25_array, J = 4, "la8")
# arg_pbs = commandArgs(trailingOnly=TRUE)
# id = as.numeric(arg_pbs[1]); ft = toString(arg_pbs[2]); cat(id, "\t",ft, "\n")
# 
# shift = sapply(1:3, function(i){
#   la8_phase(substring(ft,i,i), as.numeric(substring(ft,4,4)))
# })
# slice_30_modwt_ft_aligned = align_objects(slice_30_modwt[ft][[1]], slice_30_hm_t25_array, max_shift = shift, min_shift = shift)
# ######Lipschitz estimates#######
# assign(paste0("maxima_", ft),local_max_3d(slice_30_modwt_ft_aligned, w = 3, jump = 3, nclus = 20))

save(list = paste0("maxima_", ft,"_",it), file = paste0("/storage/home/hqd1/NeuroProject/maxima_",ft,"_",it,".rdata"))



