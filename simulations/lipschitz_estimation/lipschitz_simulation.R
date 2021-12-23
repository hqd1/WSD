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
arg_pbs = commandArgs(trailingOnly=TRUE)
it = as.numeric(arg_pbs[1]); positive = as.logical(arg_pbs[2]); cat(it, "\t",positive, "\n")
set.seed(it)
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
x4_modwt = waveslim::modwt.3d(x4, "la8", J  = 3)

nchars = nchar(names(x4_modwt))
unique_ft = unique(substring(names(x4_modwt), 1, nchars-1))
names(unique_ft) = unique_ft #so that below, maxima will retain the filter name
nscale = max(substring(names(x4_modwt), nchars, nchars))
load_ft = names(x4_modwt)
for (i in load_ft[-(length(load_ft))]){
  load(paste0("/storage/home/hqd1/NeuroProject/maxima_",i,"_",it,".rdata"))
}

maxima_list = lapply(unique_ft[-length(unique_ft)], function(ft){ #exclude LLL filter
  ft_list_temp = vector("list",nscale)
  for(j in 1:nscale){
    ft_list_temp[[j]] = get(paste0("maxima_", ft,j))
  }
  ft_list_temp
})
if (positive){
  assign(paste0("positive_wavelet_lips_w1_",it),lipschitz_estim(ts = x4, ts_modwt = x4_modwt, maxima_list, nclus = 20, positive =  TRUE,w = 1))
  save(list = paste0("positive_wavelet_lips_w1_",it), file = paste0("/storage/home/hqd1/NeuroProject/positive_wavelet_lips_w1_",it,".rdata"))
}
if (!positive){
  assign(paste0("negative_wavelet_lips_w1_",it),lipschitz_estim(ts = x4, ts_modwt = x4_modwt, maxima_list, nclus = 20, positive =  TRUE,w = 1))
  save(list = paste0("negative_wavelet_lips_w1_",it), file = paste0("/storage/home/hqd1/NeuroProject/negative_wavelet_lips_w1_",it,".rdata"))
}