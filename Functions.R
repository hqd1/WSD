######install libraries#######
if (!require(devtools)){
  install.packages("devtools")
}
#if (!require(oro.nifti)){install.packages("oro.nifti")}
#if (!require(oro.dicom)){install.packages("oro.dicom")}
#if (!require(fslr)){devtools::install_github("muschellij2/fslr")}
#if (!require(cmaker)){devtools::install_github("stnava/cmaker")}
#if (!require(ITKR)){devtools::install_github("stnava/ITKR")}
#if (!require(ANTsR)){devtools::install_github("stnava/ANTsR")}
#if (!require(extrantsr)){devtools::install_github("muschellij2/extrantsr")}
if (!require(AnalyzeFMRI)){install.packages("AnalyzeFMRI")}
#if (!require(WaveletComp)) {install.packages("WaveletComp")}
if (!require(dplyr)) {install.packages("dplyr")}
if (!require(ggplot2)) {install.packages("ggplot2")}
#if (!require(wavelets)) {install.packages("wavelets")}
#if (!require(wavethresh)){install.packages("wavethresh")}
if (!require(waveslim)){install.packages("waveslim")}
library(devtools);
#library(oro.nifti);library(oro.dicom);library(fslr); library(ITKR); library(ANTsR);library(extrantsr);library(AnalyzeFMRI);library(cmaker);library(scales); 
#library(WaveletComp)
#library(wavelets); 
#library(wavethresh); 
library(waveslim)
library(dplyr); library(ggplot2); library(tidyr)
library(RColorBrewer)
jBrewColors = brewer.pal(n = 9, name = "BuPu"); color_contrast = c(jBrewColors[2],jBrewColors[9])
shadesOfGrey <- colorRampPalette(c("grey0", "grey100"))
library(foreach)
library(doParallel) #In addition, foreach supports a parallelizable operator %dopar% from the doParallel package. This allows each iteration through the loop to use different cores or different machines in a cluster. Here, we demonstrate with using all the cores on the current machine

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
# A function that apply roll, yaw, pitch
ryp = function(xyz,alpha = pi/2,beta = pi/2,gamma = pi/2){
  # source: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
  # These matrices produce the desired effect only if they are used to premultiply column vectors, 
  # and (since in general matrix multiplication is not commutative) only if they are applied in the specified order 
  # alpha: roll angle, beta: yaw angle, gamma: pitcch angle
  # xyz: coordinates
  R_z = matrix(c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, 0, 0, 1), byrow = TRUE, nrow = 3)
  R_y = matrix(c(cos(beta), 0, sin(beta), 0, 1, 0, -sin(beta), 0, cos(beta)), byrow = TRUE, nrow = 3)
  R_x = matrix(c(1,0,0,0, cos(gamma), -sin(gamma), 0,sin(gamma),cos(gamma)), byrow = TRUE, nrow = 3)
  #return(R_z %*% R_y %*% R_x %*% xyz)
  return(xyz %*% t(R_x) %*% t(R_y) %*% t(R_z))
}
# A function to align functions, images, surfaces
# In the case of 2 functions, we need to circularly translate/shift one function by an amount h
# f_h(x) = f(x+h)
circular_shift = function(ts, shift){
  shift = shift %% length(ts)
  # sapply(1:length(ts),function(i){
  #   ifelse(((i + shift)/length(ts))%%1 == 0, ts[length(ts)],ts[max((i+shift) %% length(ts),1)]) #don't return 0 index
  #   #ifelse needed bc: cn%%n = 0
  # }) #too slow for vectors of length a million and up, which is needed for circular_shift_3d
  c(ts[(1+shift):length(ts)], ts[1:shift])
}
circular_shift_2d = function(im, shift){
  dim1_shift = t(apply(im, 1, circular_shift, shift[1]))
  dim2_dim1_shift = apply(dim1_shift, 2, circular_shift, shift[2])
  return(array(dim2_dim1_shift, dim = dim(im)))
}
circular_shift_3d = function(vol, shift){
  n1 = dim(vol)[1]; n2 = dim(vol)[2]; n3 = dim(vol)[3] 
  index = expand.grid(1:n1,1:n2,1:n3) %>% 
    tibble()
  colnames(index) = c("x","y","z")
  shift_index = index %>%
    arrange(z,y) %>%
    mutate(shift_x = circular_shift(x,shift[1])) %>%
    arrange(x,z) %>% 
    mutate(shift_y = circular_shift(y,shift[2])) %>%
    arrange(y,x) %>%
    mutate(shift_z = circular_shift(z,shift[3])) %>%
    as.matrix(ncol = 6) #as.matrix() rearrange rows, so must keep x and shift_x in the same matrix
  vol[shift_index[,1:3]] = vol[shift_index[,4:6]]
  return(vol)
}
# calculate correlation between 2d images
cor_high = function(f,g){
  cor(as.vector(f), as.vector(g))
}
# calculate la8 phase shift
la8_phase = function(lh, j){
  #percival 2000
  if (lh == "l" | lh == "L"){
    return((2^j - 1)*3)
  }
  if (lh == "h" | lh == "H"){
    return(((2^j -1)*7 + 1)/2)
  }
  stop("lh must be either 'l' for low or 'h' for high")
}

align_objects = function(o1,o2, max_shift = NA, min_shift = 0, parallel = FALSE){
  if (!is.array(o1) | !is.array(o2)) {
    warning('Objects supplied are not arrays and will be converted into arrays')
    o1 = as.array(o1); o2 = as.array(o2)
  }
  if (any(dim(o1) != dim(o2))) stop('Objects must have the same dimensions')
  if (length(dim(o1)) == 1){
    #1d array can be thought of as discretized 1d function
    n = ifelse(is.na(max_shift),dim(o1), max_shift)
    cor = sapply(min_shift:(n-1), function(shift){ #this is a vector that contains correlation of 2 functions over all possible shift
      cor(circular_shift(o1, shift),o2)
    })
    return(circular_shift(o1, which.max(cor)))
  }
  if (length(dim(o1)) == 2){
    #2d array can be thought of as discretized 2d image
    n1 = dim(o1)[1]; n2 = dim(o1)[2]
    # if (all(!is.na(max_shift))){
    # n1 = max_shift[1]; n2= max_shift[2]
    # }
    if (all(!is.na(max_shift))){
      if (length(max_shift) == 1){max_shift = rep(max_shift, 2)}
      n1 = max_shift[1] + 1; n2= max_shift[2] + 1
    }
    if(length(min_shift) == 1){min_shift = rep(min_shift, 2)}
    shift = expand.grid(min_shift[1]:(n1-1),min_shift[2]:(n2-1))
    cor = sapply(1:dim(shift)[1], function(i){
      c(cor_high(circular_shift_2d(o1, as.numeric(shift[i,])), o2),i)
    })
    best = cor[2,which.max(cor[1,])]
    return(circular_shift_2d(o1, as.numeric(shift[best,])))
  }
  if (length(dim(o1)) == 3){
    #3d array can be thought of as discretized 3d volume
    n1 = dim(o1)[1]; n2 = dim(o1)[2]; n3 = dim(o1)[3]
    # if (!is.na(max_shift)){
    #   n1 = max_shift + 1; n2= max_shift + 1; n3 = max_shift + 1
    # }
    if (all(!is.na(max_shift))){
      if (length(max_shift) == 1){max_shift = rep(max_shift, 3)}
      n1 = max_shift[1] + 1; n2= max_shift[2] + 1; n3 = max_shift[3] + 1
    }
    if(length(min_shift) == 1){min_shift = rep(min_shift, 3)}
    shift = expand.grid(min_shift[1]:(n1-1),min_shift[2]:(n2-1), min_shift[3]:(n3-1))
    if (parallel == FALSE){
      cor = sapply(1:dim(shift)[1], function(i){
        c(cor_high(circular_shift_3d(o1, as.numeric(shift[i,])), o2),i) #as.numeric turns data.frame to vector
      })
    }
    if (parallel == TRUE){
      cor = foreach(i = 1:dim(shift)[1], .combine = cbind) %dopar% {
        c(cor_high(circular_shift_3d(o1, as.numeric(shift[i,])), o2),i)
      }
    }
    best = cor[2,which.max(cor[1,])]
    return(circular_shift_3d(o1, as.numeric(shift[best,])))
  }
}
# a function that finds local maxima
local_max_3d = function(mod3d_ft, w = 1, jump = 1, nclus = 4){
  dim = dim(mod3d_ft)
  mod3d_ft_df = reshape2::melt(mod3d_ft) %>%
    mutate(x = Var1, y = Var2, z = Var3, value = value, .keep =  "none")
  
  registerDoParallel(nclus) 
  maxima_location = foreach(i = 1:floor(dim[1]/jump), .combine = rbind) %:%
    foreach(j = 1:floor(dim[2]/jump), .combine = rbind) %:%
    foreach (k = 1:floor(dim[3]/jump), .combine = rbind) %dopar% {
      mod3d_ft_df %>%
        filter(between(x,jump*i - w, jump*i + w), between(y,jump*j - w, jump*j + w), between(z, jump*k - w, jump*k + w)) %>%
        filter(abs(value) == max(abs(value)), abs(value) != min(abs(value))) %>% #in case max is same as min
        select(x,y,z, value)
    }
  stopImplicitCluster()
  return(maxima_location)
}
# a function that chains the local maxima across scales
chain_max = function(maxima, w = 1, nclus = 4){
  # maxima is a list that contains local maxima of all filters of the same type across different scales, for e.g. HHH1, HHH2, HHH3, HHH4
  # nclus is the number of cores to run parallel
  ft_num = length(maxima)
  registerDoParallel(nclus) 
  forward = vector("list", length = ft_num)
  
  forward[[ft_num]] = maxima[[ft_num]] #search forward from highest to lowest scale for neighbors
  for(j in (ft_num - 1):1) {
    forward[[j]] = (foreach(i = 1:dim(forward[[j + 1]])[1], .combine = rbind) %dopar%
                      {
                        maxima[[j]] %>%
                          filter(between(x, forward[[j + 1]][i,]$x-w, forward[[j + 1]][i,]$x+w),
                                 between(y, forward[[j + 1]][i,]$y-w, forward[[j + 1]][i,]$y+w),
                                 between(z, forward[[j + 1]][i,]$z-w, forward[[j + 1]][i,]$z+w)) %>%
                          distinct() %>%
                          mutate(x = x, y = y, z = z, value = value, scale = j, .keep = "none")
                      }) %>% distinct() %>%
      arrange(x,y,z)
  }
  
  backward = vector("list", length = ft_num) #search backward from lowest to highest scale for neighbors
  backward[[1]] = forward[[1]]
  for(j in 2:ft_num) {
    backward[[j]] = (foreach(i = 1:dim(backward[[j - 1]])[1], .combine = rbind) %dopar%
                       {
                         forward[[j]] %>%
                           filter(between(x, backward[[j - 1]][i,]$x-w, backward[[j - 1]][i,]$x+w),
                                  between(y, backward[[j - 1]][i,]$y-w, backward[[j - 1]][i,]$y+w),
                                  between(z, backward[[j - 1]][i,]$z-w, backward[[j - 1]][i,]$z+w)) %>%
                           distinct() %>%
                           mutate(x = x, y = y, z = z, value = value, scale = j, .keep = "none")
                       }) %>% 
      distinct() %>%
      arrange(x,y,z) #need to arrange in order to cluster the maxima based on row id later
  }
  
  chain_maxima_df = tibble()
  for (i in 1:length(backward)){
    chain_maxima_df = chain_maxima_df %>% rbind(backward[[i]])
  }
  
  chain_initial_group = (foreach(i = 1:dim(chain_maxima_df)[1], .combine = rbind) %dopar%
                           {
                             chain_maxima_df %>%
                               filter(between(x, chain_maxima_df[i,]$x-w, chain_maxima_df[i,]$x+w),
                                      between(y, chain_maxima_df[i,]$y-w, chain_maxima_df[i,]$y+w),
                                      between(z, chain_maxima_df[i,]$z-w, chain_maxima_df[i,]$z+w)) %>%
                               distinct() %>%
                               mutate(x = x, y = y, z = z, value = value, scale = scale, group = i)
                           }) %>% distinct()
  
  distinct_coor = chain_initial_group %>% distinct(x,y,z)
  group_extract = foreach(i = 1:dim(distinct_coor)[1]) %dopar% {
    chain_initial_group %>%
      filter(x == distinct_coor[i,]$x, y == distinct_coor[i,]$y, z == distinct_coor[i,]$z) %>%
      select(group) %>% 
      distinct() %>%
      pull(group) #turn column data frame into vector
  }
  
  #merge groups that overlap/share a common element
  for (i in 1:(length(group_extract)-1)){
    for (j in (i + 1):length(group_extract)) {
      if(any(group_extract[[j]] %in% group_extract[[i]])) {
        group_extract[[i]] = c(group_extract[[i]], group_extract[[j]]) %>%
          unique()
        group_extract[[j]] = -1
      }
    }
  }
  
  group_unique = unique(group_extract)
  
  chain_final_group = (foreach(i = 1:length(group_unique), .combine = rbind) %dopar%
                         {
                           chain_initial_group %>%
                             filter(group %in% group_unique[[i]]) %>%
                             mutate(x = x, y = y, z = z, value = value, scale = scale, group = i) %>%
                             distinct() 
                         }) %>% arrange(x,y,z)
  
  stopImplicitCluster()
  return(chain_final_group)
}
#write a function that do this simultaneously for all wavelet details and smooth
phase_1d = function(modwt, ts){
  #modwt must be a modwt object from waveslim package and the like
  modwt_matrix = matrix(unlist(modwt), ncol=length(modwt), byrow=FALSE) 
  aligned_modwt = apply(modwt_matrix, 2, align_objects, ts)
  return(aligned_modwt)
}

plot_dwt3d = function(ts, ts_modwt, ft, clus){
  shift = sapply(1:3, function(i){
    la8_phase(substring(ft,i,i), as.numeric(substring(ft,4,4)))
  })
  start_time <- Sys.time()
  aligned_ft = align_objects(ts_modwt[ft][[1]], ts, max_shift = shift, min_shift = shift)
  end_time <- Sys.time()
  cat('takes ',end_time - start_time,'seconds to align wavelet coefficients \n')
  
  aligned_ft_df = reshape2::melt(aligned_ft) %>%
    mutate(x = Var1, y = Var2, z = Var3, value = value,.keep =  "none")
  set.seed(1)
  ft_coefs_cluster = kmeans(aligned_ft_df$value, 3)
  aligned_ft_df = aligned_ft_df %>% 
    mutate(group = ft_coefs_cluster$cluster)
  p = plotly::plot_ly(aligned_ft_df[aligned_ft_df$group == clus,],
                      x = ~x,
                      y = ~y,
                      z = ~z,
                      type = "scatter3d",
                      mode = "markers",
                      size = 1)
  return(list(aligned_df = aligned_ft_df, p = p))
  # p = plotly::plot_ly(aligned_ft_df[aligned_ft_df$group == which.max(ft_coefs_cluster$centers),],
  #                     x = ~x,
  #                     y = ~y,
  #                     z = ~z,
  #                     type = "scatter3d",
  #                     mode = "markers",
  #                     size = 1); p
}
