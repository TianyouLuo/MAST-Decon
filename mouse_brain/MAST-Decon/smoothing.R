#### fitBulk step is super fast, while choose_sigma_c step is not ####
#### we start from the dataset that has run fitBulk and choose_sigma_c #####
#### generate distance matrix ####
## input: physical locations, stored in RCTD@spatialRNA@coords, spot-by-2 matrix with rownames
dist.mat = function(location){
  return(as.matrix(dist(location)))
}


## distance mask matrix: 0/1 matrix with distance matrix and r as input ##
distmask.mat = function(distmat, r){
  mask = matrix(NA, nrow = dim(distmat)[1], ncol = dim(distmat)[2])
  mask[distmat <= r] = 1
  mask[distmat > r] = 0
  return(mask)
}
################

#### find a good starting value of radius ####
find.r0 = function(){
  
}
############################



#### generate the weight matrix ######
## location weights ##
get.locw = function(distmat, r){
  locw = 1 - distmat / r
  locw[locw < 0] = 0
  return(locw)
}

## L2 distance matrix
l2.dist <- function(Z){
  ZZ = Z %*% t(Z)
  diag_ZZ = diag(ZZ)
  diag_ZZ_F = matrix(diag_ZZ, nrow = dim(Z)[1], ncol = dim(Z)[1])
  diag_ZZ_F_T = matrix(diag_ZZ, nrow = dim(Z)[1], ncol=dim(Z)[1], byrow = TRUE)
  sq_dist = diag_ZZ_F + diag_ZZ_F_T - 2*ZZ
  dist = sqrt(abs(sq_dist))
  return(dist)
}


## similarity weights ##
get.simw = function(theta, c=1, var_mat = NULL, comp_var = F){
  nspot = dim(theta)[1]
  if (comp_var == F){
    simw = l2.dist(theta)
  }
  reg = c*median(simw)
  simw = exp(-simw/reg)
  return(simw)
}


## final weight matrix ##
get.Wmat = function(distmat, r, theta, c=1, comp_var = F){
  Wmat = get.locw(distmat = distmat, r = r) * 
    get.simw(theta = theta, c=c, comp_var = comp_var)
  return(Wmat)
}

##################



#### pixel-level log-likelihood for one neighboring pixel p(Y_i(d')|theta(d)), we ignore the epsilon part for now ####
## cell_type_mean_norm is a matrix of dimension spot-by-celltype , theta is a vector
## output is a vector with each element correspond to a gene
pixel.pixlogLik = function(Y, theta, UMI, cell_type_mean_norm){
  cell_type_mean_norm = as.matrix(cell_type_mean_norm)
  pois_par = UMI * (cell_type_mean_norm %*% theta)
  return(sum(log(dpois(Y, lambda = pois_par)+1e-320),na.rm = T))
}


#### pixel-level log-likelihood for pixel i ####
### input: weight matrix: ith row is the weight between ith spot and other spots
## distance mask matrix, n_cell_types, UMI vector
## spatialRNA counts: gene-by-spots
## theta matrix: every row is a spot
## cell_type_mean_norm matrix: gene-by-celltype
pixel.logLik = function(dmaski, weight_vec, spatialRNA, thetai, UMIi, cell_type_mean_norm){
  n_cell_type = dim(cell_type_mean_norm)[2]
  #dmaski = dmask[i,]
  #weight_vec = weight[i, ]
  #thetai = theta[i,]
  #UMIi = UMI[i]
  
  # initialize the log-likelihood
  loglik = 0
  for (j in 1:length(dmaski)){
    if (dmaski[j]==1){
      loglik = loglik + weight_vec[j] * pixel.pixlogLik(spatialRNA[,j], theta = thetai, UMI = UMIi, 
                                                        cell_type_mean_norm = cell_type_mean_norm)
    }
  }
  return(unname(loglik))
}
###############



#### first gradient ####
#### pixel-level first gradient of log-likelihood for one neighboring pixel p(Y_i(d')|theta(d)), we ignore the epsilon part for now ####
## cell_type_mean_norm is a matrix of dimension spot-by-celltype, theta is a vector
## output is a vector with each element correspond to a gene
pixel.pixderiv = function(Y, theta, UMI, cell_type_mean_norm){
  cell_type_mean_norm = as.matrix(cell_type_mean_norm)
  pois_par = as.vector(UMI * (cell_type_mean_norm %*% theta))
  deriv = (Y/pois_par - 1) * UMI * cell_type_mean_norm
  return(colSums(deriv))
}

#### pixel-level first gradient for pixel i ####
pixel.deriv = function(dmaski, weight_vec, spatialRNA, thetai, UMIi, cell_type_mean_norm){
  n_cell_type = dim(cell_type_mean_norm)[2]
  #dmaski = dmask[i,]
  #weight_vec = weight[i, ]
  #thetai = theta[i,]
  #UMIi = UMI[i]
  
  # initialize the log-likelihood
  deriv = numeric(n_cell_type)
  for (j in 1:length(dmaski)){
    if (dmaski[j]==1){
      deriv = deriv + weight_vec[j] * pixel.pixderiv(spatialRNA[,j], theta = thetai, UMI = UMIi, 
                                                        cell_type_mean_norm = cell_type_mean_norm)
    }
  }
  return(unname(deriv))
}
###############





#### optimization function ####
## get the current MLE for pixel i ##
get.pixelmle = function(dmaski, weight_vec, spatialRNA, thetai, UMIi, cell_type_mean_norm, tol = 10^-6, trace = 0){
  fit = optimx(
    par = thetai, # initial values for the parameters. 
    fn = function(x, dmaski, weight_vec, spatialRNA, UMIi, cell_type_mean_norm){
      pixel.logLik(dmaski = dmaski, weight_vec = weight_vec, spatialRNA = spatialRNA, 
                   thetai = x, UMIi = UMIi, cell_type_mean_norm = cell_type_mean_norm)
      }, # log likelihood
    #gr = function(x, dmaski, weight_vec, spatialRNA, UMIi, cell_type_mean_norm){
    #  pixel.deriv(dmaski = dmaski, weight_vec = weight_vec, spatialRNA = spatialRNA, 
    #              thetai = x, UMIi = UMIi, cell_type_mean_norm = cell_type_mean_norm)
    #  }, # gradient/1st derivative
    method = "bobyqa",
    dmaski = dmaski,
    weight_vec = weight_vec,
    spatialRNA = spatialRNA,
    UMIi = UMIi,
    cell_type_mean_norm = cell_type_mean_norm,
    lower = numeric(length(thetai)),
    #itmax = maxit, # max number of iterations
    control = list(
      trace = trace, # higher number print more detailed output
      maximize = T, # default is to minimize
      abstol= tol
      # parscale = vector of scaling values to correct for scale differencs in parameters.  Should be the same length as the num of parameters
    )
  )
  theta_new = fit[1,1:dim(cell_type_mean_norm)[2]]
  return(as.matrix(unname(theta_new)))
}


######## update theta in each iteration #######
new.iter = function(distmat, spatialRNA, theta, UMI, cell_type_mean_norm, r, c, 
                    tol = 10^-6, trace = 0){
  dmask = distmask.mat(distmat = distmat, r = r)
  weight = get.Wmat(distmat = distmat, r=r, theta=theta, c=c, comp_var = F)
  theta_new = theta
  nspot = length(UMI)
  for (i in 1:nspot){
    print(paste0("spot: ", i))
    dmaski = dmask[i,]
    weight_vec = weight[i, ]
    thetai = theta[i,]
    UMIi = UMI[i]
    theta_new[i,] = get.pixelmle(dmaski = dmaski, weight_vec = weight_vec, 
                                 spatialRNA = spatialRNA, thetai = thetai, UMIi = UMIi, 
                                 cell_type_mean_norm = cell_type_mean_norm, tol = tol, trace = trace)
  }
  return(theta_new)
}


#### Function to keep only top celltypes and celltypes with proportion above certain threshold ####
keep.maxCT = function(x, keep = 4){
  x[rank(-x, ties.method = "random") > keep] = 0
  #x = x/sum(x)
  return(x)
}

#### Function to keep only celltypes with proportion above certain threshold ####
keep.threshold = function(x, threshold = 0.02){
  x = x/sum(x)
  x[x <= threshold] = 0
  x = x/sum(x)
  return(x)
}



##### write iterations ######
train.smooth = function(spatialRNA, theta0, UMI, cell_type_mean_norm, distmat, r0 = 505, trace = 0,
                        s = 1.42, c = 1, iter = 4, radius_seq = NULL, tol = 10^-6,
                        max_CT = NULL, threshold = NULL){
  theta_update = theta0
  if (!is.null(radius_seq)){
    n_iter = length(radius_seq)
    for (ite in 1:n_iter){
      print(paste0("iter: ", ite))
      print(paste0("r: ", radius_seq[ite]))
      theta_update = new.iter(distmat = distmat, spatialRNA = spatialRNA, r = radius_seq[ite],
                              theta = theta_update, UMI = UMI, c = c, tol = 10^-6, trace = trace,
                              cell_type_mean_norm = cell_type_mean_norm)
    }
  } else {
    for (ite in 1:iter){
      print(paste0("iter: ", ite))
      print(paste0("r: ",r0*s^(ite-1)))
      theta_update = new.iter(distmat = distmat, spatialRNA = spatialRNA, r = r0*s^(ite-1),
                              theta = theta_update, UMI = UMI, c = c, tol = 10^-6, trace = trace,
                              cell_type_mean_norm = cell_type_mean_norm)
    }
  }
  if (!is.null(max_CT)){
    theta_update = t(apply(theta_update, 1, keep.maxCT, keep = max_CT))
  }
  if (!is.null(threshold)){
    theta_update = t(apply(theta_update, 1, keep.threshold, threshold = threshold))
  }
  return(theta_update)
}



