#' Direct estimation of higher-level correlations
#'
#' Direct estimation of higher-level correlations from lower-level measurements without data aggregation.
#'
#' @param Z An nXq matrix of lower-level measurements; n is the sample size; q is the number of lower-level variables.
#' @param A A qXp matrix characterizing the binding relationship between higher- and lower-level variables. A should have 0-1 entries.
#' If the (j,k) entry of A is 1, then the j-th lower-level variable belongs to the k-th higher-level variable.
#' @return A list of two objects:
#'  \itemize{
#'   \item \code{cov.est}: direct covariance estimator.
#'   \item \code{cor.est}: direct correlation estimator.
#' }
#' Note that they may not be positive definite.
#' @export
#' @examples
#' data(binding_mat)
#' data(gene_expression_stage1)
#' direct_est(Z=gene_expression_stage1, A=binding_mat)
direct_est <- function (Z, A){
  Sigma <-  cov(Z) # Estimate covariance matrix of z
  p <- ncol(A)
  C2 <- matrix(rep(0, p*p),p,p)
  n <- nrow(Z)
  for (l in seq(1:p)){
      for  (k in seq(1:p)){
          Il <- which(rowSums(A)==1 & A[,l]==1)
          Ik <- which(rowSums(A)==1 & A[,k]==1)
          if(l==k){
              s <- 0
              for (i in Il){
                  for (j in Il) {
                      if (i !=j){
                          s <- s + abs(Sigma[i,j]);
                      }
                  }
              }
              C2[l,k] <- s/length(Il)/(length(Il)-1)
          }
          else{
              s <- 0
              for (i in Il){
                  for (j  in Ik){
                      s <- s + Sigma[i,j]
                  }
              }
              C2[l,k] <- s/length(Il)/length(Ik);
          }
      }
  }

  # force C2 to be positive definite
  #eigen.C2 <- eigen(C2)
  #eigen.C2.values <- eigen.C2$values
  #eigen.C2.values[eigen.C2.values <= 0] <- min(eigen.C2.values[eigen.C2.values > 0])/2

  #C2.new <- eigen.C2$vectors%*%diag(eigen.C2.values)%*%t(eigen.C2$vectors)

  #estimate_Inv <- clime_AIC_cov(C2, n = n, lambda = lambda) # The inverse of covariance matrix. PPI_Matrix estimated
  #estimate_P <- estimate_Inv$Omega
  #return(list(PPI_P = estimate_P, PPI_C = C2))
  return( list(cov.est=C2, cor.est=cov2cor(C2)) )

}


#' Make a kernel positive-definite.
#'
#' An internal function making a pre-specified kernel positive definite.
#'
#' @param K a pre-specified symmetric kernel matrix.
#' @param thres a pre-specified threshold. Any singular values of K that is below thres will be set to 0. The default is 1e-9.
#' @return A new positive-kernel kernel
KPD <- function(K, thres = 1e-9){

  values <- eigen(K)$values
  vectors <- eigen(K)$vectors
  values.new <- values
  values.new[values.new <= thres] <- (values.new[values.new > thres])[sum(values.new>thres)]/2
  K.new <- vectors%*%diag(values.new)%*%t(vectors)
  return(K.new)

}

#' Estimation of the alpha, beta, and gamma
#'
#' An internal function for generating the optimal shrinkage estimator
#'
#' @param Sigma The direct covariance estimator from direct_est().
#' @param Z An nXq matrix of lower-level measurements; n is the sample size; q is the number of lower-level variables.
#' @param A A qXp matrix characterizing the binding relationship between higher- and lower-level variables. A should have 0-1 entries.
#' If the (j,k) entry of A is 1, then the j-th lower-level variable belongs to the k-th higher-level variable.
#' @param trace Whether only the trace of the matrix is calculated. The default is TRUE
#' @param Vhat Whether a pre-specified value of Vhat is given. The default is NULL
#' @param parallel Whether mclapply() is used. The default is TRUE.
#' @param num_cores If parallel is TRUE, then how many cores are used? The default is 8.
#' @return A list of alpha, beta, and gamma estimates.
beta_thres <- function(Sigma, A, Z, trace = T, Vhat = NULL, parallel=T, num_cores=8){

  #options(bigalgebra.mixed_airthmetic_returns_R_matrix=FALSE)

  Z <- scale(Z, center = T, scale = F)
  p <- dim(Sigma)[1]
  q <- dim(A)[1]
  n <- dim(Z)[1]

  ## convert Z into a list----
  Z.list <- as.list(data.frame(t(Z)))
  Chat <- cov(Z)

  if(trace == F){

    Thetahat <- matrix(0, p^2, p^2)
    for(l1 in 1:p){

      index_l1 <- which(rowSums(A) == 1 & A[,l1] == 1)
      S_l1 <- length(index_l1)

      for(k1 in 1:p){

        index_k1 <- which(rowSums(A) == 1 & A[,k1] == 1)
        S_k1 <- length(index_k1)

        m_l1k1 <- rep(0, q^2)
        for(templ1 in 1:S_l1){
          for(tempk1 in 1:S_k1){
            index_temp = index_l1[templ1] + (index_k1[tempk1] - 1)*q
            m_l1k1[index_temp] <- as.numeric(index_l1[templ1] != index_k1[tempk1])
          }
        }

        r = l1 + (k1-1)*p

        for(l2 in 1:p){

          index_l2 <- which(rowSums(A) == 1 & A[,l2] == 1)
          S_l2 <- length(index_l2)

          for(k2 in 1:p){

            index_k2 <- which(rowSums(A) == 1 & A[,k2] == 1)
            S_k2 <- length(index_k2)

            m_l2k2 <- rep(0, q^2)
            for(templ2 in 1:S_l2){
              for(tempk2 in 1:S_k2){
                index_temp = index_l2[templ2] + (index_k2[tempk2] - 1)*q
                m_l2k2[index_temp] <- as.numeric(index_l2[templ2] != index_k2[tempk2])
              }
            }

            if(l1 == k1 & l2 == k2)
              const = 1/(S_l1*(S_l1-1)*S_l2*(S_l2-1))
            if(l1 == k1 & l2 != k2)
              const = 1/(S_l1*(S_l1-1)*S_l2*S_k2)
            if(l1 != k1 & l2 == k2)
              const = 1/(S_l1*S_k1*S_l2*(S_l2-1))
            if(l1 != k1 & l2 != k2)
              const = 1/(S_l1*S_k1*S_l2*S_k2)

            Thetahat[r,s] = const*as.numeric(t(m_l1k1)%*%Vhat%*%m_l2k2)
          }
        }
      }
    }

    return(Thetahat)
  }

  if(trace == T){


    #diag_Thetahat <- rep(0, p^2) # only consider l1=l2=l, k1=k2=k

    ### can we do mcapply here as well?

    diag_cal <- function(l1){

      index_l1 <- which(rowSums(A) == 1 & A[,l1] == 1)
      S_l1 <- length(index_l1)

      diag_Thetahat <- rep(0, p)

      for(k1 in 1:p){

        index_k1 <- which(rowSums(A) == 1 & A[,k1] == 1)
        S_k1 <- length(index_k1)

       # m_l1k1 <- rep(0, q^2)
        m_l1k1 <- matrix(0, q, q)
        for(templ1 in 1:S_l1){
          for(tempk1 in 1:S_k1){
            #index_temp = index_l1[templ1] + (index_k1[tempk1] - 1)*q
            m_l1k1[index_l1[templ1], index_k1[tempk1]] <- as.numeric(index_l1[templ1] != index_k1[tempk1])
            #m_l1k1[index_temp] <- as.numeric(index_l1[templ1] != index_k1[tempk1])
          }
        }

        #r = l1 + (k1-1)*p
        if(l1 == k1) const = 1/(S_l1*(S_l1 - 1))^2
        if(l1 != k1) const = 1/(S_l1*S_k1)^2

        #diag_Thetahat[r] = const*as.numeric(t(m_l1k1)%*%Vhat%*%m_l1k1)

        grid_indices <- expand.grid(index_l1, index_k1)

        # Create i and j vectors using rep()
        i.ind <- grid_indices[, 1]
        j.ind <- grid_indices[, 2]

        m_l1k1 <- Matrix::sparseMatrix(i = i.ind,
                             j = j.ind,
                             dims=c(q,q),
                             x = as.vector(m_l1k1[index_l1, index_k1]))
        #print(m_l1k1)
        #temp1 <- mclapply(Z.list, function(x){  (t(m_l1k1)%*%kronecker.prod(x, x))^2  }, mc.cores = 2)

        if(parallel==T) temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric((t(x)%*%m_l1k1%*%x)^2) }, mc.cores = num_cores)
        if(parallel==F) temp1 <- lapply(Z.list, function(x){  as.numeric((t(x)%*%m_l1k1%*%x)^2) })

        temp2 <- mean(unlist(temp1))
        #print(temp2)
        diag_Thetahat[k1] <- const*(temp2 - sum( as.vector(m_l1k1)*c(Chat))^2)

    }

      return(diag_Thetahat)

  }

   if(parallel==T) diag_Thetahat.list <- parallel::mclapply(1:p, diag_cal, mc.cores = num_cores)
   if(parallel==F) diag_Thetahat.list <- lapply(1:p, diag_cal)

    return(diag_Thetahat.list)

}

}


#' Estimation of the shrinkage estimator
#'
#' Given the direct covariance estimator, this function generates the optimal shrinkage which is guaranteed positive definiteness.
#'
#' @param Sigma The direct covariance estimator from direct_est().
#' @param Z An nXq matrix of lower-level measurements; n is the sample size; q is the number of lower-level variables.
#' @param A A qXp matrix characterizing the binding relationship between higher- and lower-level variables. A should have 0-1 entries.
#' If the (j,k) entry of A is 1, then the j-th lower-level variable belongs to the k-th higher-level variable.
#' @param kappa A positive shrinkage parameter or a vector of shrinkage parameters. The default is 0.1.
#' @param par1 Whether mclapply() is used. The default is TRUE.
#' @param num_cores The number of cores for parallel computing.
#' @param weights Whether the shrinkage weights are pre-specified. The default is NULL, which means that the weights need to be estimated.
#' @return A list of four objects.
#'  \itemize{
#'   \item \code{weights}: estimates of alpha, beta, and gamma.
#'   \item \code{rho}: a vector of optimal threshold for each value of kappa.
#'   \item \code{Sigma_opt.list}: a list of optimally thresholded covariance estimator for each value of kappa.
#'   \item \code{R_opt.list}: a list of optimally thresholded correlation estimator for each value of kappa.
#' }
#' @export
#' @examples
#' \dontrun{
#' data(toy_lower_dat)
#' data(toy_binding)
#' Sigma_hat <- direct_est(Z=toy_lower_dat, A=toy_binding)$cov.est
#' fit <- opt_thres(Sigma=Sigma_hat, A=toy_binding, Z=toy_lower_dat, kappa = 0.1, par1=FALSE, weights=NULL)
#' }
opt_thres <- function(Sigma, A, Z, kappa = 0.1, par1=T, weights=NULL, num_cores=8){

  p <- dim(Sigma)[1]
  q <- dim(A)[1]
  n <- dim(Z)[1]
  #mu <- sum(diag(Sigma))/p
  #alpha2 <- sum((Sigma - mu*diag(p))^2)
  #beta2 <- sum(beta_thres(Sigma, A, Z, trace = T))/n

  if(is.null(weights)){
    diag_Theta.list <- beta_thres(Sigma, A, Z, trace = T,parallel=par1, num_cores=num_cores)
    diag_Theta <- c(do.call(rbind, diag_Theta.list))
    beta2 <- sum(diag_Theta[ (p+1)*(1:p) - p ])/n
    alpha2 <- norm(Sigma, 'F')^2 - 2*sum(diag(Sigma%*%diag(diag(Sigma))))+ sum(diag(Sigma^2)) + beta2
    gamma2 <- sum(diag_Theta)/n
  }
  else{
    alpha2=weights[1]
    beta2=weights[2]
    gamma2=weights[3]
  }

  R<-cov2cor(Sigma)
  lambda_min <- eigen(R)$values[p]
  rho_LW <- (gamma2 - beta2)/(alpha2 + gamma2 - 2*beta2)
  #print(rho_LW)
  #rho2_LW <- alpha2/(alpha2 + beta2)

  n.kappa <- length(kappa)

  if(lambda_min >= 0){
    rho_opt = rho_LW
    Sigma_opt = rho_opt*diag(diag(Sigma)) + (1 - rho_opt)*Sigma
    #rho2_opt = rho2_LW
    Sigma_opt.list=Sigma_opt
  }



  if(lambda_min < 0){

    rho_opt <- rep(0, n.kappa)
    Sigma_opt.list <- list(n.kappa)

    for(s in 1:n.kappa){

      rho_UPC = (1+kappa[s])*abs(lambda_min)/(1 + (1+kappa[s])*abs(lambda_min))
      rho_UPC = max( c(rho_LW, rho_UPC) )
      #rho2_UPC = 1 - rho1_UPC/mu
      rho_opt[s] = rho_UPC
      Sigma_opt.list[[s]] = rho_opt[s]*diag(diag(Sigma)) + (1 - rho_opt[s])*Sigma

    }

    names(Sigma_opt.list)=paste0("kappa = ", kappa)

    }



  # if(lambda_min < 0){
  #   rho_opt = rho_UPC
  #   #rho2_opt = rho2_UPC
  # }

  R_opt.list <- lapply(Sigma_opt.list, FUN=cov2cor)
  return( list(weights = c(alpha2, beta2, gamma2), rho=rho_opt, Sigma_opt.list = Sigma_opt.list, R_opt.list = R_opt.list) )
}


#' Estimation of the shrinkage estimator based on cross validation
#'
#' Given the direct covariance estimator, this function generates the optimal shrinkage which is guaranteed positive definiteness, where kappa is selected by cross validation.
#'
#' @param Z An nXq matrix of lower-level measurements; n is the sample size; q is the number of lower-level variables.
#' @param A A qXp matrix characterizing the binding relationship between higher- and lower-level variables. A should have 0-1 entries.
#' If the (j,k) entry of A is 1, then the j-th lower-level variable belongs to the k-th higher-level variable.
#' @param kappa A vector of candidate values of kappa. The default is c(0.1, 0.5, 1, 5, 10, 50, 100).
#' @param par1 Whether mclapply() is used. The default is TRUE.
#' @param par2 Whether additional parallel computing is used. The default is FALSE.
#' @param weights Whether the shrinkage weights are pre-specified. The default is NULL, which means that the weights need to be estimated.
#' @param num_cores How many cores are used for parallel computing? The default is 8.
#' @param n.split How many sample splittings will be performed? The default is 50.
#' @return A list of two objects.
#'  \itemize{
#'   \item \code{cv_result}: cross-validation-based F-norm statistics for all kappas.
#'   \item \code{kappa.min}: the optimal value of kappa from the pool.
#' }
#' @export
#' @examples
#' \dontrun{
#' data(toy_lower_dat)
#' data(toy_binding)
#' fit.cv <- cv_opt_thres(A=toy_binding, Z=toy_lower_dat, n.split=20)
#' Sigma_hat <- direct_est(Z=toy_lower_dat, A=toy_binding)$cov.est
#' fit <- opt_thres(Sigma=Sigma_hat, A=toy_binding, Z=toy_lower_dat, kappa = fit.cv$kappa.min)
#' }

cv_opt_thres <- function(A, Z, kappa=c(0.1, 1, 5, 10, 50, 100), par1=T, par2=F,  weights=NULL, n.split=50, num_cores=8){

  n <- dim(Z)[1]
  p <- dim(A)[2]
  q <- dim(A)[1]
  n.kappa <- length(kappa)
  #kappa.cv <- 0

  for(i.split in 1:n.split){

    cal_cv <- function(i.split){

      set.seed(1992 + 719*i.split)
      idx <- sample.int(n, size=floor(n/2))
      Z1 <- Z[idx,] ## shrinkage ---
      Z2 <- Z[-idx,] ## original ---

      Sigma.est.Z1 <- direct_est(Z=Z1, A=A)
      Sigma.est.Z2 <- direct_est(Z=Z2, A=A)

      fit.shrink.Z1 <- opt_thres(Sigma=Sigma.est.Z1$cov.est, A=A, Z=Z1, kappa = kappa, par1=par1, weights=weights, num_cores=num_cores)
      Sigma.est.Z1.thres <- fit.shrink.Z1$Sigma_opt.list

      kappa.cv <- unlist(lapply(Sigma.est.Z1.thres, function(x){ norm(x - Sigma.est.Z2$cov.est, 'F')  } ))

    }

    ## mclappy does not really help.
    if(par2 == T) kappa_cv_list <- mclapply(1:n.split, cal_cv, mc.cores=num_cores)
    if(par2 == F) kappa_cv_list <- lapply(as.list(1:n.split), cal_cv)

    kappa_cv_mat <- do.call("rbind", kappa_cv_list)

    #}
    cv_result = colMeans(kappa_cv_mat)
    kappa.est = kappa[which.min(cv_result)]


    return( list( cv_result = cv_result,
                  kappa.min =  kappa.est) )

  }
}


#' Generating p-values for all higher-level correlations
#'
#' Given the direct covariance estimator, this function generates p-values for all higher-level correlations.
#'
#' @param Sigma The direct covariance estimator from direct_est().
#' @param Z An nXq matrix of lower-level measurements; n is the sample size; q is the number of lower-level variables.
#' @param A A qXp matrix characterizing the binding relationship between higher- and lower-level variables. A should have 0-1 entries.
#' If the (j,k) entry of A is 1, then the j-th lower-level variable belongs to the k-th higher-level variable.
#' @param Vhat An internal object. The default is NULL. Vhat is estimated in the algorithm.
#' @param xi A value or a vector of values from 0 to 1 in the null hypothesis. The default is 0.
#' @param parallel Whether mclapply() is used. The default is TRUE.
#' @param num_cores If parallel is TRUE, then how many cores are used? The default is 8.
#' @param sparseMat Whether sparse-matrix representations are used. The default is TRUE, which can accelerate the algorithm especially when q is large.
#' @return A list of p-value matrices, each item corresponding to each value of xi.
#' @export
#' @examples
#' \dontrun{
#' data(toy_lower_dat)
#' data(toy_binding)
#' Sigma_hat <- direct_est(Z=toy_lower_dat, A=toy_binding)$cov.est
#' inf <- direct_inf(Sigma=Sigma_hat, A=toy_binding, Z=toy_lower_dat, xi=c(0,0.1), parallel=FALSE)
#' }
direct_inf <- function(Sigma, Z, A, xi=0, Vhat = NULL, parallel=T, num_cores=8, sparseMat=T){

  Z <- scale(Z, center = T, scale = F)
  n <- dim(Z)[1]
  p <- dim(Sigma)[1]
  q <- dim(Z)[2]

  Z.list <- as.list(data.frame(t(Z)))

  # if(is.null(Vhat)){
  Chat <- cov(Z)
  xi.l <- length(xi)

  pvalue_cal <- function(l){


    # temp.return.pval <- rep(0, p)
    temp.return.pval <- matrix(0, xi.l, p)
    #temp.return.test.plus <- rep(0, p)
    #temp.return.test.minus <- rep(0, p)

    sigma_ll <- Sigma[l,l]
    index_l <- which(rowSums(A) == 1 & A[,l] == 1)
    S_l <- length(index_l)

    #m_ll <- rep(0, q^2) # this would be very sparse ==> may speed up
    m_ll <- matrix(0, q, q)
    #index_vec <- NULL
    for(templ in 1:S_l){
      for(tempk in 1:S_l){
        #index_temp = index_l[templ] + (index_l[tempk] - 1)*q
        #index_vec <- c(index_vec, index_temp)
        #m_ll[index_temp] <- as.numeric(index_l[templ] != index_l[tempk])
        m_ll[index_l[templ], index_l[tempk]] <- as.numeric(index_l[templ] != index_l[tempk])
      }
    }

    if(sparseMat==T){

      grid_indices <- expand.grid(index_l, index_l)

      # Create i and j vectors using rep()
      i.ind <- grid_indices[, 1]
      j.ind <- grid_indices[, 2]

      m_ll <- Matrix::sparseMatrix(i = i.ind,
                           j = j.ind,
                           dims=c(q,q),
                           x = as.vector(m_ll[index_l, index_l]))
    }

    #print(m_ll[index_l, index_l])
    #print(m_ll_s[index_l, index_l])

    #m_ll <- sparseVector(x=rep(1, length(index_vec)), i=index_vec, length=q^2)

    for(k in (l+1):p){

      #sigma_ll <- Sigma[l,l]
      sigma_kk <- Sigma[k,k]
      sigma_lk <- Sigma[l,k]
      r_lk <- sigma_lk/sqrt(sigma_ll*sigma_kk)
      f_lk <- c(1/sqrt(sigma_ll*sigma_kk), -1/2*r_lk/sigma_ll, -1/2*r_lk/sigma_kk)


      #index_l <- which(rowSums(A) == 1 & A[,l] == 1)
      #S_l <- length(index_l)
      index_k <- which(rowSums(A) == 1 & A[,k] == 1)
      S_k <- length(index_k)

      #m_lk <- rep(0, q^2)
      #index_vec <- NULL
      m_lk <- matrix(0, q, q)
      for(templ in 1:S_l){
        for(tempk in 1:S_k){
          #index_temp = index_l[templ] + (index_k[tempk] - 1)*q
          #m_lk[index_temp] <- as.numeric(index_l[templ] != index_k[tempk])
          m_lk[index_l[templ], index_k[tempk]] <- as.numeric(index_l[templ] != index_k[tempk])
          #index_vec <- c(index_vec, index_temp)
        }
      }

      if(sparseMat==T){

        grid_indices <- expand.grid(index_l, index_k)

        # Create i and j vectors using rep()
        i.ind <- grid_indices[, 1]
        j.ind <- grid_indices[, 2]

        m_lk <- Matrix::sparseMatrix(i = i.ind,
                             j = j.ind,
                             dims=c(q,q),
                             x = as.vector(m_lk[index_l, index_k]))
      }

      #m_kk <- rep(0, q^2)
      m_kk <- matrix(0, q, q)
      #index_vec <- NULL
      for(templ in 1:S_k){
        for(tempk in 1:S_k){
          #index_temp = index_k[templ] + (index_k[tempk] - 1)*q
          #m_kk[index_temp] <- as.numeric(index_k[templ] != index_k[tempk])
          #index_vec <- c(index_vec, index_temp)
          m_kk[index_k[templ], index_k[tempk]] <- as.numeric(index_k[templ] != index_k[tempk])
        }
      }
      #m_kk <- sparseVector(x=rep(1, length(index_vec)), i=index_vec, length=q^2)

      if(sparseMat==T){

        grid_indices <- expand.grid(index_k, index_k)

        # Create i and j vectors using rep()
        i.ind <- grid_indices[, 1]
        j.ind <- grid_indices[, 2]

        m_kk <- Matrix::sparseMatrix(i = i.ind,
                             j = j.ind,
                             dims=c(q,q),
                             x = as.vector(m_kk[index_k, index_k]))
      }

      const_ll = 1/(S_l*(S_l - 1))
      const_lk = 1/(S_l*S_k)
      const_kk = 1/(S_k*(S_k - 1))

      Upsilon_lk <- matrix(0, 3, 3)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_lk%*%x)  }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_lk%*%x)  })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_ll)*c(Chat))*sum(as.vector(m_lk)*c(Chat))
      Upsilon_lk[1,2] <- temp3*(const_ll*const_lk)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_lk%*%x) * as.numeric(t(x)%*%m_kk%*%x)  }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_lk%*%x) * as.numeric(t(x)%*%m_kk%*%x)  })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_lk)*c(Chat))*sum(as.vector(m_kk)*c(Chat))
      Upsilon_lk[1,3] <- temp3*(const_kk*const_lk)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_kk%*%x) }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_kk%*%x) })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_ll)*c(Chat))*sum(as.vector(m_kk)*c(Chat))
      Upsilon_lk[2,3] <- temp3*(const_kk*const_ll)

      Upsilon_lk <- Upsilon_lk + t(Upsilon_lk)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_lk%*%x) * as.numeric(t(x)%*%m_lk%*%x) }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_lk%*%x) * as.numeric(t(x)%*%m_lk%*%x) })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_lk)*c(Chat))*sum(as.vector(m_lk)*c(Chat))
      diag(Upsilon_lk)[1] <- temp3*(const_lk*const_lk)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_ll%*%x) }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_ll%*%x) * as.numeric(t(x)%*%m_ll%*%x) })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_ll)*c(Chat))*sum(as.vector(m_ll)*c(Chat))
      diag(Upsilon_lk)[2] <- temp3*(const_ll*const_ll)

      if(parallel==T)
        temp1 <- parallel::mclapply(Z.list, function(x){  as.numeric(t(x)%*%m_kk%*%x) * as.numeric(t(x)%*%m_kk%*%x)  }, mc.cores = num_cores)
      if(parallel==F)
        temp1 <- lapply(Z.list, function(x){  as.numeric(t(x)%*%m_kk%*%x) * as.numeric(t(x)%*%m_kk%*%x)  })

      temp2 <- mean(unlist(temp1))
      temp3 <- temp2 - sum(as.vector(m_kk)*c(Chat))*sum(as.vector(m_kk)*c(Chat))
      diag(Upsilon_lk)[3] <- temp3*(const_kk*const_kk)

      #var11[l,k] <- diag(Upsilon_lk)[1]

      var_lk <- t(f_lk)%*%Upsilon_lk%*%f_lk # this becomes a negative number?
      #print(var_lk)
      #print(eigen(Upsilon_lk)$values)

      for(index.xi in 1:xi.l){

        xi0 <- xi[index.xi]
        test.plus <- sqrt(n)*(r_lk - xi0)*(r_lk - xi0 > 0)/sqrt(var_lk)
        test.minus <- sqrt(n)*(r_lk + xi0)*(r_lk + xi0 < 0)/sqrt(var_lk)
        #pvalue[l,k] <- pnorm( max(c(test_stat_plus[l,k], 0)), lower.tail = F ) +pnorm( min(c(test_stat_minus[l,k], 0)), lower.tail = T )

        temp.return.pval[index.xi, k] = 2*pnorm( max( c( abs(test.plus), abs(test.minus) )), lower.tail = F)

      }
    # print(temp.return.pval)
    }

    return(temp.return.pval)
  }

  if(parallel==T) pvalue.list <- parallel::mclapply(1:(p-1), pvalue_cal, mc.cores=num_cores)
  if(parallel==F) pvalue.list <- lapply(as.list(1:(p-1)), pvalue_cal)

  pvalue.xi <- list(xi.l)
  for(index.xi in 1:xi.l){

    pvalue.list.temp <- lapply(pvalue.list, function(x){x[index.xi,]})
    pvalue <- do.call("rbind", pvalue.list.temp)
    pvalue <- rbind(pvalue, rep(0, p)) ## add a final row
    pvalue <- pvalue + t(pvalue)
    pvalue.xi[[index.xi]] <- pvalue
  }
  names(pvalue.xi) <- paste("xi0 = ", xi)
  return(pvalue = pvalue.xi)
}

#' Gene expression data
#'
#' This data set contains normalized gene-expression profiles from 278 Stage-I lung cancer patients.
#'
#' @format A data frame with 278 subjects and 2775 genes. These genes belong to 109 gene pathways.
#' @source Quantitative Biomeical Research Center. (2019). Lung Cancer Explorer; http://lce.biohpc.swmed.edu/ (15 June 2021, date last accessed).
"gene_expression_stage1"

#' Gene-pathway binding information
#'
#' This data set characterizes the binding information for the 2775 genes and 109 pathways.
#'
#' @format A binary matrix with 2775 genes and 109 gene pathways.
#' @source Information is obtained from the KEGG database.
"binding_mat"

#' Toy lower-level data
#'
#' This data set contains randomly generated lower-level data from 50 subjects.
#'
#' @format A data frame with 50 subjects and 150 lower-level variables which belong to 20 higher-level variables.
"toy_lower_dat"

#' Toy binding matrix
#'
#' This data set contains the binding information for the 150 lower-level variables and the 20 higher-level variables.
#'
#' @format A binary matrix with 150 rows and 20 columns.
"toy_binding"
