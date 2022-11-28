
#source("/home/qiyiwen/function_package/BFGAMCMC_info.R")
#source("/home/qiyiwen/function_package/GBFA.R")

source("E:/Post_doc/BFA_01/MCMC/BFGAMCMC_info.R")
source("E:/Post_doc/BFA_01/MCMC/GBFA.R")


#library(iClusterPlus)
#library(parallel)
###### What's followed is to generate Gaussian data #####



for (mm in 1:1) {
  arg = mm
  set.seed(19402+arg)
  
  
  n = 180
  p = 90
  L = 3
  T = 100
  n_nz = p/L  # number of non-zero elements of each col (W)
  
  
  
  ## Data Generating
  
  # m = runif(p,min=0,max=1) # we suppose m = 0 
  
  # W: p x L factor loading matrix
  W = matrix(0,nrow=p,ncol=L)
  # true graph consistent with W 
  true_graph = diag(0,nrow = p, ncol = p)
  
  
  ## different ways to generate W
  ## a) overlapped pathway
  # for (i in 1:L) {
  #   start_ind = sample(seq(1,p-n_nz+1,by=1),size=1)
  #   end_ind   = start_ind + n_nz-1
  #   W[start_ind:end_ind,i] = as.matrix(rnorm(n_nz,1.5,0.1)*sample(c(-1,1),n_nz,replace=T))
  # }
  
  ## b)non-overlapped pathway
  
  for (i in 1:L) {
    start_ind = i*n_nz-(n_nz-1)
    end_ind   = start_ind + n_nz-1
    true_graph[start_ind, (start_ind+1):end_ind] = rep(1,n_nz-1)
    W[start_ind:end_ind,i] = as.matrix(rnorm(n_nz,1.5,0.1)*sample(c(-1,1),n_nz,replace=T))
  }
  
  # star-like true graph
  true_graph = forceSymmetric(true_graph,uplo="U")
  
  # fully connecte true graph  
  
  #true_graph = true_graph+1
  #diag(true_graph) = 0
  
  
  # Z: L x n latent factor matrix
  Z = matrix(0,nrow=L,ncol=n)
  for (i in 1:L) {
    temp = rpois(1,30)
    ind = sort(sample(1:n,size=temp))   
    Z[i,ind] = t(rnorm(temp,1.5,0.1)) 
  }
  
  # X: p x n matrix 
  mu = W %*% Z
  X = W %*% Z + matrix(rnorm(p*n,0,1),nrow=p,ncol=n,byrow=T)
  
  trials = matrix(0,nrow = p, ncol = n)
  
  data_tbs = list('true_W'= W,'true_Z'=Z,'true_X'=X)
  
  
  # simulation dataset_1
  
  # ###################################################
  # 
  # data_icluster = iClusterPlus(dt1=X,type = 'gaussian')
  # 
  # data_tbs[['icluster']] = data_icluster
  # 
  # ###################################################
  # data_fabia = fabia(X,p=5)
  # data_tbs[['fabia']] = data_fabia
  
  
  
  #####################################################
  
  Sigma = rep(0.01,p)  # variance for each entry of m
  
  Q = diag(4,nrow = p,ncol = p) # proposal density
  
  ######## working graph####################################
  
  graph = working_graph(1)
  
  
  ############# argument for GBFA_EM #####################
  
  edge_index = which(graph==1,arr.ind = TRUE)
  edge_index= edge_index[order(edge_index[,1],edge_index[,2],decreasing=FALSE),]
  
  
  ########################################################
  grid_L  = c(3)
 # grid_L = c(3)
  w_ini_l = list()
  z_ini_l = list()
  for (l in 1:length(grid_L)) {
    LL= grid_L[l]
    w_ini_l[[l]]=matrix(rnorm(LL*p,0,1),nrow = p,ncol = LL)
  }
  
  for (l in 1:length(grid_L)) {
    LL= grid_L[l]
    z_ini_l = matrix(0,LL,n)
  }
  
  
  
  for (l in 1:length(grid_L)) {
    LL= grid_L[l]
    data_GBFAEM = GBFA_EM(X,rep(0,p),param=rep(1,p),edge_index,L=LL,v0=0.2,v1=1,fit.m=TRUE,m.init=0,W.init=w_ini_l[[l]],delta=1,eta=0.5)
    names = paste0('GBFAEM_',l)
    data_tbs[[names]] = data_GBFAEM
    
  }
  

  
  ## initialization 
  
  #grid_nu_1 = c(-1.5,-1,-0.5,0,0.5,1)   # mean 
  #grid_nu_2 = c(0.25,0.5,1) # variance 
  
  # grid_nu_1 = c(-1.5)
  # grid_nu_2 = c(0.25) 
  # 
  
  # 3/23
 grid_nu_1 = c(0)
 grid_nu_2 = c(1)
  eta = 10
  eps = 0.2
  tun_arg = 1
  
  
  
  
  for (j in 1:length(grid_nu_1)) {
    for (i in 1:length(grid_nu_2)) {
       
          nu_1 = grid_nu_1[j]
          nu_2 = grid_nu_2[i]
          rho_ini = matrix(1,nrow=p,ncol=n )
          
          tau_temp = matrix(1,nrow=p,ncol=L )
          tau_ini  = tau_temp
          
          alpha_ini = matrix(nu_1,nrow=p,ncol=L )
          
          #alpha_ini = matrix(2,nrow=p,ncol=L ) 
          #alpha_ini[which(W!=0,arr.ind = TRUE)] = -2.3 
          
          # W and Z are initialized with the true value
          #w_ini = W 
          #z_ini = Z 
          
          # with wrong value 
          w_ini = w_ini_l[[1]] 
          z_ini = z_ini_l[[1]] 
          
          ## omega should be positive definit
          omega_ini = diag(1,nrow=p,ncol=p)  # compatible with graph
          
          ind = which(graph!=0,arr.ind = T)  # indice of the nonzero elements
          for (i in 1:dim(ind)[1]) {
            if(ind[i,1]>ind[i,2]){
              r_ind = ind[i,1]
              c_ind = ind[i,2]
              omega_ini[r_ind,c_ind] = 0.05
            }
          }
          
          omega_temp = as.matrix(forceSymmetric(omega_ini,uplo = 'L'))
          
          inv_omega_temp = solve(omega_temp) 
          
          empty_chain <- function(n,p,T,L,rho_ini,w_ini,z_ini,alpha_ini){
            chain_m   = matrix(0,ncol = p, nrow = T+2)
            
            chain_rho   = matrix(0,ncol = p*n, nrow = T+2)
            chain_rho[1,]=as.vector(t(rho_ini))
            
            chain_w   = matrix(0,ncol = p*L, nrow = T+2)
            chain_w[1,]=as.vector(t(w_ini))
            
            chain_tau   = matrix(0,ncol = p*L, nrow = T+2)
            chain_tau[1,]=as.vector(t(tau_ini))
            
            chain_z   = matrix(0,ncol = L*n, nrow = T+2)
            chain_z[1,]=as.vector(t(z_ini))
            
            chain_alpha =matrix(0,ncol = p*L, nrow = T+2)
            chain_alpha[1,]=as.vector(t(alpha_ini))
            
            
            
            list(chain_tau=chain_tau,chain_alpha=chain_alpha,chain_rho=chain_rho,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z)
            
          }
          
          mcmc_box = empty_chain(n,p,T,L,rho_ini,w_ini,z_ini,alpha_ini)
          
          chain_rho = mcmc_box$chain_rho
          chain_w = mcmc_box$chain_w
          chain_z = mcmc_box$chain_z
          chain_alpha = mcmc_box$chain_alpha
          chain_m = mcmc_box$chain_m
          chain_tau = mcmc_box$chain_tau
          w_temp   = matrix(as.numeric(chain_w[1,]),nrow=p ,ncol=L,byrow=T)
          z_temp   = matrix(as.numeric(chain_z[1,]),nrow=L ,ncol=n,byrow=T)
          rho_temp   = matrix(as.numeric(chain_rho[1,]),nrow=p ,ncol=n,byrow=T)
          alpha_temp = matrix(as.numeric(chain_alpha[1,]),nrow=p ,ncol=L,byrow=T)
          m_temp = rep(0,p)
          
          
          
          trial_result = BFGA_MCMC(0,T=T,L=L,p=p,n=n,X,trials=trials,nu_1=nu_1,nu_2=nu_2,Sigma,Q,eta,eps,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,th=0.1)
          trial_name = paste0('mcmc_',tun_arg)
          data_tbs[[trial_name]] = trial_result
          tun_arg = tun_arg+1 
        }
      }
    
  
  saveRDS(data_tbs,paste0('tuning_alpha_info_90180_',arg,'.rds'))  
  
}






#########################

