## load data ##

X = as.matrix(TCGA_data[3:116,])
X = matrix(as.numeric(X),ncol=ncol(X))
X = rbind(X[39:114,],X[1:38,])
graph = TCGA_graph
n = 240
p = 114
loc = 1
#########################################################################

grid_nu_1 = c(-2,-1,0,1,2)
grid_nu_2 = c(0.5,1)
grid_L  = rep(c(6,7,8,9,10,11,12),4)
grid_eta = c(rep(c(5,10),each=7),rep(c(5,10),each=7))
grid_eps = rep(c(0.1,0.2),each = 14)


w_ini_l = list()
z_ini_l = list()

L_seed = 1113
for (l in 1:length(grid_L)) {
  set.seed(l%%7+L_seed)
  LL= grid_L[l]
  w_ini_l[[l]]=matrix(rnorm(LL*p,0,1),nrow = p,ncol = LL)
}

for (l in 1:length(grid_L)) {
  set.seed(l%%7+L_seed)
  LL= grid_L[l]
  z_ini_l[[l]] = matrix(rnorm(LL*n,0,1),nrow=LL,ncol=n)
}


data_tbs = list()

################################ MCMC ###################################
T = 5000
Sigma = rep(0.1,p)  # variance for each entry of m
Q = diag(4,nrow = p,ncol = p) # proposal density
tun_arg = 1

L   =  grid_L[loc]
eta =  grid_eta[loc]
eps =  grid_eps[loc]
rho_ini = matrix(1,nrow=p,ncol=n )
tau_temp = matrix(1,nrow=p,ncol=L )
tau_ini  = tau_temp

w_ini = w_ini_l[[loc]] 
z_ini = z_ini_l[[loc]] 

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
start = Sys.time()
for (j in 1:length(grid_nu_1)) {
  for (i in 1:length(grid_nu_2)) {
    nu_1 = grid_nu_1[j]
    nu_2 = grid_nu_2[i]
    alpha_ini = matrix(nu_1,nrow=p,ncol=L )
    mcmc_box = empty_chain(n,p,T,L=L,rho_ini,w_ini,z_ini,alpha_ini)
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
    trials = matrix(0,p,n)
    trial_result = BFGA_MCMC(3,T=T,L=L,p=p,n=n,X=X,trials=trials,nu_1=nu_1,nu_2=nu_2,Sigma,Q,eta,eps,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,th=0.1)
    trial_name = paste0('mcmc_',tun_arg)
    data_tbs[[trial_name]] = trial_result
    tun_arg = tun_arg+1 
  }
}
end = Sys.time()
start-end


saveRDS(data_tbs,paste0('TCGA_',loc,'.rds'))  

