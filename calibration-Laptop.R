library(rstan)
library(shinystan)
library(akima)
library(stats)
library(RColorBrewer)
library(latex2exp)
library(R.matlab)
library(HDInterval)
library(progress)

ss = function(d,p,a){
  
  strain = (4/3/pi)*d[-1]/a[-1]
  stress = p[-1]/pi/a[-1]**2*1e-6
  return(list(s=stress,e=strain))
}

##################################################
#        DATA PROCESSING
##################################################

#set dir
#this_dir <- function(directory)
#setwd( file.path(getwd(), directory) )
if (interactive()){
  setwd('C:\\Users\\Patxi\\OneDrive\\MELKOTE\\Experimental\\OFHC\\nanoindentation\\rstan_single\\Inverse-Spherical-Indentation')
  setwd('D:\\OneDrive\\MELKOTE\\Experimental\\OFHC\\nanoindentation\\rstan_single\\Inverse-Spherical-Indentation')
} else{
  this_dir <- function(directory)
    setwd( file.path(getwd(), directory) )
}





#####################################################################
#          load and manipulate the computer code output
#####################################################################

Mexp=100

files = list.files('LH_u_ofhc_val','[0-9]_output.txt')

files = files[order(as.numeric(gsub("[^0-9]+","",files)))] # get files ordered relative to design pnts
f = files[1]
x = as.matrix(read.csv(paste('LH_u_ofhc_val\\',f,sep=''),header=FALSE,sep=","));

out = ss(x[,3],x[,2],x[,4])
E = matrix(0,length(files),length(out$e))
S = matrix(0,length(files),length(out$e))
E[1,]=out$e
S[1,]=out$s


count=1
for (f in files[-1]){
  count=count+1
  x = as.matrix(read.csv(paste('LH_u_ofhc_val\\',f,sep=''),header=FALSE,sep=","));

  out = ss(x[,3],x[,2],x[,4])
  E[count,]=out$e
  S[count,]=out$s
}
e1 = max(E[,1])
e2 = min(E[,dim(E)[2]])


escale = e2-e1
ecenter = e1
e = seq(e1,e2,length.out=15)
scaling = matrix(c(50,200,100,1000,10,180,e1,e2),nrow=2)
S = t(sapply(seq(1,dim(S)[1]),function(i) {approx(E[i,],S[i,],e)$y}))
e = (e-ecenter)/escale



yscale = var(as.vector(S))**0.5
ycenter = mean(S)
S = (S-ycenter)/yscale

#################################################################
#          get design and compute the distance matrix
#################################################################
# load data and scale the tunning-params
design  = as.matrix(read.csv("LH_u_ofhc_val\\n30d3.txt",header=FALSE,sep=","));





# map the design to 0-1
design_norm = sweep(sweep(design,2,scaling[1,seq(1,3)],'-'),2,
                    scaling[2,seq(1,3)]-scaling[1,seq(1,3)],'/')


distNx = array(0,dim=c(dim(design)[1],dim(design)[1],dim(design)[2]))
for (i in 1:dim(design)[2]){
  distNx[,,i] = as.matrix(dist(design_norm[,i])^2)
}
distNx <- aperm(distNx, c(3,1,2))

distMx = array(0,dim=c(length(e),length(e),1))
distMx[,,1] = as.matrix(dist(e)^2)
distMx <- aperm(distMx, c(3,1,2))


###############################################################
#        RSTAN PART - STAGE 1 get the computer model
###############################################################
#
# phi_E, phi_sig0, phi_k 
# phi_eps
# 
# 5 regressors [((np.exp(eps)*E*1e3)**-b5 + ((b1*Y+b4*K*1e3*(b2+b3*np.exp(eps))))**-b5)**(-(1/b5))] 
#
# sig2
#


betaprior = matrix(c(c(log(2),log(1e-4),log(1),log(2),log(2)),
                     c(0.1,0.1,0.1,0.1,0.1)),nrow=5);
sigprior = c(4,2)

phiprior = c(50,10,50,10,50,10,1,1)

y       = as.vector(S)
inputdata1 <- list(
  N     = dim(S)[1],
  M     = dim(S)[2],
  sigprior= sigprior,
  phiprior= phiprior,
  betaprior=betaprior,

  Nnugg    = 1e-6,
  Mnugg    = 1e-4,
  y       = as.vector(S),
  distN    = distNx,
  distM    = distMx,
  design   = design_norm,
  eps      = e,
  scaling = scaling,
  yscaling = c(ycenter,yscale)
)
init = list('r_betax'=c(log(2),log(0.0001),log(1),log(2.0),log(2)),
            'sig2'=0.1,'phi'=c(rep(3,3),0.1))

m1 = stan_model(file='stage1.rstan')
map1 <- optimizing(
  object=m1,                              # Stan program
  data = inputdata1,                      # named list of data
  verbose = T,
  as_vector=F,
  iter=10000,
  #algorithm='Newton',
  init=init)

betax   = map1$par$betax
lambda  = map1$par$lambda
phi     = map1$par$phi
mu      = map1$par$mu
sig2      = map1$par$sig2
invCn      = map1$par$inv_Cn
invCm      = map1$par$inv_Cm
logdetCn   = 2*(sum(log(diag(map1$par$LN))))
logdetCm   = 2*(sum(log(diag(map1$par$LM))))
save(list=c('logdetCn','logdetCm','betax','sig2','phi','invCn','invCm','y','mu','S','e','scaling','yscale','ycenter'),
     file='stage1_params.RData')

########################################################################
#           prediction function
########################################################################

N = dim(S)[1]
M = dim(S)[2]


#-------------------------------------------
#    mean func
#-------------------------------------------
mean_func = function(x, eps){
  
  # choose to scale things so that the derived betas' have some physical meaning
  # and we can interpret them / evaluate them against our intuition
  Ehat   = x[1]*(scaling[2,1]-scaling[1,1]) + scaling[1,1];
  Yhat   = x[2]*(scaling[2,2]-scaling[1,2]) + scaling[1,2];
  Khat   = x[3]*(scaling[2,3]-scaling[1,3]) + scaling[1,3];
  epshat   = eps*(scaling[2,4]-scaling[1,4]) + scaling[1,4];
  
  
  mu = (((epshat*Ehat*1e3/(1-0.3**2))**(-betax[5])) + 
          (betax[1]*Yhat + betax[4]*Khat*1e3*(betax[2]+betax[3]*epshat))**(-betax[5])
  )**(-(1/betax[5]))
  
  
  return((mu-ycenter)/yscale)
}

#-------------------------------------------
#   map 2&3 computer code predictor
#-------------------------------------------
y_pred =  function(x,eexp){
  
  dstN   = sweep(design_norm,2,x,'-')^2
  dstN = dstN %*% diag(phi[1:3])
  rtN = exp(-apply(dstN,1,sum))
  dstM   = sweep(matrix(rep(e,length(eexp)),nrow=length(e)),2,eexp,'-')^2
  rtM = exp(-dstM*phi[4])
  
  
  Sighat = sig2*rtM - rtM%x%(t(rtN)%*%invCn%*%rtN)
  
  
  muhat = c(t(rtN)%*%invCn%*%(matrix(map1$par$err,
                                     nrow=dim(invCn)[1])%*%invCm)%*%rtM)
  muhat = muhat+mean_func(x,eexp)
  
  return(list('muhat'=muhat,'Sighat'=0))
}


########################################################################
#         cross validation of model with values in other design
########################################################################

files = list.files('LH_u_ofhc_val','[0-9]_output.txt')
files = files[order(as.numeric(gsub("[^0-9]+","",files)))] # get files ordered relative to design pnts
f = files[1]
x = as.matrix(read.csv(paste('LH_u_ofhc_val\\',f,sep=''),header=FALSE,sep=","));
out = ss(x[,3],x[,2],x[,4])
E2 = matrix(0,length(files),length(out$e))
S2 = matrix(0,length(files),length(out$e))
E2[1,]=out$e
S2[1,]=out$s


count=1
for (f in files[-1]){
  count=count+1
  x = as.matrix(read.csv(paste('LH_u_ofhc_val\\',f,sep=''),header=FALSE,sep=","));
  out = ss(x[,3],x[,2],x[,4])
  E2[count,]=out$e
  S2[count,]=out$s
}

inds = as.integer(gsub("[^0-9]+","",files))
design2  = as.matrix(read.csv("LH_u_ofhc_val\\n30d3.txt",header=FALSE,sep=","))[inds,];
scaling2 = scaling

# map the design to 0-1

inds = which(design2[,1]>design2[,3])
design2 = sweep(sweep(design2,2,scaling[1,seq(1,3)],'-'),2,
                scaling[2,seq(1,3)]-scaling[1,seq(1,3)],'/')
design2=design2[inds,]
E2=E2[inds,]
S2=S2[inds,]


inds = which(apply(design2,1,min)>0.0 & apply(design2,1,max)<1.0)
design2=design2[inds,]
E2=(E2[inds,]-ecenter)/escale
S2=S2[inds,]

E2=E2[,seq(2,58)]
S2=S2[,seq(2,58)]

Sp = S2
for (i in 1:nrow(design2)){
  Sp[i,] = y_pred(design2[i,],E2[i,])$muhat*yscale + ycenter
}

png("val.png", width=8, height=5, units="in", res=300)
par(mfrow=c(1,2),pty="s")
R2 = 1-(sum((S2-Sp)**2)/sum((S2-mean(S2))**2))
RMSE = mean((S2-Sp)**2)**0.5
plot(Sp,S2,ylim=c(0,4000),xlim=c(0,4000),
     main=TeX(sprintf('Ypred: RMSPE=%1.1f (MPa)',RMSE)),
     ylab=bquote(paste(sigma[ind] (MPa))),
     xlab=bquote(hat(sigma)[ind] (MPa)))
lines(c(0,4000),c(0,4000),col='red')


Sp = S2
for (i in 1:nrow(design2)){
  Sp[i,] = mean_func(design2[i,],E2[i,])*yscale + ycenter
}


R2 = 1-(sum((S2-Sp)**2)/sum((S2-mean(S2))**2))
RMSE = mean((S2-Sp)**2)**0.5
plot(Sp,S2,ylim=c(0,4000),xlim=c(0,4000),
     main=TeX(sprintf('Mean func: RMSPE=%1.1f (MPa)',RMSE)),
     ylab=bquote(paste(sigma[ind] (MPa))),
     xlab=bquote(hat(sigma)[ind] (MPa)))
lines(c(0,4000),c(0,4000),col='red')
dev.off()


#------------------------------------------------
#    plot NL-regression bit
#------------------------------------------------
rbPal <- colorRampPalette(brewer.pal(dim(design)[1],'YlGnBu'))

sX = sweep(sweep(design_norm,2,scaling[2,1:3]-scaling[1,1:3],'*'),2,scaling[1,1:3],'+')
se = e*(scaling[2,4]-scaling[1,4]) + scaling[1,4]

muhat = (
  (sX[,1]%*%t(se)*1e3/(1-0.3^2))**-betax[5]+
    (sX[,2]*betax[1]+1e3*betax[4]*sX[,3]%*%t(betax[2]+betax[3]*se))**-betax[5]
)**(-1/betax[5])


png("data_nonlinear_means.png", width=5, height=5, units="in", res=300)
par(mfrow=c(1,1),family='serif')
sD = S*yscale + ycenter
plot(se,sD[1,],ylim=c(0,4000),
     xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
     family="serif",pch=16,col=rbPal(dim(design)[1])[1])
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey", border = "black") 
points(se,sD[1,],ylim=c(0,4000),
       xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
       family="serif",pch=16,col=rbPal(dim(design)[1])[1])
lines(se,muhat[1,],col=rbPal(dim(design)[1])[1])
for (j in 2:dim(design)[1]){
  points(se,sD[j,],ylim=c(0,4000),
         xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
         family="serif",pch=16,col=rbPal(dim(design)[1])[j])
  lines(se,muhat[j,],col=rbPal(dim(design)[1])[j])
}
title('FEA data and Mean function prediction',outer=TRUE,line=-3)
dev.off()
###############################################################
#      LEAVE ONE OUT CROSS VALIDATION USING SHORTCUT FORMULA
###############################################################
cv =  function(){
  err = invCn %*% (matrix(y,nrow=N,ncol=M)-matrix(mu,nrow=N,ncol=M))
  cv = sweep(err,1,diag(invCn),'/')
  return(as.vector(cv))
}


Y = y*yscale+ycenter
Y_m_Yhat = cv()
pct = abs(Y_m_Yhat)/Y*100

png("loo-crossval.png", width=8, height=5, units="in", res=300)
par(mfrow=c(1,1),pty="s")
R2 = 1-(sum(Y_m_Yhat**2)/sum((Y-mean(Y))**2))
RMSE = mean(Y_m_Yhat**2)**0.5
plot(Y-Y_m_Yhat,Y,ylim=c(0,4000),xlim=c(0,4000),
     main=sprintf('R2=%1.3f RMSPE=%1.1f (MPa)',R2,RMSE),
     ylab=bquote(paste(sigma[ind] (MPa))),
     xlab=bquote(hat(sigma)[ind] (MPa)))
lines(c(0,4000),c(0,4000),col='red')
dev.off()

##############################################################################
#      Stage 2 - "fit" the model given the experimentally estimates hat(mu)
##############################################################################
#
#----------------------------------------------------------------------------
#                load the experimental data
#----------------------------------------------------------------------------

f = '002.txt'
expdata  = read.csv(file.path(getwd(),'experimental',f),header=FALSE,sep=",");
inds=sort(expdata[,1],decreasing=F,index.return=T)$ix
expdata[,1]=expdata[inds,1]
expdata[,2]=expdata[inds,2]
plot(expdata[seq(1,nrow(expdata),length.out=50),1],expdata[seq(1,nrow(expdata),length.out=50),2],
       col='green',pch=16)
lines(expdata[,1],expdata[,2],col='green')


expdata = expdata[which(expdata[,1]>e1)[1]:(min(which(expdata[,1]>e2)[1],nrow(expdata),na.rm=TRUE)-1),]
expdata[,1]=(expdata[,1]-ecenter)/escale
expdata[,2]=(expdata[,2]-ycenter)/yscale
eexp=seq(min(expdata[,1]),max(expdata[,1]),length.out = Mexp)
yexp = approx(expdata[,1],expdata[,2],eexp)$y


see=expdata[,1]*escale+ecenter
ss=expdata[,2]*yscale+ycenter


io=which(see>0.002)[1]
X=cbind(rep(1,io),see[1:io])
sig_int = (solve(t(X)%*%X,t(X)%*%ss[1:io])[1])/yscale
sig_int_sig = (sum(err**2)/(length(err)-1)*(solve(t(X)%*%X,diag(rep(1,2)))[1,1]))**0.5/yscale
Eprior=solve(t(X)%*%X,t(X)%*%ss[1:io])[2]*(1-0.3**2)
estar = see[which((ss-(see-0.001)*Eprior) < 0)[1]]
err = (ss[1:io]-solve(t(X)%*%X,t(X)%*%ss[1:io])[1]-see[1:io]*Eprior)
Esig = (sum(err**2)/(length(err)-1)*(solve(t(X)%*%X,diag(rep(1,2)))[2,2]))**0.5/1000/(scaling[2,1]-scaling[1,1])*(1-0.3**2)
Yindprior=which((ss-Eprior*(see-0.0002))<0)[1]
Yindprior=ss[Yindprior]/2
Eprior=((Eprior/1000)-scaling[1,1])/(scaling[2,1]-scaling[1,1])
Yindprior=(Yindprior-scaling[1,2])/(scaling[2,2]-scaling[1,2])


sigprior = c(5,1)
xprior=c(Eprior,Esig,Yindprior,0.5,0.5,1)


inputdata2     = list(
  N = dim(invCn)[1],
  M = dim(invCm)[1],
  xprior=xprior,
  Me=Mexp,
  yComp=as.vector(S),
  Me = Mexp,
  epsexp = eexp,
  yObs =yexp-sig_int,
  xprior=xprior,
  design=design_norm,
  eps=e,
  scaling=scaling,
  yscaling=c(ycenter,yscale),
  Rinv=invCn,
  logdetR=logdetCn,
  logdetSig=logdetCm,
  Siginv=invCm,
  phi=phi,
  betax=betax,
  mu=mu,
  err=matrix(map1$par$err,nrow=nrow(S)),
  lambdaprior=c(2,4),
  sig_int_prior=c(sig_int,sig_int_sig),
  estar=estar
)

m2 = stan_model(file='stage2.rstan')
map2 <- optimizing(
  object=m2,                               # Stan program
  data = inputdata2,                       # named list of data
  verbose = T,
  as_vector=F,
  #  algorithm='Newton',
  #  tol_obj = 1e-4,
  init = list('sig2eps'=0.005,'X'=c(0.5,0.5,0.1),'epsbeta'=rep(0,3)),
  iter=1000)

sprintf('E=%1.1f GPa, sig0=%1.1f MPa, K=%1.1f GPa',map2$par$sX[1],map2$par$sX[2],map2$par$sX[3])
X_map=map2$par$X
savefile=sprintf('%s/experimental/%s_map_params.RData',getwd(),strsplit(f,'.txt')[[1]])
save(list=c('X_map','sig2eps_map'),file=savefile)



#------------------------------------------------
#    plot the fit
#------------------------------------------------
png(file.path(getwd(),'experimental',"fit_map.png"), width=5, height=5, units="in", res=300)
x=map2$par$X
par(mfrow=c(1,1),family='serif')
plot(expdata[,1]*escale+ecenter,expdata[,2]*yscale+ycenter,
     ylim=c(0,max(expdata[,2]*yscale+ycenter)),xlim=c(0,1.1*max(expdata[,1]*escale+ecenter)),
     xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
     family="serif",pch=16)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey", border = "black") 
points(expdata[,1]*escale+ecenter,expdata[,2]*yscale+ycenter,
       ylim=c(0,max(expdata[,2]*yscale+ycenter)),
       xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
       family="serif",pch=1,col='black')
lines(eexp*escale+ecenter,(sig_int+y_pred(x,eexp)$muhat)*yscale+ycenter,col='red',lwd=1)



dev.off()

#------------------------------------------------
#    plot residuals
#------------------------------------------------
png(file.path(getwd(),'experimental',"resid_map.png"), width=5, height=5, units="in", res=300)
x=map2$par$X
ub = 1.1*max(c(map2$par$sig2epsv**0.5*2*yscale,map2$par$err2*yscale))
lb = -ub

par(mfrow=c(1,1),family='serif')
plot(eexp*escale+ecenter,map2$par$sig2epsv**0.5*yscale,
     ylim=c(lb,ub),xlim=c(0,max(eexp*escale+ecenter)),
     xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
     family="serif",pch=16)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey", border = "black") 
lines(eexp*escale+ecenter,2*map2$par$sig2epsv**0.5*yscale,
       xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
       family="serif",pch=16)
lines(eexp*escale+ecenter,-2*map2$par$sig2epsv**0.5*yscale,
      xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
      family="serif",pch=16)
points(eexp*escale+ecenter,map2$par$err2*yscale,col='black')
title('Residuals and 95% Observation error bounds')
dev.off()


########################################################################
#         MCMC SAMPLING ON THETA/FRIC
########################################################################


mcmc <- stan(
  file='stage3.rstan',         # Stan program
  data = c(inputdata2,list('epsbeta'=map2$par$epsbeta)),            # named list of data
  chains = 1,                # number of Markov chains
  warmup=1000,               # warmup
  iter = 1000+10*1000,             # total number of iterations per chain
  cores = 4,                 # number of cores (using 2 just for the vignette)
  refresh = 10*10,            # show progress every 'refresh' iterations
  thin =10 ,              # thining
  verbose=FALSE,
  init = list(list('sig2eps'=sig2eps_map,'X'=X_map)))


savefile=file.path(getwd(),'experimental',sprintf('%s_stats.txt',strsplit(f,'.txt')[[1]]))
s=summary(mcmc,pars=c('lp__','sX','sig2eps'))$summary
capture.output(s,file=savefile)


savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_chain.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
traceplot(mcmc,par=c('lp__','sX','sig2eps'), inc_warmup = F, nrow = 5)
dev.off()

savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_params.RData',strsplit(f,'.txt')[[1]]))
mcmc = extract(mcmc,permuted=TRUE)
save(list=c('s','mcmc'),file=savefile)


savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_pairs.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
pairs(mcmc$sX)
dev.off()

savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_post_sX.txt',strsplit(f,'.txt')[[1]]))
write.table(mcmc$sX, file=savefile, row.names=FALSE, col.names=FALSE)

savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_hist.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
par(mfrow=c(1,3),family='serif')
hist(mcmc$sX[,1],breaks=20)
hist(mcmc$sX[,2],breaks=20)
hist(mcmc$sX[,3],breaks=20)
dev.off()
########################################################################
#         MCMC POSTERIOR PREDICTIONS
########################################################################


savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_posterior_mean.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
par(mfrow=c(1,1),family='serif')
plot(expdata[,1]*escale+ecenter,expdata[,2]*yscale+ycenter,
     ylim=c(0,max(expdata[,2]*yscale+ycenter)),xlim=c(0,1.1*max(expdata[,1]*escale+ecenter)),
     xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
     family="serif",pch=16)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey", border = "black") 
points(expdata[,1]*escale+ecenter,expdata[,2]*yscale+ycenter,
       ylim=c(0,max(expdata[,2]*yscale+ycenter)),
       xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
       family="serif",pch=1,col='black')

YP = t(sweep(apply(mcmc$X,1,function(x) y_pred(x,eexp)$muhat),1,sig_int,'+')*yscale+ycenter)
lines(eexp*escale+ecenter,apply(YP,2,mean),col='red',lwd=1)
bnds = apply(YP,2,function(x) hdi(x,credMass=0.99))
polygon(c(eexp*escale+ecenter,rev(eexp*escale+ecenter)),
      c(bnds[1,],rev(bnds[2,])),col=rgb(1,0,0,0.5),border=NA)
dev.off()

########################################################################
#         residuals
########################################################################


savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_residuals.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
par(mfrow=c(1,1),family='serif')
x=map2$par$X
ub = 1.1*max(c(mcmc$sig2epsv**0.5*2*yscale,mcmc$err2*yscale))
lb = -ub

par(mfrow=c(1,1),family='serif')
plot(eexp*escale+ecenter,apply(mcmc$sig2epsv**0.5*2*yscale,2,mean),
     ylim=c(lb,ub),xlim=c(0,max(eexp*escale+ecenter)),
     xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
     family="serif",pch=16)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey", border = "black") 
lines(eexp*escale+ecenter,apply(mcmc$sig2epsv**0.5*2*yscale,2,mean),
      xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
      family="serif",pch=16)
lines(eexp*escale+ecenter,apply(mcmc$sig2epsv**0.5*-2*yscale,2,mean),
      xlab=TeX('$\\epsilon$'),ylab=TeX('$\\sigma\\,\\,(MPa)$'),
      family="serif",pch=16)
points(eexp*escale+ecenter,apply(mcmc$err2*yscale,2,mean),col='black')

bnds = apply(mcmc$sig2epsv**0.5*2*yscale,2,function(x) hdi(x,credMass=0.99))
polygon(c(eexp*escale+ecenter,rev(eexp*escale+ecenter)),
        c(bnds[1,],rev(bnds[2,])),col=rgb(0.3,0.3,0.3,0.4),border=NA)
polygon(c(eexp*escale+ecenter,rev(eexp*escale+ecenter)),
        c(-bnds[1,],rev(-bnds[2,])),col=rgb(0.3,0.3,0.3,0.4),border=NA)

title('Mean residuals and 95% Observation error bounds')
dev.off()


########################################################################
#         ACF
########################################################################

savefile=file.path(getwd(),'experimental',sprintf('%s_mcmc_acf.png',strsplit(f,'.txt')[[1]]))
png(savefile, width=5, height=5, units="in", res=300)
acf(apply(mcmc$err2,2,mean),main="")
title('ACF on the mean residuals')
dev.off()

