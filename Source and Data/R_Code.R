library(xtable)
library(plyr)
library(stargazer)
library(lubridate)
library(ggplot2)


rm(list = ls())

# set output directory
out.dir <- getwd()

# industry data using Fama French Data
file.dir <- "Industry Portfolios Monthly"

# stock data
snp.file <- "SNP_50_A.csv"
snp.file2 <- "SNP_50_B.csv"

file.names <-  paste(file.dir,list.files(file.dir),sep = "/")

{
  # snp data
  snp.ds <- read.csv(snp.file,stringsAsFactors = F)
  snp.ds <- snp.ds[snp.ds[,1] < 201600,]
  
  snp.ds2 <- read.csv(snp.file2,stringsAsFactors = F)
  snp.ds2 <- snp.ds2[snp.ds2[,1] < 201600,]
  
  #### read data ########
  ds.list <- numeric()
  d.seq <- c(5,10,17,30,48)
  for(i.d in 1:length(d.seq)) {
    d <- d.seq[i.d]
    file.i <- grep(d,file.names)
    ds <- tolower(scan(file.names[file.i],what = character()))
    i <- grep("monthly",ds)[1]
    j <- grep("average",ds)[2]
    i.seq <- (i+1):(j-1)
    ds <- c("date",ds[i.seq])
    ds <- t(matrix(ds,d+1))
    var.n <- ds[1,]
    ds <- data.frame(apply(ds[-1,],2,as.numeric))
    names(ds) <- var.n
    ds.list[i.d] <- list(ds)
  }
  
  ds.f <- function(ds) {
    ds <- ds[order(ds$date),]
    ds[,-1][ds[,-1] == -99.99 ] <- NA
    ds[,-1] <- ds[,-1]/100
    ds <- na.omit(ds)
    ds <- ds[as.numeric(substr(ds$date,1,4)) >=1970,]
    ds <- ds[as.numeric(substr(ds$date,1,4)) <= 2015,]
    return(ds)
  }
  
  ds.list <- lapply(ds.list,ds.f)
  ## ------------> drop the 5-FF data set <---------------------
  ds.list[[1]] <- NULL
  
  # add the S&P 50 A and B data
  ds.list[length(ds.list)+ 1] <- list(snp.ds)
  ds.list[length(ds.list)+ 1] <- list(snp.ds2)
  
 }

# I have have 6 datasets : 4 for FF and 2 for industry
sapply(ds.list, ncol)
######################################################################

#### MOTIVATION ##### 

# EXPECTED UTILITY
EU <- function(X,M,S,k) c(t(X)%*%M) - 0.5*k*c(t(X)%*%S%*%X)

MV_portfolio <- function(M,S,k,short) {
  eps <- 10^-3
  # CONSTRAINTS FOR SUM TO ONE
  d <- nrow(S)
  A <- matrix(1,1,d)
  A <- rbind(A,-A)
  B <- c(0.999,-1.001)
  
  # in case short constrainst are added
  if(short) { 
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A <- rbind(A,A2)
    B <- c(B,B2)
   }
  
  # set initial guess
  X0 <- rep(1/d,d)
  
  f <- function(x) -EU(x,M,S,k)
  g <- function(X) -M + k*S%*%X
  X1 <- constrOptim(X0,f,grad = g,ui = A,ci = B)$par
  
  return(X1)
}

GMV_portfolio <- function(S,short) {
  d <- nrow(S)
  A <- matrix(1,1,d)
  A <- rbind(A,-A)
  B <- c(0.999,-1.001)
  
  # in case short constrainst are added
  if(short) { 
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A <- rbind(A,A2)
    B <- c(B,B2)
  }
  
  X0 <- rep(1/d,d)
  
  portfolio_variance <- function(X) c(t(X)%*%S%*%X)
  g <- function(X) 2*c(S%*%X)
  X0 <- constrOptim(X0,portfolio_variance,grad = NULL,ui = A,ci = B)$par
  return(X0)
  }

# FUNCTION TO PLOT THE MVE FRONTIER FOR A GIVEN BC
MVE_function <- function(ds.i) {  
  
  ds <- ds.list[[ds.i]]
  R <- ds[,-1]
  
  i.index <- 1:nrow(R)
  set.seed(13)
  i.in <- sample(i.index,floor(0.5*nrow(R))  )
  i.out <- i.index[!i.index %in% i.in]
  length(intersect(i.in,i.out) ) == 0 # should be true
  
  R.in <- R[i.in,] 
  R.out <- R[i.out,]
  
  M.in <- apply(R.in, 2, mean)
  S.in <- var(R.in)
  
  M.out <- apply(R.out, 2, mean)
  S.out <- var(R.out)
  
  # FUNCTION TO PRODUCE A LIST OF PORTFOLIOS FOR DIFFERENT k VALUES
  MV_portfolios <- function(M,S,short) {
    d <- ncol(S) # number of assets
    e <- as.matrix(rep(1,d)) # vector of ones
    MV_portfolio_k <- function(k) MV_portfolio(M,S,k,short)
    k.seq <- c(seq(2.5,5,length = 20),seq(5,10,length = 20),seq(10,100,length = 60))
    X.list <- lapply(k.seq,MV_portfolio_k)
    return(X.list)
  }

  # OUT-OF-SAMPLE CONSTRUCTED PORTFOLIOS
  M1 <- M.in
  M2 <- M.out
  
  S1 <- S.in
  S2 <- S.out
  
  X.list <- MV_portfolios(M2,S2,short)
  # IN-SAMPLE CONSTRUCTED PORTFOLIOS
  X2.list <- MV_portfolios(M1,S1,short)
  
  # OUT-OF-SAMPLE FRONTIER (HYPOTHETICAL CASE)
  M_p <- sapply(X.list,function(X) t(X)%*%M2 )
  V_p <- sqrt(sapply(X.list,function(X) t(X)%*%S2%*%X))
  MVE <- data.frame(M = M_p,V = V_p)
  
  # IN-SAMPLE FRONTIER (REALISTIC CASE)
  M2_p <- sapply(X2.list,function(X) t(X)%*%M2 )
  V2_p <- sqrt(sapply(X2.list,function(X) t(X)%*%S2%*%X))
  MVE2 <- data.frame(M = M2_p,V = V2_p)
  
  X_gmv1 <- GMV_portfolio(S2,short)
  X_gmv2 <- GMV_portfolio(S1,short)
  
  # HIGHLIGHT THE GMV POINT
  M_0 <- t(X_gmv1)%*%M2
  V_0 <- sqrt(t(X_gmv1)%*%S2%*%X_gmv1)
  MVE <- rbind(MVE,c(M_0,V_0))
  
  M_02 <- t(X_gmv2)%*%M2
  V_02 <- sqrt(t(X_gmv2)%*%S2%*%X_gmv2)
  MVE2 <- rbind(MVE2,c(M_02,V_02))
  
  # ORDER THE FRONTIER
  MVE <- MVE[order(MVE$M),]
  MVE2 <- MVE2[order(MVE2$M),]
  
  # ADD THE NAIVE PORTFOLIO
  X_N <- rep(1/ncol(R),ncol(R))
  M_N <- t(X_N)%*%M2
  V_N <- sqrt(t(X_N)%*%S2%*%X_N)
  
  list(MVE2 = MVE2,MVE = MVE,NAIVE = c(M_N,V_N), GMV = c(M_02,V_02) )
}

short <- F
for(ds.i in 1:length(ds.list)){
  MVE.i <- MVE_function(ds.i)
  
  MVE.out <- MVE.i$MVE
  MVE.in <- MVE.i$MVE2
  
  MVE.N <- data.frame(t(MVE.i$NAIVE),Type = "Naive")
  names(MVE.N)[1:2] <- c("M","V")
  
  mve.ds <- data.frame(MVE.out,Type = "Out")
  mve.ds <- rbind(mve.ds,data.frame(MVE.in,Type = "In"))
  mve.ds <- rbind(mve.ds,MVE.N)
  summary(mve.ds)
  
  p <- ggplot() 
  p <- p + geom_point(data = mve.ds[mve.ds$Type %in% c("In","Out"),],aes(x = V, y = M, colour = Type), size = 0.5)
  p <- p + stat_smooth(data = mve.ds[mve.ds$Type %in% c("In","Out"),],aes(x = V, y = M, colour = Type),size = 0.5,fill = NA, method = "loess")
  p <- p + ylab(expression(mu[p]))   + xlab(expression(sigma[p]))
  p <- p + geom_point(data = mve.ds[mve.ds$Type %in% c("Naive"),],aes(x = V, y = M, colour = Type), size = 3)
  
  
  out.file <- paste(out.dir,"/MVE_in_out_",ds.i,"_",short*1,".pdf",sep = "")
  pdf(out.file)
  print(p)
  dev.off()
  }


##############################################################
####### FUNCTION TO CONSTUCT THE EFFICIENT FRONTIER ##########
##############################################################

### this the MVEF constructed using our model

plot_mve <- function(ds,n,k,add.plot = F) { 
  
  # take the following as the true parameters
  RET <- ds[,-1]
  
  d <- ncol(RET)
  S <- var(RET)
  M <- apply(RET,2,mean)
  V <- solve(S)
  e <- as.matrix(rep(1,d))
  alpha0 <- V%*%e/sum(V)
  B <- V%*%(diag(1,d) - e%*%t(alpha0))
  alpha1 <- B%*%M
  
  eta0 <- c(t(alpha0)%*%M)
  eta1 <- c(t(alpha1)%*%M)
  sigma0 <- c(t(alpha0)%*%S%*%alpha0)
  
  ##### LEMMA 1 ####
  I. <- diag(rep(1,d))
  A <- I. - e%*%t(alpha0)
  var_Z <- V%*%A%*%S%*%t(A)%*%V
  # check if true
  norm(var_Z - B); round(var_Z,2) == round(B,2)
  
  
  #### PROPOSITION 1 ####
  #  look at the frontier under estimation riskand full information
  mve1 <- function(s) eta0 + sqrt(eta1*(s-sigma0))
  mve2 <- function(s) {
    a <- n*(eta1^2)/(d-1+(n+1)*eta1)
    b <- eta0 + sqrt(a*(s-sigma0))
    return(b)
  }
  
  theta <- ((d-1)+eta1)/(2*n*k)
  
  p <- numeric()
  
  if(add.plot) {
    
    x.range <- quantile(seq(sigma0,max(diag(S)),length = 100))
    x <- seq(x.range[[1]],x.range[[2]],length = 40)
    x <- c(x,seq(x.range[[2]],x.range[[4]],length = 20))
    x <- c(x,seq(x.range[[4]],x.range[[5]],length = 10))
    
    y1 <- mve1(x) 
    lo1 <- predict(loess(y1~x))
    y2 <- mve2(x)
    lo2 <- predict(loess(y2~x))
    x <- sqrt(x)
    
    mve.ds <- data.frame(x = x, y = y2,Type = "Estimation")
    mve.ds <- rbind(mve.ds,data.frame(x = x,y = y1,Type = "Full"))

    p <- ggplot() 
    p <- p + geom_point(data = mve.ds,aes(x = x, y = y, colour = Type), size = 0.5)
    p <- p + stat_smooth(data = mve.ds,aes(x = x, y = y, colour = Type),size = 0.5,fill = NA, method = "loess")
    p <- p + ylab(expression(mu[p]))   + xlab(expression(sigma[p]))
    
  }
  
  results <- c(d,n,eta1,theta,eta1 < (d-1)/(n-1))
  list(results,p)
  
}

# CRETAE PLOTS FOR MVE USING PROPOSITION 1
plot_n_48 <- function(n) print(plot_mve(ds.list[[4]],n,100,TRUE)[[2]] )

n.seq <- c(120,240,360,480,1200)
for(n in n.seq) {
  file.n <- paste(out.dir,"/mve_",n,".pdf",sep = "")
  pdf(file.n)
  plot_n_48(n)
  dev.off()
}

#######################################################################################################
#######################################################################################################

###################################
#### CONSTRUCT PORTFOLIOS #########
###################################

BC.d <- function(d,alph = 0.05) {
  A <- matrix(1,1,d)
  A <- rbind(A,-A)
  B <- c(0.99999,-1.00001)

  # SHORT-SALES CONSTRAINTS
  if(short) { 
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A <- rbind(A,A2)
    B <- c(B,B2)
    }
  
  # upper limits as in Jaganathan and Ma 2003
  A3 <- -diag(rep(1,d))
  B3 <- -rep(1/d + alph,d)
  
  A <- rbind(A,A3)
  B <- c(B,B3)
  
  A4 <- diag(rep(1,d))
  B4 <- rep(1/d - alph,d)
  A <- rbind(A,A4)
  B <- c(B,B4)
  
  BC <- list(A,B)
  return(BC)
  }  


MV_portfolio <- function(M,S) {
  #f <- function(X) -c(t(X)%*%M)/sqrt(c(t(X)%*%S%*%X))
  f <- function(X) -c(t(X)%*%M)/(c(t(X)%*%S%*%X))
  g <- function(X) -((1/c(t(X)%*%S%*%X)^2)*( c(t(X)%*%S%*%X)*M - 2*c(t(X)%*%M)*S%*%X  ))
  d <- ncol(S)
  BC <- BC.d(d)
  A <- BC[[1]]
  B <- BC[[2]]
  XN <- rep(1/d,d)
  X1 <- constrOptim(XN,f,grad = g,ui = A,ci = B)$par
  return(X1)
}

# FUNCTION FOR GMV PORTFOLIO
gmv_portfolio <- function(S) {
  d <- ncol(S)
  portfolio_variance <- function(X) c(t(X)%*%S%*%X)
  g <- function(X) 2*S%*%X
  XN <- rep(1/d,d)
  BC <- BC.d(d)
  A <- BC[[1]]
  B <- BC[[2]]
  # solving this should give the minimum variance portfolio (GMV)
  X_gmv <- constrOptim(XN,portfolio_variance,grad = g,ui = A,ci = B)$par
  return(X_gmv)
}

############################################################
###### PERFORMANCE OVER TIME TEP VERSUS MVP  ###############
############################################################

portfolio.rolling <- function(ds.list,ds.i,T.) {
  # Pi <- prod.cond1(ds.i,T.)[[1]]
  t.seq <- function(t.) (t. - T.+1):t.
  ds2 <- ds.list[[ds.i]]
  n <- T.
  d <- ncol(ds2)-1
  
  port.t <- function(t.) {
    cat(t.,"\n")
    ds3 <- ds2[t.seq(t.),]
    RET1 <- ds3[,-1]
    RET_out <- ds2[t.seq(t.)+1,-1]
    
    d <- ncol(RET1)
    M <- apply(RET1,2,mean)
    S <- var(RET1[1:n,])
    V <- solve(S)
    e <- as.matrix(rep(1,d))
    X0 <- V%*%e/sum(V)
    B <- V%*%(diag(1,d) - e%*%t(X0))
    eta1 <- c(t(M)%*%B%*%M)
    
    A <- c(t(e)%*%V%*%M)
    X <- X0 + (1/A)*B%*%M
    
    if(constrain) {
      X0 <- gmv_portfolio(S)
      X <- MV_portfolio(M,S)
      }
    
    # a practical implementation
    muBmu <- eta1
    mBm.sd <- sqrt(2*(d-1 + 2*n*muBmu))
    sigma_m <- sd(apply(RET1,2,mean))
    sig <- mean(sqrt(diag(S)))
    R_cor <- cor(RET1)
    rho <- mean(R_cor[upper.tri(R_cor)])
    sig_cond1 <- (sigma_m/sig)
    
    # the last 12 monthly returns
    RET12 <- RET1[(nrow(RET1)-12+1):nrow(RET1),] # last 12 months returns
    
    X_N <- rep(1/d,d)
    
   
    f.star <- function(PI) {
      Pi <- PI[1]
      Pi2 <- PI[2]
      #Pi <- 1 - pchisq(n*(d-1)/(n-1) + K*mBm.sd,d-1,n*muBmu)
      X_pi <- Pi*X +
        (1-Pi)*( X0*Pi2   + (1-Pi2)*X_N  )
      ret.pi <- as.matrix(RET12) %*% X_pi
      return(-mean(ret.pi)/sd(ret.pi))
      }
    
    # if needed to optimize
    if(opt)
      Pi.star <- nlminb(runif(2), f.star,lower = c(0,0),upper = c(1,1))[[1]]
    
    else 
      Pi.star <- c(0,1)
    
    
    Pi <- Pi.star[1]
    Pi2 <- Pi.star[2]
    
    
    X_pi <- (1-Pi)*(X0*Pi2 + X_N*(1-Pi2) ) + Pi*X
    
    list(X_0 = X0, X_p = X,  X_pi =  X_pi, Pi = Pi,Pi2 = Pi2,X_N = X_N,muBmu = muBmu, sigma_m = sigma_m, sig_cond1 = sig_cond1,A = A)
  }
  
  t.seq2 <- t.start:(nrow(ds2)-1) # in
  t.seq3 <- (t.start+1):nrow(ds2) # out
  
  port.l <- lapply(t.seq2, port.t )
  names(port.l) <-  ds2[t.seq2,"date"]
  
  Pi <- sapply(port.l,function(x) x$Pi) 
  Pi2 <- sapply(port.l,function(x) x$Pi2) 
  
  X.list <- list(X_0 =  t(sapply(port.l, function(x) x$X_0 )) , X_p = t(sapply(port.l, function(x) x$X_p)), X_pi = t(sapply(port.l, function(x) x$X_pi)), X_N = t(sapply(port.l, function(x) x$X_N))   )
  
  R.list <- lapply(X.list, function(x)  apply(x* ds2[t.seq3,-1],1,sum)    )
  
  SR <- sqrt(12)*sapply(R.list, mean)/sapply(R.list, sd)
  TO <- sapply(X.list, function(X) mean(apply(abs(X[-1,] - X[-nrow(X),] )  ,1,sum) ) ) 
  
  muBmu <- sapply(port.l,function(x) x$muBmu) 
  sigma_m <- sapply(port.l,function(x) x$sigma_m) 
  sig_cond1 <- sapply(port.l,function(x) x$sig_cond1) 
  
  # compute the next period variance
  sigma_0 <- var(R.list$X_0)
  sigma_N <- var(R.list$X_N)
  sig_cond2 <- sqrt(sigma_N/sigma_0 )
  
  
  list(SR = SR,TO = TO,Pi = Pi,Pi2 = Pi2, muBmu = muBmu, sigma_m = sigma_m, sig_cond1 = sig_cond1 , sig_cond2 = sig_cond2)
}

# set the start date
t.start <- 200
constrain <- FALSE
opt <- TRUE # 
d.seq <- 1:length(ds.list)
n.seq <- seq(60,120, by = 30)
par.seq <- expand.grid(n.seq,d.seq)

SR <- numeric()
TO <- numeric()
Pi.list <- numeric()
Pi2.list <- numeric()
muBmu.list <- numeric()
sigma_m.list <- numeric()
sig_cond1.list <- sig_cond2.list <- numeric()

for(i in 1:nrow(par.seq)) {
  cat(i," out of ", nrow(par.seq), "\n")
  n.i <- par.seq[i,1]
  d.i <- par.seq[i,2]
  
  result.i <- portfolio.rolling(ds.list,d.i,n.i) 
  SR <- rbind(SR,result.i$SR)
  TO <- rbind(TO,result.i$TO)
  Pi.list[i] <- list(result.i$Pi)
  Pi2.list[i] <- list(result.i$Pi2)
  muBmu.list[i] <- list(result.i$muBmu)
  sigma_m.list[i] <- list(result.i$sigma_m)
  sig_cond1.list[i] <- list(result.i$sig_cond1)
  sig_cond2.list[i] <- list(result.i$sig_cond2)
}


names(par.seq) <- c("n","d")

ds.SR <- data.frame(par.seq,SR)
ds.SR$pi <- sapply(Pi.list, mean)
ds.SR$pi2 <- sapply(Pi2.list, mean)
ds.SR$muBmu <- sapply(muBmu.list, mean)
ds.SR$sigma_m <- sapply(sigma_m.list, mean)
ds.SR$sig_cond1 <- sapply(sig_cond1.list, mean)
ds.SR$sig_cond2 <- sapply(sig_cond2.list, mean)

ds.SR$d <- mapvalues(ds.SR$d, from = d.seq, to = sapply(ds.list,ncol)-1 )

# add B and S to stock label
ds.SR$d <- as.character(ds.SR$d)
ds.SR$d[13:15] <- paste(50,"B",sep = "")
ds.SR$d[16:18] <- paste(50,"S",sep = "")

sd.SR.plot <- ds.SR
sd.SR.plot$d <- as.factor(sd.SR.plot$d)
sd.SR.plot$n <- as.factor(sd.SR.plot$n)

order.names <- c("n","d","X_p","X_0","X_N","X_pi","pi","pi2","sig_cond1","sig_cond2")
sd.SR.plot <- sd.SR.plot[,order.names]
print(xtable(sd.SR.plot),include.rownames = F)

########################

# MV versus GMV
file.i <- paste(out.dir,"/SR1_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(sd.SR.plot, aes(x = X_p, y = X_0, colour = n,shape = d)) + geom_point(size = 3) 
p <- p +  geom_point(aes(group=interaction(d,n)))
p <- p + geom_abline(intercept = 0, slope = 1,linetype = "dashed")
p <- p + labs(x = expression(X), y = expression(X[0])) 
print(p)
dev.off()

# GMV versus Naive
file.i <- paste(out.dir,"/SR2_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(sd.SR.plot, aes(x = X_N, y = X_0, colour = n,shape = d)) + geom_point(size = 3) 
p <- p +  geom_point(aes(group=interaction(d,n)))
p <- p + geom_abline(intercept = 0, slope = 1,linetype = "dashed")
p <- p + labs(x = expression(X[N]), y = expression(X[0])) 
print(p)
dev.off()

# MV versus Naive
file.i <- paste(out.dir,"/SR3_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(sd.SR.plot, aes(x = X_N, y = X_p, colour = n,shape = d)) + geom_point(size = 3) 
p <- p +  geom_point(aes(group=interaction(d,n)))
p <- p + geom_abline(intercept = 0, slope = 1,linetype = "dashed")
p <- p + labs(x = expression(X[N]), y = expression(X)) 
print(p)
dev.off()


######### LINK DIFFERENCE IN SR W.R.T CONDITIONS

file.i <- paste(out.dir,"/C1_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(sd.SR.plot[sd.SR.plot$n == "60",], aes(x = sig_cond1, y = X_p - X_0)) + geom_smooth(method=lm,col = 1,lwd = 0.5) 
p <- p + geom_point(aes(x = sig_cond1, y = X_p - X_0, colour = d),size = 3)  +
  labs(y = expression(SR[p]-SR[0]), x = "Condition 1" ) 
print(p)
dev.off()


for(n.i in c("60","90","120")) {
  file.i <- paste(out.dir,"/C1_",constrain,"_",n.i,".pdf",sep = "")
  pdf(file.i)
  p <- ggplot(sd.SR.plot[sd.SR.plot$n == n.i,], aes(x = sig_cond1, y = X_p - X_0)) + geom_smooth(method=lm,col = 1,lwd = 0.5) 
  p <- p + geom_point(aes(x = sig_cond1, y = X_p - X_0, colour = d),size = 3)  +
  labs(y = expression(SR[p]-SR[0]), x = "Condition 1" ) 
  print(p)
  dev.off()
  }

file.i <- paste(out.dir,"/C2_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(sd.SR.plot, aes(x = sig_cond2, y = X_0 - X_N)) + geom_smooth(method=lm,col = 1,lwd = 0.5) 
p <- p + geom_point(aes(x = sig_cond2, y = X_0 - X_N, colour = d),size = 3)  +
  labs(y = expression(SR[0]-SR[N]), x = "Condition 2" ) 
print(p)
dev.off()


for(n.i in c("60","90","120")) {
  file.i <- paste(out.dir,"/C2_",constrain,"_",n.i,".pdf",sep = "")
  pdf(file.i)
  p <- ggplot(sd.SR.plot[sd.SR.plot$n == n.i,], aes(x = sig_cond2, y = X_0 - X_N)) + geom_smooth(method=lm,col = 1,lwd = 0.5) 
  p <- p + geom_point(aes(x = sig_cond2, y = X_0 - X_N, colour = d),size = 3)  +
    labs(y = expression(SR[0]-SR[N]), x = "Condition 2" ) 
  print(p)
  dev.off()
}


############### MIXED STRATEGY ########################

# summarize the mixed strategy over each singe portfolio
sd.SR.plot <- ds.SR


ds.plot <- with(sd.SR.plot, data.frame(SR_diff = X_pi - X_p, Type = "MV" , d = d, n = n )  )
ds.plot <- rbind(ds.plot,with(sd.SR.plot, data.frame(SR_diff = X_pi - X_0, Type = "GMV", d = d, n = n )  ) )
ds.plot <- rbind(ds.plot,with(sd.SR.plot, data.frame(SR_diff = X_pi - X_N, Type = "Naive", d = d, n = n )  ) )
ds.plot$Type <- as.factor(ds.plot$Type)
ds.plot$n <- as.factor(ds.plot$n)


file.i <- paste(out.dir,"/SR_diff1_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(ds.plot,aes(x = factor(Type), y = SR_diff)) 
p <- p +  geom_violin() 
p <- p + geom_jitter(aes(colour = d),size = 3,height = 0, width = 0.1) 
p <- p + geom_abline(intercept = 0, slope = 0, linetype ="dashed")
p <- p +  labs(x = "Portfolio", y = "Difference in SR")
print(p)
dev.off()

file.i <- paste(out.dir,"/SR_diff2_",constrain,".pdf",sep = "")
pdf(file.i)
p <- ggplot(subset(ds.plot,!d  %in% c("50B","50S") ),aes(x = factor(Type), y = SR_diff)) 
p <- p +  geom_violin() 
p <- p + geom_jitter(aes(colour = d),size = 3,height = 0, width = 0.1) 
p <- p + geom_abline(intercept = 0, slope = 0, linetype ="dashed")
p <- p +  labs(x = "Portfolio", y = "Difference in SR")
print(p)
dev.off()


for (n.i in c("60","90","120")) {
  file.i <- paste(out.dir,"/SR_diff_FF",constrain,"_",n.i,".pdf",sep = "")
  pdf(file.i)
  p <- ggplot(ds.plot[(ds.plot$n == n.i) & (!ds.plot$d  %in% c("50B","50S")),],aes(x = factor(Type), y = SR_diff)) 
  p <- p +  geom_violin() 
  p <- p + geom_jitter(aes(colour = d),size = 3,height = 0, width = 0.1) 
  p <- p + geom_abline(intercept = 0, slope = 0, linetype ="dashed")
  p <- p +  labs(x = "Portfolio", y = "Difference in SR")
  print(p)
  dev.off()
}



ds.SR <- dlply(ds.SR,"d",data.frame )

output.table <- function(i) {
  v <- print(xtable(ds.SR[[i]][,-2],digits = c(0,0,3,3,3,3)), include.rownames = F)
  v <- sapply(strsplit(v,"\n"),function(x)  x[9:(9+nrow(ds.SR[[i]])-1)] )
  cat(v, file = paste(output.dir,"/SR_",i,"_",short*1,".txt",sep = ""))
}

lapply(1:length(ds.SR), output.table)
