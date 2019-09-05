## ----setup, set.seed(2), include=FALSE-----------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(2)

## ------------------------------------------------------------------------
set.seed(2)
library(corrcoverage)

## ------------------------------------------------------------------------
# simulate reference haplotypes

nsnps <- 200
nhaps <- 1000
lag <- 6 # genotypes are correlated between neighbouring variants
maf <- runif(nsnps+lag,0.05,0.5) # common SNPs

laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
haps <- laghaps[,1:nsnps]
for(j in 1:lag) 
    haps <- haps + laghaps[,(1:nsnps)+j]
haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
snps <- colnames(haps) <- paste0("s",1:nsnps)
freq <- as.data.frame(haps+1)
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)
MAF <- colMeans(freq[,snps]-1)

# SNP correlation matrix

LD <- cor2(haps)

## ------------------------------------------------------------------------
CV <- sample(snps[which(colMeans(haps)>0.1)],1)
iCV <- sub("s", "", CV) # index of cv
OR <- 1.1
N0 <- 5000 # cases
N1 <- 5000 # controls

## ------------------------------------------------------------------------
z0_file <- system.file('extdata', 'z-scores.RDS', package='corrcoverage')
z0 <- readRDS(z0_file)

length(z0)
z0[1:5]

## ----eval = FALSE--------------------------------------------------------
#  library(simGWAS)
#  z0 <- simulated_z_score(N0=5000, # number of controls
#                N1=5000, # number of cases
#                snps=snps, # column names in freq of SNPs
#                W=CV, # causal variants, subset of snps
#                gamma.W=log(OR), # log odds ratios
#                freq=freq # reference haplotypes
#                )

## ----warning=FALSE, error=FALSE, fig.width=4-----------------------------
varbeta <- Var.data.cc(f = MAF, N = N1+N0, s = N1/(N0+N1)) # variance of 
                                                         # estimated effect size

postprobs <- ppfunc(z = z0, V = varbeta)

plot(postprobs, xlab = "SNP index", ylab = "Posterior probability")
abline(v = iCV, col = 2)

## ----warning=FALSE, error=FALSE------------------------------------------
thresh <- 0.9

credibleset <- credset(pp = postprobs, thr = thresh)

credibleset

## ----warning=FALSE, error=FALSE------------------------------------------
# is the true CV contained within the credible set?

thresh <- 0.9

credibleset <- credset(pp = postprobs, CV = iCV, thr = thresh)

credibleset

## ------------------------------------------------------------------------
corrected_cov_estimate <- corrcov(z = z0, f = MAF, N0, N1, Sigma = LD, thr = thresh)

## ----echo = FALSE--------------------------------------------------------
df <- data.frame("Corrected Coverage" = corrected_cov_estimate, 
                 "Claimed Coverage" = credibleset$claimed.cov, "Threshold" = thresh)
print(df, row.names = FALSE)

## ------------------------------------------------------------------------
z0_rep_file <- system.file('extdata', 'rep-z-scores.RDS', package='corrcoverage')
z0.rep <- readRDS(z0_rep_file)

## ----eval = FALSE--------------------------------------------------------
#  z0.rep <- simulated_z_score(N0=5000, # number of controls
#                              N1=5000, # number of cases
#                              snps=snps, # column names in freq
#                              W=CV, # causal variants, subset of snps
#                              gamma.W=log(OR), # log odds ratios
#                              freq=freq, # reference haplotypes
#                              nrep = 1000 # 1000 simulations
#  )

## ------------------------------------------------------------------------
pps <- ppfunc.mat(zstar = z0.rep, V = varbeta) # find pps
cs <- apply(pps, 1, function(x) credset(x, CV = iCV, thr = thresh)$cov)
true.cov.est <- mean(cs)
true.cov.est

## ----echo = FALSE--------------------------------------------------------
df.new <- data.frame("Empirical Coverage" = true.cov.est, 
                     "Corrected Coverage" = corrected_cov_estimate, "Claimed Coverage" = credibleset$claimed.cov, "Threshold" = thresh)[1,]
print(df.new, row.names = FALSE)

