# Revised two-component residence time fitting with robust convergence
# Verified for syntax and several logic fixes applied.
# Requires: minpack.lm, ggplot2, dplyr

suppressPackageStartupMessages({
  library(minpack.lm)
  library(ggplot2)
  library(dplyr)
})

# ================= CONFIG =================
FILE_A <- "2-12-26- SRCAP HALO - 25PM.csv"
FILE_B <- "4-20-26- NLS SNAP - 25PM.csv"

COL_TRAJECTORY <- "Trajectory"
COL_FRAME <- "Frame"

FRAME_INTERVAL <- 0.2   # verify this matches acquisition interval
MIN_TRAJ_FRAMES <- 2

LABEL_A <- "SRCAP-HALO"
LABEL_B <- "NLS-SNAP"

B <- 2000
T_MAX <- 100

# ================= DATA =================
load_times <- function(path){

  df <- read.csv(path, stringsAsFactors=FALSE, row.names=1)

  frames <- as.numeric(df[[COL_FRAME]])
  traj <- df[[COL_TRAJECTORY]]

  lens <- tapply(
    frames,
    traj,
    function(x) max(x)-min(x)+1
  )

  lens <- lens[lens >= MIN_TRAJ_FRAMES]
  times <- as.numeric(lens)*FRAME_INTERVAL
  times <- times[times>0]

  cat("\n",path,"\n")
  cat("n trajectories:",length(times),"\n")
  cat("unique dwell times:",length(unique(times)),"\n")
  cat("range:",range(times),"\n\n")

  times
}

make_survival <- function(times){
  t <- sort(times)
  n <- length(t)
  data.frame(
    t=t,
    surv=1-(seq_len(n)/n)
  )
}

# ================= MODELS =================
pred2 <- function(t,a,ks,kns){
  a*exp(-ks*t)+(1-a)*exp(-kns*t)
}

fit_single <- function(surv_df){

  starts <- c(.001,.005,.01,.03,.1)
  bestfit <- NULL
  bestrss <- Inf

  for(k0 in starts){

    fit <- try(
      nlsLM(
        surv ~ exp(-k*t),
        data=surv_df,
        start=list(k=k0),
        lower=.000001,
        control=nls.lm.control(maxiter=500)
      ),
      silent=TRUE
    )

    if(inherits(fit,"try-error")) next

    rss <- sum(residuals(fit)^2)

    if(rss < bestrss){
      bestrss <- rss
      bestfit <- fit
    }
  }

  if(is.null(bestfit)) return(NULL)

  list(
    model="single",
    coef=c(
      alpha=1,
      ks=unname(coef(bestfit)["k"]),
      kns=unname(coef(bestfit)["k"])
    ),
    rss=bestrss
  )
}

estimate_starts <- function(surv_df){

  t <- surv_df$t
  s <- surv_df$surv

  alpha_grid <- seq(.1,.9,.2)
  ks_grid <- c(.001,.003,.005,.01,.03,.05,.1,.2)
  kns_grid <- c(.03,.05,.1,.2,.5,1,2)

  best <- NULL
  best_rss <- Inf

  for(a in alpha_grid){
    for(ks in ks_grid){
      for(kns in kns_grid){

        if(ks >= kns) next

        rss <- sum((s-pred2(t,a,ks,kns))^2)

        if(rss < best_rss){
          best_rss <- rss
          best <- c(alpha=a,ks=ks,kns=kns)
        }
      }
    }
  }

  best
}

fit_two_component <- function(surv_df){

  starts <- estimate_starts(surv_df)

  candidate_starts <- list(
    starts,
    c(alpha=.3,ks=.005,kns=.1),
    c(alpha=.5,ks=.01,kns=.3),
    c(alpha=.7,ks=.001,kns=.05)
  )

  bestfit <- NULL
  bestrss <- Inf

  for(st in candidate_starts){

    fit <- try(
      nlsLM(
        surv ~ alpha*exp(-ks*t)+(1-alpha)*exp(-kns*t),
        data=surv_df,
        start=as.list(st),
        lower=c(alpha=.001,ks=.000001,kns=.000001),
        upper=c(alpha=.999,ks=5,kns=10),
        control=nls.lm.control(maxiter=1000)
      ),
      silent=TRUE
    )

    if(inherits(fit,"try-error")) next

    cf <- coef(fit)

    if(abs(cf["ks"]-cf["kns"]) < .005) next

    rss <- sum(residuals(fit)^2)

    if(rss < bestrss){
      bestrss <- rss
      bestfit <- fit
    }
  }

  singlefit <- fit_single(surv_df)

  if(is.null(bestfit) && !is.null(singlefit))
    return(singlefit)

  if(is.null(bestfit)){
    return(list(
      model="failed",
      coef=c(alpha=NA,ks=NA,kns=NA),
      rss=NA
    ))
  }

  two <- list(
    model="two",
    coef=coef(bestfit),
    rss=bestrss
  )

  if(!is.null(singlefit)){
    n <- nrow(surv_df)

    # protect against log(0)
    rss2 <- max(two$rss,1e-12)
    rss1 <- max(singlefit$rss,1e-12)

    bic_two <- n*log(rss2/n)+3*log(n)
    bic_one <- n*log(rss1/n)+1*log(n)

    if(bic_one < bic_two)
      return(singlefit)
  }

  two
}

# ================= BOOTSTRAP =================
bootstrap_fit <- function(times,B=2000){

  n <- length(times)

  out <- matrix(
    NA,
    B,
    3,
    dimnames=list(NULL,c("alpha","ks","kns"))
  )

  cat("Running",B,"bootstrap replicates\n")
  pb <- txtProgressBar(min=0,max=B,style=3)

  fail <- 0

  for(i in seq_len(B)){

    samp <- sample(times,n,replace=TRUE)

    fit <- fit_two_component(
      make_survival(samp)
    )

    out[i,] <- fit$coef

    if(any(is.na(fit$coef))) fail <- fail+1

    setTxtProgressBar(pb,i)
  }

  close(pb)
  cat("\nFailures:",fail,"/",B,"\n")

  out
}

curve_ci <- function(boot,tgrid){

  boot <- boot[
    complete.cases(boot),
    ,
    drop=FALSE
  ]

  if(nrow(boot)==0) stop("No successful bootstrap fits.")

  curves <- sapply(
    seq_len(nrow(boot)),
    function(i){
      p <- boot[i,]
      pred2(tgrid,p[1],p[2],p[3])
    }
  )

  if(is.vector(curves)) curves <- matrix(curves,nrow=length(tgrid))

  data.frame(
    t=tgrid,
    lo=apply(curves,1,quantile,.025),
    hi=apply(curves,1,quantile,.975),
    med=apply(curves,1,median)
  )
}

# ================= RUN =================
times_A <- load_times(FILE_A)
times_B <- load_times(FILE_B)

fit_A <- fit_two_component(make_survival(times_A))
fit_B <- fit_two_component(make_survival(times_B))

cat("\nObserved fits\n")
print(fit_A$model)
print(fit_A$coef)

print(fit_B$model)
print(fit_B$coef)

set.seed(42)
boot_A <- bootstrap_fit(times_A,B)
boot_B <- bootstrap_fit(times_B,B)

cat("\n95% CIs\n")
for(p in c("alpha","ks","kns")){

  c1 <- quantile(
    boot_A[,p],
    c(.025,.975),
    na.rm=TRUE
  )

  c2 <- quantile(
    boot_B[,p],
    c(.025,.975),
    na.rm=TRUE
  )

  cat(
    p,"\n",
    LABEL_A,c1,"\n",
    LABEL_B,c2,"\n\n"
  )
}

for(p in c("alpha","ks","kns")){

  d <- boot_A[,p]-boot_B[,p]

  ci <- quantile(
    d,
    c(.025,.975),
    na.rm=TRUE
  )

  pval <- 2*min(
    mean(d<0,na.rm=TRUE),
    mean(d>0,na.rm=TRUE)
  )

  cat(
    p,
    "delta CI:",
    ci,
    "p=",
    pval,
    "\n"
  )
}

grid <- seq(0,T_MAX,length.out=500)

ciA <- curve_ci(boot_A,grid)
ciB <- curve_ci(boot_B,grid)

obsA <- data.frame(
 t=grid,
 y=pred2(
   grid,
   fit_A$coef[1],
   fit_A$coef[2],
   fit_A$coef[3]
 ),
 condition=LABEL_A
)

obsB <- data.frame(
 t=grid,
 y=pred2(
   grid,
   fit_B$coef[1],
   fit_B$coef[2],
   fit_B$coef[3]
 ),
 condition=LABEL_B
)

ggplot()+
 geom_ribbon(
   data=ciA,
   aes(t,ymin=lo,ymax=hi),
   alpha=.2
 )+
 geom_ribbon(
   data=ciB,
   aes(t,ymin=lo,ymax=hi),
   alpha=.2
 )+
 geom_line(
   data=obsA,
   aes(t,y)
 )+
 geom_line(
   data=obsB,
   aes(t,y)
 )+
 theme_classic()
