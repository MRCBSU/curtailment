#' @import gridExtra
#' @import ggplot2
#' @import grid
#' @import ggthemes

####
#### These functions are used to run the code that finds the designs used.
#### The main function, used to find designs, is findSCdes.
####

# rmDominatedDesigns <- function(df, essh0="EssH0", essh1="Ess", n="n"){
#     discard <- rep(NA, nrow(df))
#     if("tbl_df" %in% class(df)){
#         essh0.vec <- df[[essh0]]
#         essh1.vec <- df[[essh1]]
#         n.vec <- df[[n]]
#       for(i in 1:nrow(df)){
#         discard[i] <-  any(essh0.vec[i] > essh0.vec & essh1.vec[i] > essh1.vec & n.vec[i] >= n.vec)
#       }
#     } else {
#         essh0.vec <- df[, essh0]
#         essh1.vec <- df[, essh1]
#         n.vec <- df[, n]
#         for(i in 1:nrow(df)){
#           discard[i] <- any(essh0.vec[i] > essh0.vec & essh1.vec[i] > essh1.vec & n.vec[i] >= n.vec)
#         }
#       }
#     newdf <- df[discard==FALSE,,drop=FALSE]
#     newdf
# }


find2armBlockOCs <- function(n,r, Bsize, mat, thetaF, thetaE, power, alpha, pat.cols, prob.vec, prob.vec.p0, blank.mat, zero.mat){
######################## UPDATE CP MATRIX USING thetaF/1 VALUES:
for(i in (n+r+1):1){
  for(j in pat.cols){  # Only look every Bsize patients (no need to look at final col)
    if(i-1<=j){ # Condition: Sm<=m
      newcp <- sum(prob.vec*mat[i:(i+Bsize), j+Bsize])
      if(newcp > thetaE) mat[i,j] <- 1
      if(newcp < thetaF) mat[i,j] <- 0
      if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
    }
  }
}


###### STOP if design is pointless, i.e either failure or success is not possible:


# IF DESIGN GUARANTEES FAILURE (==0) or SUCCESS (==2) at n=C:
first.cohort <- sum(mat[,Bsize], na.rm = T)

if(first.cohort==Bsize+1){
  return(c(n, r, Bsize, 1, 1,  NA, NA, thetaF, thetaE, NA))
}

if(first.cohort==0){
  return(c(n, r, Bsize, 0, 0,  NA, NA, thetaF, thetaE, NA))
}




########################### FIND PROB. OF REACHING EACH POINT:
################# START WITH AN INDICATOR MATRIX OF NON-TERMINAL POINTS:

tp.mat <- blank.mat
tp.mat[which(mat==0 | mat==1)] <- 0
tp.mat[which(mat>0 & mat<1)] <- 1


############ CREATE MATRIX OF "POSSIBLE POINTS" -- IE, LIKE MATRIX OF NON-TERMINAL POINTS ***PLUS*** THE TERMINAL POINTS:
##### THIS WILL BE USED AS AN INDICATOR MATRIX FOR WHICH POINTS TO CALCULATE THE PROB'Y OF REACHING.

# Start with non-terminal points and add the terminal points
poss.mat <- tp.mat
poss.mat[1:(Bsize+1), Bsize] <- 1 # It is of course possible to reach 0,1,...Bsize after Bsize patients. This line is included in case CP=0 or CP=1 at first check (i.e. m=Bsize)


# Failures first:
fail.mat <- zero.mat

rows.with.cp0 <- which(apply(mat, 1, function(x) {any(x==0, na.rm = T)}))
fail.n <- apply(mat[rows.with.cp0,], 1, which.min)

for(i in 1:length(fail.n)){
  poss.mat[names(fail.n)[i],fail.n[i]] <- 1
  fail.mat[names(fail.n)[i],fail.n[i]] <- 1
}

# Now successes: what are the successful terminal points?
# Points with mat[i,j]==1 AND (0>mat[i,j-2]>1 OR 0>mat[i-1,j-2]>1 OR 0>mat[i-2,j-2]>1)
success.mat <- zero.mat
rows.with.cp1 <- which(apply(mat, 1, function(x) {any(x==1, na.rm = T)}))

# browser()


for(i in rows.with.cp1){
  for(j in seq(Bsize, 2*n, by=Bsize)){
    if(i-1<=j & mat[i,j]==1 & (j==Bsize | any(tp.mat[i:max(1,(i-Bsize)), j-Bsize]==1, na.rm=TRUE))){ # max() condition to take care of cases where
        # Conditions ensure CP=1 and that it is possible to actually reach the point (tp.mat==1 indicates non terminal point)
        poss.mat[i, j] <- 1
        success.mat[i,j] <- 1
    }
  }
}


######################## PROBABILITY OF REACHING EACH POINT
############################## FIRSTLY, UNDER PT=PT

final.probs.mat <- poss.mat

# First n=Bsize rows (0,1,... Bsize-1) are special cases:
# fill in first column:
final.probs.mat[1:(Bsize+1), Bsize] <- prob.vec


for(i in 1:Bsize){
  row.index <- which(poss.mat[i,]==1)[-1] # First entry (ie first column) has been inputted already, directly above, hence [-1]
  for(j in row.index){
    final.probs.mat[i, j] <- sum(prob.vec[1:i]*final.probs.mat[i:1, j-Bsize]*tp.mat[i:1, j-Bsize])
  }
}


#if(n==16 & r==8) browser()

# For the remaining rows:
for(i in (Bsize+1):nrow(final.probs.mat)){ # Skip first Bsize rows; they have been taken care of above.
  for(j in seq(2*Bsize, 2*n, by=Bsize)){ # skip first column of patients (n=2) -- again, they have been taken care of above.
    if(i-1<=j & poss.mat[i,j]==1){
      final.probs.mat[i,j] <- sum(prob.vec*final.probs.mat[i:(i-Bsize), j-Bsize]*tp.mat[i:(i-Bsize), j-Bsize], na.rm = TRUE)
    }
  }
}


# IMPORTANT: end early if pwr < power.
# Note 2: Could comment this out to ensure that the final test in undertaken for all runs, not just feasible ones.
prob.success <- final.probs.mat[success.mat==1]
pwr <- sum(prob.success)

if(pwr < power) # | pwr < power+tol )
{
  return(c(n, r, Bsize, NA, pwr,  NA, NA, thetaF, thetaE, NA))
}


############################## SECONDLY, UNDER PT=PC

final.probs.mat.p0 <- poss.mat

# First n=Bsize rows (0,1,... Bsize-1) are special cases:
# fill in first column:
final.probs.mat.p0[1:(Bsize+1), Bsize] <- prob.vec.p0

for(i in 1:Bsize){
  row.index <- which(poss.mat[i,]==1)[-1] # First entry (ie first column) has been inputted already, directly above, hence [-1]
  for(j in row.index){
    final.probs.mat.p0[i, j] <- sum(prob.vec.p0[1:i]*final.probs.mat.p0[i:1, j-Bsize]*tp.mat[i:1, j-Bsize])
  }
}

# For the remaining rows:
for(i in (Bsize+1):nrow(final.probs.mat.p0)){ # Skip first Bsize rows; they have been taken care of above.
  for(j in seq(2*Bsize, 2*n, by=Bsize)){ # skip first column of patients (n=2) -- again, they have beene taken care of above.
    if(i-1<=j & poss.mat[i,j]==1){
      final.probs.mat.p0[i,j] <- sum(prob.vec.p0*final.probs.mat.p0[i:(i-Bsize), j-Bsize]*tp.mat[i:(i-Bsize), j-Bsize], na.rm = TRUE)
    }
  }
}


prob.success.p0 <- final.probs.mat.p0[success.mat==1]
typeIerr <- sum(prob.success.p0)

# IMPORTANT: end early if type I error > alpha
# Note 2: Could comment this out to ensure that the final test in undertaken for all runs, not just feasible ones.
if( typeIerr > alpha) # | alpha < alpha-tol)
{
  return(c(n, r, Bsize, typeIerr, pwr,  NA, NA, thetaF, thetaE, NA))
}



########################## ESS FOR SUCCESS POINTS
success.n <- which(success.mat==1, arr.ind = T)[,"col"] # Note: this is potentially dangerous if this and final.probs.mat[success.mat==1] are found in a different order, though this shouldn't be the case.

success.df <- data.frame(prob=prob.success, prob.p0=prob.success.p0, ess=prob.success*success.n, essH0=prob.success.p0*success.n)


################## ESS FOR FAILURE POINTS
fail.n <- which(fail.mat==1, arr.ind = T)[,"col"] # Note: this is potentially dangerous if this and final.probs.mat[success.mat==1] are found in a different order, though this shouldn't be the case.

prob.fail <- final.probs.mat[fail.mat==1]
prob.fail.p0 <- final.probs.mat.p0[fail.mat==1]

fail.df <- data.frame(prob=prob.fail, prob.p0=prob.fail.p0, ess=prob.fail*fail.n, essH0=prob.fail.p0*fail.n)

all.df <- rbind(fail.df, success.df)

###### CHECK PROBS ALL SUM TO 1
# sum(all.df$prob)

if(sum(all.df$prob)+sum(all.df$prob.p0)-2 > 1e-8)  stop("Total probability of failure + success =/= 1. Something has gone wrong." , call. = FALSE)

ess <- sum(all.df$ess)
essH0 <- sum(all.df$essH0)



############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops

cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
possible.cps <- apply(mat, 2, function(x) {sum(!is.na(x))})

effective.n <- min(which(cp.colsums==possible.cps))
return(data.frame(n=n, r=r, Bsize=Bsize, typeIerr=typeIerr, pwr=pwr, EssH0=essH0, Ess=ess, thetaF=thetaF, thetaE=thetaE, eff.n=effective.n))
}


################ Function for finding the Prob(reponses on treatment + non-responses on control)=0, 1, 2,... Bsize:
findProbVec <- function(Bsize, pt=pt, qt=qt, pc=pc, qc=qc){
  prob.vec <- rep(NA, Bsize+1)
  for(i in 1:(Bsize+1)){
    positives <- i-1
    full.vec <- expand.grid(rep(list(0:1), Bsize))
    positive.mat <- full.vec[rowSums(full.vec) == positives,]
    negative.mat <- -1*(positive.mat-1)

    positive.vec <- rep(c(pt,qc), each=Bsize/2)
    negative.vec <- rep(c(qt,pc), each=Bsize/2)

    posneg.mat <- t(t(positive.mat)*positive.vec) + t(t(negative.mat)*negative.vec)
    prob.vec[i] <- sum(apply(posneg.mat, 1, prod))
  }
  if(sum(prob.vec)-1 > 1e-8) stop("Probabilities do not sum to 1.")
  prob.vec
}

################ Function for finding the uncurtailed CP matrix:
findBlock2armUncurtailedMatrix <- function(n, r, Bsize, pat.cols, prob.vec){

  cpmat <- matrix(3, ncol=2*n, nrow=min(n+r+Bsize+2, 2*n+1))
  rownames(cpmat) <- 0:(nrow(cpmat)-1)
  cpmat[(n+r+2):nrow(cpmat),] <- 1
  cpmat[1:(n+r+1),2*n] <- 0 # Fail at end

  for(i in (n+r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        cpmat[i,j] <- ifelse(test=j-(i-1) >= n-r+1, yes=0, no=sum(prob.vec*cpmat[i:(i+Bsize), j+Bsize]))
        # IF success is not possible (i.e. [total no. of pats-Xa+Ya-Xb] >= n-r+1), THEN set CP to zero. Otherwise, calculate it based on "future" CPs.
      }
    }
  }

  for(i in 3:nrow(cpmat)){
    cpmat[i, 1:(i-2)] <- NA
  }
  cpmat
}

#' Find two-arm trial designs that use stochastic curtailment
#'
#' This function finds admissible design realisations for two-arm binary outcome trials, using stochastic curtailment.
#' The output can be used as the sole argument in the function 'drawDiagram', which will return the stopping boundaries for the
#' admissible design of your choice. Monitoring frequency can set in terms of block size.
#' @param nmin.arm Minimum permitted sample size *per arm*. Should be a multiple of block size.
#' @param nmax.arm Maximum permitted sample size *per arm*. Should be a multiple of block size.
#' @param block.size Block size.
#' @param pc Anticipated response rate on the control arm.
#' @param pt Anticipated response rate on the treatment arm.
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max.
#' @param minthetaE Minimum value of upper threshold theta_E_min.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param rm.dominated.designs Logical. If TRUE, dominated designs will be
#' removed from final output. Defaults to TRUE.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param fast.method Logical. If FALSE, design search is conducted over all combinations of
#' (theta_F, theta_E). If TRUE, a much faster, though less thorough, design search is undertaken.
#' Defaults to FALSE.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return Output is a list of two dataframes. The first, $input, is a one-row data frame that contains all the arguments used in the call.
#' The second, $all.des, contains the operating characteristics of all admissible designs found.
#' @examples
#' \donttest{
#' des <- twoarmDesign(nmin.arm=20,
#' nmax.arm=24,
#' block.size=8,
#' pc=0.1,
#' pt=0.4,
#' alpha=0.1,
#' power=0.8,
#' maxthetaF=0.4,
#' minthetaE=0.7,
#' max.combns=1e4)
#' }
#' @export
twoarmDesign <- function(nmin.arm,
                            nmax.arm,
                            block.size,
                            pc,
                            pt,
                            alpha,
                            power,
                            maxthetaF=NULL,
                            minthetaE=0.7,
                            bounds="ahern",
                            fixed.r=NULL,
                            max.combns=1e6,
                            rm.dominated.designs=TRUE,
                            exact.thetaF=NULL,
                            exact.thetaE=NULL,
                            fast.method=FALSE)
{
  Bsize <- block.size

  if(Bsize%%2!=0) stop("Block size must be an even number")

  if((2*nmin.arm)%%Bsize!=0) stop("2*nmin.arm must be a multiple of block size")
  if((2*nmax.arm)%%Bsize!=0) stop("2*nmax.arm must be a multiple of block size")

  nposs <- seq(from=nmin.arm, to=nmax.arm, by=Bsize/2)

  qc <- 1-pc
  qt <- 1-pt





  prob.vec <- findProbVec(Bsize=Bsize, pt=pt, qt=qt, pc=pc, qc=qc)
  prob.vec.p0 <- findProbVec(Bsize=Bsize, pt=pc, qt=qc, pc=pc, qc=qc)

  pat.cols.list <- lapply(nposs, function(x) seq(from=2*x, to=Bsize, by=-Bsize)[-1])
  names(pat.cols.list) <- nposs

  if(is.null(maxthetaF)){
    maxthetaF <- pt
  }

  r.list <- list()
  for(i in 1:length(nposs))
  {
    r.list[[i]] <- 0:(nposs[i]-2) # r values: 0 to nposs[i]-2
  }

  ns <- NULL

  for(i in 1:length(nposs)){
    ns <- c(ns, rep(nposs[i], length(r.list[[i]])))
  }

  sc.subset <- data.frame(n=ns, r=unlist(r.list))

  if(!is.null(bounds)){
    # Incorporate A'Hern's bounds:
    if(bounds=="ahern")  {
      #sc.subset <- sc.subset[sc.subset$r >= pc*sc.subset$n & sc.subset$r <= pt*sc.subset$n, ] # One-arm case
      sc.subset <- sc.subset[sc.subset$r >= 1 & sc.subset$r <= pt*sc.subset$n, ] # Try this for two-arm case -- interval [1, pt*Narm]
    }

    if(bounds=="wald"){
      # Even better to incorporate Wald's bounds:
      denom <- log(pt/pc) - log((1-pt)/(1-pc))
      accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-pc)/(1-pt))/denom
      accept.null <- floor(accept.null)

      reject.null <- log((power)/alpha) / denom  + nposs * log((1-pc)/(1-pt))/denom
      reject.null <- ceiling(reject.null)

      r.wald <- NULL
      ns.wald <- NULL
      for(i in 1:length(nposs)){
        r.wald <- c(r.wald, accept.null[i]:reject.null[i])
        ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
      }
      sc.subset <- data.frame(n=ns.wald, r=r.wald)
      sc.subset <- sc.subset[sc.subset$n - sc.subset$r >=2, ]
    }
  }

  # In case you want to specify values for r:
  if(!is.null(fixed.r))  {
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }

  ###### Find thetas for each possible {r, N} combn:


  mat.list <- vector("list", nrow(sc.subset))
  for(i in 1:nrow(sc.subset)){
    mat.list[[i]] <- findBlock2armUncurtailedMatrix(n=sc.subset[i,"n"], r=sc.subset[i,"r"], Bsize=Bsize, pat.cols=pat.cols.list[[paste(sc.subset$n[i])]], prob.vec=prob.vec)
  }

  store.all.thetas <- lapply(mat.list, function(x) {sort(unique(c(x))[unique(c(x)) <= 1])})


  ##### To cut down on computation, try cutting down the number of thetas used:
  ##### max.combns:=max. number of (thetaF, thetaE) combinations.
  ##### n.thetas*(n.thetas-1)/2 = n.combns, so if n.thetas > sqrt(2*max.combns), take out every other value, excluding 0 and 1.
  ##### Note: further below, more combns are removed if constraints on maxthetaF and minthetaE are specified.
  # check ####
  if(max.combns!=Inf){
    maxthetas <- sqrt(2*max.combns)
    for(i in 1:nrow(sc.subset))
    {
      while(length(store.all.thetas[[i]]) > maxthetas)
      {
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
  }

  if(!is.null(exact.thetaF) & !is.null(exact.thetaE)){ # if exact thetas are given (to speed up a result check):
    for(i in 1:length(store.all.thetas)){
      keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-3 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-3
      store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
    }
  }

  h.results.list <- vector("list", nrow(sc.subset)) #

  pb <- txtProgressBar(min = 0, max = nrow(sc.subset), style = 3)

  # Now, find the designs, looping over each possible {r, N} combination, and within each {r, N} combination, loop over all combns of {thetaF, thetaE}:
  for(h in 1:nrow(sc.subset)){
    #k <- 1
    blank.mat <- matrix(NA, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(blank.mat) <- 0:(nrow(blank.mat)-1)
    zero.mat <- matrix(0, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(zero.mat) <- rownames(blank.mat)
    pat.cols.single <- pat.cols.list[[paste(sc.subset$n[h])]]

    if(fast.method==TRUE){
      h.results <- fastSearch(thetas=store.all.thetas[[h]],
                              maxthetaF=maxthetaF,
                              minthetaE=minthetaE,
                              sc.h=sc.subset[h,],
                              Bsize=Bsize,
                              mat.h=mat.list[[h]],
                              blank.mat=blank.mat,
                              zero.mat=zero.mat,
                              power=power,
                              alpha=alpha,
                              pat.cols.single=pat.cols.single,
                              prob.vec=prob.vec,
                              prob.vec.p0=prob.vec.p0)
    }else{
      h.results <- slowSearch(thetas=store.all.thetas[[h]],
                              maxthetaF=maxthetaF,
                              minthetaE=minthetaE,
                              sc.h=sc.subset[h,],
                              Bsize=Bsize,
                              mat.h=mat.list[[h]],
                              blank.mat=blank.mat,
                              zero.mat=zero.mat,
                              power=power,
                              alpha=alpha,
                              pat.cols.single=pat.cols.single,
                              prob.vec=prob.vec,
                              prob.vec.p0=prob.vec.p0)
    }


    setTxtProgressBar(pb, h)
    h.results.df <- do.call(rbind, h.results)

    if(!is.null(h.results.df)){
      # Remove all "skipped" results:
      colnames(h.results.df) <- c("n", "r", "block", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n")
      h.results.df <- h.results.df[!is.na(h.results.df[, "Ess"]),]
      if(nrow(h.results.df)>0){
        # Remove dominated designs:
        if(rm.dominated.designs==TRUE){
          discard <- rep(NA, nrow(h.results.df))
          for(i in 1:nrow(h.results.df)){
            discard[i] <- sum(h.results.df[i, "EssH0"] > h.results.df[, "EssH0"] & h.results.df[i, "Ess"] > h.results.df[, "Ess"] & h.results.df[i, "n"] >= h.results.df[, "n"])
            }
          h.results.df <- h.results.df[discard==0,, drop=FALSE]
        }
        # Remove duplicated designs:
        if(is.matrix(h.results.df)){ # i.e. if there is more than one design (if not, h.results.df is a vector)
        duplicates <- duplicated(h.results.df[, c("n", "Ess", "EssH0"), drop=FALSE])
        h.results.df <- h.results.df[!duplicates,, drop=FALSE]
        }
        h.results.list[[h]] <- h.results.df
      }
    }
  } # End of "h" loop


  full.results <- do.call(rbind, h.results.list)
  #if(length(full.results)==0) stop("There are no feasible designs for this combination of design parameters" , call. = FALSE)
  if(length(full.results)>0){
    # Discard all "inferior" designs:
    discard <- rep(NA, nrow(full.results))
    for(i in 1:nrow(full.results)){
      discard[i] <- sum(full.results[i, "EssH0"] > full.results[, "EssH0"] & full.results[i, "Ess"] > full.results[, "Ess"] & full.results[i, "n"] >= full.results[, "n"])
      #print(i)
    }
    subset.results <- full.results[discard==0,,drop=FALSE]


    # Remove duplicates:
    duplicates <- duplicated(subset.results[, c("n", "EssH0", "Ess"), drop=FALSE])
    all.des <- subset.results[!duplicates,,drop=FALSE]
    all.des$stage <- all.des[,"eff.n"]/all.des[,"block"]
    names(all.des)[names(all.des)=="n"] <- "n.arm"
    names(all.des)[names(all.des)=="eff.n"] <- "n"
    input <- data.frame(nmin.arm=nmin.arm,
                        nmax.arm=nmax.arm,
                        block=block.size,
                        pc=pc,
                        pt=pt,
                        alpha=alpha,
                        power=power,
                        maxthetaF=maxthetaF,
                        minthetaE=minthetaE,
                        bounds=bounds,
                        max.combns=max.combns,
                        fast.method=fast.method)
    final.output <- list(input=input,
                         all.des=all.des)
    class(final.output) <- append(class(final.output), "curtailment_twoarm")
    return(final.output)
  }
}


findN1N2R1R2twoarm <- function(nmin, nmax, e1=FALSE){
    nposs <- nmin:nmax
    n1.list <- list()
    n2.list <- list()

    for(i in 1:length(nposs)){
      n1.list[[i]] <- 1:(nposs[i]-1)
      n2.list[[i]] <- nposs[i]-n1.list[[i]]
    }

    # All possibilities together:
    n1 <- rev(unlist(n1.list))
    n2 <- rev(unlist(n2.list))
    n <- n1 + n2
    ns <- cbind(n1, n2, n)

    ################################ FIND COMBNS OF R1 AND R ###############################

    r1.list <- vector("list")
    ns.list <- vector("list")
    for(i in 1:nrow(ns)){
      r1.list[[i]] <- -n1[i]:n1[i] # r1 values: -n1 to n1, for each possible n1
      #ns.list[[i]] <-
    }

    rownames(ns) <- 1:nrow(ns)
    ns <- ns[rep(row.names(ns), sapply(r1.list, length)), ] # duplicate each row so that there are sufficient rows for each r1 value
    ns <- cbind(ns, unlist(r1.list))
    colnames(ns)[4] <- "r1"

    ######### Add possible r values:
    r.list1 <- apply(ns, 1, function(x) {(x["r1"]-x["n2"]):x["n"]})  # r must satisfy r > r1 and r < n. Also, number of responses required in stage 2 (r2-r1) must be at most n2

    how.many.rs <- sapply(r.list1, length)

    row.names(ns) <- 1:nrow(ns)
    ns <- ns[rep(row.names(ns), how.many.rs), ] # duplicate each row a certain number of times
    ns <- cbind(ns, unlist(r.list1))
    colnames(ns)[5] <- "r2"

    ### Finally, add e1 for stopping for benefit:

    if(e1==TRUE)
    {

    } else {
      rownames(ns) <- 1:nrow(ns)
      ns <- data.frame(ns)
    }

    return(ns)
}


slowSearch <- function(thetas,
                       maxthetaF,
                       minthetaE,
                       sc.h,
                       Bsize,
                       mat.h,
                       blank.mat,
                       zero.mat,
                       power,
                       alpha,
                       pat.cols.single,
                       prob.vec,
                       prob.vec.p0
                       ){
  k <- 1
  h.results <- vector("list", ceiling(0.5*length(thetas)^2))
  all.thetas <- rev(thetas)[-length(thetas)] # All thetaE values, decreasing, not including the final value, thetaE=0.
  all.thetas <- all.thetas[all.thetas>=minthetaE]
  for(i in 1:length(all.thetas)){ # For each thetaE,
    thetaFs.current <- thetas[thetas<all.thetas[i] & thetas<=maxthetaF]
    for(m in 1:length(thetaFs.current)){
      h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=as.numeric(thetaFs.current[m]), thetaE=all.thetas[i], mat=mat.h,
                                         power=power, alpha=alpha, pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
      k <- k+1
      # Add lines here: if power decreases below desired value, break:
      if(h.results[[k-1]][5] < power){
        break
      }
    }
  } # end of "i" loop
  return(h.results)
}

fastSearch <- function(thetas,
                       maxthetaF,
                       minthetaE,
                       sc.h,
                       Bsize,
                       mat.h,
                       blank.mat,
                       zero.mat,
                       power,
                       alpha,
                       pat.cols.single,
                       prob.vec,
                       prob.vec.p0){
  k <- 1
  ########### START 2D BISECTION
  thetaF.vec <- thetas[thetas<=maxthetaF]
  thetaE.vec <- thetas[thetas>=minthetaE]
  h.results <- vector("list", length(thetaF.vec)*length(thetaE.vec))
  # Bounds for thetaF:
  a0 <- 1
  b0 <- length(thetaF.vec)
  d0 <- ceiling((b0-a0)/2)
  # Bounds for thetaE:
  a1 <- 1
  b1 <- length(thetaE.vec)
  d1 <- ceiling((b1-a1)/2)
  minthetaF <- NA
  maxthetaE <- NA
  while(min((b0-a0),(b1-a1))>1 & is.na(minthetaF)){ # Break/move on when bisection method fails to find anything OR when final feasible design is found.
    output <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[d0], thetaE=thetaE.vec[d1], mat=mat.h,  power=power, alpha=alpha,
                               pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
    if(!is.na(output[6])){ # If ESS is not NA, then design IS feasible, and do:
      feasible <- TRUE
      maxthetaE <- thetaE.vec[d1]
      while((feasible==TRUE) & d0<length(thetaF.vec)){
        d0 <- d0+1
        h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[d0], thetaE=maxthetaE, mat=mat.h,  power=power, alpha=alpha,
                                           pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
        feasible <- !is.na(h.results[[k]][6])
        k <- k+1
      } # Once the final feasible design for the given thetaF/1 is found (or we reach the largest thetaF), record thetaF and make it a limit:
      minthetaF <- thetaF.vec[d0-1]
    } else { # If design isn't feasible, decrease thetaF, increase thetaE and test again:
      b0 <- d0
      a1 <- d1
      d0 <- a0 + floor((b0-a0)/2)
      d1 <- a1 + floor((b1-a1)/2)
    }
  }
  if(!is.na(minthetaF)){ # If at least one feasible design was found, then minthetaF exists, and we search over all thetaF/1 combinations subject to the new limits we have just created:
    thetaF.vec <- thetaF.vec[thetaF.vec>=minthetaF]
    thetaE.vec <- thetaE.vec[thetaE.vec<=maxthetaE]
    for(i in 1:length(thetaE.vec)){
      for(j in 1:length(thetaF.vec)){
        #  print(paste(thetaF.vec[i], thetaE.vec[j]))
        h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[j], thetaE=thetaE.vec[i], mat=mat.h,
                                           power=power, alpha=alpha, pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
        k <- k+1
      }
    }
  } # if no feasible designs found, do nothing and let loop end.
  return(h.results)
}


