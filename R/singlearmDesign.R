#' Find single-arm trials using stochastic curtailment
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
#' The output can be used as the sole argument in the function 'drawDiagram', which will return the stopping boundaries for the
#' admissible design of your choice. Monitoring frequency can set in terms of block(/cohort) size ("C") or number of stages ("stages").
#' @param n.max Maximum sample size. Can be a single value or a vector of length 2, indicating the range of maximum sample sizes to use. Value(s) should be a multiple of block size (C) or number of stages.
#' @param C Possible block size(s) to use. Vectors are permitted. Default=1, i.e. sequential monitoring. Overridden by specifying a "stages" argument.
#' @param stages Number of interim analyses or "stages" to use. Overrides block size C. Vectors are permitted.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param minstop Minimum permitted sample size at the first interim analysis
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max. Defaults to power.
#' @param minthetaE Minimum value of upper CP threshold theta_E_min. Defaults to power.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param return.only.admissible Logical. Returns only admissible design realisations if TRUE, otherwise returns all feasible designs. Defaults to TRUE.
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e4.
#' @param maxthetas Provide a maximum number of CP values used to create ordered pairs (theta_F, theta_E). Can be used instead of max.combns. Defaults to NA.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param progressBar Logical. If TRUE, shows progress bar. Defaults to FALSE.
#' @export
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return Output is a list of two dataframes. The first, $input, is a one-row data frame that contains all the arguments used in the call. The second, $all.des, contains the operating characteristics of all admissible designs found.
#' @examples output <- singlearmDesign(n.max = c(25, 30),
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
singlearmDesign <- function(
                        n.max,
                        p0,
                        p1,
                        alpha,
                        power,
                        C=1,
                        stages=NULL,
                        minstop=NULL,
                        maxthetaF=power,
                        minthetaE=power,
                        bounds="wald",
                        return.only.admissible=TRUE,
                        max.combns=1e4,
                        maxthetas=NA,
                        fixed.r=NA,
                        exact.thetaF=NA,
                        exact.thetaE=NA,
                        progressBar=TRUE){
  use.stages <- !is.null(stages)
  if(any(C!=1) & !is.null(stages)) stop("Values given for both cohort/block size C and number of stages. Please choose one only.")

if(is.null(minstop)){
  minstop <- 1
}

if(length(n.max)==2){
  nmin <- n.max[1]
  nmax <- n.max[2]
}
if(length(n.max)==1){
    nmax <- n.max
    nmin <- ifelse(use.stages, yes=floor(n.max/stages), no=max(C, minstop))
}

  if(minstop>nmin) stop("earliest stopping point (minstop) must not be greater than any potential maximum sample sizes.")



  if(use.stages==TRUE){
    intermediate.output <- lapply(stages,
                                  FUN = function(x) findDesignsGivenCohortStage(stages=x,
                                                         C=NA,
                                                         nmin=nmin,
                                                         nmax=nmax,
                                                         p0=p0,
                                                         p1=p1,
                                                         alpha=alpha,
                                                         power=power,
                                                         minstop=minstop,
                                                         maxthetaF=maxthetaF,
                                                         minthetaE=minthetaE,
                                                         bounds=bounds,
                                                         return.only.admissible=return.only.admissible,
                                                         max.combns=max.combns,
                                                         maxthetas=maxthetas,
                                                         fixed.r=fixed.r,
                                                         exact.thetaF=exact.thetaF,
                                                         exact.thetaE=exact.thetaE,
                                                         progressBar=progressBar)
           )
  }else{
    intermediate.output <- lapply(C,
                                  FUN = function(x) findDesignsGivenCohortStage(stages=NA,
                                                         C=x,
                                                         nmin=nmin,
                                                         nmax=nmax,
                                                         p0=p0,
                                                         p1=p1,
                                                         alpha=alpha,
                                                         power=power,
                                                         minstop=minstop,
                                                         maxthetaF=maxthetaF,
                                                         minthetaE=minthetaE,
                                                         bounds=bounds,
                                                         return.only.admissible=return.only.admissible,
                                                         max.combns=max.combns,
                                                         maxthetas=maxthetas,
                                                         fixed.r=fixed.r,
                                                         exact.thetaF=exact.thetaF,
                                                         exact.thetaE=exact.thetaE,
                                                         progressBar=progressBar)
    )
  }
  input <- data.frame(nmin=nmin,
                      nmax=nmax,
                      Cmin=ifelse(!use.stages, C[1], NA),
                      Cmax=ifelse(!use.stages, C[length(C)], NA),
                      stagesmin=ifelse(use.stages, stages[1], NA),
                      stagesmax=ifelse(use.stages, stages[length(stages)], NA),
                      minstop=minstop,
                      p0=p0,
                      p1=p1,
                      alpha=alpha,
                      power=power,
                      maxthetaF=maxthetaF,
                      minthetaE=minthetaE,
                      bounds=bounds,
                      fixed.r=fixed.r,
                      return.only.admissible=return.only.admissible,
                      max.combns=max.combns,
                      maxthetas=maxthetas,
                      exact.thetaF=exact.thetaF,
                      exact.thetaE=exact.thetaE)
  intermediate.output.combined <- do.call(rbind, intermediate.output)
  intermediate.output.combined[, c("alpha", "power", "EssH0", "Ess", "thetaF", "thetaE")] <- signif(intermediate.output.combined[, c("alpha", "power", "EssH0", "Ess", "thetaF", "thetaE")], digits=4)
  final.output <- list(input=input,
                       all.des=intermediate.output.combined)
  class(final.output) <- "curtailment_single"
  return(final.output)
}
