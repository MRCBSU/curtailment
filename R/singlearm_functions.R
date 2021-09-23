#' @import data.table
#' @import ggplot2
#' @import ggthemes
#' @import pkgcond
#' @import stats
#' @importFrom grDevices dev.off pdf
#' @importFrom utils combn head setTxtProgressBar txtProgressBar
#requireNamespace("data.table")

#### EXPLANATION OF CODE
# singlearmDesign: used to search for m-stage designs
# findDesignOCs: called within singlearmDesign
# findCPmatrix: called within singlearmDesign
# findSingleSimonDesign: finds operating characteristics of a single Simon design
# findWald: finds ESS under H0 and H1 for Wald's design







#' Find two-stage designs
#'
#' This function finds two-stage designs for a given set of design parameters, allowing
#' stopping for benefit at the interim (Mander and Thompson's design) or no stopping
#' for benefit at the interim (Simon's design). It returns not
#' only the optimal and minimax design realisations, but all design realisations that could
#' be considered "best" in terms of expected sample size under p=p0 (EssH0), expected
#' sample size under p=p1 (Ess), maximum sample size (n) or any weighted combination of these
#' three optimality criteria.
#'
#' @param nmin Minimum permitted sample size. Should be a multiple of block size or number of stages.
#' @param nmax Maximum permitted sample size. Should be a multiple of block size or number of stages.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta)
#' @param benefit Allow the trial to end for a go decision and reject the null hypothesis at the interim analysis (i.e., the design of Mander and Thompson)
#' @return A list of class "curtailment_simon" containing two data frames. The first data frame, $input,
#' has a single row and contains all the inputted values. The second data frame, $all.des, contains one
#' row for each design realisation, and contains the details of each design, including sample size,
#' stopping boundaries and operating characteristics. To see a diagram of any obtained design realisation
#' and its corresponding stopping boundaries, simply call the function drawDiagram with this output as the only argument.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @examples
#' \donttest{
#' find2stageDesigns(nmin=23,
#'  nmax=27,
#'  p0=0.75,
#'  p1=0.92,
#'  alpha=0.22,
#'  power=0.95,
#'  benefit=TRUE)
#'  }
#' @references
#' \doi{10.1016/j.cct.2010.07.008}{A.P. Mander, S.G. Thompson,
#' Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials,
#' Contemporary Clinical Trials,
#' Volume 31, Issue 6,
#' 2010,
#' Pages 572-578}
#'
#' \doi{10.1016/0197-2456(89)90015-9}{Richard Simon,
#' Optimal two-stage designs for phase II clinical trials,
#' Controlled Clinical Trials,
#' Volume 10, Issue 1,
#' 1989,
#' Pages 1-10}
#' @export
find2stageDesigns <- function(nmin, nmax, p0, p1, alpha, power, benefit=FALSE)
{

  if(benefit==FALSE){
    nr.lists <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=FALSE)
    simon.df <- apply(nr.lists, 1, function(x) {findSingleSimonDesign(n1=x[1], n2=x[2], r1=x[4], r=x[5], p0=p0, p1=p1)})
    simon.df <- t(simon.df)
    simon.df <- as.data.frame(simon.df)
    names(simon.df) <- c("n1", "n2", "n", "r1", "r", "alpha", "power", "EssH0", "Ess")
  } else{
    nr.lists <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=TRUE)
    simon.df <- apply(nr.lists, 1, function(x) {simonEfficacy(n1=x[1], n2=x[2], r1=x[4], r=x[5], p0=p0, p1=p1, e1=x[6])})
    simon.df <- t(simon.df)
    simon.df <- as.data.frame(simon.df)
    names(simon.df) <- c("n1", "n2", "n", "r1", "r", "alpha", "power", "EssH0", "Ess", "e1")
  }
  correct.alpha.power <- simon.df$alpha < alpha & simon.df$power > power
  simon.df <- simon.df[correct.alpha.power, ]

  if(nrow(simon.df)==0){
    stop("No suitable designs exist for these design parameters.")
  }
  # Discard all dominated designs. Note strict inequalities as EssH0 and Ess will (almost?) never be equal for two designs:
  discard <- rep(NA, nrow(simon.df))
  for(i in 1:nrow(simon.df))
  {
    discard[i] <- sum(simon.df$EssH0[i] > simon.df$EssH0 & simon.df$Ess[i] > simon.df$Ess & simon.df$n[i] >= simon.df$n)
  }

  simon.df <- simon.df[discard==0,]
  simon.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power)
  simon.output <- list(input=simon.input,
                       all.des=simon.df)
  class(simon.output) <- append(class(simon.output), "curtailment_simon")
  return(simon.output)
}


findCoeffs <- function(n, p0, p1){
  ##### Matrix of coefficients: Probability of a SINGLE path leading to each point:
  coeffs <- matrix(NA, ncol=n, nrow=n+1)
  coeffs.p0 <- coeffs
  q0 <- 1-p0
  q1 <- 1-p1
  for(i in 1:(n+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+n-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+n-1))
  }
  coeffs.list <- list(coeffs.p0=coeffs.p0,
                      coeffs=coeffs)
}







#' drawDiagram
#'
#' This function produces both a data frame and a diagram of stopping boundaries.
#' The function takes a single argument: the output from the function singlearmDesign.
#' If the supplied argument contains more than one admissible designs, the user is offered a choice of which design to use.
#' @param  findDesign.output Output from either the function singlearmDesign or find2stageDesigns
#' @param  print.row Choose a row number to directly obtain a plot and stopping boundaries for a particular design realisation. Default is NULL.
#' @param  xmax,ymax Choose values for the upper limits of the x- and y-axes respectively. Helpful for comparing two design realisations. Default is NULL.
#' @export
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return The output is a list of two elements. The first, $diagram, is a ggplot2 object showing how the trial should proceed: when to to undertake an interim analysis, that is, when to check if a stopping boundary has been reached (solid colours) and what decision to make at each possible point (continue / go decision / no go decision). The second list element, $bounds.mat, is a data frame containing three columns: the number of participants at which to undertake an interim analysis (m), and the number of responses at which the trial should stop for a go decision (success) or a no go decision (fail).
#' @examples output <- singlearmDesign(nmin = 30,
#'  nmax = 30,
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
#'  dig <- drawDiagram(output, print.row = 2)
drawDiagram <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  UseMethod("drawDiagram")
}

#' @export
drawDiagram.curtailment_single <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
    des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
} # end of drawDiagram()

#' @export
drawDiagram.curtailment_simon <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
}


#function generator
defunct <- function(msg = "This function is depreciated") function(...) return(stop(msg))
# @export
#' SCfn = defunct("SCfn changed name to findSCdesigns")
