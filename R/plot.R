
#' plot
#'
#' This function produces both a data frame and a diagram of stopping boundaries.
#' The function takes two arguments: the output from the function singlearmDesign or find2stageDesigns and the row number of the design to plot.
#' @param  findDesign.output Output from either the function singlearmDesign or find2stageDesigns
#' @param  print.row Choose a row number to directly obtain a plot and stopping boundaries for a particular design realisation.
#' @param  xmax,ymax Choose values for the upper limits of the x- and y-axes respectively. Helpful for comparing two design realisations. Default is NULL.
#' @param ... Any further arguments to be passed
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
#'  dig <- plot(output, print.row=2)
#' @export
plot.curtailment_single <- function(findDesign.output, print.row, xmax=NULL, ymax=NULL, ...){
  des <- findDesign.output$all.des
  des <- des[print.row, , drop=FALSE]
  des.input <- findDesign.output$input
  plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax, ...)
  return(plot.and.bounds)
}

#' @export
plot.curtailment_2stage <- function(findDesign.output, print.row, xmax=NULL, ymax=NULL){
  des <- findDesign.output$all.des
  des <- des[print.row, , drop=FALSE]
  des.input <- findDesign.output$input
  plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
  return(plot.and.bounds)
}
