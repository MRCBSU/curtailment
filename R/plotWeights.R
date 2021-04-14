findLossSingleDesign <- function(single.des, weight.mat){
    current.des.ess0 <- single.des["EssH0"]
    current.des.ess1 <- single.des["Ess"]
    current.des.n <- single.des["n"]
    vector.of.loss <- apply(weight.mat, 1, function(x) x[1]*current.des.ess0 + x[2]*current.des.ess1 + x[3]*current.des.n)
    vector.of.loss
}

plotWeightsOneCohortSize <- function(des){
  w <- seq(0, 1, by=0.01)
  weights.matrix <- expand.grid(w, w)
  weights.matrix <- weights.matrix[rowSums(weights.matrix)<=1, ]
  weights.matrix <- cbind(weights.matrix, 1-rowSums(weights.matrix))
  weights.matrix <- weights.matrix[-1, ]

  loss.matrix <- apply(des, 1, findLossSingleDesign, weights.matrix)
  # What is the best design for each weight?
  best.des <- apply(loss.matrix, 1, which.min)
  weights.df <- data.frame(weights.matrix, factor(best.des))
  names(weights.df) <- c("w_0", "w_1", "1-w_0-w_1", "Design")
  design.details <- apply(des, 1, function(x) paste("{", x["r"], "/", x["n"], ", ", format(round(as.numeric(x["thetaF"]),2), nsmall=2), "/", format(round(as.numeric(x["thetaE"]),2), nsmall=2), "}", sep=""))

  output <- ggplot2::ggplot(weights.df, ggplot2::aes(x=w_1, y=w_0)) +
    ggplot2::geom_raster(ggplot2::aes(fill = Design)) +
    ggplot2::scale_fill_discrete(name="Design",labels=as.vector(design.details)) +
    ggplot2::xlab(expression(w[1])) +
    ggplot2::ylab(expression(w[0])) +
    ggplot2::theme_bw()+
    ggplot2::theme(#axis.text.x=element_text(),
      legend.text=ggplot2::element_text(size=ggplot2::rel(1.5)),
      legend.title=ggplot2::element_text(size=ggplot2::rel(1.5)),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.4), hjust=0.5),
      axis.ticks=ggplot2::element_blank(),
      axis.line=ggplot2::element_blank(),
      axis.title=ggplot2::element_text(size=ggplot2::rel(1.5)),
      axis.title.y=ggplot2::element_text(angle = 0, vjust=0.5),
      axis.text=ggplot2::element_text(size=ggplot2::rel(1.5)),
      legend.key.size = ggplot2::unit(1, "cm"),
      legend.position = c(0.8,0.75),
      panel.border=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(color='#eeeeee'))+
    ggplot2::ggtitle(paste("Cohort size=", des[1, "C"], ", No. of stages=", des[1, "stage"]))

output
}

#' @export
plotWeights <- function(main.output, split.by.cohort.size=TRUE){
  if(var(main.output$all.des[, "C"])==0 | split.by.cohort.size==FALSE){
      output <- plotWeightsOneCohortSize(main.output$all.des)
  }else{
    des.df <- as.data.frame(main.output$all.des)
    all.des.list <- split(des.df, des.df$C)
    output <- lapply(all.des.list, plotWeightsOneCohortSize)
  }
  output
}
