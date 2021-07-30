#' @export
drawDiagram.power.htest <- function(findDesign.output, print.row=NULL, save.plot=FALSE, xmax=NULL, ymax=NULL){
  m <- Sm <- decision <- analysis <- NULL
  des <- findDesign.output
  # des is output from power.prop.test()
  n <- 2*ceiling(des$n)
  coords <- expand.grid(0:n, 1:n)
  diag.df <- data.frame(Sm=as.numeric(coords[,1]),
                        m=as.numeric(coords[,2]),
                        decision=rep("Continue", nrow(coords)))
  diag.df$decision <- as.character(diag.df$decision)
  diag.df$decision[coords[,1]>coords[,2]] <- NA
  diag.df$decision[diag.df$m==n] <- "Make decision"
  # Add shading:
  diag.df.subset <- diag.df[!is.na(diag.df$decision),]
  diag.df.subset$analysis <- "No"
  diag.df.subset$analysis[diag.df.subset$m==n] <- "Yes"

  plot.title <- "Single-stage two-arm design"
  sub.text1 <- paste("Max no. of analyses: 1. Max(N): ", n, sep="")
  plot.subtitle2 <- bquote(.(sub.text1))
  diagram <- pkgcond::suppress_warnings(ggplot2::ggplot(data=diag.df.subset, mapping = aes(x=m, y=Sm, fill=decision, alpha=analysis))+
                                          scale_alpha_discrete(range=c(0.5, 1)),
                                        "Using alpha for a discrete variable is not advised")
  diagram <- diagram +
    geom_tile(color="white")+
    labs(fill="Decision",
         alpha="Analysis",
         x="Number of participants",
         y="Number of responses on treatment + non-responses on control",
         title=plot.title,
         subtitle = plot.subtitle2
         )+
    coord_cartesian(expand = 0)+
    theme_minimal()
  if(!is.null(xmax)){
    diagram <- diagram +
      expand_limits(x=xmax)
  }
  if(!is.null(ymax)){
    diagram <- diagram +
      expand_limits(y=ymax)
  }
  print(diagram)
  if(save.plot==TRUE){
    save.plot.yn <- readline("Save plot? (y/n): ")
    if(save.plot.yn=="y"){
      plot.filename <- readline("enter filename (no special characters or inverted commas): ")
      plot.filename <- paste(plot.filename, ".pdf", sep="")
      plot.width <- 8
      scaling <- 1.25
      plot.height <- 8*scaling*(des$r+1)/des$n
      pdf(plot.filename, width = plot.width, height = plot.height)
      print(diagram)
      dev.off()
    }
  }
  return(diagram)
}
