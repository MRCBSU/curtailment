#'@export
findLoss <- function(main.output, w0, w1){
  weight.vec <- data.frame(w0, w1, 1-w0-w1)
  loss <- apply(main.output$all.des, 1, findLossSingleDesign, weight.vec)
  loss.ranked <- rank(loss)
  loss.rounded <- round(loss, 1)
  main.output$all.des <- cbind(main.output$all.des, loss.rounded, loss.ranked)
  main.output
}
