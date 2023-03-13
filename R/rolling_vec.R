rolling_vec <- function(x,i){
  len <- length(x)
  out <-
    sapply(seq_along(x),function(x1){
      i2 <- i + x1
      if(i2 <= len) return(i2)
      if(i2 - len > len) stop("index is greater than vector length", call. = FALSE)
      if(i2 > len) return(i2 - len)

  })
  return(x[out])

}
