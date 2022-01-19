
# weight_to_cumul : function converting weight loss to cumulated CO2 production with the 
# assumption that the lost weight is only/mostly CO2
# Author : Hugo Devillers

weight_to_cumul <- function(w, start = 2, t.ref = 2, vol = 15.2) {
  
  # Compute the weight of reference
  w.ref <- NULL
  w.len <- length(w)
  if( length(t.ref) == 1 ) {
    w.ref <- rep(w[t.ref], w.len)
  } else {
    w.ref <- rep(mean(w[t.ref]), w.len)
  }
  
  # Compute the cumulated values
  out <- (w.ref - w) / vol * 1000
  
  # Reset to 0 first values if required
  if( start > 1) {
    out[seq(1, start-1, by=1)] <- 0
  }
  
  return(out)
}


# moving_diff : function computing a moving difference between points separated 
# of l-1 points to estimate the local slope (i.e. CO2 flow rate, here) in the
# middle of these two points.
# Author : Hugo Devillers

moving_diff <- function(x, ti, l=3) {
  
  # Extract parameters
  r <- floor(l/2)
  n <- length(x)
  
  # Extend time
  dt1 <- ti[2] -ti[1]
  dtn <- ti[n] - ti[n-1]
  exti1 <- seq(ti[1] - n * dt1, ti[1] - dt1, by = dt1)
  extin <- seq(ti[n] + dtn, ti[n] + n * dtn, by= dtn)
  ti <- c(exti1, ti, extin)
  
  # Extend x (values)
  x <- c(rep(x[1], r), x, rep(x[n], r))
  
  # New parameters
  n <- length(x)
  from <- 1:(n-l+1)
  to <- l:n
  
  # Compute sliding weighted diff
  tmp <- mapply(function(fr, to, r, x, ti){
    num <- mean(x[(fr+r):to]) - mean(x[fr:(fr+r)])
    #den <- mean(ti[(fr+r):to]) - mean(ti[fr:(fr+r)])
    den <- mean(ti[c((fr+r),to)]) - mean(ti[c(fr,(fr+r))])
    return( num / den )
  }, fr = from, to = to, MoreArgs = list(x = x, r=r, ti=ti), SIMPLIFY = TRUE)
  
  # Return the output
  return(tmp) 
}