
moyenne_glissante <- function(x, l=3) {
  # Extend boundarie
  r <- floor(l/2)
  n <- length(x)
  x <- c(rep(x[1], r), x, rep(x[n], r))
  n <- length(x)
  from <- 1:(n-l+1)
  to <- l:n
  tmp <- mapply(function(fr, to, x){
    return(mean(x[fr:to]))
  }, fr = from, to = to, MoreArgs = list(x = x), SIMPLIFY = TRUE)
  return(tmp)
}

#Dérivée de la moyenne glissante
diff_glissante <- function(x, ti, l=3) {
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

weight_to_cumul <- function(w, start = 2, t.ref = 2, vol = 15.4) { #15.4 to change to 15.2
  # Compute the weight of reference
  w.ref <- NULL
  w.len <- length(w)
  if( length(t.ref) == 1 ) {
    w.ref <- rep(w[t.ref], w.len)
  } else {
    w.ref <- rep(mean(w[t.ref]), w.len)
  }
  
  # Compute the cumuled values
  out <- (w.ref - w) / vol * 1000
  
  # Reset to 0 first values if required
  if( start > 1) {
    out[seq(1, start-1, by=1)] <- 0
  }
  
  return(out)
}

get_raw_max_cumul <- function(m) {
  return(max(m[,6]))
}

get_stat_from_cumul <-function(cum, l=5) {
  c.max <- max(cum[,2])
  c.lat <- cum[min(which(cum[,2] > 1)), 1]
  d.val <- diff_glissante(cum[,2], cum[,1], l = l)
  d.top <- which.max(d.val)
  d.max <- d.val[d.top]
  d.tma <- cum[d.top, 1]
  return(c(
    latency = c.lat,
    co2.max = c.max,
    v.max = d.max,
    tv.max = d.tma
  ))
}
