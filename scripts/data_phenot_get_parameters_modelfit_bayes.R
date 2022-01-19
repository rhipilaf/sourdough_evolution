
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(growthcurver) #contains the growthdata dataset
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

long_growthdata <- growthdata %>% gather(well,absorbance,-time)
glimpse(long_growthdata)

ggplot(long_growthdata,aes(time,absorbance,group=well)) +
  geom_line() + 
  theme_bw()

logisticgrowth_model <- "
functions {
  real[] logisticgrowth(real t,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i
                  ) {
    real dydt[x_i[1]];
    for (i in 1:x_i[1]){
      dydt[i] = theta[1] * y[i] * (1-y[i]/theta[2]);
    }
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> n_wells;
  real y0[n_wells];
  real z[T,n_wells];
  real t0;
  real ts[T];
}
transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1] = n_wells;
}
parameters {
  real<lower=0> theta[2];
  real<lower=0> sigma;
}
model {
  real y_hat[T,n_wells];
  theta ~ cauchy(0,2.5);
  sigma ~ normal(0,0.01);
  y_hat = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    for (i in 1:n_wells) {
      z[t,i] ~ normal(y_hat[t,i], sigma);
    }
  }
}
generated quantities{
  real y_pred[T,n_wells];
  real z_pred[T,n_wells];
  y_pred = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i );
  for (t in 1:T) {
    for(i in 1:n_wells){
      z_pred[t,i] = y_pred[t,i] + normal_rng(0,sigma);
    }
  }
}
"

nSamples = nrow(growthdata) - 1 #use time=0 as initial condition, take this as fixed
y0 = filter(growthdata,time==0) %>% select(-time) %>% unlist #initial condition
t0 = 0.0
ts = filter(growthdata,time>0) %>% select(time) %>% unlist
z = filter(growthdata,time>0) %>% select(-time)
n_wells = 9 #running on all wells can be slow
estimates <- sampling(model_code = ,
                      data = list (
                        T  = nSamples,
                        n_wells = n_wells,
                        y0 = y0[1:n_wells],
                        z  = z[,1:n_wells],
                        t0 = t0,
                        ts = ts
                      ),
                      seed = 123,
                      chains = 4,
                      iter = 1000,
                      warmup = 500
)

parametersToPlot = c("theta","sigma","lp__")
print(estimates, pars = parametersToPlot)