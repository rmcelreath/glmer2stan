glmer2stan
==========

Define Stan models using glmer-style (lme4) formulas.

Stan (mc-stan.org) is a Hamiltonian Monte Carlo engine for fitting Bayesian models to data.

glmer2stan compiles design formulas, such as y ~ (1|id) + x, into Stan model code. This allows the specification of simple multilevel models, using familiar formula syntax of the kind many people have learned from popular R packages like lme4.

The Stan code that glmer2stan produces is didactic, not efficient. So instead of using design matrices for varying and fixed effects, it builds out explicit additive expressions. This makes the code easier for novices to understand. But it also reduces speed a little. Since I originally conceived of glmer2stan as a teaching tool, the Stan code will probably remain didactic. But I might eventually make an option to use faster-but-opaque techniques.

Status
==========

Current state: All models pass tests, DIC computations validated. Basic prediction function added, but still needs work. WAIC is in and working for single-formula models, but still needs some more thorough numerical testing, to find corner cases.

Horizon: I'm not planning on developing glmer2stan much further. The problem is that I find it unsatisfying both as a teaching tool and as a multilevel regression tool. 

It is unsatisfying for teaching, because students don't tend use it as a way to learn the connections between Stan code and glmer formulas, but rather as a way to avoid learning Stan code. For those who are using glmer2stan as a teaching tool, it remains and will remain useful, as it exists now. I don't think everyone needs to learn Stan code, mind you. But I would like a tool that helps more in that regard.

It is unsatisfying as a multilevel regression tool, because there's no natural way to allow prior specifications. And the natural default priors are the same flat (or nearly flat) priors as glmer, which leads to lots of problems with Stan's HMC. Also, complex issues like specifying covariance priors are very hard to approach at all with glmer's syntax. This isn't a criticism of glmer (which remains a great tool), but rather just a comment on how far one can reasonably push its syntax.

So caught between bad defaults and no easy input syntax (other than editing Stan code directly), I'd rather give my time to another project that addresses all these issues at once. That project (map2stan) is in my other package, rethinking, and it will be described in a book I'm writing. map2stan uses templating and all the other rational things I didn't use when hacking together glmer2stan. So it'll also be able to grow as Stan grows.

I think it might be worth using some glmer2stan code to translate glmer formulas into map2stan input formulas, which do lay bare all the model assumptions. But I'm still figuring out how to proceed. My primary interests are pedagogical here, so I'm trying to tune the tools so they balance convenience and understanding.

Examples
==========

(1) Gaussian with varying intercepts and slopes
-------

```
library(lme4)
library(glmer2stan)

data(sleepstudy) # built into lme4

# fit with lme4
m1_lme4 <- lmer( Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE )

# construct subject index --- glmer2stan forces you to manage your own cluster indices
sleepstudy$subject_index <- as.integer(as.factor(sleepstudy$Subject))

# fit with glmer2stan
m1_g2s <- lmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy )

# [timings for 3GHz i5, no memory limitations]
#Elapsed Time: 46.1606 seconds (Warm-up)
#              32.3429 seconds (Sampling)
#              78.5035 seconds (Total)

summary(m1_lme4)

#Linear mixed model fit by maximum likelihood ['lmerMod']
#Formula: Reaction ~ Days + (Days | Subject) 
#   Data: sleepstudy
#
#      AIC       BIC    logLik  deviance 
#1763.9393 1783.0971 -875.9697 1751.9393 
#
#Random effects:
# Groups   Name        Variance Std.Dev. Corr
# Subject  (Intercept) 565.52   23.781       
#          Days         32.68    5.717   0.08
# Residual             654.94   25.592       
#Number of obs: 180, groups: Subject, 18
#
#Fixed effects:
#            Estimate Std. Error t value
#(Intercept)  251.405      6.632   37.91
#Days          10.467      1.502    6.97
#
#Correlation of Fixed Effects:
#     (Intr)
#Days -0.138

stanmer(m1_g2s)

#glmer2stan model: Reaction ~ Days + (Days | subject_index) [gaussian]
#
#Level 1 estimates:
#            Expectation StdDev   2.5%  97.5%
#(Intercept)      249.73   7.70 234.51 265.10
#Days              10.40   1.67   7.08  13.73
#sigma             25.76   1.50  22.99  28.80
#
#Level 2 estimates:
#(Std.dev. and correlations)
#
#Group: subject_index (18 groups / imbalance: 0)
#                 
#                    (1)  (2)
#  (1) (Intercept) 32.33   NA
#  (2) Days        -0.03 7.56
#
#DIC: 1711   pDIC: 31.9   Deviance: 1647.2

# compare varying effect estimates with:
ranef(m1_lme4)
varef(m1_g2s)$expectation

# extract posterior samples
posterior <- extract(m1_g2s)
str(posterior)

# compute posterior predictions
pp <- stanpredict(m1_g2s,data=sleepstudy)
str(pp)

#List of 1
# $ Reaction:List of 3
#  ..$ mu    : num [1:180] 252 272 292 312 332 ...
#  ..$ mu.ci : num [1:2, 1:180] 225 278 249 294 273 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:2] "2.5%" "97.5%"
#  .. .. ..$ : NULL
#  ..$ obs.ci: num [1:2, 1:180] 195 308 216 325 238 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:2] "2.5%" "97.5%"
#  .. .. ..$ : NULL

# compare to lme4 MAP/MLE predictions
predict(m1_lme4)

# plot comparison
plot( predict(m1_lme4) , pp$Reaction$mu )
```

You can expose the Stan model code by pulling out m1_g2s@stanmodel:
```
data{
    int N;
    real Reaction[N];
    real Days[N];
    int subject_index[N];
    int N_subject_index;
}

transformed data{
    vector[2] zeros_subject_index;
    for ( i in 1:2 ) zeros_subject_index[i] <- 0;
}

parameters{
    real Intercept;
    real beta_Days;
    real<lower=0> sigma;
    vector[2] vary_subject_index[N_subject_index];
    cov_matrix[2] Sigma_subject_index;
}

model{
    real vary[N];
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    beta_Days ~ normal( 0 , 100 );
    sigma ~ uniform( 0 , 100 );
    // Varying effects
    for ( j in 1:N_subject_index ) vary_subject_index[j] ~ multi_normal( zeros_subject_index , Sigma_subject_index );
    // Fixed effects
    for ( i in 1:N ) {
        vary[i] <- vary_subject_index[subject_index[i],1]
                + vary_subject_index[subject_index[i],2] * Days[i];
        glm[i] <- vary[i] + Intercept
                + beta_Days * Days[i];
    }
    Reaction ~ normal( glm , sigma );
}

generated quantities{
    real dev;
    real vary[N];
    real glm[N];
    dev <- 0;
    for ( i in 1:N ) {
        vary[i] <- vary_subject_index[subject_index[i],1]
                + vary_subject_index[subject_index[i],2] * Days[i];
        glm[i] <- vary[i] + Intercept
                + beta_Days * Days[i];
        dev <- dev + (-2) * normal_log( Reaction[i] , glm[i] , sigma );
    }
}
```

(2) Binomial with varying intercepts
--------

```
data(cbpp) # built into lme4

m2_lme4 <- glmer( cbind(incidence,size-incidence) ~ period + (1|herd) , data=cbpp , family=binomial )

cbpp$herd_index <- as.integer(cbpp$herd)
m2_g2s <- glmer2stan( cbind(incidence,size-incidence) ~ period + (1|herd_index) , data=cbpp , family="binomial" )

summary(m2_lme4)
stanmer(m2_g2s)
```
The Stan model code, accessed by m2_g2s@stanmodel:
```
data{
    int N;
    int incidence[N];
    real period2[N];
    real period3[N];
    real period4[N];
    int herd_index[N];
    int bin_total[N];
    int N_herd_index;
}

parameters{
    real Intercept;
    real beta_period2;
    real beta_period3;
    real beta_period4;
    real vary_herd_index[N_herd_index];
    real<lower=0> sigma_herd_index;
}

model{
    real vary[N];
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    beta_period2 ~ normal( 0 , 100 );
    beta_period3 ~ normal( 0 , 100 );
    beta_period4 ~ normal( 0 , 100 );
    sigma_herd_index ~ uniform( 0 , 100 );
    // Varying effects
    for ( j in 1:N_herd_index ) vary_herd_index[j] ~ normal( 0 , sigma_herd_index );
    // Fixed effects
    for ( i in 1:N ) {
        vary[i] <- vary_herd_index[herd_index[i]];
        glm[i] <- vary[i] + Intercept
                + beta_period2 * period2[i]
                + beta_period3 * period3[i]
                + beta_period4 * period4[i];
        glm[i] <- inv_logit( glm[i] );
    }
    incidence ~ binomial( bin_total , glm );
}

generated quantities{
    real dev;
    real vary[N];
    real glm[N];
    dev <- 0;
    for ( i in 1:N ) {
        vary[i] <- vary_herd_index[herd_index[i]];
        glm[i] <- vary[i] + Intercept
                + beta_period2 * period2[i]
                + beta_period3 * period3[i]
                + beta_period4 * period4[i];
        dev <- dev + (-2) * binomial_log( incidence[i] , bin_total[i] , inv_logit(glm[i]) );
    }
}
```

(3) Zero-augmented (zero-inflated) gamma model
-------
This is really a two-outcome non-linear factor analysis (or Gaussian process) of a sort, using varying intercepts to relate outcomes from the same individuals. It demonstrates how to specify models with more than one formula and use varying effects to link them together. This model cannot be fit with lme4, but the formula syntax still follows lme4 conventions (mostly).
```
data(Ache) # built into glmer2stan

# prepare two outcome formula list
the_formula <- list(
		iszero ~ (1|hunter.id) ,
		nonzero ~ (1|hunter.id)
	)

# note the list of family names in call to glmer2stan
m3 <- glmer2stan( the_formula , data=Ache , family=list("binomial","gamma") )

stanmer(m3)

# plot varying intercepts across outcomes, showing correlation
v <- varef(m3)
plot( v$expectation$hunter_id[,1] , v$expectation$hunter_id[,2] , xlab="hunter effect (failures)" , ylab="hunter effect (harvests)" )
```
As before, you can view the Stan code by extracting m3@stanmodel.
