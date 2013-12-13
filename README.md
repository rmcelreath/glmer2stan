glmer2stan
==========

Define Stan models using glmer-style (lme4) formulas.

Stan (mc-stan.org) is a Hamiltonian Monte Carlo engine for fitting Bayesian models to data.

glmer2stan compiles design formulas, such as y ~ (1|id) + x, into Stan model code. This allows the specification of simple multilevel models, using familiar formula syntax of the kind many people have learned from popular R packages like lme4.

Status
==========

Current state: Production usable. All models pass tests, DIC computations validated. Basic prediction function added, but still needs work. WAIC is in and working for single-formula models, but still needs some more thorough numerical testing, to find corner cases.

Horizon: Code will eventually need to be refactored, to reduce redundancy and clean up some of the compilation algorithm. Speed of Stan code could be improved during that refactoring.

Examples
==========

(1) Gaussian with varying intercepts and slopes

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
varef(m1_g2s)$expection
```
