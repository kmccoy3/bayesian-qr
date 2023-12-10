# bayesian-qr
Bayesian quantile regression project completed for STAT 525 at Rice University. The final report is included as [report.pdf](./report.pdf). Additionally, we have included an associated slideshow presentation in [slides.pdf](./slides.pdf). This report specifically discusses standard Bayesian quantile regression and longitudinal Bayesian quantile regression. We have also included code that implements the methods discussed in the report. 

# Repository Contents
- [OBQR_syntheticdata.R](./OBQR_syntheticdata.R): Ordinary Bayesian quantile regression applied to synthetic data set. This code generates Figures 1-4 and 9-12.
- [OBQR_realdata.R](./OBQR_realdata.R): Ordinary Bayesian quantile regression applied to real data set. This code generates Figures 5.
- [BQRcode.R](./BQRcode.R): Ordinary Bayesian quantile regression applied on a sythentic data generated from a linear mixed-effects model (uses the bayesQR R package).
- [BQRGScode.R](./BQRGScode.R): Bayesian quantile regression Gibbs sampler for data generated from a linear mixed-effects model with 2 fixed effects and 2 random effects.
- [BQRMHcode.R](./BQRMHcode.R): Bayesian quantile regression Metropolis-Hastings sampler for data generated from a linear mixed-effects model with 2 fixed effects and 2 random effects.

