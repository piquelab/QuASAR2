# QuASAR2 Refactor: Mathematical Details

## 1. Model Definition

We model the number of reference reads $R_{ij}$ for SNP $i$ in observation $j$ using a beta-binomial distribution:

$$ R_{ij} \sim \text{BetaBinomial}(n_{ij}, \pi_{ij}, M_i) $$

where $n_{ij} = R_{ij} + A_{ij}$ is the total coverage, $\pi_{ij}$ is the probability of a reference read, and $M_i$ is the concentration parameter.

### 1.1. Parameterization

We reparameterize the concentration parameter $M_i$ using the dispersion parameter $\phi_i$:

$$ \phi_i = \log\left(\frac{1}{M_i}\right) \implies M_i = e^{-\phi_i} $$

Larger values of $\phi_i$ correspond to higher dispersion (smaller $M_i$).

The logistic model for $\pi_{ij}$ is:

$$ \text{logit}(\pi_{ij}^0) = X_j \beta_i + \text{offset}_{ij} $$

where $\pi_{ij}^0$ is the true underlying proportion. Accounting for a sequencing error rate $\epsilon$:

$$ \pi_{ij} = \pi_{ij}^0 (1 - \epsilon) + (1 - \pi_{ij}^0) \epsilon $$

## 2. Likelihood Functions

The log-likelihood for a single observation $(R, n, \pi, \phi)$ is:

$$ l(R; n, \pi, \phi) = \log \Gamma(e^{-\phi}) + \log \Gamma(R + \pi e^{-\phi}) + \log \Gamma(n - R + (1 - \pi) e^{-\phi}) - \log \Gamma(n + e^{-\phi}) - \log \Gamma(\pi e^{-\phi}) - \log \Gamma((1 - \pi) e^{-\phi}) $$

## 3. Cox-Reid Adjusted Profile Likelihood

To estimate the dispersion $\phi_i$ while accounting for the uncertainty in the estimated coefficients $\beta_i$, we use the Cox-Reid adjusted profile likelihood (APL):

$$ l_{APL}(\phi_i) = \sum_j l(R_{ij}; n_{ij}, \hat{\pi}_{ij}(\phi_i), \phi_i) - \frac{1}{2} \log \det(X^T W_i X) $$

where $\hat{\pi}_{ij}(\phi_i)$ are the probabilities evaluated at the MLE $\hat{\beta}_i$ for a given $\phi_i$. The weight matrix $W_i = \text{diag}(w_{ij})$ is defined by the expected information:

$$ w_{ij} = \frac{(\partial \mu_{ij} / \partial \eta_{ij})^2}{\text{Var}(R_{ij})} $$

With $\mu_{ij} = n_{ij} \pi_{ij}$ and $\text{Var}(R_{ij}) = n_{ij} \pi_{ij} (1 - \pi_{ij}) [1 + (n_{ij}-1)\theta_i]$, where $\theta_i = \frac{1}{M_i+1} = \frac{1}{e^{-\phi_i}+1}$.

## 4. Dispersion Moderation

Following the empirical Bayes approach (similar to edgeR and DESeq2), we moderate the SNP-specific dispersions $\phi_i$ towards a trended value $\phi_{trend}(A_i)$, where $A_i$ is the average log-abundance of SNP $i$.

The moderated dispersion $\hat{\phi}_i$ is obtained by maximizing the penalized APL:

$$ l_{mod}(\phi_i) = l_{APL}(\phi_i) - \frac{1}{2s^2}(\phi_i - \phi_{trend}(A_i))^2 $$

where $s^2$ represents the prior variance of the dispersions across SNPs.

## 5. Statistical Inference

### 5.1. Likelihood Ratio Test (LRT)

To test a null hypothesis $H_0: L\beta_i = 0$ (where $L$ is a contrast matrix), we fit both the full model and the reduced model, and calculate the test statistic:

$$ \Lambda_i = 2 [l(\hat{\beta}_{full}, \hat{\phi}_i) - l(\hat{\beta}_{reduced}, \hat{\phi}_i)] $$

Under the null hypothesis, $\Lambda_i$ follows a $\chi^2$ distribution with degrees of freedom equal to the difference in the number of parameters.
