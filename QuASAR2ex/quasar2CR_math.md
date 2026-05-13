# QuASAR2 with Cox–Reid Profile Likelihood and Empirical-Bayes Dispersion Moderation

This note documents the statistical model and algorithms for `fitQuasar2CR`, a rewrite of `fitQuasar2` that

1.  parameterizes the beta–binomial precision $M$ as $\psi = \log(1/M)$ so that the dispersion lives on the whole real line,
2.  estimates $\psi$ per SNP using the **Cox–Reid adjusted profile likelihood** (rather than a coverage-binned grid),
3.  **moderates** the per-SNP $\psi$ across SNPs using an `edgeR`/`DESeq2`-style trend + empirical-Bayes shrinkage, and
4.  exposes a design-matrix interface so that **likelihood ratio tests** and **contrasts** become first-class operations.

Throughout we keep the existing QuASAR2 sequencing-error mixture for $\varepsilon$, which is shared across SNPs and profiled out jointly.

------------------------------------------------------------------------

## 1. Notation and model

For SNP $g = 1,\dots,G$ we observe $i = 1,\dots,n_g$ replicates / conditions. Let

-   $R_{gi}$ = reads on the reference allele,
-   $A_{gi}$ = reads on the alternate allele,
-   $N_{gi} = R_{gi} + A_{gi}$,
-   $\mathbf{x}_{gi} \in \mathbb{R}^{p_g}$ = covariate row from `model.matrix(design, data_g)`,
-   $o_{gi}$ = optional offset (e.g. $\text{qlogis}(\text{DNA proportion})$ for MPRA).

Let $\boldsymbol{\beta}_g \in \mathbb{R}^{p_g}$ be the per-SNP regression coefficients, $M_g > 0$ the beta–binomial precision (sum of the Beta shape parameters), and $\varepsilon \in [0, 1/2]$ a single shared sequencing error.

### 1.1 Linear predictor and error-corrected mean

$$
\eta_{gi} = \mathbf{x}_{gi}^\top \boldsymbol{\beta}_g + o_{gi},
\qquad
\mu_{gi} = \mathrm{logit}^{-1}(\eta_{gi}) = \frac{1}{1 + e^{-\eta_{gi}}},
$$

$$
\tilde\mu_{gi} = (1-\varepsilon)\mu_{gi} + \varepsilon(1 - \mu_{gi})
              = \varepsilon + (1 - 2\varepsilon)\mu_{gi}.
$$

The mixture $\tilde\mu_{gi}$ is the probability a single read is called as the reference allele after a Bernoulli$(\varepsilon)$ allele swap due to sequencing error.

### 1.2 Beta–binomial likelihood with $\psi = \log(1/M)$

$$
R_{gi} \mid \boldsymbol{\beta}_g, M_g, \varepsilon
\sim \mathrm{BetaBinom}\!\left(N_{gi},\; \alpha_{gi} = M_g \tilde\mu_{gi},\;
                                       \beta^*_{gi} = M_g (1-\tilde\mu_{gi})\right).
$$

Reparametrizing:

$$
\boxed{\psi_g := \log(1/M_g) = -\log M_g, \qquad M_g = e^{-\psi_g}.}
$$

This is the natural log-dispersion: it is **unbounded**, so unconstrained optimization (BFGS, Brent) is well-posed, and the prior in §4 is naturally placed on $\psi_g$. The mapping to the more familiar intra-class correlation is $\rho_g = 1/(M_g + 1)$, so $\psi_g = \log(\rho_g / (1 - \rho_g))$ is the logit of $\rho_g$ asymptotically; small $\psi_g$ means a tight beta–binomial (close to binomial) and large $\psi_g$ means heavy over-dispersion.

The log-likelihood for one observation is

$$
\ell_{gi}(\boldsymbol{\beta}_g, \psi_g, \varepsilon)
= \log\binom{N_{gi}}{R_{gi}}
+ \log \mathrm{B}(R_{gi} + \alpha_{gi}, A_{gi} + \beta^*_{gi})
- \log \mathrm{B}(\alpha_{gi}, \beta^*_{gi}),
$$

with $\mathrm{B}(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)$. Expanding into log-gammas:

$$
\ell_{gi} =
\log\binom{N_{gi}}{R_{gi}}
+ \log\Gamma(R_{gi} + M_g \tilde\mu_{gi})
+ \log\Gamma(A_{gi} + M_g (1-\tilde\mu_{gi}))
- \log\Gamma(N_{gi} + M_g)
- \log\Gamma(M_g \tilde\mu_{gi})
- \log\Gamma(M_g (1-\tilde\mu_{gi}))
+ \log\Gamma(M_g).
$$

The per-SNP log-likelihood is $\ell_g(\boldsymbol{\beta}_g, \psi_g, \varepsilon) = \sum_{i=1}^{n_g} \ell_{gi}$.

The full log-likelihood is $\ell(\boldsymbol{\theta}, \varepsilon) = \sum_g \ell_g$ with $\boldsymbol{\theta} = (\boldsymbol{\beta}_1, \psi_1, \dots, \boldsymbol{\beta}_G, \psi_G)$.

------------------------------------------------------------------------

## 2. Score and information for $\boldsymbol{\beta}_g$

We need analytic derivatives both for fast BFGS fitting and for the Cox–Reid adjustment.

### 2.1 Score with respect to $\tilde\mu$

Differentiating $\ell_{gi}$ with respect to $\tilde\mu_{gi}$ and using $\frac{d}{dx}\log\Gamma(x) = \psi^{(0)}(x)$ (digamma):

$$
\frac{\partial \ell_{gi}}{\partial \tilde\mu_{gi}}
= M_g\Big[
\psi^{(0)}(R_{gi} + M_g \tilde\mu_{gi}) - \psi^{(0)}(A_{gi} + M_g (1-\tilde\mu_{gi}))
- \psi^{(0)}(M_g \tilde\mu_{gi}) + \psi^{(0)}(M_g (1-\tilde\mu_{gi}))
\Big].
$$

### 2.2 Chain rule to $\boldsymbol{\beta}_g$

Since $\partial \tilde\mu / \partial \mu = 1 - 2\varepsilon$, $\partial \mu / \partial \eta = \mu(1-\mu)$, and $\partial \eta / \partial \boldsymbol{\beta} = \mathbf{x}^\top$, the score is

$$
\mathbf{s}_g(\boldsymbol{\beta}_g) := \frac{\partial \ell_g}{\partial \boldsymbol{\beta}_g}
= \sum_{i=1}^{n_g}
\underbrace{(1 - 2\varepsilon)\mu_{gi}(1-\mu_{gi})
\frac{\partial \ell_{gi}}{\partial \tilde\mu_{gi}}}_{=: u_{gi}}
\mathbf{x}_{gi}
= \mathbf{X}_g^\top \mathbf{u}_g.
$$

### 2.3 Observed information

Differentiating once more, and writing $\psi^{(1)} = \text{trigamma}$:

$$
\frac{\partial^2 \ell_{gi}}{\partial \tilde\mu_{gi}^2}
= M_g^2\Big[
\psi^{(1)}(R_{gi} + M_g\tilde\mu_{gi}) + \psi^{(1)}(A_{gi} + M_g(1-\tilde\mu_{gi}))
- \psi^{(1)}(M_g\tilde\mu_{gi}) - \psi^{(1)}(M_g(1-\tilde\mu_{gi}))
\Big] \leq 0.
$$

Applying the chain rule to $\boldsymbol{\beta}_g$ and dropping the score-times-curvature cross-term (which vanishes in expectation, and in practice at the MLE),

$$
\mathcal{I}_g(\boldsymbol{\beta}_g \mid \psi_g, \varepsilon)
= -\frac{\partial^2 \ell_g}{\partial \boldsymbol{\beta}_g \partial \boldsymbol{\beta}_g^\top}
\approx \mathbf{X}_g^\top \mathbf{W}_g \mathbf{X}_g,
$$

with diagonal working weights

$$
(W_g)_{ii} = -\frac{\partial^2 \ell_{gi}}{\partial \tilde\mu_{gi}^2}
             (1-2\varepsilon)^2 \mu_{gi}^2(1-\mu_{gi})^2 \geq 0.
$$

This is the form used as the Cox–Reid correction term below.

For inner-loop $\boldsymbol{\beta}_g$ optimization we use BFGS with the analytic score; $\mathcal{I}_g^{-1}$ at convergence gives the Wald covariance matrix used for contrasts.

------------------------------------------------------------------------

## 3. Cox–Reid adjusted profile likelihood for $\psi_g$

The naive profile log-likelihood

$$
\ell_g^{\text{prof}}(\psi_g) := \ell_g\big(\hat{\boldsymbol{\beta}}_g(\psi_g), \psi_g\big),
$$

with $\hat{\boldsymbol{\beta}}_g(\psi_g) = \arg\max_{\boldsymbol{\beta}} \ell_g(\boldsymbol{\beta}, \psi_g)$, is downward-biased in $\psi_g$ when $p_g$ is non-negligible relative to $n_g$, because it treats $\hat{\boldsymbol{\beta}}_g$ as exact.

Cox & Reid (1987) propose the **adjusted profile log-likelihood**

$$
\boxed{
\mathrm{APL}_g(\psi_g)
= \ell_g\big(\hat{\boldsymbol{\beta}}_g(\psi_g), \psi_g\big)
- \tfrac{1}{2}\log\det \mathcal{I}_g\big(\hat{\boldsymbol{\beta}}_g(\psi_g) \mid \psi_g\big),
}
$$

which is, equivalently, the Laplace approximation to the marginal log-likelihood of $\psi_g$ under a flat prior on $\boldsymbol{\beta}_g$:

$$
\log \int \exp\lbrace\ell_g(\boldsymbol{\beta}_g, \psi_g)\rbrace\, d\boldsymbol{\beta}_g
\approx \ell_g(\hat{\boldsymbol{\beta}}_g, \psi_g) + \tfrac{p_g}{2}\log(2\pi)
       - \tfrac{1}{2}\log\det \mathcal{I}_g(\hat{\boldsymbol{\beta}}_g \mid \psi_g).
$$

This is exactly the adjustment used in `edgeR` (McCarthy, Chen & Smyth 2012) for negative-binomial dispersion estimation, and in `DESeq2` (Love, Huber & Anders 2014) for the gene-wise dispersion step. It generalizes REML to the non-Gaussian case.

### 3.1 Per-SNP gene-wise estimate

$$
\hat\psi_g^{\text{gw}} := \arg\max_{\psi_g} \mathrm{APL}_g(\psi_g),
$$

obtained via a 1-D Brent search over a wide interval (we use $[-15, 15]$ in $\psi$-space, i.e. $M \in [e^{-15}, e^{15}]$). For each candidate $\psi_g$ the inner optimization $\hat{\boldsymbol{\beta}}_g(\psi_g)$ is warm-started from the previous $\hat{\boldsymbol{\beta}}_g$.

------------------------------------------------------------------------

## 4. Empirical-Bayes moderation across SNPs

The $\hat\psi_g^{\text{gw}}$ are noisy when $n_g$ is small. We shrink them towards a smooth trend, exactly as in `DESeq2`.

### 4.1 Trend fit

Let $\bar N_g = \frac{1}{n_g}\sum_i N_{gi}$ be the per-SNP mean coverage. We fit a smooth function

$$
\psi^{\text{trend}}(\bar N_g) = f\big(\log \bar N_g\big)
$$

by local polynomial regression on $\lbrace (\log \bar N_g, \hat\psi_g^{\text{gw}}) \rbrace_g$. We use `locfit` with the default tricube kernel when available, falling back to a natural cubic spline. The trend is monotonically non-increasing in $\bar N_g$ in typical ASE / MPRA data: higher coverage yields cleaner estimates.

### 4.2 Prior variance

Define residuals from the trend:

$$
e_g = \hat\psi_g^{\text{gw}} - \psi^{\text{trend}}(\bar N_g).
$$

Their variance decomposes (approximately) as

$$
\mathrm{Var}(e_g) \approx \sigma_\psi^2 + v_g^{\text{samp}},
$$

where $\sigma_\psi^2$ is the prior variance to estimate and $v_g^{\text{samp}}$ is the sampling variance of $\hat\psi_g^{\text{gw}}$. The latter is approximated as the negative inverse of the APL second derivative at the gene-wise estimate:

$$
v_g^{\text{samp}} \approx \left(
-\frac{d^2 \mathrm{APL}_g}{d\psi^2} \Bigg|_{\psi = \hat\psi_g^{\text{gw}}}
\right)^{-1},
$$

computed by central difference. Following `DESeq2` we estimate $\sigma_\psi^2$ robustly as the robust variance of residuals minus the median sampling variance, truncated at zero:

$$
\hat\sigma_\psi^2 = \max\left(0,\;
\frac{\mathrm{median}_g\lbrace(e_g - \mathrm{median}(e))^2\rbrace}{c_{4}^{2}}
- \mathrm{median}_g\lbrace v_g^{\text{samp}}\rbrace
\right),
$$

where $c_4 \approx 0.6745$ is the normal-consistency constant ($\mathrm{MAD}^2 / c_4^2$ is a robust variance estimate).

### 4.3 MAP estimate of $\psi_g$

We **shrink** by placing a normal prior on $\psi_g$:

$$
\psi_g \sim \mathcal{N}\left(\psi^{\text{trend}}(\bar N_g),\; \hat\sigma_\psi^2\right),
$$

and computing the maximum a-posteriori estimate:

$$
\boxed{
\hat\psi_g^{\text{MAP}} = \arg\max_{\psi_g}
\left[
\mathrm{APL}_g(\psi_g)
- \frac{\big(\psi_g - \psi^{\text{trend}}(\bar N_g)\big)^2}{2\hat\sigma_\psi^2}
\right].
}
$$

When $\hat\sigma_\psi^2 = 0$ all SNPs collapse onto the trend (full pooling); when $\hat\sigma_\psi^2 \to \infty$ the prior is uninformative and the MAP returns $\hat\psi_g^{\text{gw}}$. The realized weight of the prior is data-driven through the empirical-Bayes step.

This is structurally the same as `DESeq2`'s "shrunken dispersion" with the prior placed on $\log \alpha$ (their dispersion lives on the log scale, matching our $\psi$).

------------------------------------------------------------------------

## 5. Sequencing error $\varepsilon$

The error parameter $\varepsilon$ is shared across all SNPs, exactly as in the existing `fitQuasar2`. After updating $\boldsymbol{\beta}_g$ and $\psi_g$ we update $\varepsilon$ by maximizing the pooled log-likelihood:

$$
\hat\varepsilon = \arg\max_{\varepsilon \in (0,\, 0.1]}
\sum_{g=1}^{G}\sum_{i=1}^{n_g} \ell_{gi}\big(\hat{\boldsymbol{\beta}}_g, \hat\psi_g, \varepsilon\big),
$$

using a 1-D Brent search.

------------------------------------------------------------------------

## 6. Algorithm summary

Let $\hat\psi_g^{(t)}, \hat{\boldsymbol{\beta}}_g^{(t)}, \hat\varepsilon^{(t)}$ be the iterates.

1.  **Initialize.** Set $\hat\psi_g^{(0)} = \log(1/20) = -\log 20$ for all SNPs, $\hat\varepsilon^{(0)} = 10^{-3}$, and fit $\hat{\boldsymbol{\beta}}_g^{(0)}$ by BFGS at $(\hat\psi^{(0)}, \hat\varepsilon^{(0)})$.
2.  **Repeat** until convergence:
    -   **(Dispersion, gene-wise)** For each $g$, recompute $\hat\psi_g^{\text{gw}}$ by Brent maximization of $\mathrm{APL}_g$ (warm-start the inner $\boldsymbol{\beta}_g$ fit).
    -   **(Trend)** Smooth $\hat\psi_g^{\text{gw}}$ against $\log \bar N_g$ via locfit / spline.
    -   **(Prior variance)** Compute $\hat\sigma_\psi^2$ from residuals.
    -   **(MAP)** Compute $\hat\psi_g^{\text{MAP}}$ for all $g$.
    -   **(Error)** Update $\hat\varepsilon$ by pooled maximization.
    -   **(Coefficients)** Refit $\hat{\boldsymbol{\beta}}_g$ for all $g$ at the new $(\hat\psi_g^{\text{MAP}}, \hat\varepsilon)$.
    -   **(Convergence)** Stop when $\max_g \lvert \Delta\hat\psi_g \rvert < \tau_\psi$, $\max_g \lVert \Delta\hat{\boldsymbol{\beta}}_g \rVert_\infty < \tau_\beta$, and $\lvert \Delta\hat\varepsilon \rvert < \tau_\varepsilon$.

A single pass is usually enough in practice because the inner per-SNP BFGS already gives a high-quality $\hat{\boldsymbol{\beta}}_g$ for the current dispersion; we still iterate a few times to let $\varepsilon$ and the shrinkage stabilize.

------------------------------------------------------------------------

## 7. Inference: contrasts and the LRT

After fitting we have, for each SNP $g$:

-   $\hat{\boldsymbol{\beta}}_g$ (length $p_g$),
-   $\widehat{\mathrm{Cov}}(\hat{\boldsymbol{\beta}}_g) = \mathcal{I}_g^{-1}$ evaluated at $(\hat{\boldsymbol{\beta}}_g, \hat\psi_g^{\text{MAP}}, \hat\varepsilon)$.

### 7.1 Wald test for an arbitrary contrast

For a contrast matrix $\mathbf{L} \in \mathbb{R}^{k \times p_g}$ of rank $k$ (e.g. "Treatment_B $-$ Treatment_C", or several contrasts stacked), the Wald statistic is

$$
W_g = \big(\mathbf{L}\hat{\boldsymbol{\beta}}_g\big)^\top
       \left(\mathbf{L}\mathcal{I}_g^{-1}\mathbf{L}^\top\right)^{-1}
       \big(\mathbf{L}\hat{\boldsymbol{\beta}}_g\big)
\;\overset{H_0}{\sim}\; \chi^2_{k}.
$$

The single-coefficient case ($\mathbf{L} = \mathbf{e}_j^\top$) gives the usual $z = \hat\beta_{g,j} / \widehat{\mathrm{se}}(\hat\beta_{g,j})$.

### 7.2 Likelihood ratio test

For a reduced design $\mathbf{X}_g^{(r)}$ (a subset of columns of $\mathbf{X}_g$, or a formula with fewer terms), refit $\hat{\boldsymbol{\beta}}_g^{(r)}$ **at the same** $\hat\psi_g^{\text{MAP}}$ and $\hat\varepsilon$ as the full model. Then

$$
\Lambda_g = 2\Big[\ell_g(\hat{\boldsymbol{\beta}}_g, \hat\psi_g, \hat\varepsilon)
                - \ell_g(\hat{\boldsymbol{\beta}}_g^{(r)}, \hat\psi_g, \hat\varepsilon)\Big]
\;\overset{H_0}{\sim}\; \chi^2_{p_g - p_g^{(r)}}.
$$

Holding $\psi_g$ fixed at the full-model MAP is the standard `edgeR`/`DESeq2` choice and avoids spurious LRT inflation due to dispersion re-estimation under the null.

### 7.3 Multiple testing

Per-SNP $p$-values are adjusted with `p.adjust(..., method = "BH")` by default.

------------------------------------------------------------------------

## 8. What this changes vs. `fitQuasar2`

| Aspect | `fitQuasar2` (current) | `fitQuasar2CR` (new) |
|----|----|----|
| Dispersion parameter | $M$, bounded $[0.1, 10^4]$ | $\psi = \log(1/M)$, unbounded |
| Estimation | Per coverage-bin profile MLE | Per-SNP Cox–Reid APL |
| Sharing | Constant within ad-hoc quantile bin | Smooth trend + EB shrinkage |
| $\hat{\boldsymbol{\beta}}_g$ uncertainty in $\hat\psi_g$ | Ignored | CR adjustment |
| Inference | Single-coefficient $z$ from inner Hessian | Wald (any $\mathbf{L}$) + LRT |
| Contrasts | None native | `makeContrastsQuasar()` |
| Formula | `~ x` | `~ x` (unchanged) |

------------------------------------------------------------------------

## 9. References

-   Cox, D. R., & Reid, N. (1987). *Parameter orthogonality and approximate conditional inference.* Journal of the Royal Statistical Society, Series B, 49(1), 1–39.
-   McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). *Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.* Nucleic Acids Research, 40(10), 4288–4297.
-   Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 15:550.
-   Harvey, C. T., Moyerbrailean, G. A., Davis, G. O., Wen, X., Luca, F., & Pique-Regi, R. (2015). *QuASAR: quantitative allele-specific analysis of reads.* Bioinformatics, 31(8), 1235–1242.
-   Kalita, C. A., Moyerbrailean, G. A., Brown, C., Wen, X., Luca, F., & Pique-Regi, R. (2018). *QuASAR-MPRA: Accurate allele-specific analysis for massively parallel reporter assays.* Bioinformatics, 34(5), 787–794.
