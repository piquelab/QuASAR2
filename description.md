---
editor_options:
  markdown:
    wrap: 72
output:
  pdf_document: default
  html_document:
    df_print: paged
---

QuASAR by default assumes that under the null hypothesis of no allelic
imbalance the reference and alternate allele read counts should be at
1:1 ratio. However, in MPRA, the proportion $r_l$ of the reference reads
is not necessarily $0.5$ across all the $l$ genetic variants, due to
differences in PCR amplification, as well as cloning and transformation
efficiencies.

Here, we have extended QuASAR to test for differences between the
proportion of reference reads in DNA $r_l$ and the proportion obtained
from RNA reads $\rho_l$. To reject the null hypothesis $\rho_l=r_l$, we
extend QuASAR's beta-binomial model. The observed reference $R_l$ and
alternate $A_l$ allele read counts at a given $l$ are modeled as:

$$\Pr\left(R_{l} | N_{l}, \psi_{l}, M_l\right) =
{N_{l} \choose R_{l}} \frac{\Gamma \left( M_l \right)\Gamma \left( R_{l}+\psi_{l} M_l \right)\Gamma \left( A_{l} + \left( 1-\psi_{l} \right) M_l \right)}{\Gamma\left( N_{l} + M_l \right) \Gamma \left( \psi_{l} M_l \right) \Gamma \left(\left(1-\psi_{l}\right) M_l\right)} 
$$

$$
\psi_{l} = \left[\rho_l(1-\epsilon)+(1-\rho_l)\epsilon\right]
$$

where $N_l=R_l+A_l$ is the total read count at $l$, and $M_l$ is the
concentration parameter that controls over-dispersion of the mean
proportion centered around $\psi_{l}$, which also incorporates in the
model a base-calling error $\epsilon$ and the underlying allelic ratio
$\rho_l$. We can estimate $\epsilon$ using an EM procedure
\citep{Harvey2014}, but here for MPRA we fixed $\hat \epsilon=0.001$ as
a conservative estimate of the true base-calling error rate.

We have found previously for ASE that overdispersion decreases with
greater depth of coverage (Figure S9 in \citet{Moyerbrailean2016a}).
Therefore here, as compared to our previous implementation of QuASAR, we
use different $M_l$ parameters depending on the sequencing depth $N_l$.
We bin $N_l$ into different quantiles (here deciles) and we estimate
$M_l=M_b$ for each bin separately using a grid search:

$$
{\hat M_b}
    =\arg\max_{M_b}\left( \prod_{l=1}^L\Pr\left( R_{l} | N_{l}, \hat \epsilon ,\rho_{l}=r_l, M_b\right) \right) 
$$

This should work well when the number of sites (i.e., SNPs tested) is
relatively large so each bin $b$ has $>$ 200 observations to estimate
$M_b$. In our experience sequencing depth is a major determinant for M,
and because we estimate M under the null, we tend to be conservative
(i.e., M is the worst case scenario for all the constructs that belong
to the same group). As a consequence, the QuASAR-MPRA $p$-values remain
well calibrated (or in the worst case scenario they will tend to be
slightly conservative).

We estimate ${\hat \rho}_{l}$ using with $M_b=\hat M_b$ from
\eqref{eq:Mhat} and a standard gradient method (L-BFGS-B) to maximize
the log-likelihood function:

$$
l(\rho_{l};\hat M_b, \hat \epsilon) 
    = \log \Pr\left(R_{l} | N_{l}, \psi_{l}=\psi(\rho_{l},\hat \epsilon_s), \hat M_b\right) 
$$

Finally, all parameters are used to calculate the LRT statistic,
contrasting $H_1: \rho_{l}=\hat \rho _{l}$ to $H_0: \rho_{l} = r_{l}$
and the resulting $p$-value.

## Extension to multiple replicates

Now we want to extend this if we have multiple replicates across
conditions. Let's assume we have $K$ replicates indexed with $k$, and
that we measure the same $l$ SNP genetic variants. In this case, we will
model the observed reference $R_{lk}$ and alternate $A_{lk}$ allele read
counts at a given $l$ for the $k$ samples again with a beta-binomial:

$$\Pr\left(R_{lk} | N_{lk}, \psi_{lk}, M_l\right) =
{N_{lk} \choose R_{lk}} \frac{\Gamma \left( M_l \right)\Gamma \left( R_{lk}+\psi_{lk} M_l \right)\Gamma \left( A_{lk} + \left( 1-\psi_{lk} \right) M_l \right)}{\Gamma\left( N_{lk} + M_l \right) \Gamma \left( \psi_{lk} M_l \right) \Gamma \left(\left(1-\psi_{lk}\right) M_l\right)} 
$$

$$
\psi_{lk} = \left[\rho_{lk}(1-\epsilon)+(1-\rho_{lk})\epsilon\right]
$$

In this case $M_l$ is the concentration parameter that controls
over-dispersion of the mean proportion centered around $\psi_{lk}$,
which also incorporates in the model a base-calling error $\epsilon$ and
the underlying allelic ratio $\rho_{lk}$. The concentration parameter
$M_l$ is assumed to perhaps change by SNP but be the same across
replicates, and the \$\\epsilon\$ parameter will be fixed across $l$ and
$k$. Now we can also also model the allelic ratio using a logistic
model:
$\log(\rho_{lk}/(1-\rho_{lk}))= \log(r_{l}/(1-r_{l})) + \sum_j \beta_{jl} X_{kj}$

Where $r_{l}$ is the allelic ratio from the DNA proportion, which is
what is expected under the null hypothesis. For regular ASE in diploid
heterozygote SNP $r_l=0.5$ $\log(r_{l}/(1-r_{l}))=0$. If we are just
testing for the same amount of ASE across all $k$ replicates then
$\log(\rho_{lk}/(1-\rho_{lk}))= \log(r_{l}/(1-r_{l}))+ \beta_l$ and we
would do inference on $\beta_l$. In general $X_{kj}$ could be any
reasonable design matrix for a logistic model, and the parameters can be
used to test for treatment specific effects, etc.

Again we can fix $\epsilon$ and we need to estimate the $M_l$ and the
$\beta_{jl}$. We can also use this model to make a Likelihood Ratio Test
to test hypotheses about the $\beta_{jl}$ parameters. The $M_l$
parameters also need to be estimated. We can use the same approach that
estimates this parameters under the null (assuming $\rho{lk}=r_{lk}$,
and share them across $k$:

$$
{\hat M_b}
    =\arg\max_{M_b}\left( \prod_{k=1}^K \prod_{l=1}^L\Pr\left( R_{lk} | N_{lk}, \hat \epsilon ,\rho_{lk}=r_{lk}, M_b\right) \right) 
$$

\subsection{QuASAR meta-analysis}

Using the QuASAR approach, we can generate summary statistics of the
allelic imbalance that can be used for downstream analyses. For example,
to compare DNA to RNA, or between RNA of different cell-types, or to
perform meta-analysis of multiple MPRA libraries. Instead of using an
estimate of the allelic proportion $\rho_l$, in the QuASAR approach we
report the estimate of $\beta_l = \log(\rho_l/(1-\rho_l))$ and its
standard error $\hat \sigma_{l}$ using the second derivative (i.e.
Hessian) of the log-likelihood function in \eqref{eq:lrtrho}. We prefer
the logistic transformed parameter $\beta_l$ as it provides a more
robust fit and the second derivative is better behaved than that of
$\rho_l$ on the edges.

To illustrate this for the \citeauthor{Tewhey2016} data, we combined the
summary statistics for the two LCL individuals using standard fixed
effects meta-analysis. The effect size $\beta_{l,n}$ of each replicate
$n$ is weighted by $w_{n,l} = 1/{\hat \sigma_{n,l}}^2$, to calculate the
overall effect size and standard error: \begin{align}
    \beta_{l}^* = \frac{1}{w_l^*} \sum_n{\beta_{n,l}\,w_{n,l}} & & \sigma_l^*=\sqrt{1/w_l^*}  \label{eq:metabeta}
\end{align} where $w_l^*= \sum_n{w_n,l}$. We can then calculate the
$Z$-score and $p$-value to test for an overall change between all the
RNA replicates combined with respect to the original DNA proportion
$\beta_0$: \begin{align}
    Z_l = \frac{  \beta_{l}^* - \beta_0}{\sigma_l^*} \,,\, & \beta_0=\log{\frac{r_l}{1-r_l}}  \,, & p = 2 \Phi(-|Z_l|)
\end{align}

Across all the paper, $p$-values were corrected for multiple testing
using the Benjamini-Hochberg's (BH) method \citep{Benjamini1995}. To
compare the different approaches we quantified the genomic inflation
parameter, $\lambda$, for a set of $p$-values \citep{Yang2011}. For this
we calculated the ratio of the median of the $p$-value distribution to
the expected median, thus quantifying the extent of the bulk inflation
and the excess false positive rate. We also use a rank sum paired test
to assess statistical significance in the p-value inflation between
QuASAR-MPRA and other methods with similar performance.
