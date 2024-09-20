# ASNA Reliability

## 1. Objectives

The objective of this project is to assess the reliability of animal social network analysis (ASNA) methods, specifically focusing on two types of sampling biases:

  1.  Sampling effort bias: Variability in observations between individuals.
  2.  Interaction bias: Unobserved interactions between individuals during sampling.

## 2. Methods

### 2.1. Identifying issues with current simulations for generating sampling effort and interaction biases

Since 2000, the main standard for ASNA has been permutations. However, recent simulations have highlighted problems related to the rates of type I and type II errors. These simulations follow a scenario outlined by Farine in 2017, but even this scenario has some issues.
The main problem with the simulation is that when links to certain individuals are removed, entire observations are often lost, resulting in a new bias in sampling effort between different phenotypes in the population. Additionally, these simulations assume that females can be equally social or have higher sociality, but never lower sociality. Refer to the R script ["Appendix 1.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/1.%20Codes/1.%20current%20simulation%20issue.R) for more details.

### 2.2. Introducing a new simulation approach
To address these concerns, I have developed a new simulation that allows for the specification of differences between two categorical phenotypes of individuals, while independently controlling the degree of bias introduced in sampling effort and censoring. For more information, please see the R script  ["Support_Functions.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/1.Codes/Support_Functions.R).

With these new simulation approaches, we can now test various scenarios while controlling for sociality differences between phenotypes. The scenarios include:

  1.  No sampling effort bias and no censoring bias.
  2.  Transition from negative to positive sampling effort bias.
  3.  Transition from negative to positive censoring bias.
  4.  Transition from negative to positive sampling effort bias and censoring bias.

### 2.3. Reliability estimation

We employ different ASNA methods aimed at correcting for sampling biases, see ["Support_Functions.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/1.Codes/Support_Functions.R).

  1.  Rates of interactions per unit of time (commonly used by primatologists).
  2.  Proportion of time a dyad has been observed together, known as SRI (typically used by behavioral ecologists).
  3.  Bayesian approach : STRAND, BISONR, AMEN.
  4.  An extension of STRAND that we developped : STRAND ME.

### 2.4. Bayesian approach developed by our team
We use [dyadic connections model](https://www.youtube.com/watch?v=XDoAglqd7ss&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=15&pp=iAQB&ab_channel=RichardMcElreath)  combine with [measurement error approach](https://www.youtube.com/watch?v=PIuqxOBJqLU&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=17&ab_channel=RichardMcElreath). The model can be summarized as follows, considering that predictors are centered:
<p align="center">
$X_{ij} \sim Binomial(x_{ij}, N_{ij})$
<p align="center">
$x_{ij} = logistic(\alpha + \beta_{i1}P_{i1} + \beta_{j1}P_{j1} +... + \beta_{in}P_{in} + \beta_{jn}P_{jn})$
<p align="center">
$N_{ij} = S_{i} + S_{j} +  ...$

Where $X_{ij}$ represents "true network", $x_{ij}$ represents the the observed network,  $N_{ij}$ represents the exposure for individuals $i$ and $j$ and $P_{1}$ to $P_{n}$ are the predictors of research interest for individuals $i$ and $j$. This is the 'measurement error approach' where we model the "true network" as a function of the observed network and a variance based on the error measurement.

We can use 'dyadic connections model' to estimate ($\sigma_{ij}$) as the outcome of the sampling effort of individual $i$ ($S_{j}$) and individual $j$ ($S_{j}$). Additionally, we can include any additional phenothypic traits ($...$) that might influence interaction observations (e.g., sex, age, hierarchical rank of $i$ and $j$) while examining their effects on $x_{ij}$.

### 2.5 Censoring bias
In contrast to sampling effort bias, interaction bias is not known to the observer. Thus, we expect that neither approach can adequately tackle this bias. To address this, we can add an additional layer to the Bayesian approach outlined above:

<p align="center">
$\hat{X_{ij}} \sim Binomial(X_{ij}, \hat{\theta_{ij}})$
  <p align="center">
$\hat{\theta_{ij}} = \theta{ij}(1- \phi_{i} \phi_{j}) $

Where $\hat{X_{ij}}$ represents the "true network" defined by the previous model $X_{ij}$, including a variance denoted by $\hat{\theta_{ij}}$ that indicates the probability of missing interactions between individuals $i$ and $j$, which can be expressed as $\phi_{i} \phi_{j}$.

However, as previously stated, this probability of missing interactions is unknown. Researchers have several potential methods to estimate $\phi_{i} \phi_{j}$:

1. Conduct telemetric experiments concurrently with observations to estimate $\phi_{i} \phi_{j}$.
2. Utilize mark-recapture information to estimate $\phi_{i} \phi_{j}$.
3. Calculate the mean duration of individuals' focal points to estimate $\phi_{i} \phi_{j}$.

If none of this information is available, researchers can conduct a sensitivity analysis to assess the extent of bias in interactions needed to unvalidate the hypothesis.

