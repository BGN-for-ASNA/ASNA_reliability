# ASNA Reliability

## 1. Objectives

The objective of this project is to assess the reliability of animal social network analysis (ASNA) methods, specifically focusing on two types of sampling biases:

  1.  Sampling effort bias: Variability in observations between individuals.
  2.  Interaction bias: Unobserved interactions between individuals during sampling.

## 2. Methods

### 2.1. Identifying issues with current simulations for generating sampling effort and interaction biases

Since 2000, the main standard for ASNA has been permutations. However, recent simulations have highlighted problems related to the rates of type I and type II errors. These simulations follow a scenario outlined by Farine in 2017, but even this scenario has some issues.
The main problem with the simulation is that when links to certain individuals are removed, entire observations are often lost, resulting in a new bias in sampling effort between different phenotypes in the population. Additionally, these simulations assume that females can be equally social or have higher sociality, but never lower sociality. Refer to the R script ["1. current simulation issue.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/1.%20current%20simulation%20issue.R) for more details.

### 2.2. Introducing a new simulation approach
To address these concerns, I have developed a new simulation that allows for the specification of differences between two categorical phenotypes of individuals, while independently controlling the degree of bias introduced in sampling effort and interactions. For more information, please see the R script ["2. generate_biased_network.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/2.%20generate_biased_network.R).

Cody has also created a simulation that provides more refined details and can generate both types of biases. This simulation is documented in the R script ["3. Data simulation.R"](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/3.%20Data%20simulation.R), which we will use for our publications.

With these new simulation approaches, we can now test various scenarios while controlling for sociality differences between phenotypes. The scenarios include:

  1.  No sampling effort bias and no interaction bias.
  2.  Transition from negative to positive sampling effort bias.
  3.  Transition from negative to positive interaction bias.
  4.  Transition from negative to positive sampling effort bias and interaction bias.

### 2.3. Reliability estimation

We employ different ASNA methods aimed at correcting for sampling biases, see ["4.simulation.R](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/4.%20Simulation.R). These methods include:

  1.  Rates of interactions per unit of time (commonly used by primatologists).
  2.  Proportion of time a dyad has been observed together, known as SRI (typically used by behavioral ecologists).
  3.  Bayesian approach developed by our team.

### 2.4. Bayesian approach developed by our team
Using rethinking approaches, we use [model dyadic connections](https://www.youtube.com/watch?v=XDoAglqd7ss&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=15&pp=iAQB&ab_channel=RichardMcElreath)  combine with [measurement error approach](https://www.youtube.com/watch?v=PIuqxOBJqLU&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=17&ab_channel=RichardMcElreath). The model can be summarized as follows, considering that predictors are centered:
<p align="center">
$X_{ij} \sim binomial(x_{ij}, \sigma_{ij})$
<p align="center">
$x_{ij} = logistic(\alpha + \beta_{i1}P_{i1} + \beta_{j1}P_{j1} +... + \beta_{in}P_{in} + \beta_{jn}P_{jn})$
<p align="center">
$\sigma_{ij} = \alpha_{2} + \beta_{i}S_{i} + \beta_{j}S_{j} +  ...$

Where $X_{ij}$ represents "true network", $x_{ij}$ represents the the observed network,  $\sigma_{ij}$ represents the measurement error and $P_{1}$ to $P_{n}$ are the predictors of research interest for individuals $i$ and $j$. This is the 'measurement error approach' where we model the "true network" as a function of the observed network and a variance based on the error measurement.

We can use 'dyadic connections model' to estimate ($\sigma_{ij}$) as the outcome of the sampling effort of individual $i$ ($S_{j}$) and individual $j$ ($S_{j}$). Additionally, we can include any additional phenothypic traits ($...$) that might influence link observations (e.g., sex, age, hierarchical rank of $i$ and $j$) while examining their effects on $x_{ij}$.

### 2.5 Interaction bias
We could introduce an extra layer to the Bayesian approach we've developed to account for this bias. However, unlike sampling effort bias, interaction bias remains unknown to the observer and can't be defined within the model. Therefore, we anticipate that neither approach can effectively address this bias.

Currently, two potential methods may estimate interaction bias:

1. Running telemetric experiments concurrently with observations.
2. Use $\sigma_{ij}$ results as a proxy for interaction bias. For instance, if gender affects sampling effort, we could assume it also influences interaction bias to the same degree. Yet, it's important to note that this approach may yield unreliable results.

## 3. Results

Script ["4.simulation.R](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/4.%20Simulation.R) provide two types of plots:
1.  Estimated effect size vs. true effect size for each approach.
2.  Difference between true effect size and estimated effect size.
![alt text](https://github.com/BGN-for-ASNA/ASNA_reliability/blob/main/Simulation%20results.png)

# 4.To-Do
1.  Simulation doesn't have interaction bias.
2.  Run simulations.
3.  Manuscript
