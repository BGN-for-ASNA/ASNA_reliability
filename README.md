# ASNA Reliability

## 1. Objectives

The objective of this project is to assess the reliability of animal social network analysis (ASNA) methods, specifically focusing on two types of sampling biases:

  1.  Sampling effort bias: Variability in observations between individuals.
  2.  Interaction bias: Unobserved interactions between individuals during sampling.

### 1.1. Identifying issues with current simulations for generating sampling effort and interaction biases

Since 2000, the main standard for ASNA has been permutations. However, recent simulations have highlighted problems related to the rates of type I and type II errors. These simulations follow a scenario outlined by Farine in 2017, but even this scenario has some issues.
The main problem with the simulation is that when links to certain individuals are removed, entire observations are often lost, resulting in a new bias in sampling effort between different phenotypes in the population. Additionally, these simulations assume that females can be equally social or have higher sociality, but never lower sociality. Refer to the R script "1. current simulation issue.R" for more details.

###1.2. Introducing a new simulation approach
To address these concerns, I have developed a new simulation that allows for the specification of differences between two categorical phenotypes of individuals, while independently controlling the degree of bias introduced in sampling effort and interactions. For more information, please see the R script "2. generate_biased_network.R".

Cody has also created a simulation that provides more refined details and can generate both types of biases. This simulation is documented in the R script "3. cody_simulation.R," which we will use for our publications.

With these new simulation approaches, we can now test various scenarios while controlling for sociality differences between phenotypes. The scenarios include:

  1.  No sampling effort bias and no interaction bias.
  2.  Transition from negative to positive sampling effort bias.
  3.  Transition from negative to positive interaction bias.
  4.  Transition from negative to positive sampling effort bias and interaction bias.

### 1.3. Reliability estimation

We employ different ASNA methods aimed at correcting for sampling biases. These methods include:

  1.  Rates of interactions per unit of time (commonly used by primatologists).
  2.  Proportion of time a dyad has been observed together, known as SRI (typically used by behavioral ecologists).
  3.  Bayesian approach developed by our team.

### 1.4. Bayesian approach developed by our team
Using rethinking approaches, we model dyadic connections (as discussed in video 15) utilizing a measurement error approach (as explained in video 17). The model can be summarized as follows, considering that predictors are centered:
