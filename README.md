# Markov Chain Monte Carlo methods in Bayesian Machine Learning

## Abstract

Markov chain Monte Carlo (MCMC) is a commonly used method for inference in Bayesian
statistics and machine learning. In this project, we compare the performance of popular
MCMC algorithms on different families of Bayesian machine learning tasks involving both
simulated and real datasets. We first present results on linear regression tasks, a simple
baseline setting where all algorithms successfully converge to the stationary distribution.
We then present results on more difficult inference tasks involving neural network models,
which require MCMC algorithms that exploit gradient information. These results show that
Hamiltonian Monte Carlo (HMC) gives the best and most reliable performance on these tasks.
Simpler and less compute-intensive gradient based algorithms can sometimes be successfully
be applied, but they are however more limited in scalability and reliability.


## Introduction

This project conducts a study and evaluation among different MCMC methods namely, the
Metropolis Hastings (MH), Gibbs sampling, Halmitonian Monte Carlo (HMC) (Liu et al., 2021),
the Metropolis-Adjusted Langevin algorithm (MALA) (Roberts and Tweedie, 1996), and the unadjusted Langevin algorithm (ULA) (Parisi, 1981). The MCMC methods are used for obtaining
a sequence of random samples to a probability distribution from which direct sampling is impossible. The algorithms generate samples by constructing a Markov Chain that has the desired
distributions as its stationary distribution (Neal, 2012).

## Objective

The goal of this project is to derive, implement and compare different MCMC methods in
Bayesian machine learning tasks such as nonlinear regression with Bayesian neural networks
(BNNs) and time series prediction with BNNs. Such neural networks do not only enjoy the
expressive power of traditional neural networks but also benefit from the advantages of the
Bayesian approach for inference such as flexible uncertainty quantification. The project will also
evaluate both quantitatively and qualitatively how well samples from MCMC chain represent
samples from posterior distribution of interest. The algorithms are also compared by measuring
the resulting models accuracy on validation data
