# Introduction to Recursive Bayesian Inference

R repository I made on the study of sequential Monte Carlo (filtering). It's intended for educational purposes.

## Usage

### Folders

* **IRP**. This folder contains a project where I provide an overview of two-stage MCMC and why it would be advantageous to single-stage MCMC if we're doing Bayesian inference of large spatio-temporal ordinal data. It closely follows the comparisons study in "Two-stage MCMC for Fast Bayesian Inference of Large Spatio-temporal Ordinal Data, with Application to US Drought" by Hepler et al. (2025) [1].
* **RecursiveBayes** This folder is a prototype R package on recursive Bayesian inference. It allows the user to run prior-recursive and prior-proposal sequential Monte Carlo in a couple lines of code. Test the package by downloading the folder and installing the package "RecursiveBayes" in R.

<br>

### Notebooks

* **0_MCMC.Rmd** Introduces Markov Chain Monte Carlo (MCMC) and the adaptive Metropolis-Hastings algorithm. Written from scratch.
* **1_Prior_Recursive.Rmd** A recurisve take on MCMC where I consider batches of data sequentially. Follows the prior-recursive sequential Monte Carlo algorithm in Hooten et al. [2] whereby the prior for the next batch in the sequence uses a KDE defined on the posterior from the previous batch. Written from scratch.
* **2_PriorProposal.Rmd** Follows the prior-proposal sequential Monte Carlo (PP-RB) algorithm in Hooten et al. [2] whereby the prior is fixed and the current proposal is the samples from the previous posterior in the sequence. Written from scratch.
* **3_HMC.Rmd** Implements Hamiltonian Monte Carlo (HMC). Uses NIMBLE package [3].
* **3_HMC_notes.ipynb** Python jupyter notebook explaining the mathematics motivating HMC.
* **4_RP_nimble.Rmd** Same as 1, but uses HMC instead and the prior is a normal distribution defined by a mean and standard deviation of the posterior from the previous batch in the sequence.
* **5_PR_kde.Rmd** Same as 1, but uses HMC.

<br>

### References
1) [Hepler et al.](https://arxiv.org/pdf/2505.24594v1)

2) [Hooten et al.](https://www.tandfonline.com/doi/pdf/10.1080/00031305.2019.1665584)

3) [NIMBLE](https://r-nimble.org/)


## License


[MIT](https://choosealicense.com/licenses/mit/)

