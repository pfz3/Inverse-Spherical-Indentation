# Inverse Spherical Indentation

Methods employed in this repo are explained in full detail in our recent publication [Fernandez-Zelaia et al. 2018](https://www.sciencedirect.com/science/article/pii/S0264127518302168). Please appropriately reference that work and this repo if you decide to use these codes.

The experimental data utilized in these methods corresponds to indentation stress-strain curves as defined by Pathak and Kalidindi's original works [[1](https://www.sciencedirect.com/science/article/pii/S1359645408002413),[2](https://www.sciencedirect.com/science/article/pii/S1359646208008221)]. There have been a number of studies published since those early works and an excellent summary can be found in [[3](https://www.sciencedirect.com/science/article/pii/S0927796X15000157)].

## Methodology

The code in this repo solves the inverse indentation problem: given indentation data we wish to extract the uniaxial equivalent mechanical constitute properties. Since there are no closed form solutions which describe the post-elastic indentation response a complex finite element (FE) model may be used. FE models however can be costly - in our simulations an isotropic linear-hardening asymmetric model required 10-25 hrs per simulation to run. Therefore using the FE model directly in unfeasible. 

The solution methodology follows the strategy layed out in the seminal work by Kennedy and O'Hagan [[4](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00294)]. 

1. Build a Gaussian process (GP) model of the complex computer code from a reasonable number of runs
2. use a Bayesian strategy to perform the inversion. 

Our exact implementation is a bit different and closer to a standard regression analysis; Kennedy and O'Hagan sought to build a *predictive* model by *fusing* experiments with simulations. In the indentation setting we don't care about predicting future indentation responses - the experiments are a means to an ends which enable inference of the underlying intrinsic mechanical behavior. 

For additional detail please refer to our paper provided above.

## Surrogate Model

The surrogate model is built by evaluating the surrogate model at specified settings (the experimental design). A  constrained Maximum Projection design was utilized [[5](https://academic.oup.com/biomet/article-abstract/102/2/371/246859)]. Note a very nice R-package exists for MaxPro [[6](https://cran.r-project.org/web/packages/MaxPro/index.html)]. Constraint is imposed by using rejection sampling and generating the MaxPro design conditional on some sample population of candidate points satisfying the constraints - see the rejoinder corresponding to [[7](https://www.tandfonline.com/doi/abs/10.1080/08982112.2015.1100447)]. The constraint says that the hardening slope cannot be larger than the modulus (no doy). 

The surrogate was then programmed in the statistical programming language RStan [[8](http://mc-stan.org/users/interfaces/rstan)]. Everything is built from scratch because off the shelf GP packages couldn't do exactly what was needed. The model can be considered a multiple-output function - one finite element simulation produces many values of stress-strain pairs. Things get a bit burdensome if you don't find and exploit some structure so we utilized a Kronecker product formulation to make things efficient - see the paper. Further the inference part was a bit messy as well with some heteroskedasticity present in the data. Therefore the problem was sufficiently specific that investing in a tailored solution seemed appropriate.

The model was built using 30 FE evaluations and validated on another 50-run MaxPro design. A "fancy" mean function was utilized that gets very "close" to emulating the FE model itself. Adding the GP on top of that makes things even better. The validation is exceptionally good [see pic](https://github.com/pfz3/Inverse-Spherical-Indentation/blob/master/val.png). The RMSPE is ~6 GPa... Try to beat that! I double-dog-dare you!

## Example

The FE summary outputs are provided as well as some experimental data. The FE simulations were run over a design that is considered to be appropriate for some of my other research work. Therefore do not expect that you can use this for steels, Ti-alloys, etc.. It MIGHT work for some aluminum alloys. Just pay attention to the bounds/constraints used and think if it will work for you. If not, run your own FE simulations. All the scripts used to automate FE model submissions can be found in the folder "FE_scripting". This should get you started on manipulating .inp files, post processing results, submitting jobs in a HPC environment, etc..  

Data from a single indentation experiment is provided. The material indented is a copper sample subject to severe plastic deformation. The [posterior predicted curve](https://github.com/pfz3/Inverse-Spherical-Indentation/blob/master/experimental/002_mcmc_posterior_mean.png) fits the data exceptionally well. Further, note that in this code some additional complexity is included to account for the apparent [heteroskedastic noise behavior](https://github.com/pfz3/Inverse-Spherical-Indentation/blob/master/experimental/002_mcmc_residuals.png). Attempts to use common transformation (Box-Cox) and also a more modern transformation ([Yeo-Johnson](https://academic.oup.com/biomet/article-abstract/87/4/954/232908)) were made yet neither was able to mitigate the observed heteroskedasticity.