# Simulate Data for Classification with Parallel Factor Analysis

Simulates a two-way matrix, three-way array, or four-way array; and
simulates a set of class labels that are related to the simulated matrix
or array through one mode of that matrix or array. Data matrix or array
is simulated using one of the following component models without
constraints: Parafac, Parafac2, or PCA. Weights for mode weight matrices
can be drawn from 12 common probability distributions. Alternatively,
custom weights can be provided for any mode.

## Usage

``` r
simcpfa(arraydim = NULL, model = "parafac", nfac = 2, nclass = 2, 
        smethod = "logistic", nreps = 100, onreps = 10, props = NULL,
        corresp = NULL, meanpred = NULL, modes = 3, corrpred = NULL, 
        pf2num = NULL, Amat = NULL, Bmat = NULL, Cmat = NULL, Dmat = NULL, 
        Gmat = NULL, Emat = NULL, technical = list())
```

## Arguments

- arraydim:

  Numeric vector containing the number of dimensions for each mode of
  the simulated data matrix or array. Must contain integers greater than
  or equal to 2. Defaults to `c(10, 100)` for `modes = 2`; to
  `c(10, 10, 100)` for `modes = 3`; and to `c(10, 10, 10, 100)` for
  `modes = 4`.

- model:

  Character specifying the model to use for simulating the data array.
  Must be 'parafac', 'parafac2', or 'pca'.

- nfac:

  Number of components in the component model. Must be an integer
  greater than or equal to 1.

- nclass:

  Number of classes in simulated class labels. Must be an integer
  greater than or equal to 2.

- smethod:

  Simulation method, either "logistic" or "eigende". The former
  implements an iterative Monte Carlo rejection sampling method based on
  the generalized linear model (slower). The latter uses a multivariate
  normal distribution and constructs a joint covariance matrix,
  employing a eigendecomposition and assuming class labels arise from
  discretizing continuous latent variables (faster).

- nreps:

  Number of replications for simulating class labels for a given set of
  classification mode component weights.

- onreps:

  Number of replications for simulating a set of classification mode
  component weights.

- props:

  Target proportions for simulated class labels in output 'y'. Defaults
  to equal proportions across classes.

- corresp:

  Numeric vector of target correlations between simulated class labels
  and columns of the classification mode component weight matrix. Must
  have length equal to 'nfac'. Defaults to 0.5 for all components.

- meanpred:

  Numeric vector of means used to generate the classification mode
  component weights. Must be real numbers. Operates as the mean vector
  parameterizing a multivariate normal distribution from which
  classification mode component weights are generated. Length must be
  equal to input `nfac`. Defaults to a vector of zeros.

- modes:

  Single integer of 2, 3, or 4, indicating whether to simulate a two-way
  matrix, three-way array, or four-way array, respectively.

- corrpred:

  A positive definite correlation matrix containing the target
  correlations for the classification mode component weights. Must have
  number of rows and columns equal to input 'nfac'. Operates as the
  covariance matrix parameterizing a multivariate normal distribution
  from which classification mode component weights are generated.
  Defaults to a correlation matrix with 1 on the diagonal and 0.2 on the
  off-diagonals.

- pf2num:

  When `model = 'parafac2'`, number of rows for each simulated matrix in
  the list of matrices `Amat`. Replaces the first element of input
  `arraydim` because, for the Parafac2 model, the number of rows in each
  simulated matrix can vary. If not specified when `model = 'parafac2'`,
  defaults to
  `rep(c((nfac + 1), (nfac + 2), (nfac + 3)), length.out = arraydim[modes])`.

- Amat:

  When `model = 'parafac'`, a matrix of A mode weights with number of
  rows equal to the first element of input 'arraydim' and with number of
  columns equal to input 'nfac'. When `model = 'parafac2'`, a list with
  length equal to the last element of input 'arraydim', where each list
  element contains a matrix with number of rows of at least 2 and with
  number of columns equal to input 'nfac'. When provided, replaces a
  simulated `Amat`. When `model = 'pca'`, a matrix of loadings with
  number of rows equal to the second element of input 'arraydim' and
  with number of columns equal to input 'nfac'.

- Bmat:

  A matrix of B mode weights with number of rows equal to the second
  element of input 'arraydim' and with number of columns equal to the
  input 'nfac'. When provided, replaces a simulated `Bmat`. When
  `model = 'pca'`, a matrix of scores with number of rows equal to the
  first element of input 'arraydim' and with number of columns equal to
  input 'nfac'. If provided when `modes = 2`, `onreps` is reduced to one
  when `smethod` is `logistic`.

- Cmat:

  A matrix of C mode weights with number of rows equal to the third
  element of input 'arraydim' and with number of columns equal to the
  input 'nfac'. When provided, replaces a simulated `Cmat` when
  `modes = 4`. When `modes = 3`, replaces the simulated classification
  mode weight matrix. If provided when `modes = 3`, `onreps` is reduced
  to one when `smethod` is `logistic`. When `modes` is 2, this argument
  is ignored.

- Dmat:

  A matrix of D mode weights with number of rows equal to the fourth
  element of input 'arraydim' and with number of columns equal to the
  input 'nfac'. When `modes = 4`, replaces the simulated classification
  mode weight matrix. When `modes` is 2 or 3, this argument is ignored.
  If provided when `modes = 4`, `onreps` is reduced to one when
  `smethod` is `logistic`.

- Gmat:

  When `model = 'parafac2'`, a matrix of G mode weights with number of
  rows equal to input 'nfac' and with number of columns equal to input
  'nfac'. When provided, replaces a simulated `Gmat`.

- Emat:

  When `model = 'parafac'`, a 3-way or 4-way array containing noise to
  be added to the corresponding elements in the simulated data array.
  Error array dimensions must be equal to the values contained in
  `arraydim`. When `model = 'parafac2'`, a list containing either
  matrices (i.e., when `modes = 3`) or three-way arrays (i.e., when
  `modes = 4`) whose elements contain noise to be added to corresponding
  elements in the simulated data array. When provided, replaces a
  simulated `Emat`. When `model = 'pca'`, a matrix whose elements
  contain noise to be added to corresponding elements in the simulated
  data matrix. When provided, replaces a simulated `Emat`.

- technical:

  List containing arguments related to distributions from which to
  simulate data. When specified, must contain one or more of the
  following:

  distA

  :   List containing arguments specifying the distribution from which
      deviates are drawn for A mode weights contained in `Amat`.
      Defaults to standard normal distribution when not specified. See
      Details section for additional information on acceptable
      arguments.

  distB

  :   List containing arguments specifying the distribution from which
      deviates are drawn for B mode weights contained in `Bmat`.
      Defaults to standard normal distribution when not specified. See
      Details section for additional information on acceptable
      arguments. Ignored if `model = 'pca'`.

  distC

  :   For when `modes = '4'`, list containing arguments specifying the
      distribution from which deviates are drawn for C mode weights
      contained in `Cmat`. Defaults to standard normal distribution when
      not specified. See Details section for additional information on
      acceptable arguments. Ignored when `modes = 3`.

  distG

  :   For when `model = 'parafac2'`, list containing arguments
      specifying the distribution from which deviates are drawn for G
      weights contained in `Gmat`. Defaults to standard normal
      distribution when not specified. See Details section for
      additional information on acceptable arguments.

  distE

  :   List containing arguments specifying the distribution from which
      deviates are drawn for error contained in `Emat`. Defaults to
      standard normal distribution when not specified. See Details
      section for additional information on acceptable arguments.

## Details

By selecting `smethod = "logistic"`, the data array simulation consists
of two steps. First, a Monte Carlo simulation is conducted to simulate
class labels using a binomial logistic (i.e., in the binary case) or
multinomial logistic (i.e., in the multiclass case) regression model.
Specifically, columns of the classification mode weights matrix (e.g.,
`Cmat` when `modes = 3`) are generated from a multivariate normal
distribution with mean vector `meanpred` and with covariance matrix
`corrpred`. Values are then drawn randomly from a uniform or a normal
distribution and serve as beta coefficients. A linear combination of
these beta coefficients and the generated classification weights
produces a linear systematic part, which is passed through the logistic
function (i.e., the sigmoid) in the binary case or through the softmax
function in the multiclass case. Resulting probabilities are used to
assign class labels. The simulation repeats classification weights
generation `onreps` times and repeats class label generation, within
each `onreps` iteration, a total of `nreps` times. The generated class
labels that correlate best with the generated classification weights
(i.e., with correlations closest to `corresp`) are retained as the final
class labels with corresponding final classification weights. An
adaptive sampling technique is used during the simulation such that
optimal beta coefficients from previous iterations are used to
parameterize a normal distribution, from which new coefficients are
drawn in subsequent iterations. Note that, if any simulation replicate
produces a set of class labels where all labels are the same (i.e., have
no variance), that replicate is discarded. Note also that `onreps` is
ignored when the classification mode weight matrix (i.e., `Bmat` when
`modes = 2`, `Cmat` when `modes = 3`, or `Dmat` when `modes = 4`) is
provided; in this case, class labels are simulated with respect to the
provided classification mode weight matrix.

Second, depending on the chosen model (i.e., Parafac, Parafac2, or PCA)
specified via `model`, and depending on the number of modes specified
via `modes`, component matrices are randomly generated for each mode of
the data matrix or array. A data matrix or array is then constructed
using a Parafac, Parafac2, or PCA structure from these weight matrices,
including the generated classification mode weight matrix (i.e., `Bmat`,
`Cmat`, or `Dmat`) from the first step. Alternatively, weight matrices
can be provided to override random generation for any weight matrix with
the exception of the classification mode. When provided, weight matrices
are used to form the final data matrix or array. Finally, random noise
is added to each value in the matrix or array. The resulting output is a
synthetic data matrix or array paired, through one mode of that matrix
or array, with a simulated binary or multiclass response.

Alternatively, by selecting `smethod = "eigende"`, the function
simulates component weights and a latent response variable
simultaneously from a joint multivariate normal distribution using an
eigendecomposition. The covariance structure is defined by `corrpred`
and by a corrected version of `corresp` that accounts for the
attenuation in correlation caused by discretizing the latent response.
This continuous latent response vector is then discretized into class
labels using quantile cuts defined by the cumulative sum of `props`.
This method offers an efficient, non-iterative alternative that
satisfies the target covariance structure without rejection sampling.

The `technical` argument controls the probability distributions used to
simulate weights for different modes. Currently, `technical` is highly
structured. In particular, `technical` must be provided as a named list
whose elements must be one of 'distA', 'distB', 'distC', 'distG', or
'distE', with the last letter of each name designating a mode or, in the
case of 'distE', designating error. Each element provided must itself be
a list where the first inner list element is named 'dname', specifying
the distribution to be used to generate weights for a given mode or for
error. There are 12 'dname' options: 'normal', 'uniform', 'gamma',
'beta', 'binomial', 'poisson', 'exponential', 'geometric',
'negbinomial', 'hypergeo', 'lognormal', and 'cauchy'. Additional
arguments can be added to each inner list to parameterize the
probability distribution being used. These arguments can be one of the
following, for each distribution allowed:

For `dname = 'normal'`, allowed arguments are `mean` or `sd` (i.e.,
function `rnorm` is called).

For `dname = 'uniform'`, allowed arguments are `min` or `max` (i.e.,
function `runif` is called).

For `dname = 'gamma'`, allowed arguments are `shape` or `scale` (i.e.,
function `rgamma` is called).

For `dname = 'beta'`, allowed arguments are `shape1` or `shape2` (i.e.,
function `rbeta` is called).

For `dname = 'binomial'`, allowed arguments are `size` or `prob` (i.e.,
function `rbinom` is called).

For `dname = 'poisson'`, allowed argument is `lambda` (i.e., function
`rpois` is called).

For `dname = 'exponential'`, allowed argument is `rate` (i.e., function
`rexp` is called).

For `dname = 'geometric'`, allowed argument is `prob` (i.e., function
`rgeom` is called).

For `dname = 'negbinomial'`, allowed arguments are `size` or `prob`
(i.e., function `rnbinom` is called).

For `dname = 'hypergeo'`, allowed arguments are `m`, `n`, or `k` (i.e.,
function `rhyper` is called).

For `dname = 'lognormal'`, allowed arguments are `meanlog` or `sdlog`
(i.e., function `rlnorm` is called).

For `dname = 'cauchy'`, allowed arguments are `location` or `scale`
(i.e., function `rcauchy` is called).

Note that if a weight matrix and technical information are both provided
for a given mode (or for error), the weight matrix is used while
technical information is ignored. See Examples below for an example of
how to set up `technical`.

## Value

- X:

  Simulated data matrix or array with dimensions specified by `arraydim`
  and, when `model = 'parafac2'`, also by `pf2num`. When
  `model = 'parafac'`, `X` is an object of class 'array'. When
  `model = 'parafac2'`, `X` is an object of class 'list'. When
  `model = 'pca'`, `X` is an object of class 'matrix'.

- y:

  Simulated class labels provided as an object of class 'factor'. When
  `model = 'parafac'` or `model = 'parafac2'`, the number of labels is
  equal to the last element of `arraydim`. When `model = 'pca'`, the
  number of labels is equal to the first element of `arraydim`.

- model:

  Character value indicating the component model that was used to
  simulate the data matrix or array.

- Amat:

  Simulated A mode weights. When `model = 'parafac'`, output is a matrix
  with number of rows equal to the first element of `arraydim` and with
  number of columns equal to the number of components `nfac`. When
  `model = 'parafac2'`, output is a list of matrices with number of rows
  for each matrix equal to those specified by `pf2num` and with number
  of columns equal to `nfac`. When `model = 'pca'`, output is a matrix
  with number of rows equal to the second element of `arraydim`. If
  `Amat` was supplied, returns original `Amat` instead of a simulated
  `Amat`.

- Bmat:

  Simulated B mode weights provided as a matrix with number of rows
  equal to the second element of `arraydim` and with number of columns
  equal to the number of components `nfac`. When `model = 'pca'`, output
  is a matrix with number of rows equal to the first element of
  `arraydim`. If `Bmat` was supplied, returns original `Bmat` instead of
  a simulated `Bmat`.

- Cmat:

  Simulated C mode weights provided as a matrix with number of rows
  equal to the third element of `arraydim` and with number of columns
  equal to the number of components `nfac`. If `Cmat` was supplied when
  `modes = 4`, returns original `Cmat` instead of a simulated `Cmat`.
  Not provided when `modes = 2`.

- Dmat:

  Simulated D mode weights provided when `modes = 4`. Output is a matrix
  with number of rows equal to the fourth element of `arraydim` and with
  number of columns equal to the number of components `nfac`. Not
  provided when `modes = 2` or when `modes = 3`.

- Gmat:

  Simulated G weights provided when `model = 'parafac2'`. Provided as a
  matrix with number of rows and columns equal to `nfac`. If `Gmat` was
  supplied, returns original `Gmat` instead of a simulated `Gmat`.

- Emat:

  Error matrix, array, or list containing noise added to corresponding
  elements of simulated data matrix or array. Output has dimensions
  specified by `arraydim` and, when `model = 'parafac2'`, also by
  `pf2num`. When `model = 'parafac'`, `Emat` is an object of class
  'array'. When `model = 'parafac2'`, `Emat` is an object of class
  'list'. When `model = 'pca'`, `Emat` is an object of class 'matrix'.

## Author

Matthew Asisgress \<mattgress@protonmail.ch\>

## References

See help file for function `cpfa` for a list of references.

## Examples

``` r
########## Parafac2 example with 4-way array and multiclass response ##########
if (FALSE) { # \dontrun{
# set seed for reproducibility
set.seed(5)

# define list of arguments specifying distributions for A and G weights
techlist <- list(distA = list(dname = "poisson", 
                              lambda = 3),                 # for A weights
                 distG = list(dname = "gamma", shape = 2, 
                              scale = 4))                  # for G weights

# define target correlation matrix for columns of D mode weights matrix
cormat <- matrix(c(1, .6, .6, .6, 1, .6, .6, .6, 1), nrow = 3, ncol = 3)

# simulate a four-way ragged array connected to a response
data <- simcpfa(arraydim = c(10, 11, 12, 100), model = "parafac2", nfac = 3, 
                nclass = 3, nreps = 1e2, onreps = 10, corresp = rep(.6, 3), 
                meanpred = rep(2, 3), modes = 4, corrpred = cormat,
                technical = techlist, smethod = "eigende")
                
# examine correlations among columns of classification mode matrix Dmat
cor(data$Dmat)

# examine correlations between columns of classification mode matrix Dmat and
# simulated class labels
yproc <- as.numeric(data$y) - 1
cor(data$Dmat, yproc)
} # }
```
