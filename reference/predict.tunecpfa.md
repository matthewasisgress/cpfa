# Predict Method for Tuning for Classification with Parallel Factor Analysis

Obtains predicted class labels from a 'tunecpfa' model object generated
by function `tunecpfa`.

## Usage

``` r
# S3 method for class 'tunecpfa'
predict(object, newdata = NULL, newdata.z = NULL, 
        method = NULL, type = c("response", "prob", "classify.weights"), 
        threshold = NULL, ...)
```

## Arguments

- object:

  A fit object of class 'tunecpfa' produced by function `tunecpfa`.

- newdata:

  An optional two-way matrix, three-way array, or four-way data array
  used to predict Parafac, Parafac2, or PCA component weights using
  estimated Parafac, Parafac2, or PCA model component weights from the
  input object. For Parafac2, can be a list of length `K` where the
  `k`-th element is a matrix or three-way array associated with the
  `k`-th observation. Matrix, array, or list must contain only real
  numbers. Dimensions must match dimensions of original data for all
  modes except the classification mode. If omitted, the original data
  are used.

- newdata.z:

  An optional two-way matrix containing additional features used to
  predict class labels. If omitted, and original `z` was NULL, defaults
  to NULL. Otherwise, if provided, this argument must have number of
  features equal to the number of columns of the original `z`. As such,
  the original `z` must not have been NULL.

- method:

  Character vector indicating classification methods to use. Possible
  methods include penalized logistic regression (PLR); support vector
  machine (SVM); random forest (RF); feed-forward neural network (NN);
  regularized discriminant analysis (RDA); and gradient boosting machine
  (GBM). If none selected, default is to use methods used in original
  `tunecpfa` run.

- type:

  Character vector indicating type of prediction to return. Possible
  values include: (1) `"response"`, returning predicted class
  labels; (2) `"prob"`, returning predicted class probabilities; or (3)
  `"classify.weights"`, returning predicted component weights used for
  classification in the specified component models. Must be specified.

- threshold:

  For binary classification, value indicating prediction threshold over
  which observations are classified as the positive class. If not
  provided, calculates threshold using class proportions in original
  data. For multiclass classification, `threshold` is not currently
  implemented.

- ...:

  Additional predict arguments. Currently ignored.

## Details

Predicts class labels for a binary or a multiclass outcome.
Specifically, predicts component weights for one mode of a Parallel
Factor Analysis-1 (Parafac) model, one mode of a Parallel Factor
Analysis-2 (Parafac2) model, or scores from a PCA model, using new data
and previously estimated mode weights from original data. Passes
predicted component weights (or scores for PCA) to one or several
classification methods as new data for predicting class labels.

Tuning parameters optimized by k-fold cross-validation are used for each
classification method (see help for `tunecpfa`). If not supplied in
argument `threshold`, prediction threshold for all classification
methods is calculated using proportions of class labels for original
data in the binary case (and the positive class proportion is set as the
threshold). For multiclass case, class with highest probability is
chosen.

## Value

Returns one of the following, depending on the choice for argument
`type`:

- type = "response":

  A data frame containing predicted class labels for each component
  model and classification method selected (see argument `type`). Number
  of columns is equal to number of methods times number of component
  models. Number of rows is equal to number of predicted observations.

- type = "prob":

  A list containing predicted probabilities for each component model and
  classification method selected (see argument `type`). The number of
  list elements is equal to the number of methods times the number of
  component models.

- type = "classify.weights":

  List containing predicted component weights for each component model.
  Length is equal to number of component models that were fit.

## Note

Applications of this function to real datasets can be explored at the
following repository:
<https://github.com/matthewasisgress/multiway-classification/>.

## Author

Matthew Asisgress \<mattgress@protonmail.ch\>

## References

See help file for function `tunecpfa` for a list of references.

## Examples

``` r
########## Parafac example with 3-way array and binary response ##########
if (FALSE) { # \dontrun{
# set seed and simulate a three-way array related to a binary response
set.seed(5)

# define list of arguments specifying distributions for A and G weights
techlist <- list(distA = list(dname = "poisson", 
                              lambda = 3),                 # for A weights
                 distG = list(dname = "gamma", shape = 2, 
                              scale = 4))                  # for G weights

# define target correlation matrix for columns of C mode weights matrix
cormat <- matrix(c(1, .6, .6, .6, 1, .6, .6, .6, 1), nrow = 3, ncol = 3)

# simulate a three-way array connected to a response
data <- simcpfa(arraydim = c(11, 12, 100), model = "parafac", nfac = 3, 
                nclass = 2, nreps = 1e2, onreps = 10, corresp = rep(.6, 3), 
                meanpred = rep(2, 3), modes = 3, corrpred = cormat,
                technical = techlist, smethod = "eigende")

# initialize
alpha <- seq(0, 1, length = 2)
gamma <- c(0, 0.01)
cost <- c(1, 2)
method <- c("PLR", "SVM")
family <- "binomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost)
model <- "parafac"
nfolds <- 3
nstart <- 3

# constrain first mode weights to be orthogonal
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tunecpfa(x = data$X[, , 1:80], y = data$y[1:80], 
                        model = model, nfac = 3, nfolds = nfolds, 
                        method = method, family = family, 
                        parameters = parameters, parallel = FALSE, 
                        const = const, nstart = nstart)
                    
# predict class labels
predict.labels <- predict(object = tune.object, newdata = data$X[, , 81:100], 
                          type = "response")

# print predicted labels
predict.labels
} # }
```
