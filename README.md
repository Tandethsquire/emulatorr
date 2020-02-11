# emulatorr
Bayes Linear Emulation Package (R)

## A package for Bayes Linear Emulation and History Matching ##

### Background ###
In many fields, it is necessary to build a complex model to represent real-life systems. Such simulations can encompass a large, or high-dimensional, input space, and the simulator can be slow to evaluate over points in the space. Examples include [simulation of galaxy formation](http://dro.dur.ac.uk/15070/), [HIV transmission models](http://dro.dur.ac.uk/22952/), or [hormonal crosstalk in biological systems](https://arxiv.org/abs/1801.01538). In these, and other, situations, it is useful to be able to generate large portions of the simulator space without explicit recourse to the simulator.

An emulator is a statistical model of the simulator, using a relatively small number of sample points taken from simulator runs. The broad, simplistic form of the emulator is \(f(x) = g(x)\beta + u(x)\), where the \(g(x)\) and \(\beta\) create a regression surface, and the \(u(x)\) encodes correlation between points in the input space. In this way, the emulator represents the simulator results, but also has a built-in understanding of the uncertainties of both the emulation and of the simulator uncertainties. The value of the emulator at a point x is given by \(\text{E}[f(x)]\), with associated variance \(\text{Var}[f(x)]\).

Given such an emulator structure, outputs from the simulator can be used to 'train' the emulator using Bayesian methods. In this package, we forego a full Bayesian approach and instead focus on Bayes linear updates (for details, see e.g. [this book](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470065662)). This has the advantage of being quick to evaluate, and does not require us to supply full probabilistic specifications for all parameters of the emulator.

Often, we don't care about the structure of the emulator (or the simulator), but instead the following question: **what input values would give us agreement with observations?** To answer this question, we can use History Matching. We define an *implausibility measure* I(x), which in some sense indicates how likely a set of input values are to give a particular output. By applying a suitable cut-off to this measure, we can generate a non-implausible space for a desired output, from which we can sample more points, from which we can train a better emulator. The process can be repeated until either:
- We have enough non-implausible points to make an informed statement about the structure of the model;
- The emulator uncertainty is smaller than the uncertainties inherent in the system (e.g. observation error, model discrepancy)

The method of Bayes Linear emulation and history matching is a powerful tool to efficiently and robustly explore large, complex, time-intensive models, and find informative regions of input space that could potentially be overlooked by recourse only to the simulator.

### Package Details ###
**Note: the package is a work in progress! Use with care.**

The package consists of a number of objects (written as R6 classes), and functions that use these objects. For details of usage, see the associated help files and examples therein. The following is a summary.

 The R6 classes are
- ``Correlator``: Encodes the correlation structure for an emulator. A common correlation structure is \(\sigma^2 c(x,x')\), where \(c(x,x')\) is an exponential-squared function \(\exp(-(x-x')^2/\theta^2)\) for some correlation length \(\theta\). This is the default correlation structure in this package (implemented as ``exp_sq(x,xp,theta)``); more will be forthcoming, or user-defined functions can be added.
- ``Emulator``: This is the key object in the package. Given the second-order specifications for \(\beta\) and \(u(x)\) (namely their expectations, variances, and covariances; some of this will come from the ``Correlator`` object), and the basis functions \(g(x)\), an emulator is generated. From this we can find expectation and variance of the emulator at input points; we can use data points from a simulator to train new emulators based on these prior specifications; and we can determine implausibility for a given input against a desired output.

Associated functions (not an exhaustive list):
- ``emulator_from_data``: Generates basic univariate emulators from simulator data, along with some simple prior information
- diagnostic functions (``standard_error``, ``comparison_diagnostics``, ...): Shows diagnostic tests on an emulator given a validation set of data
- ``nth_implausible``: For a collection of outputs, and univariate emulators trained to each, finds the \(n\)-th maximum implausibility for an input point against a desired set of outputs.

### Installation ###
Not sure! At present, using RStudio with ``devtools`` enabled, the command ``devtools::install_github("Tandethsquire/emulatorr")`` should install it. But this is currently a private repository, so it may not work.

### Comments and Suggestions ###
Suggestions for improvement, and comments on functionality, are welcome. For preference these should be addressed via pull requests and/or issue logging, but if in doubt email andrew.iskauskas@durham.ac.uk.
