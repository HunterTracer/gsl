GSL-GM —— GNU Scientific Library with support for Gaussian Mixture distribution
============================

This is GSL-GM, which extends the functionality of GNU Scientific Library (GSL) to enable convenient modeling of Gaussian Mixture distributions.

About this repository
=====================

Many modifications have been made on GSL for the modeling of Gaussian Mixture distribution. The additions and modifications in directories include:

- `cdf`: `gaussmixture.c`, `gaussmixtureinv.c`, `gsl_cdf.h`, `Makefile.am`, `Makefile.in`
- `randist`: `gsl_randist.h`, `gaussmixture.c`, `Makefile.am`, `Makefile.in`
- `ampl/src`: `amplgsl.cc`

How to use GSL-GM for optimization modeling or calculation?
========================

Build with AMPL bindings
------------------------

To build the AMPL bindings: 

1. Make sure that the ASL submodule is initialized (in ```/ampl/thirdparty/asl```). Using the latest version of ASL is recommended. For obtaining latest ASL, just clone this repository and then run `git submodule update --init --recursive` and `git submodule update --remote --merge`.

2. Create a build directory and move there:

   ```
   mkdir build
   cd build
   ```

3. Initialize the build files with your desired [generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), for example:

   `cmake .. -G"Unix Makefiles"` or `cmake .. -G"Visual Studio 17 2022"`

4. Build the resulting build files accordingly, for example:

   `make .` or `cmake --build . --target install --config "Release"`


C Function List
----------------------

In the following part, the cumulative distribution function for the lower tail $P(x)$ is defined by the integral $P(x) = \int_{-\infty}^{x} dx' p(x')$, and gives the probability of a variate taking a value less than $x$.

The cumulative distribution function for the upper tail $Q(x)$ is defined by the integral, $Q(x) = \int_{x}^{\infty} dx' p(x')$, and gives the probability of a variate taking a value greater than $x$.

Here are the list of C functions.

`double gsl_ran_gaussian_mixture (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])`

- This function return a Gaussian Mixture random variate, which has `K` mixture components with mixture weights `w`, means `mu` and standard deviations `sigma`. The probability distribution for Gaussian random variates is $p(x)dx = \sum\limits_{k=1}^K {w_k \frac{1}{\sqrt{2 \pi \sigma_k^2}}} \exp (-\frac{(x -\mu_k)^2} {2\sigma_k^2}) dx$.

`double gsl_ran_gaussian_mixture_pdf (const double x, const size_t K, const double w[], const double mu[], const double sigma[])`

- This function computes the probability density $p(x)$ at $x$ for a Gaussian Mixture distribution, which has `K` mixture components with mixture weights `w`, means `mu` and standard deviations `sigma`, using the formula given above.

`double gsl_ran_gaussian_mixture_ziggurat (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])`

`double gsl_ran_gaussian_mixture_ratio_method (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])`

- This function computes a Gaussian Mixture random variate using the alternative Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods. The Ziggurat algorithm is the fastest available algorithm in most cases.

`double gsl_cdf_gaussian_mixture_P(const double x, const size_t K, const double w[], const double mu[], const double sigma[])`

`double gsl_cdf_gaussian_mixture_Q(const double x, const size_t K, const double w[], const double mu[], const double sigma[])`

`double gsl_cdf_gaussian_mixture_Pinv(const double P, const size_t K, const double w[], const double mu[], const double sigma[])`

`double gsl_cdf_gaussian_mixture_Qinv(const double Q, const size_t K, const double w[], const double mu[], const double sigma[])`

- These functions compute the cumulative distribution functions $P(x), Q(x)$ and their inverses for the Gaussian Mixture distribution, which has `K` mixture components with mixture weights `w`, means `mu` and standard deviations `sigma`.

Use GSL-GM in Python
----------------------

1. First, run `pip install pyomo` to install pyomo. Then, add the folder where `amplgsl.dll` is located to Environment Variable `PATH`.

2. For calculating the value, first-order derivatives, second-order derivatives of a given Gaussian Mixture CDF or PDF, or the value of inverse Gaussian Mixture CDF, the following code can be used as an example:

   ```
   import pyomo.core.base.external as external
   import pyomo.common.gsl as gsl
   
   dll = gsl.find_GSL()
   gm_pdf = external.ExternalFunction(library=dll, function="gsl_ran_gaussian_mixture_pdf")
   gm_cdf = external.ExternalFunction(library=dll, function="gsl_cdf_gaussian_mixture_P")
   gm_cdf_inv = external.ExternalFunction(library=dll, function="gsl_cdf_gaussian_mixture_Pinv")
   
   K = 3
   w = [0.3, 0.5, 0.2]
   mu = [-3., 1., 5.]
   sigma = [0.5, 1.0, 1.5]
   p = 0.9
   x = 2
   cdf_val, cdf_derivs, cdf_hes = gm_cdf.evaluate_fgh(args=(x, *w, *mu, *sigma), fgh=2)
   pdf_val, pdf_derivs, pdf_hes = gm_pdf.evaluate_fgh(args=(x, *w, *mu, *sigma), fgh=2)
   cdf_val, cdf_derivs_x, cdf_hes_x = gm_cdf.evaluate_fgh(args=(x, *w, *mu, *sigma), fgh=2, fixed=[False] + [True] * (3 * K))
   pdf_val, pdf_derivs_x, pdf_hes_x = gm_pdf.evaluate_fgh(args=(x, *w, *mu, *sigma), fgh=2, fixed=[False] + [True] * (3 * K))
   cdf_inv_val, _, _ = gm_cdf_inv.evaluate_fgh(args=(p, *w, *mu, *sigma), fgh=0)
   ```

   Here, `w`, `mu`, `sigma` can be also 1-D numpy arrays. `fgh=0` means that only function value is returned, `fgh=2` means that function value, first derivatives, and hessian matrix are returned, `fixed` indicates if the corresponding argument value is fixed.

3. For modeling Gaussian Mixture CDF or PDF in optimization problems and solve it via NLP or MINLP solvers, the following code can be used as an example:

   ```
   import pyomo.environ as pyo
   import pyomo.core.base.external as external
   import pyomo.common.gsl as gsl
   
   dll = gsl.find_GSL()
   model = pyo.ConcreteModel()
   model.gm_pdf = external.ExternalFunction(library=dll, function="gsl_ran_gaussian_mixture_pdf")
   model.gm_cdf = external.ExternalFunction(library=dll, function="gsl_cdf_gaussian_mixture_P")
   model.gm_cdf_inv = external.ExternalFunction(library=dll, function="gsl_cdf_gaussian_mixture_Pinv")
   
   w = [0.3, 0.5, 0.2]
   mu = [-3., 1., 5.]
   sigma = [0.5, 1.0, 1.5]
   p = 0.9
   
   model.x = pyo.Var(domain=pyo.Reals)
   model.cons = pyo.ConstraintList()
   model.cons.add(model.gm_cdf(model.x, *w, *mu, *sigma) >= p)
   model.obj = pyo.Objective(rule=lambda m: m.x, sense=pyo.minimize)
   
   # change executable to directory where ipopt.exe is
   optimizer = pyo.SolverFactory("ipopt", executable=r"C:\Program Files (x86)\IpOpt\bin\ipopt.exe")
   optimizer.options["linear_solver"] = "ma97"
   # change hsllib to directory where coinhsl.dll is
   optimizer.options["hsllib"] = r"C:\coinhsl\bin\coinhsl.dll"
   results = optimizer.solve(model, tee=True)
   assert (results.solver.status == pyo.SolverStatus.ok) and (
   	results.solver.termination_condition == pyo.TerminationCondition.optimal)
   print(pyo.value(model.x), model.gm_cdf_inv.evaluate_fgh(args=(p, *w, *mu, *sigma), fgh=0)[0])
   ```

   
