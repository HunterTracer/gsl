/* cdf/gaussmixtureinv.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the inverse guass mixture cumulative distribution function 
 * according to Newton method
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_test.h>

 /* stopping parameters */
const double EPSREL = (10 * GSL_DBL_EPSILON);
const double EPSABS = (10 * GSL_DBL_EPSILON);

typedef struct                /* struct for gaussian mixture cdf inverse */
{
  double P;
  size_t K;
  const double* w;
  const double* mu;
  const double* sigma;
} gsl_gaussian_mixture_params_t;

double
gsl_cdf_gaussian_mixture_P_f(const double x, void* params)
{
  gsl_gaussian_mixture_params_t* gm_params = (gsl_gaussian_mixture_params_t*)(params);
  return gsl_cdf_gaussian_mixture_P(x, gm_params->K, gm_params->w, gm_params->mu, gm_params->sigma) - gm_params->P;
}

double
gsl_cdf_gaussian_mixture_P_df(const double x, void* params)
{
  gsl_gaussian_mixture_params_t* gm_params = (gsl_gaussian_mixture_params_t*)(params);
  return gsl_ran_gaussian_mixture_pdf(x, gm_params->K, gm_params->w, gm_params->mu, gm_params->sigma);
}

void
gsl_cdf_gaussian_mixture_P_fdf(const double x, void* params, double* f, double* df)
{
  *f = gsl_cdf_gaussian_mixture_P_f(x, params);
  *df = gsl_cdf_gaussian_mixture_P_df(x, params);
}

double
gsl_cdf_gaussian_mixture_Pinv(const double P, const size_t K, const double w[], const double mu[], const double sigma[])
{
  int status;
  double prev = 0.0;
  gsl_root_fsolver* fsolver;
  gsl_root_fdfsolver* fdfsolver;
  gsl_function f;
  gsl_function_fdf fdf;
  gsl_gaussian_mixture_params_t params;
  double x_lower = mu[0];
  double x_upper = mu[0];
  double cdf_x_lower;
  double cdf_x_upper;
  int i;

  if (P == 1.0)
  {
    return GSL_POSINF;
  }
  else if (P == 0.0)
  {
    return GSL_NEGINF;
  }

  for (i = 1; i < K; ++i) {
    if (mu[i] < x_lower)
      x_lower = mu[i];
    else if (mu[i] > x_upper)
      x_upper = mu[i];
  }
  cdf_x_lower = gsl_cdf_gaussian_mixture_P(x_lower, K, w, mu, sigma);
  cdf_x_upper = gsl_cdf_gaussian_mixture_P(x_upper, K, w, mu, sigma);

  if (P <= cdf_x_lower)
    prev = x_lower;
  else if (P >= cdf_x_upper)
    prev = x_upper;

  params.P = P;
  params.K = K;
  params.w = w;
  params.mu = mu;
  params.sigma = sigma;

  if (P <= cdf_x_lower || P >= cdf_x_upper) {
    fdf.f = &gsl_cdf_gaussian_mixture_P_f;
    fdf.df = &gsl_cdf_gaussian_mixture_P_df;
    fdf.fdf = &gsl_cdf_gaussian_mixture_P_fdf;
    fdf.params = &params;

    fdfsolver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);
    gsl_root_fdfsolver_set(fdfsolver, &fdf, prev);

    do
    {
      gsl_root_fdfsolver_iterate(fdfsolver);
      status = gsl_root_test_delta(gsl_root_fdfsolver_root(fdfsolver), prev, EPSABS, EPSREL);
      prev = gsl_root_fdfsolver_root(fdfsolver);
    } while (status == GSL_CONTINUE);
    gsl_root_fdfsolver_free(fdfsolver);

    return prev;
  }
  else {
    f.function = &gsl_cdf_gaussian_mixture_P_f;
    f.params = &params;

    fsolver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(fsolver, &f, x_lower, x_upper);

    do
    {
      gsl_root_fsolver_iterate(fsolver);
      x_lower = gsl_root_fsolver_x_lower(fsolver);
      x_upper = gsl_root_fsolver_x_upper(fsolver);
      status = gsl_root_test_interval(x_lower, x_upper, EPSABS, EPSREL);
    } while (status == GSL_CONTINUE);
    prev = gsl_root_fsolver_root(fsolver);
    gsl_root_fsolver_free(fsolver);

    return prev;
  }
}

double
gsl_cdf_gaussian_mixture_Qinv(const double Q, const size_t K, const double w[], const double mu[], const double sigma[])
{
  return gsl_cdf_gaussian_mixture_Pinv(1.0 - Q, K, w, mu, sigma);
}
