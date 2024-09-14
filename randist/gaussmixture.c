/* randist/gaussmixture.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006, 2007 James Theiler, Brian Gough
 * Copyright (C) 2006 Charles Karney
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double
gsl_ran_gaussian_mixture (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])
{
  gsl_ran_discrete_t *g = NULL;
  size_t v;
  
  g = gsl_ran_discrete_preproc (K, w);
  v = gsl_ran_discrete (r, g);
  gsl_ran_discrete_free (g);
  
  return mu[v] + gsl_ran_gaussian (r, sigma[v]);
}

double
gsl_ran_gaussian_mixture_ratio_method (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])
{
  gsl_ran_discrete_t *g = NULL;
  size_t v;
  
  g = gsl_ran_discrete_preproc (K, w);
  v = gsl_ran_discrete (r, g);
  gsl_ran_discrete_free (g);
  
  return mu[v] + gsl_ran_gaussian_ratio_method (r, sigma[v]);
}

double
gsl_ran_gaussian_mixture_ziggurat (const gsl_rng * r, const size_t K, const double w[], const double mu[], const double sigma[])
{
  gsl_ran_discrete_t *g = NULL;
  size_t v;
  
  g = gsl_ran_discrete_preproc (K, w);
  v = gsl_ran_discrete (r, g);
  gsl_ran_discrete_free (g);
  
  return mu[v] + gsl_ran_gaussian_ziggurat (r, sigma[v]);
}

double
gsl_ran_gaussian_mixture_pdf (const double x, const size_t K, const double w[], const double mu[], const double sigma[])
{
  double p = 0.0;
  double u;
  int i;

  for (i = 0; i < K; ++i) {
    u = (x - mu[i]) / sigma[i];
    p += w[i] / sigma[i] * exp(-u * u / 2.0);
  }
  p /= sqrt(2.0 * M_PI);

  return p;
}

