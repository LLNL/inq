/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include <stdlib.h>


/* ------------------------------------------------------ */
gsl_complex gsl_complex_step_real (gsl_complex a)
{        
  gsl_complex z;
	
  if (GSL_REAL(a) < 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 1, 0);
    }
  
  return z;
}


/* ------------------------------------------------------ */
gsl_complex gsl_complex_min_real (gsl_complex a, gsl_complex b)
{
  gsl_complex z;
  double min;
	
  /* just consider real parts */
  min = GSL_REAL(a) < GSL_REAL(b) ? GSL_REAL(a) : GSL_REAL(b);
  GSL_SET_COMPLEX (&z, min, 0);
	
  return z;
}


/* ------------------------------------------------------ */
gsl_complex gsl_complex_max_real (gsl_complex a, gsl_complex b)
{
  gsl_complex z;
  double max;
	
  /* just consider real parts */
  max = GSL_REAL(a) > GSL_REAL(b) ? GSL_REAL(a) : GSL_REAL(b);
  GSL_SET_COMPLEX (&z, max, 0);
	
  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_carg (gsl_complex a)
{        
  gsl_complex z;

  GSL_SET_COMPLEX (&z, gsl_complex_arg(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_cabs (gsl_complex a)
{        
  gsl_complex z;

  GSL_SET_COMPLEX (&z, gsl_complex_abs(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_cabs2 (gsl_complex a)
{        
  gsl_complex z;

  GSL_SET_COMPLEX (&z, gsl_complex_abs2(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_clogabs (gsl_complex a)
{        
  gsl_complex z;

  GSL_SET_COMPLEX (&z, gsl_complex_logabs(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_erf (gsl_complex a)
{        
  gsl_complex z;

  GSL_SET_COMPLEX (&z, gsl_sf_erf(GSL_REAL(a)), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_arctan2 (gsl_complex a, gsl_complex b)
{        
  gsl_complex z, p;

  if(GSL_REAL(b) != 0.0)
    {
      z = gsl_complex_arctan(gsl_complex_div(a, b));
      if(GSL_REAL(b) < 0.0){
	GSL_SET_COMPLEX (&p, M_PI, 0);
	if(GSL_REAL(a) >= 0.0)
	  z = gsl_complex_add(z, p);
	else
	  z = gsl_complex_sub(z, p);
      }
    }
  else
    {
      if(GSL_REAL(a) >= 0.0)
	{
	  GSL_SET_COMPLEX (&z, M_PI/2.0, 0.0);
	}
      else
	{
	  GSL_SET_COMPLEX (&z, -M_PI/2.0, 0.0);
	}
    }

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_realpart (gsl_complex a)
{
  gsl_complex z;

  GSL_SET_COMPLEX (&z, GSL_REAL(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_imagpart (gsl_complex a)
{
  gsl_complex z;

  GSL_SET_COMPLEX (&z, GSL_IMAG(a), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_round (gsl_complex a)
{
  gsl_complex z;

  GSL_SET_COMPLEX (&z, trunc(GSL_REAL(a)), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_ceiling (gsl_complex a)
{
  gsl_complex z;

  GSL_SET_COMPLEX (&z, ceil(GSL_REAL(a)), 0);

  return z;
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_floor (gsl_complex a)
{
  gsl_complex z;

  GSL_SET_COMPLEX (&z, floor(GSL_REAL(a)), 0);

  return z;
}

/* ------------------------------------------------------ */
void gsl_complex_rand_seed(long a)
{
  srandom(a);
}

/* ------------------------------------------------------ */
gsl_complex gsl_complex_rand (gsl_complex a)
{

  double r = random()/((double) RAND_MAX);
  return gsl_complex_rect(r, 0.0);
}

