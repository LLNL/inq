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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "liboct_parser.h"
#include "symbols.h"
#include "string_f.h"

#ifndef FC_FUNC
#error "Unknown Fortran name-mangling. Check configuration."
#endif

/* --------------------- Interface to the parsing routines ---------------------- */

/* --------------------------------------------------------- */
/* Initialization of the library                             */
/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_init, OCT_PARSE_INIT)
  (STR_F_TYPE s, int *dont_write STR_ARG1)
{
  int r;
  char *s_c;
  
  TO_C_STR1(s, s_c);
  r = parse_init(s_c, dont_write); 
  free(s_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_int, OCT_PARSE_PUTSYM_INT)
	(STR_F_TYPE s, int *i  STR_ARG1)
{
  char *s_c;

  TO_C_STR1(s, s_c);
  parse_putsym_int(s_c, *i);
  free(s_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_double, OCT_PARSE_PUTSYM_DOUBLE)
	(STR_F_TYPE s, double *d STR_ARG1)
{
  char *s_c;

  TO_C_STR1(s, s_c);
  parse_putsym_double(s_c, *d);
  free(s_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_complex, OCT_PARSE_PUTSYM_COMPLEX)
	(STR_F_TYPE s, gsl_complex *c STR_ARG1)
{
  char *s_c;

  TO_C_STR1(s, s_c);
  parse_putsym_complex(s_c, *c);
  free(s_c);
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_input, OCT_PARSE_INPUT)
       (STR_F_TYPE s, int *set_used STR_ARG1)
{
  int r;
  char *s_c;
  
  TO_C_STR1(s, s_c);
  r = parse_input(s_c, *set_used);
  free(s_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_environment, OCT_PARSE_ENVIRONMENT)
	(STR_F_TYPE s STR_ARG1)
{
  char *s_c;
  
  TO_C_STR1(s, s_c);
  parse_environment(s_c); 
  free(s_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_end, OCT_PARSE_END)
	()
{
  parse_end(); 
}

/* --------------------------------------------------------- */
void FC_FUNC_(oct_sym_output_table, OCT_SYM_OUTPUT_TABLE)
        (int *only_unused, int *mpiv_node)
{
  sym_output_table(*only_unused, *mpiv_node); 
}


/* --------------------------------------------------------- */
/* Parser functions                                          */
/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_isdef, OCT_PARSE_ISDEF)
	(STR_F_TYPE name STR_ARG1)
{ 
  int r;
  char *name_c;
  
  TO_C_STR1(name, name_c);
  r = parse_isdef(name_c); 
  free(name_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_int, OCT_PARSE_INT)
	(STR_F_TYPE name, int64_t *def, int64_t *res STR_ARG1)
{ 
  char *name_c;

  TO_C_STR1(name, name_c);
  *res = parse_int(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_double, OCT_PARSE_DOUBLE)
	(STR_F_TYPE name, double *def, double *res STR_ARG1)
{
  char *name_c;

  TO_C_STR1(name, name_c);
  *res = parse_double(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_complex, OCT_PARSE_COMPLEX)
	(STR_F_TYPE name, gsl_complex *def, gsl_complex *res STR_ARG1)
{
  char *name_c;

  TO_C_STR1(name, name_c);
  *res = parse_complex(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_string, OCT_PARSE_STRING)
	(STR_F_TYPE name, STR_F_TYPE def, STR_F_TYPE res STR_ARG3)
{
  char *c, *name_c, *def_c;
  
  TO_C_STR1(name, name_c);
  TO_C_STR2(def, def_c);
  c = parse_string(name_c, def_c); 
  TO_F_STR3(c, res);             /* convert string to Fortran */
  free(name_c); free(def_c);     /* this has to be *after* the to_f_str or we will have memory problems */
  free(c);
}


/* --------------------------------------------------------- */
static void parse_block_error(const char *type, const char *name, int l, int c){
  fprintf(stderr, "Parser error: block \"%s\" does not contain a %s in line %d and col %d.\n",
					name, type, l, c);
  exit(1);
}

/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block, OCT_PARSE_BLOCK)
	(STR_F_TYPE name, sym_block **blk STR_ARG1)
{
  int r;
  char *block_name;

  TO_C_STR1(name, block_name);
  r = parse_block(block_name, blk);
  free(block_name);
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_end, OCT_PARSE_BLOCK_END)
	(sym_block **blk)
{
  parse_block_end(blk);
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block_n, OCT_PARSE_BLOCK_N)
	(sym_block **blk)
{
  return parse_block_n(*blk);
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block_cols, OCT_PARSE_BLOCK_COLS)
	(sym_block **blk, int *l)
{
  return parse_block_cols(*blk, *l);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_int, OCT_PARSE_BLOCK_INT)
	(sym_block **blk, int *l, int *c, int *res)
{
  if(parse_block_int(*blk, *l, *c, res) != 0)
    parse_block_error("int", (*blk)->name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_double, OCT_PARSE_BLOCK_DOUBLE)
	(sym_block **blk, int *l, int *c, double *res)
{
  if(parse_block_double(*blk, *l, *c, res) != 0)
    parse_block_error("double", (*blk)->name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_complex, OCT_PARSE_BLOCK_COMPLEX)
	(sym_block **blk, int *l, int *c, gsl_complex *res)
{
  if(parse_block_complex(*blk, *l, *c, res) != 0)
    parse_block_error("complex", (*blk)->name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_string, OCT_PARSE_BLOCK_STRING)
	(sym_block **blk, int *l, int *c, STR_F_TYPE res STR_ARG1)
{
  char *s;
  
  if(parse_block_string(*blk, *l, *c, &s) != 0)
    parse_block_error("string", (*blk)->name, *l, *c);
  else{
    TO_F_STR1(s, res);
    free(s);
  }
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_expression, OCT_PARSE_EXPRESSION)
     (double *re, double *im, const int *ndim, const double *x, const double *r, const double *t, STR_F_TYPE pot STR_ARG1)
{
  symrec *rec;
  parse_result c;
  char *s_c;

  TO_C_STR1(pot, s_c);
  
  rec = putsym("x", S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, x[0], 0);
  rec->def = 1;

  if(*ndim > 1){
    rec = putsym("y", S_CMPLX);
    GSL_SET_COMPLEX(&rec->value.c, x[1], 0);
    rec->def = 1;
  }
  
  if(*ndim > 2){
    rec = putsym("z", S_CMPLX);
    GSL_SET_COMPLEX(&rec->value.c, x[2], 0);
    rec->def = 1;
  }

  if(*ndim > 3){
    rec = putsym("w", S_CMPLX);
    GSL_SET_COMPLEX(&rec->value.c, x[3], 0);
    rec->def = 1;
  }

  rec = putsym("r", S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, *r, 0);
  rec->def = 1;

  rec = putsym("t", S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, *t, 0);
  rec->def = 1;

  parse_exp(s_c, &c);

  /* clean up */
  rmsym("x");
  if(*ndim > 1) rmsym("y");
  if(*ndim > 2) rmsym("z");
  if(*ndim > 3) rmsym("w");
  rmsym("r");
  rmsym("t");

  *re = GSL_REAL(c.value.c);
  *im = GSL_IMAG(c.value.c);

  free(s_c);
}

/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_expression1, OCT_PARSE_EXPRESSION1)
  (double *re, double *im, STR_F_TYPE variable, double *val, STR_F_TYPE string STR_ARG2)
{

  symrec *rec;
  parse_result c;
  char *s_c, *var_c;

  TO_C_STR1(variable, var_c);
  TO_C_STR2(string, s_c);

  rec = putsym(var_c, S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, *val, 0);
  rec->def = 1;

  parse_exp(s_c, &c);

  rmsym(var_c);

  *re = GSL_REAL(c.value.c);
  *im = GSL_IMAG(c.value.c);

  free(s_c); free(var_c);
}
