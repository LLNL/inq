/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch, D. Strubbe

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>

#include "symbols.h"
#include "gsl_userdef.h"

/* The symbol table: a chain of `struct symrec'.  */
symrec *sym_table = (symrec *)0;

void str_tolower(char *in)
{
  for(; *in; in++)
    *in = (char)tolower(*in);
}

void sym_mark_table_used ()
{
  symrec *ptr;

  for (ptr = sym_table; ptr != (symrec *) 0;
       ptr = (symrec *)ptr->next)
  {
    ptr->used = 1;
  }
}

symrec *putsym (const char *sym_name, symrec_type sym_type)  
{
  symrec *ptr;
  ptr = (symrec *)malloc(sizeof(symrec));

  /* names are always lowercase */
  ptr->name = strdup(sym_name);
  str_tolower(ptr->name);
  
  ptr->def  = 0;
  ptr->used = 0;
  ptr->type = sym_type;
  GSL_SET_COMPLEX(&ptr->value.c, 0, 0); /* set value to 0 even if fctn.  */
  ptr->next = (struct symrec *)sym_table;
  sym_table = ptr;
  return ptr;
}

symrec *getsym (const char *sym_name)
{
  symrec *ptr;
  for (ptr = sym_table; ptr != (symrec *) 0;
       ptr = (symrec *)ptr->next)
    if (strcasecmp(ptr->name,sym_name) == 0){
      ptr->used = 1;
      return ptr;
    }
  return (symrec *) 0;
}

int rmsym (const char *sym_name)
{
  symrec *ptr, *prev;
  for (prev = (symrec *) 0, ptr = sym_table; ptr != (symrec *) 0;
       prev = ptr, ptr = ptr->next)
    if (strcasecmp(ptr->name,sym_name) == 0){
      if(prev == (symrec *) 0)
	sym_table = ptr->next;
      else
	prev->next = ptr->next;
      free(ptr->name);
      free(ptr);
      
      return 1;
    }
  
  return 0;
}

struct init_fntc{
  char *fname;
  int  nargs;
  gsl_complex (*fnctptr)();
};

void sym_notdef (symrec *sym)
{
  fprintf(stderr, "Parser error: symbol '%s' used before being defined.\n", sym->name);
  exit(1);
}

void sym_redef (symrec *sym)
{
  fprintf(stderr, "Parser warning: redefining symbol, previous value ");
  sym_print(stderr, sym);
  fprintf(stderr, "\n");
}

void sym_wrong_arg (symrec *sym)
{
  if(sym->type == S_BLOCK) {
    fprintf(stderr, "Parser error: block name '%s' used in variable context.\n", sym->name);
  } else if(sym->type == S_STR) {
    fprintf(stderr, "Parser error: string variable '%s' used in expression context.\n", sym->name);
  } else {
    fprintf(stderr, "Parser error: function '%s' requires %d argument(s).\n", sym->name, sym->nargs);
  }
  exit(1);
}

static struct init_fntc arith_fncts[] = {
  {"sqrt",   1, (gsl_complex (*)()) &gsl_complex_sqrt},
  {"exp",    1, (gsl_complex (*)()) &gsl_complex_exp},
  {"ln",     1, (gsl_complex (*)()) &gsl_complex_log},
  {"log",    1, (gsl_complex (*)()) &gsl_complex_log},
  {"log10",  1, (gsl_complex (*)()) &gsl_complex_log10},
  {"logb",   2, (gsl_complex (*)()) &gsl_complex_log_b}, /* takes two arguments logb(z, b) = log_b(z) */

  {"arg",    1, (gsl_complex (*)()) &gsl_complex_carg},
  {"abs",    1, (gsl_complex (*)()) &gsl_complex_cabs},
  {"abs2",   1, (gsl_complex (*)()) &gsl_complex_cabs2},
  {"logabs", 1, (gsl_complex (*)()) &gsl_complex_clogabs},

  {"conjg",  1, (gsl_complex (*)()) &gsl_complex_conjugate},
  {"inv",    1, (gsl_complex (*)()) &gsl_complex_inverse},

  {"sin",    1, (gsl_complex (*)()) &gsl_complex_sin},
  {"cos",    1, (gsl_complex (*)()) &gsl_complex_cos},
  {"tan",    1, (gsl_complex (*)()) &gsl_complex_tan},
  {"sec",    1, (gsl_complex (*)()) &gsl_complex_sec},
  {"csc",    1, (gsl_complex (*)()) &gsl_complex_csc},
  {"cot",    1, (gsl_complex (*)()) &gsl_complex_cot},

  {"asin",   1, (gsl_complex (*)()) &gsl_complex_arcsin},
  {"acos",   1, (gsl_complex (*)()) &gsl_complex_arccos},
  {"atan",   1, (gsl_complex (*)()) &gsl_complex_arctan},
  {"atan2",  2, (gsl_complex (*)()) &gsl_complex_arctan2}, /* takes two arguments atan2(y,x) = atan(y/x) */
  {"asec",   1, (gsl_complex (*)()) &gsl_complex_arcsec},
  {"acsc",   1, (gsl_complex (*)()) &gsl_complex_arccsc},
  {"acot",   1, (gsl_complex (*)()) &gsl_complex_arccot},

  {"sinh",   1, (gsl_complex (*)()) &gsl_complex_sinh},
  {"cosh",   1, (gsl_complex (*)()) &gsl_complex_cosh},
  {"tanh",   1, (gsl_complex (*)()) &gsl_complex_tanh},
  {"sech",   1, (gsl_complex (*)()) &gsl_complex_sech},
  {"csch",   1, (gsl_complex (*)()) &gsl_complex_csch},
  {"coth",   1, (gsl_complex (*)()) &gsl_complex_coth},

  {"asinh",  1, (gsl_complex (*)()) &gsl_complex_arcsinh},
  {"acosh",  1, (gsl_complex (*)()) &gsl_complex_arccosh},
  {"atanh",  1, (gsl_complex (*)()) &gsl_complex_arctanh},
  {"asech",  1, (gsl_complex (*)()) &gsl_complex_arcsech},
  {"acsch",  1, (gsl_complex (*)()) &gsl_complex_arccsch},
  {"acoth",  1, (gsl_complex (*)()) &gsl_complex_arccoth},	
 
/* user-defined step function. this is not available in GSL, 
   but we use GSL namespacing and macros here. */
  {"step",   1, (gsl_complex (*)()) &gsl_complex_step_real},

/* Minimum and maximum of two arguments (comparing real parts) */  
  {"min",    2, (gsl_complex (*)()) &gsl_complex_min_real},
  {"max",    2, (gsl_complex (*)()) &gsl_complex_max_real},

  {"erf",    1, (gsl_complex (*)()) &gsl_complex_erf},

  {"realpart", 1, (gsl_complex (*)()) &gsl_complex_realpart},
  {"imagpart", 1, (gsl_complex (*)()) &gsl_complex_imagpart},
  {"round",   1, (gsl_complex (*)()) &gsl_complex_round},
  {"floor",   1, (gsl_complex (*)()) &gsl_complex_floor},
  {"ceiling", 1, (gsl_complex (*)()) &gsl_complex_ceiling},

  {"rand",    0, (gsl_complex (*)()) &gsl_complex_rand},

  {0, 0, 0}
};

struct init_cnst{
	char *fname;
	double re;
	double im;
};

static struct init_cnst arith_cnts[] = {
	{"pi",    M_PI, 0}, 
	{"e",      M_E, 0},
	{"i",        0, 1},
	{"true",     1, 0}, 
	{"yes",      1, 0},
	{"false",    0, 0}, 
	{"no",       0, 0},
	{0,          0, 0}
};

char *reserved_symbols[] = {
  "x", "y", "z", "r", "w", "t", 0
};

void sym_init_table ()  /* puts arithmetic functions in table. */
{
  int i;
  symrec *ptr;
  for (i = 0; arith_fncts[i].fname != 0; i++){
    ptr = putsym (arith_fncts[i].fname, S_FNCT);
    ptr->def = 1;
    ptr->used = 1;
    ptr->nargs = arith_fncts[i].nargs;
    ptr->value.fnctptr = arith_fncts[i].fnctptr;
  }

  /* now the constants */
  for (i = 0; arith_cnts[i].fname != 0; i++){
    ptr = putsym(arith_cnts[i].fname, S_CMPLX);
    ptr->def = 1;
    ptr->used = 1;
    GSL_SET_COMPLEX(&ptr->value.c, arith_cnts[i].re, arith_cnts[i].im);
  }
}

void sym_end_table()
{
  symrec *ptr, *ptr2;
  int l, col;

  for (ptr = sym_table; ptr != NULL;){
    free(ptr->name);
    switch(ptr->type){
    case S_STR:
      free(ptr->value.str);
      break;
    case S_BLOCK:
      if(ptr->value.block->n > 0){
	free(ptr->value.block->lines);
      }
      free(ptr->value.block);
      break;
    case S_CMPLX:
    case S_FNCT:
      break;
    }
    ptr2 = ptr->next;
    free(ptr);
    ptr = ptr2;
  }
  
  sym_table = NULL;
}

/* this function is defined in src/basic/varinfo_low.c */
int varinfo_variable_exists(const char * var_name);

void sym_output_table(int only_unused, int mpiv_node)
{
  FILE *f;
  symrec *ptr;
  int any_unused = 0;

  if(mpiv_node != 0) {
    return;
  }
  
  if(only_unused) {
    f = stderr;
  } else {
    f = stdout;
  }
  
  for(ptr = sym_table; ptr != NULL; ptr = ptr->next){
    if(only_unused && ptr->used == 1) continue;
    if(only_unused && varinfo_variable_exists(ptr->name)) continue;
    if(any_unused == 0) {
      fprintf(f, "\nParser warning: possible mistakes in input file.\n");
      fprintf(f, "List of variable assignments not used by parser:\n");
      any_unused = 1;
    }

    sym_print(f, ptr);
  }
  if(any_unused == 1) {
    fprintf(f, "\n");
  }
}

void sym_print(FILE *f, const symrec *ptr)
{
  fprintf(f, "%s", ptr->name);
  switch(ptr->type){
  case S_CMPLX:
    if(fabs(GSL_IMAG(ptr->value.c)) < 1.0e-14){
      fprintf(f, " = %f\n", GSL_REAL(ptr->value.c));
    } else {
      fprintf(f, " = (%f,%f)\n", GSL_REAL(ptr->value.c), GSL_IMAG(ptr->value.c));
    }
    break;
  case S_STR:
    fprintf(f, " = \"%s\"\n", ptr->value.str);
    break;
  case S_BLOCK:
    fprintf(f, "%s\n", " <= BLOCK");
    break;
  case S_FNCT:
    fprintf(f, "%s\n", " <= FUNCTION");
    break;
  }
}
