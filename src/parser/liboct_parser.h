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

#ifndef _LIB_OCT_H
#define _LIB_OCT_H

#include <gsl/gsl_complex.h>
#include <stdint.h>

int         parse_init   (const char *file_out, const int *mpiv_node);
int         parse_input  (const char *file_in, int set_used);
void        parse_environment  (const char *file_in);
void        parse_end    (void);

int         parse_isdef  (const char *name);
int64_t     parse_int    (const char *name, int64_t def);
double      parse_double (const char *name, double def);
gsl_complex parse_complex(const char *name, gsl_complex def);
char       *parse_string (const char *name, char *def);

/* Now comes stuff for the blocks */
typedef struct sym_block_line{
  int n;
} sym_block_line;

typedef struct sym_block{
  int n;
  sym_block_line *lines;
  char* name;
} sym_block;

int parse_block        (const char *name, sym_block **blk);
int parse_block_end    (sym_block **blk);
int parse_block_n      (const sym_block *blk);
int parse_block_cols   (const sym_block *blk, int l);
int parse_block_int    (const sym_block *blk, int l, int col, int *r);
int parse_block_double (const sym_block *blk, int l, int col, double *r);
int parse_block_complex(const sym_block *blk, int l, int col, gsl_complex *r);
int parse_block_string (const sym_block *blk, int l, int col, char **r);

/* from parse_exp.c */
enum pr_type {PR_NONE,PR_CMPLX, PR_STR};
typedef struct parse_result{
  union {
    gsl_complex c;
    char *s;
  } value;
  enum pr_type type;
} parse_result;

void parse_result_free(parse_result *t);

int parse_exp(char *exp, parse_result *t);

void parse_putsym_int(const char *s, int i);
void parse_putsym_double(const char *s, double d);
void parse_putsym_complex(const char *s, gsl_complex c);

#endif
