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

#ifndef _SYMBOLS_H
#define _SYMBOLS_H

#include <gsl/gsl_complex.h>
#include "liboct_parser.h"

typedef enum{
  S_CMPLX, S_STR, S_BLOCK, S_FNCT
} symrec_type;

/* Data type for links in the chain of symbols. */
typedef struct symrec{
  char *name;                  /* name of symbol */
  symrec_type type;            /* type of symbol: complex, string, block, or function */
  int def;                     /* has this symbol been defined */
  int used;                    /* this symbol has been used before */

  int nargs;                   /* if type==S_FNCT contains the number of arguments of the function */

  union {
    gsl_complex c;             /* value of a S_CMPLX */
    char *str;                 /* value of a S_STR */
    sym_block *block;          /* to store blocks */
    gsl_complex (*fnctptr)();  /* value of a S_FNCT */
  } value;

  struct symrec *next;         /* link field */
} symrec;

/* The symbol table: a chain of struct symrec. */
extern symrec *sym_table;
extern char *reserved_symbols[];

symrec *putsym (const char *sym_name, symrec_type sym_type);
symrec *getsym (const char *sym_name);
int      rmsym (const char *sym_name);

void sym_notdef(symrec *sym);
void sym_redef(symrec *sym);
void sym_wrong_arg(symrec *sym);
void sym_init_table(void);
void sym_end_table(void);
void sym_output_table(int only_unused, int mpiv_node);
void str_tolower(char *in);
void sym_mark_table_used();
void sym_print(FILE *f, const symrec *ptr);

#endif
