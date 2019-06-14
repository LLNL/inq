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

/* This is essentially the example from bison */
/* https://www.gnu.org/software/bison/manual/html_node/Infix-Calc.html */
%{
%}

%union {
  gsl_complex val;  /* For returning numbers.               */
  char *str;        /* For strings                          */
  symrec  *tptr;    /* For returning symbol-table pointers  */
}

%token <val>  NUM        /* Simple complex number   */
%token <str>  STR        /* For strings             */
%token <tptr> VAR FNCT   /* Variable and Function   */
%type  <val>  exp
%type  <str>  string

%left ','
%left ')'
%right '='
%left '<' '>' LE GE EQUAL
%left '-' '+'
%left LAND LOR
%left '*' '/'
%left NEG NOT  /* Negation--unary minus, and logical operator */
%right '^'     /* Exponentiation        */

/* Grammar follows */

%%

input:   /* empty */
| input line
;

line:
'\n'
| exp '\n'   { par_res.value.c = $1; par_res.type = PR_CMPLX; YYACCEPT;}
| string '\n'{ par_res.value.s = $1; par_res.type = PR_STR;   YYACCEPT;}
| error '\n' { yyerrok; YYABORT;}
;
     
exp: NUM                   { $$ = $1;                           }
| exp ',' exp              { fprintf(stderr, "Parser error: comma is not valid operator\n"); exit(1); }
| VAR                      { if(!$1->def) sym_notdef($1); $$ = $1->value.c; }
| VAR '=' exp              { if($1->def && (gsl_complex_abs(gsl_complex_sub($1->value.c, $3)) > 1e-9)) sym_redef($1);
                             $$ = $3; $1->value.c = $3; $1->def = 1; $1->type = S_CMPLX;}
| FNCT '(' ')'             { $$ = (*($1->value.fnctptr))();   }
| FNCT '(' exp ')'         { if($1->nargs != 1) sym_wrong_arg($1); $$ = (*($1->value.fnctptr))($3);   }
| FNCT '(' exp ',' exp ')' { if($1->nargs != 2) sym_wrong_arg($1); $$ = (*($1->value.fnctptr))($3, $5); }
| exp '+' exp              { $$ = gsl_complex_add($1, $3);      }
| exp '-' exp              { $$ = gsl_complex_sub($1, $3);      }
| exp '*' exp              { $$ = gsl_complex_mul($1, $3);      }
| exp '/' exp              { $$ = gsl_complex_div($1, $3);      }
| '-' exp  %prec NEG       { $$ = gsl_complex_negative($2);     }
| exp '^' exp              { $$ = gsl_complex_pow($1, $3);      }
| exp '<' exp              { GSL_SET_COMPLEX (&$$, GSL_REAL($1) <  GSL_REAL($3), 0); } /* Boolean comparisons use only the real part */
| exp '>' exp              { GSL_SET_COMPLEX (&$$, GSL_REAL($1) >  GSL_REAL($3), 0); } /* with the exception of '==' */
| exp LE  exp              { GSL_SET_COMPLEX (&$$, GSL_REAL($1) <= GSL_REAL($3), 0); }
| exp GE  exp              { GSL_SET_COMPLEX (&$$, GSL_REAL($1) >= GSL_REAL($3), 0); }
| exp EQUAL exp            { GSL_SET_COMPLEX (&$$, (GSL_REAL($1) == GSL_REAL($3)) && (GSL_IMAG($1) == GSL_IMAG($3)), 0); }
| exp LAND  exp            { GSL_SET_COMPLEX (&$$, GSL_REAL($1) && GSL_REAL($3), 0); }
| exp LOR   exp            { GSL_SET_COMPLEX (&$$, GSL_REAL($1) || GSL_REAL($3), 0); }
| '!' exp  %prec NOT       { GSL_SET_COMPLEX (&$$, !GSL_REAL($2), 0); }
| '{' exp ',' exp '}'      { GSL_SET_COMPLEX (&$$, GSL_REAL($2), GSL_REAL($4)); }
| '(' exp ')'              { $$ = $2;                           }

string: STR                { $$ = $1; }
| VAR '=' STR              { if($1->def) sym_redef($1); $$ = $3; $1->value.str = $3; $1->def = 1; $1->type = S_STR; }
;
%%
