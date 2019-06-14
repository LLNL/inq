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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "liboct_parser.h"
#include "symbols.h"

static char *par_string;
static int par_pos;
parse_result par_res;

static int oct_parser_lex(void);
static int oct_parser_error (char *s)  /* Called by oct_parser_parse on error */
{
  /* Do nothing */
  /* printf("%s\n", s); */
  return 0;
}

/* include the parser */
#include "grammar.c"

int parse_exp(char *exp, parse_result *r)
{
  int o;

  par_string = exp;
  par_pos = 0;

  o = oct_parser_parse();
  if(o == 0){
    r->type = par_res.type;
    if(r->type == PR_CMPLX)
      r->value.c = par_res.value.c;
    else
      r->value.s = par_res.value.s;
  }
  return o;
}

int get_real(char *s, double *d)
{
  int n=0;
  sscanf(s, "%lg", d);
  while(*s && (isdigit(*s) || *s=='.' || *s=='e' || *s=='E')){
    if((*s=='e' || *s=='E') && (*(s+1)=='+' || *(s+1)=='-')) {s++; n++;}
    s++; n++;
  }
  return n;
}

static int oct_parser_lex (){
  int c;
  char *symbuf = 0, *symbuf2 = 0;
  int length = 0;
  
  /* Ignore whitespace, get first nonwhite character.  */
  while ((c = par_string[par_pos++]) == ' ' || c == '\t');
	
  if (c == '\0')
    return '\n';
	
  /* Char starts a number => parse the number.         */
  if (c == '.' || isdigit (c)){
    par_pos--;
    par_pos += get_real(&par_string[par_pos], &GSL_REAL(yylval.val));
    return NUM;
  }

  /* get the logical operators */
  if (c == '<' && par_string[par_pos] == '='){
    par_pos++;
    return LE;
  }

  if (c == '>' && par_string[par_pos] == '='){
    par_pos++;
    return GE;
  }

  if (c == '=' && par_string[par_pos] == '='){
    par_pos++;
    return EQUAL;
  }

  /*
  if (c == ':' && par_string[par_pos] == '='){
    par_pos++;
    return SET;
  }
  */

  if (c == '&' && par_string[par_pos] == '&'){
    par_pos++;
    return LAND;
  }

  if (c == '|' && par_string[par_pos] == '|'){
    par_pos++;
    return LOR;
  }
     
  /* Char starts an identifier => read the name.       */
  if (isalpha (c) || c == '\'' || c == '\"'){
    symrec *s;
    char startc = (char)c;
    int i;
		
    /* Initially make the buffer long enough
       for a 40-character symbol name.  */
    if (length == 0)
      length = 40, symbuf = (char *)malloc (length + 1);
    
    if(startc == '\'' || startc == '\"')
      c = par_string[par_pos++];
    else
      startc = 0; /* false */
    
    i = 0;
    do{
      /* If buffer is full, make it larger.        */
      if (i == length){
	length *= 2;
	symbuf2 = (char *)realloc (symbuf, length + 1);
	if (symbuf2 == NULL) {
	  fprintf(stderr, "Parser error: failed to reallocate buffer to size %i.\n", length);
	  exit(1);
	} else {
	  symbuf = symbuf2;
	}
      }
      /* Add this character to the buffer.         */
      symbuf[i++] = (char)c;
      /* Get another character.                    */
      c = par_string[par_pos++];
    }while (c != '\0' && 
	    ((startc && c!=startc) || (!startc && (isalnum(c) || c == '_' || c == '[' || c == ']'))));
    
    if(!startc) par_pos--;
    symbuf[i] = '\0';

    if(!startc){
      s = getsym (symbuf);
      if (s == 0){
	int jj;
	for (jj = 0; reserved_symbols[jj] != 0; jj++){
	  if(strcasecmp(symbuf, reserved_symbols[jj]) == 0){
	    fprintf(stderr, "Parser error: trying to redefine reserved symbol '%s'.\n", symbuf);
	    exit(1);
	  }
	}
	
	s = putsym (symbuf, S_CMPLX);
      }
      yylval.tptr = s;

      free(symbuf);
      if(s->type == S_CMPLX)
	return VAR;
      else
	return FNCT;
    }else{
      yylval.str = strdup(symbuf);

      free(symbuf);
      return STR;
    }
  }
	
  /* Any other character is a token by itself.        */
  return c;
}
