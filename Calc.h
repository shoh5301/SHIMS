#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Utility.h"
#define MAX 500

#define _USE_MATH_DEFINE

typedef enum {lparen, rparen, plus, minus, times, divide, mod, square, eos, Sin, Cos, Tan, ln, space, operand} precedence;

int inner_product(int v1[],int v2[]);
double eval(char* post);
void change_postfix(char*, char*);
precedence sin_cos_tan(char* ptr);
precedence is_ln(char*);
precedence get_token(char*);
char print_token(precedence);
void init_stack(void);
void push(double);
double pop(void);
void Etodec(char* eq);
void time_ins(char* Gline,double const time);

