#!/bin/bash

mpicc -lm main.c MC.c Single.c MPI.c SampleMaker.c Visual.c Utility.c Calc.c Database.c -o SHIMS.out

