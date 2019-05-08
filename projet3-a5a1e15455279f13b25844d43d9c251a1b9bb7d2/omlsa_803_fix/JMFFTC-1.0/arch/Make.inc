############################ -*- Mode: Makefile -*- ###########################
## Make.inc --- Options de compilation
## 
## Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Créé le         : Mon Jun 19 14:47:22 2000
## Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Dern. mod. le   : Mon Jun 19 14:47:50 2000
###############################################################################

CC = cc
OPT = -h 2
PROF = -p
DEBUG = -g2
CFLAGS = $(OPT) -c
AR = ar
ARFLAGS = rv
RM = rm
RMFLAGS = -f
