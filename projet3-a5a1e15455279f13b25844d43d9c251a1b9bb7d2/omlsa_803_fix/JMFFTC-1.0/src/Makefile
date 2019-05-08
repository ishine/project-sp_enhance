############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Construction de libjmfft.a
## 
## Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Créé le         : Mon Jun 19 14:42:55 2000
## Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Dern. mod. le   : Mon May 24 11:28:43 2004
###############################################################################

include ../arch/Make.inc

LIBNAME = jmfft
LIB	= lib$(LIBNAME).a

OBJS = \
  $(LIB)(ipow.o) \
  $(LIB)(jmfftf2c.o) \
  $(LIB)(jmccfft.o) \
  $(LIB)(jmccfft2d.o) \
  $(LIB)(jmccfft3d.o) \
  $(LIB)(jmccfftm.o) \
  $(LIB)(jmccm1d.o) \
  $(LIB)(jmccm1d2.o) \
  $(LIB)(jmccm1d3.o) \
  $(LIB)(jmccm1d4.o) \
  $(LIB)(jmccm1d5.o) \
  $(LIB)(jmccm1dp.o) \
  $(LIB)(jmcctranspcs.o) \
  $(LIB)(jmcfftmlt.o) \
  $(LIB)(jmcftfax.o) \
  $(LIB)(jmcsfft.o) \
  $(LIB)(jmcsfft2d.o) \
  $(LIB)(jmcsfft3d.o) \
  $(LIB)(jmcsfftm.o) \
  $(LIB)(jmcsm1d.o) \
  $(LIB)(jmcsm1dxy.o) \
  $(LIB)(jmdecoup.o) \
  $(LIB)(jmdecoup3.o) \
  $(LIB)(jmerreur1.o) \
  $(LIB)(jmerreur2.o) \
  $(LIB)(jmerreur3.o) \
  $(LIB)(jmerreur4.o) \
  $(LIB)(jmfact.o) \
  $(LIB)(jmfftfax.o) \
  $(LIB)(jmgetcode.o) \
  $(LIB)(jmgeterreur.o) \
  $(LIB)(jmgetstop.o) \
  $(LIB)(jmgetmessage.o) \
  $(LIB)(jmgetnwork.o) \
  $(LIB)(jmgetsetcode.o) \
  $(LIB)(jmgetseterreur.o) \
  $(LIB)(jmgetsetstop.o) \
  $(LIB)(jmgetsetnwork.o) \
  $(LIB)(jmrfftmlt.o) \
  $(LIB)(jmscfft.o) \
  $(LIB)(jmscfft2d.o) \
  $(LIB)(jmscfft3d.o) \
  $(LIB)(jmscfftm.o) \
  $(LIB)(jmscm1d.o) \
  $(LIB)(jmscm1dxy.o) \
  $(LIB)(jmsetcode.o) \
  $(LIB)(jmseterreur.o) \
  $(LIB)(jmsetstop.o) \
  $(LIB)(jmsetnwork.o) \
  $(LIB)(jmtable.o) \
  $(LIB)(jmtransp.o)

all: $(OBJS)

.c.a:
	$(CC) $(CFLAGS) $(*).c
	$(AR) $(ARFLAGS) $(LIB) $(*).o
	$(RM) $(RMFLAGS) $(*).o

clean:
	$(RM) $(RMFLAGS) *.o *.a
