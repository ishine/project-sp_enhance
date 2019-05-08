############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Construction de la bibliotheque JMFFT
##
## Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Créé le         : Tue Dec  8 18:18:21 1998
## Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Dern. mod. le   : Mon May 24 11:29:11 2004
###############################################################################

SHELL = /bin/ksh

all     : jmfft
clean	:
	@echo
	@echo "Destruction de la bibliotheque libjmfft.a ..."
	rm -f libjmfft.a
	(cd src ; make clean)
jmfft     :
	@echo
	@echo "Construction de la bibliotheque libjmfft.a ..."
	(cd src ; make && mv libjmfft*.a ../lib/.)
