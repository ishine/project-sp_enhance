/* $Header: /home/teuler/cvsroot/lib/jmseterreur.f90,v 6.3 2000/02/18 17:42:04 teuler Exp $ */
/* JMFFTLIB : A library of portable fourier transform subroutines */
/* emulating Cray SciLib */
/* Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr) */
/*  */
/* Permission is granted to copy and distribute this file or modified */
/* versions of this file for no fee, provided the copyright notice and */
/* this permission notice are preserved on all copies. */

#include "jmfft.h"
#include "jmfft2.h"

void jmseterreur(arret)


  /* Arguments */
  INTEGER arret;

{

 /* Variables locales */

  INTEGER arret2;

  arret2 = arret;
  jmgetseterreur(&arret2,"s");

}
