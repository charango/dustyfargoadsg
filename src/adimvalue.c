#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double Msol    = 2.0e30;     // masse du soleil
const double AU      = 1.5e11;     // dist Terre-Soleil en metres
const double annee   = 3.1558e7;   // nb de secondes dans un an
const double tempref = 2.5797e5;   // T tq R/mu = 1

/* Les dimensions sont notees 
   M (masse), 
   L (longueur), 
   T (temps) et 
   U (temperature) */

int main() {
  double nl, nm, nt, nu;
  char charquoi;
  double dimvalue, adimvalue;
  double computeadimvalue ();
  double computedimvalue ();

  printf ("Voulez-vous calculer une grandeur adimensionnee ou dimensionnee? Taper A ou D...\n");
  scanf ("%s", &charquoi);
  if ( (charquoi != 'A') && (charquoi != 'D') ){
    printf ("Abruti, tu dois taper un A ou un D!! Try again :) \n");
    return 0;
  }
  if ( charquoi == 'A') {
    printf ("Taper la valeur de la grandeur que vous souhaiter calculer adimensionnée\n");
    scanf ("%lg", &dimvalue);
    printf ("Il me faut connaitre la dimension de la grandeur a adimensionner.\nNotons-la pow(M,nm) x pow(L,nl) x pow(T,nt) x pow(U,nu),\nou M=masse, L=longueur, T=temps et U=temperature\n");
    printf ("valeur de nm?\n");
    scanf ("%lg", &nm);
    printf ("valeur de nl?\n");
    scanf ("%lg", &nl);
    printf ("valeur de nt?\n");
    scanf ("%lg", &nt);
    printf ("valeur de nu?\n");
    scanf ("%lg", &nu);
    adimvalue = computeadimvalue (dimvalue, nm, nl, nt, nu);
    printf ("La valeur adimensionnee recherchee est:\n");
    printf ("adimvalue = %lg\n",adimvalue);
  }
  if ( charquoi == 'D') {
    printf ("Taper la valeur de la grandeur que vous souhaiter calculer dimensionnée\n");
    scanf ("%lg", &adimvalue);
    printf ("Il me faut connaitre la dimension de la grandeur a dimensionner.\nNotons-la pow(M,nm) x pow(L,nl) x pow(T,nt) x pow(U,nu),\nou M=masse, L=longueur, T=temps et U=temperature\n");
    printf ("valeur de nm?\n");
    scanf ("%lg", &nm);
    printf ("valeur de nl?\n");
    scanf ("%lg", &nl);
    printf ("valeur de nt?\n");
    scanf ("%lg", &nt);
    printf ("valeur de nu?\n");
    scanf ("%lg", &nu);
    dimvalue = computedimvalue (adimvalue, nm, nl, nt, nu);
    printf ("La valeur dimensionnee recherchee est:\n");
    printf ("dimvalue = %.2le en kg^%.1g m^%.1g s^%.1g kelvin^%.1g\n",dimvalue,nm,nl,nt,nu);
  }
  return 0;
}


double computeadimvalue (dimvalue, nm, nl, nt, nu)
     double dimvalue;
     double nm, nl, nt, nu;
{
  double den, adimvalue;
  den = pow(Msol,nm) * pow(AU,nl) * pow(annee/2.0/M_PI,nt) * pow(tempref,nu);
  adimvalue = dimvalue / den;
  return adimvalue;
}

double computedimvalue (adimvalue, nm, nl, nt, nu)
     double adimvalue;
     double nm, nl, nt, nu;
{
  double num, dimvalue;
  num = pow(Msol,nm) * pow(AU,nl) * pow(annee/2.0/M_PI,nt) * pow(tempref,nu);
  dimvalue = adimvalue * num;
  return dimvalue;
}
