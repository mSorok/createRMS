/**********************************************                                                    
  Data structures for precision

 Jean-loup Faulon May 91 Albuquerque N.M.
                  Sep 91 PSU PA.
**********************************************/

#include <math.h>

/* dimension */
#define  D1      1
#define  D2      2
#define  D3      3
#define  D4          4

/* Fortran */
#define  NE      (long)0
#define  GT      (long)1
#define  LT      (long)2
#define  EQ      (long)3

/* physical and math constants */
#define INFINI      1.e+16
#define ZERO        (double)0
#define UN          (double)1
#define ONE         (double)1
#define PI          (double)3.14159265359
#define AVOGADRO    (double)6.0221367e+23
#define ROUND(a)    ((a > 0) ? ((long)(a + 0.5)) & 0x00000000FFFFFFFF : -((long)(0.5 - a)) & 0x00000000FFFFFFFF)

/*  odd even */
#define EVEN(a) ((2 *((long)(a)/2) ==  (a)) ? VRAI : FAUX)
#define ODD(a)  ((2 *((long)(a)/2) !=  (a)) ? VRAI : FAUX)

/* precision for 3D construction */
#define EPS         1.e-3  
#define EPS_EQ(a,b) ((((b) - EPS <= (a)) && ((a) <= (b) + EPS)) ? VRAI : FAUX)
#define EPS_LT(a,b) (((a) <  (b) - EPS) ? VRAI : FAUX)
#define EPS_GT(a,b) (((a) >  (b) + EPS) ? VRAI : FAUX)

/* low precision */
#define EPS1        1.e-1  
#define EPS1_EQ(a,b) ((((b) - EPS1 <= (a)) && ((a) <= (b) + EPS1)) ? VRAI : FAUX)
#define EPS1_LT(a,b) (((a) <  (b) - EPS1) ? VRAI : FAUX)
#define EPS1_GT(a,b) (((a) >  (b) + EPS1) ? VRAI : FAUX)

/* high precision */
#define EPS5        1.e-5
#define EPS5_EQ(a,b) ((((b) - EPS5 <= (a)) && ((a) <= (b) + EPS5)) ? VRAI : FAUX)
#define EPS5_LT(a,b) (((a) <  (b) - EPS5) ? VRAI : FAUX)
#define EPS5_GT(a,b) (((a) >  (b) + EPS5) ? VRAI : FAUX)

#define EPS6        1.e-5
#define EPS6_EQ(a,b) ((((b) - EPS6 <= (a)) && ((a) <= (b) + EPS6)) ? VRAI : FAUX)
#define EPS6_LT(a,b) (((a) <  (b) - EPS6) ? VRAI : FAUX)
#define EPS6_GT(a,b) (((a) >  (b) + EPS6) ? VRAI : FAUX)

/* much more high precision */
#define EPS8        1.e-8
#define EPS8_EQ(a,b) ((((b) - EPS8 <= (a)) && ((a) <= (b) + EPS8)) ? VRAI : FAUX)
#define EPS8_LT(a,b) (((a) <  (b) - EPS8) ? VRAI : FAUX)
#define EPS8_GT(a,b) (((a) >  (b) + EPS8) ? VRAI : FAUX)

/* much much more high precision */
#define EPS12        1.e-12
#define EPS12_EQ(a,b) ((((b) - EPS12 <= (a)) && ((a) <= (b) + EPS12)) ? VRAI : FAUX)
#define EPS12_LT(a,b) (((a) <  (b) - EPS12) ? VRAI : FAUX)
#define EPS12_GT(a,b) (((a) >  (b) + EPS12) ? VRAI : FAUX)

/* much much more high precision (close to float precision) */
#define EPS30        1.e-30
#define EPS30_EQ(a,b) ((((b) - EPS30 <= (a)) && ((a) <= (b) + EPS30)) ? VRAI : FAUX)
#define EPS30_LT(a,b) (((a) <  (b) - EPS30) ? VRAI : FAUX)
#define EPS30_GT(a,b) (((a) >  (b) + EPS30) ? VRAI : FAUX)

/*  mathematical functions */
double sin();
double cos();
double atan();
double fabs();
double acos();
double sqrt();
double pow();



