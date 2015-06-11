/***********************************************************
  File for random request.
  Jean-loup Faulon May 91.
  WARNING : The rand() function is different on HP-Apollo
  Sun and Silicon.
  The function used is now drand48() which return a double
  number between 0 and 1.
  Modified Pennstate July 91.
***********************************************************/
#include <general.h>
#include <eps.h>



DEFINE 	t_real   	dara_random_real(n) 
/***********************************************************
 Return random value beetween 0 and n 
 **********************************************************/
t_real	     n;
{
EXTERN	double	drand48();
       	t_real  d;

        d = drand48(); d = d * n;
       	return(d);
       	}

DEFINE 	t_real   	dara_standard_gaussian_distribution() 
/***********************************************************
 Return random value following standard gaussian distribution
 (mean = 0, sdtdev = 1)
 **********************************************************/
{
LOCAL	t_integer	iset = 0;
LOCAL	t_real		gset;
	t_real 		fac,r,v1,v2;


	if (iset == 0) {
	  do {
	    v1 = 2*dara_random_real((t_real)1)-1;
 	    v2 = 2*dara_random_real((t_real)1)-1;
	    r = v1*v1+v2*v2;
	    } while (r >= 1.0);
	  fac= sqrt(-2.0*log(r)/r);
	  gset = v1*fac;
	  iset = 1;
	  return(v2*fac);
	  }
        else { iset = 0; return(gset); }
	}



DEFINE 	t_real   	dara_gaussian_distribution(m,s) 
/***********************************************************
 Return random value following gaussian distribution
 (mean = m, sdtdev = s) 
 **********************************************************/
t_real	     m,s;
{
       	t_real  	z;

	z = dara_standard_gaussian_distribution();
	return(s*z+m);
       	}

DEFINE 	t_err   	dara_random_integer(n) 
/***********************************************************
 Return random value beetween 0 and n-1 
 Pennstate July 91.
 **********************************************************/
t_integer     n;
{
EXTERN	double	drand48();
       	t_real  d;

       	d = n-1;
	d = dara_random_real(d);
	if (d < 0) d = 0;
       	return(ROUND(d));
       	}


DEFINE 	t_err   	dara_random_atom(a) 
/***********************************************************
 Give a random ID for atom a
 **********************************************************/
atom_p	a;
{
	t_err	i;

	if (a == NIL) return(ERROR);
	i = dara_random_integer(((molecule_p)((group_p)a->group)->molecule)->size);
	i = datw_set_twin(i,a);
	a->ID = i;
	return(OK);
	}
 
DEFINE 	t_err   	dara_random_molecule(m) 
/***********************************************************
 Give a random ID for all atom of molecule m.
 **********************************************************/
molecule_p	m;
{
EXTERN	t_err		dast_set_ID_atom();
	t_integer	ID = 0;

	if (m == NIL) return(ERROR);
	if (m->size < 1) return(ERROR);
	dast_all_atom(m,NIL,dast_set_ID_atom,FALSE,&ID); m->size = ID;
	datw_init_twin(m->size);
	dast_all_atom(m,NIL,dara_random_atom,FALSE);
	return(OK);
	}  

