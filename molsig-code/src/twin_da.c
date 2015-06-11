/***********************************************************
  File of duplication data.
  Jean-loup Faulon May 91 Albuquerque
***********************************************************/
#include <general.h>

LOCAL	t_pointer	TWIN[MAXATOM];
LOCAL	t_integer	MAXTWIN;

DEFINE	t_pointer	datw_inq_twin(i0,f,arg)
/**********************************************************
 Inquire the nearest twin from i0 who respect the function f
***********************************************************/
t_integer	i0;
f_traitement	f;
t_pointer	arg;
{
	t_integer	i;

	if (f == NIL) return(TWIN[i0]); i = i0;
	LOOP {
	  if (f(TWIN[i],arg) != ERROR) return(TWIN[i]);
	  i = ((i==MAXTWIN-1) ? 0 : ++i); if (i == i0) break;
	  }
	return((t_pointer)ERROR);
	}

DEFINE	t_err	datw_set_twin(i0,a)
/**********************************************************
  Set the value of the nearest TWIN from i0.
***********************************************************/
t_integer	i0;
t_pointer	a;
{
	t_integer	i;
	
	if (TWIN[i0] == NIL) {
	    TWIN[i0] = a; return(i0);
	    }
	i = i0;
	while (TWIN[i]) {
          i = ((i==MAXTWIN-1) ? 0 : ++i); if (i == i0) return(ERROR);
	  }
	TWIN[i] = a; return(i);
	}
	  

DEFINE	t_err	datw_init_twin(size)
/**********************************************************
  Initialisation of the array TWIN.
***********************************************************/
t_integer	size;
{
	t_integer	i;

	MAXTWIN = size;
	if (MAXTWIN > MAXATOM) xerr_c("cannot copy a molecule containing more than 100,000 atoms");
	for (i = 0; i <= MAXTWIN; i++)  TWIN[i] = NIL;
	return(OK);
	}

LOCAL	t_err	atom_twin(a)
/**********************************************************
***********************************************************/
atom_p	a;
{ 	if (a == NIL) return(ERROR);
        return(datw_set_twin(a->ID,a));
	}
DEFINE	t_err	datw_molecule_twin(molecule)
/**********************************************************
  Store all the atom of the molecule if the array TWIN.
***********************************************************/
molecule_p	molecule;
{
	if (molecule == NIL) return(ERROR);
	datw_init_twin(molecule->size);
	dast_all_atom(molecule,NIL,atom_twin,FALSE);
	return(OK);
	}
