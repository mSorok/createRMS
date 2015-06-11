/***********************************************************
 FILE UTILITY for geometrical purpose
       Jean-loup Faulon - May 1991 Albuquerque N.M.
***********************************************************/
#include <general.h>
#include <eps.h>

/* DIMENSION */
LOCAL	t_err	GEOMDIM = 3;
DEFINE	t_err	geut_inq_dimension() { return(GEOMDIM);}
DEFINE	t_err	geut_set_dimension(v) 
t_integer v;
{ GEOMDIM = v; return(OK);}

typedef	struct	ARG {
	t_coord	min,max;
	} arg_t, *arg_p;

LOCAL	t_err	bondary_box(a,b)
/******************************************************************************
 Compute the hule of the molecule.
******************************************************************************/
atom_p		a;
arg_p		b;
{
	if (a->geom.x < b->min.x) b->min.x = a->geom.x;
	if (a->geom.y < b->min.y) b->min.y = a->geom.y;
	if (a->geom.z < b->min.z) b->min.z = a->geom.z;
	if (a->geom.x > b->max.x) b->max.x = a->geom.x;
	if (a->geom.y > b->max.y) b->max.y = a->geom.y;
	if (a->geom.z > b->max.z) b->max.z = a->geom.z;
        return(OK);
	}

DEFINE	t_err	geut_bondary_box(molecule,min,max,mark)
/******************************************************************************
 Compute the hule of the molecule.
******************************************************************************/
molecule_p	molecule;
p_coord		min,max;
t_bool		mark;
{
	arg_t	box;
	box.min.x = box.min.y = box.min.z = INFINI;
	box.max.x = box.max.y = box.max.z = -INFINI;
        dast_all_atom(molecule,NIL,bondary_box,mark,&box);
        *min = box.min; *max = box.max;
	return(OK);
	}

DEFINE	t_err	geut_translation_atom(a,t)
/******************************************************************************
 Applied translation t on atom a;
******************************************************************************/
atom_p		a;
p_coord		t;
{
        if ((a == NIL) || (t == NIL)) return(ERROR);
	a->geom.x += t->x; a->geom.y += t->y; a->geom.z += t->z;
        return(OK);
	}


