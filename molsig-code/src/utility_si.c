/***********************************************************
  File of utility for the signature 
  Jean-loup Faulon Aug. 91 Pennstate
  Modified Sept. 93.
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"
#define	 MAXSCHEDULE	100

/* local varialble--------------------------------------- */	
LOCAL	t_integer	NSCHEDULE = 0;
LOCAL	t_integer	MINR = MINRING, MAXR = MAXRING;
	                /* number of atoms per rings */
LOCAL	t_integer	MINNR = 0, MAXNR = MAXATOM;
	                /* number of rings */
LOCAL	t_integer	MINNF = 0, MAXNF = MAXATOM;
	                /* number of fragments */
LOCAL	t_real		MINGR = -1.0, MAXGR = -1.0;
LOCAL	t_real		MINCONS = 0.005; /* 0.5 percent */
LOCAL	t_bool		SETMINGR = FALSE, SETMAXGR = FALSE;
	                /* radius of gyration */
LOCAL	t_real		TEMPI[MAXSCHEDULE],TEMPF[MAXSCHEDULE],TEMPT[MAXSCHEDULE];
LOCAL	t_real		TIMEINC[MAXSCHEDULE],TIMEDELTA[MAXSCHEDULE];
LOCAL	t_real		MINE = -INFINI,MAXE = INFINI,MINF = 0,MAXF = MAXATOM;
LOCAL	t_real		RADIUS = 1;
LOCAL	t_bool		PLANAR = FALSE;         /* for classes algo, aut code */
LOCAL	t_integer	MULTIPLICITY = 100;
LOCAL	t_real		FSCL = 1.0;
LOCAL	t_real		CINF = 1.6;
LOCAL	t_err		DANGLING_END = 0;
LOCAL	t_err		DANGLING_LOOP = 0;
LOCAL	t_bool		OVERLAP  = FALSE;
LOCAL	t_real		BLENGTH = VDW; /* default is C-C bond */
LOCAL	t_real		RGRATIO = 6.0;
LOCAL	t_integer	KNOTSIZE = 1;
LOCAL	t_bufstring	RGSTAT = "none";
LOCAL	t_bufstring	NETYPE = "direct";
LOCAL	t_integer	DUERINGv  = 4.0;
LOCAL	t_real		DUERINGf  = 4.0;
LOCAL	t_real		DUERINGNc = 50.0;
LOCAL	t_bufstring	DUERINGx = "intra";
LOCAL	t_bool		DIRECTED = FALSE;
LOCAL	t_bool		MULTIGRAPH = FALSE;

/* Inquire and seting routines--------------------------- */
DEFINE	t_bool	siut_inq_directed() { return(DIRECTED); }
DEFINE	t_err	siut_set_directed(d) t_bool d; {DIRECTED = d; return(OK); }

DEFINE	t_bool	siut_inq_multigraph() { return(MULTIGRAPH); }
DEFINE	t_err	siut_set_multigraph(d) t_bool d; {MULTIGRAPH = d; return(OK); }

DEFINE	t_err	siut_inq_min_size_ring() { return(MINR); }
DEFINE	t_err	siut_inq_min_Nc() { return(MINR); }
DEFINE	t_err	siut_inq_min_number_ring() { return(MINNR); }

DEFINE	t_err	siut_inq_max_size_ring() { return(MAXR); }
DEFINE	t_real	siut_inq_min_gr() { return(MINGR); }
DEFINE	t_real	siut_inq_max_gr() { return(MAXGR); }
DEFINE	t_bool	siut_inq_set_min_gr() { return(SETMINGR); }
DEFINE	t_bool	siut_inq_set_max_gr() { return(SETMAXGR); }
        /* Lattice code max radius for g(r) calculations */
DEFINE	t_err	siut_inq_max_number_ring() { return(MAXNR); }

DEFINE	t_err	siut_inq_schedule_temp() { return(NSCHEDULE); }

DEFINE	t_real	siut_inq_min_size_fragment() { return(MINF); }
DEFINE	t_real	siut_inq_min_DP() { return(MINF); }
DEFINE	t_real	siut_inq_max_size_fragment() {	return(MAXF); }
DEFINE	t_real	siut_inq_max_DP() {	return(MAXF); }
DEFINE	t_err	siut_inq_min_number_fragment() { return(MINNF); }
DEFINE	t_err	siut_inq_max_number_fragment() { return(MAXNF); }

DEFINE	t_real	siut_inq_init_temp(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPI[n]); }
DEFINE	t_real	siut_inq_ratio_mute(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPI[n]); }

DEFINE	t_real	siut_inq_min_energy() { return(MINE); }

DEFINE	t_real	siut_inq_max_energy() { return(MAXE); }

DEFINE	t_real	siut_inq_final_temp(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPF[n]); }
DEFINE	t_real	siut_inq_ratio_mate(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPF[n]); }

DEFINE	t_bool	siut_inq_planar() { return(PLANAR); }
DEFINE	t_err	siut_set_planar(v) t_bool v; {PLANAR = v; return(OK); }

DEFINE	t_real	siut_inq_kinetics_radius() {return(RADIUS); }
DEFINE	t_err	siut_set_kinetics_radius(v) t_real v; {RADIUS = v; return(OK); }

DEFINE	t_real	siut_inq_incr_temp(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPT[n]); }
DEFINE	t_real	siut_inq_ratio_select(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TEMPT[n]); }

DEFINE	t_real	siut_inq_incr_time(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TIMEINC[n]); }

DEFINE	t_real	siut_inq_delta_time(n) t_integer n; {if (n > MAXSCHEDULE) return((t_real)ERROR); return(TIMEDELTA[n]); }

LOCAL	t_name		SOFT_IN,SOFT_OUT;
DEFINE	t_err	siut_set_soft_in (name) t_name name; {strcpy(SOFT_IN,name);}
DEFINE	t_err	siut_inq_soft_in (name) t_name name; {if (name) {strcpy(name,SOFT_IN); return(OK);} return(ERROR);}
DEFINE	t_err	siut_set_soft_out(name) t_name name; {strcpy(SOFT_OUT,name);}
DEFINE	t_err	siut_inq_soft_out(name) t_name name; {if (name) {strcpy(name,SOFT_OUT); return(OK);} return(ERROR);}

LOCAL	t_bufstring	FF;
DEFINE	t_err	siut_set_force_field(name) t_bufstring name; {strcpy(FF,name);}
DEFINE	t_err	siut_inq_force_field(name) t_bufstring name; {if (name) {strcpy(name,FF); return(OK);} return(ERROR);}

/* Kinetics */
DEFINE	t_err	siut_inq_multiplicity() { return(MULTIPLICITY); }
DEFINE	t_err	siut_set_multiplicity(n)
t_integer n; { MULTIPLICITY = n; return(OK); }
DEFINE	t_real	siut_inq_mincons() { return(MINCONS); }
DEFINE	t_err	siut_set_mincons(v)
t_real v; { MINCONS = v; return(OK); }
LOCAL	t_bufstring	RATE;
DEFINE	t_err	siut_set_rate_technique(name) t_bufstring name; {strcpy(RATE,name);}
DEFINE	t_err	siut_inq_rate_technique(name) t_bufstring name; {if (name) {strcpy(name,RATE); return(OK);} return(ERROR);}
DEFINE	t_bool	siut_rate_technique(name) t_bufstring name; {if (strcmp(name,RATE) == OK) return(TRUE); return(FALSE);}

LOCAL	molecule_p	MOLECULE = NIL;
DEFINE	t_err		siut_set_init_molecule(m) molecule_p m; {MOLECULE = m; return(OK);}
DEFINE	molecule_p	siut_inq_init_molecule() {return(MOLECULE);}

/* Lattice */

DEFINE	t_bool	siut_inq_net_type(netype)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	netype;
{
	if (strcmp(netype,NETYPE) == OK) return(TRUE);
	return(FALSE);
	}
DEFINE	t_err	siut_set_net_type(netype)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	netype;
{	strcpy(NETYPE,netype); return(OK); }

DEFINE	t_real	siut_inq_duering_f()
/******************************************************************************
 no comment.
******************************************************************************/
{
	return(DUERINGf);
	}
DEFINE	t_err	siut_set_duering_f(v)
/******************************************************************************
 no comment.
******************************************************************************/
t_real	v;
{	DUERINGf = v; return(OK); }

DEFINE	t_real	siut_inq_duering_Nc()
/******************************************************************************
 no comment.
******************************************************************************/
{
	return(DUERINGNc);
	}
DEFINE	t_err	siut_set_duering_Nc(v)
/******************************************************************************
 no comment.
******************************************************************************/
t_real	v;
{	DUERINGNc = v; return(OK); }

DEFINE	t_bool	siut_inq_duering_intra_inter(s)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	s;
{
	if (strcmp(s,DUERINGx) == OK) return(TRUE);
	return(FALSE);
	}
DEFINE	t_err	siut_set_duering_intra_inter(s)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	s;
{	strcpy(DUERINGx,s); return(OK); }

DEFINE	t_err		siut_set_knot_size(v)
/****************************************************************************** 
******************************************************************************/
t_integer	v;
{	KNOTSIZE = v; return(OK);
	}
DEFINE	t_integer	siut_inq_knot_size()
/****************************************************************************** 
******************************************************************************/
{
	return(KNOTSIZE);
	}

DEFINE	t_real	siut_inq_bond_length()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(BLENGTH); }

DEFINE	t_err	siut_set_bond_length(blength)
/******************************************************************************
 no comment.
******************************************************************************/
t_real blength;
{	BLENGTH = blength; return(OK); }

DEFINE	t_bool	siut_inq_rg_stat(calc)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	calc;
{
	if (strcmp(calc,RGSTAT) == OK) return(TRUE);
	return(FALSE);
	}
DEFINE	t_err	siut_set_rg_stat(calc)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring	calc;
{	strcpy(RGSTAT,calc); return(OK); }

DEFINE	t_real	siut_inq_rg_ratio()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(RGRATIO); }

DEFINE	t_err	siut_set_rg_ratio(rgratio)
/******************************************************************************
 no comment.
******************************************************************************/
t_real rgratio;
{	RGRATIO = rgratio; return(OK); }


DEFINE	t_bool	siut_inq_overlap()
/******************************************************************************
 Inquire DANGLING
******************************************************************************/
{
	return(OVERLAP);
	return(OK);
	}

DEFINE	t_err	siut_set_overlap(v)
/******************************************************************************
 Set DANGLING ends to v
******************************************************************************/
t_bool	v;
{
	OVERLAP = v;
	return(OK);
	}

DEFINE	t_err	siut_inq_dangling() {return(DANGLING_END+DANGLING_LOOP); }
DEFINE	t_err	siut_inq_dangling_end()  {return(DANGLING_END); }
DEFINE	t_err	siut_inq_dangling_loop() {return(DANGLING_LOOP); }

DEFINE	t_err	siut_set_dangling_end(v)
/******************************************************************************
 Set DANGLING_END to v
******************************************************************************/
t_err	v;
{
	DANGLING_END = v;
	return(OK);
	}

DEFINE	t_err	siut_set_dangling_loop(v)
/******************************************************************************
 Set DANGLING_LOOP to v
******************************************************************************/
t_err	v;
{
	DANGLING_LOOP = v;
	return(OK);
	}

DEFINE	t_real	siut_inq_scaling_factor() {return(FSCL);}
DEFINE	t_err	siut_set_scaling_factor(v)
/******************************************************************************
 Set scaling factor in end to end distance calculations
******************************************************************************/
t_real	v;
{	FSCL = v; return(OK); }

DEFINE	t_real	siut_inq_cinf() {return(CINF); }
DEFINE	t_err	siut_set_cinf(v)
/******************************************************************************
 Set Cinfinity in end to end distance calculations
******************************************************************************/
t_real	v;
{	CINF = v;return(OK); }


LOCAL	t_real	DENSITY = ERROR;
DEFINE	t_real	siut_inq_density()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(DENSITY); }

DEFINE	t_err	siut_set_density(density)
/******************************************************************************
 no comment.
******************************************************************************/
t_real density;
{	DENSITY = density; return(OK); }

LOCAL	t_real		MINW = -1, MAXW = -1;
DEFINE	t_real	siut_inq_min_weight()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(MINW); }
DEFINE	t_real	siut_inq_max_weight()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(MAXW); }
DEFINE	t_err	siut_set_weight(min,max)
/******************************************************************************
 no comment.
******************************************************************************/
t_real min,max;
{	MINW = min; MAXW = max; return(OK); }

LOCAL	t_integer	SCALE = ERROR;
DEFINE	t_err		siut_inq_scale()
/******************************************************************************
 no comment.
******************************************************************************/
{	return(SCALE); }

DEFINE	t_err	siut_set_scale(scale)
/******************************************************************************
 no comment.
******************************************************************************/
t_bufstring scale;
{
	if (scale == NIL) return(ERROR);

	     if (strcmp(scale,"micro") == OK) SCALE = MICRO;
	else if (strcmp(scale,"MICRO") == OK) SCALE = MICRO;
	else if (strcmp(scale,"meso") == OK)  SCALE = MESO;
	else if (strcmp(scale,"MESO") == OK)  SCALE = MESO;
	else if (strcmp(scale,"macro") == OK) SCALE = MACRO;
	else if (strcmp(scale,"MACRO") == OK) SCALE = MACRO;
	else if (strcmp(scale,"all") == OK)   SCALE = MICROANDMESO;
	else if (strcmp(scale,"ALL") == OK)   SCALE = MICROANDMESO;
	else if (strcmp(scale,"microtomeso") == OK)   SCALE = MICROTOMESO;
	else if (strcmp(scale,"MICROTOMESO") == OK)   SCALE = MICROTOMESO;
	else if (strcmp(scale,"mesotomicro") == OK)   SCALE = MESOTOMICRO;
	else if (strcmp(scale,"MESOTOMICRO") == OK)   SCALE = MESOTOMICRO;
	else if (strcmp(scale,"microandmeso") == OK)  SCALE = MICROANDMESO;
	else if (strcmp(scale,"MICROANDMESO") == OK)  SCALE = MICROANDMESO;
	else if (strcmp(scale,"mesoandmicro") == OK)  SCALE = MICROANDMESO;
	else if (strcmp(scale,"MESOANDMICRO") == OK)  SCALE = MICROANDMESO;
	else SCALE = ERROR;
	return(OK);
	}

DEFINE	t_err	siut_set_size_fragment(min,max)
/******************************************************************************
 no comment
******************************************************************************/
t_real min,max;
{	MINF = min; MAXF = max; 
	if ((MINF > MAXATOM) || (MAXF > MAXATOM)) xerr_c("too many atoms per fragment");
	if ((MINF < 0) || (MAXF < 0)) xerr_c("negative number of atom per fragment");
	return(OK); 
	}

DEFINE	t_err	siut_set_number_fragment(min,max)
/******************************************************************************
 no comment
******************************************************************************/
t_integer min,max;
{	if (min > max+1) return(ERROR);
	MINNF = min; MAXNF = max; 
	return(OK); 
	}

DEFINE	t_err	siut_set_size_ring(min,max)
/******************************************************************************
 no comment
******************************************************************************/
t_integer min,max;
{	MINR = min; MAXR = max; 
	if ((MINR > MAXRING) || (MAXR > MAXRING)) printf("Warning: MINR and/or MAXR are too large");
	if ((MINR < MINRING) || (MAXR < MINRING)) printf("Warning: MINR and/or MAXR are too small");
	return(OK); 
	}

DEFINE	t_err	siut_set_min_gr(min)
/******************************************************************************
 no comment
******************************************************************************/
t_real min;
{	MINGR = min; SETMINGR = TRUE;
	return(OK); 
	}

DEFINE	t_err	siut_set_max_gr(max)
/******************************************************************************
 no comment
******************************************************************************/
t_real max;
{	MAXGR = max;  SETMAXGR = TRUE;
	return(OK); 
	}

DEFINE	t_err	siut_set_number_ring(min,max)
/******************************************************************************
 no comment
******************************************************************************/
t_integer min,max;
{	if (min > max+1) return(ERROR);
	MINNR = min; MAXNR = max; 
	return(OK); 
	}

DEFINE	t_err	siut_set_energy(min,max)
/******************************************************************************
 no comment
******************************************************************************/
t_real min,max;
{	MINE = min; MAXE = max; return(OK); 
	}

DEFINE	t_err	siut_set_schedule_temp(n)
/******************************************************************************
 Set the number of temperature schedule
******************************************************************************/
t_integer n;
{	if (n > MAXSCHEDULE) xerr_c("too many temperature schedules");
	NSCHEDULE = n; return(OK); 
	}

DEFINE	t_err	siut_set_init_temp(v,n)
/******************************************************************************
 Set the Temperature
******************************************************************************/
t_real	v;
t_integer	n;
{	TEMPI[n] = v; return(OK); }

DEFINE	t_err	siut_set_final_temp(v,n)
/******************************************************************************
 Set the Temperature
******************************************************************************/
t_real	v;
t_integer	n;
{	TEMPF[n] = v; return(OK); }

DEFINE	t_err	siut_set_incr_temp(v,n)
/******************************************************************************
 Set the Temperature
******************************************************************************/
t_real	v;
t_integer	n;
{	TEMPT[n] = v; return(OK); }

DEFINE	t_err	siut_set_incr_time(v,n)
/******************************************************************************
 Set the time increment.
******************************************************************************/
t_real	v;
t_integer	n;
{	TIMEINC[n] = v; return(OK); }


DEFINE	t_err	siut_set_delta_time(v,n)
/******************************************************************************
 Set the time increment.
******************************************************************************/
t_real	v;
t_integer	n;
{	TIMEDELTA[n] = v; return(OK); }


DEFINE	t_err siut_inq_atom_signature(as,dim,dis,s,min,max)
/**********************************************************
 Return as's records
***********************************************************/
atom_signature_p	as;
t_integer		*dim,*dis;
p_string		*s;
t_real			*min,*max;
{

	if ((as == NIL) || (dim == NIL) || (dis == NIL) ||
             (s == NIL) || (min == NIL) || (max == NIL)) return(ERROR);	
	*dim = as->dim; *dis = as->dis; *s = as->s; *min = as->min; *max = as->max;
	return(OK);
	}


LOCAL	p_string neighbour(dim,hight,h,bond,b,s,found)
/**********************************************************
 Recursive function to find the neighbour.
***********************************************************/
t_integer	dim,hight,h,bond,b;
p_string	s;
t_bool		*found;
{
	t_integer	B = ((h == 0) ? VALENCE : VALENCE-1);

	if ((h == hight) && (b == bond)) *found = TRUE;
        if (*found) return(&s[0]); /* this one */ 
	if ((s[0] == '*') && (s[1] == '_')) return(&s[0]); /* next one */
	if ((s[0] == '.') && (s[1] == '_')) return(&s[0]); /* next one */
	if (dim == h) return(&s[0]); /* next one */

        for (b = 0; b < B; b++) {
	    s = neighbour(dim,hight,h+1,bond,b,&s[2],found);
	    if (*found) return(&s[0]);
            }
	return(&s[0]);
	}

LOCAL	p_string built_atom_signature(dim,hight,bond,s)
/**********************************************************
 Return the  signature of an atom defined by
 its hight and its bond.
 - s is the signature of the root,
 - hight is the hight of the atom in the tree,
 - bond is the number of the bond between the atom and 
   the root.
***********************************************************/
t_integer	dim,hight,bond;
p_string	s;
{
	t_integer	size = 0;
        t_bool		found = FALSE;

	if (s == NIL) return(NIL);
	if (bond > VALENCE) return(NIL);
	if ((s[0] == '*') && (s[1] == '_')) return(&s[0]);
	if ((s[0] == '.') && (s[1] == '_')) return(&s[0]);
	if (hight == 0) return(&s[0]);
        return(neighbour(dim,hight,hight - 1,bond,0,s,&found));
	}

LOCAL	t_bool equal(value,R,i,N,I)
/**********************************************************
 Search a configuration equal to r. 
 TRUE if found FALSE if not found.
 R is the matrix of relation between the two signatures,
 i is the line,
 value is the minimum value expected : EQ > LT > GT > NE 
 N is the last line,
 I is the selected arrow for the lines 0,1,..,i-1,
 WARNING algorithm not simple.
***********************************************************/
t_err	value,R[],N,i,I[];
{
	t_integer	j,k;
	t_bool		next = FALSE;
		
        if (i == N) return(TRUE);	
	for (j = 0; j < N; j++)
            if (R[i * N + j] >= value) {
	       next = FALSE; for (k = 0; k < i; k++) if (I[k] == j) next = TRUE;
               if (next) continue;
               I[i] = j;
               if (equal(value,R,i+1,N,I)) return(TRUE);
               }
        return(FALSE);
	}

DEFINE	t_err siut_minor_major(R,N)
/**********************************************************
 Compute the best relation between two signatures.
 R[i][j] contain the relation between the branch i and j.
 The method can be describe as follow :
 The values EQ,LT,GT are successively researched.
 For a given value all the permutations in R[i][j] are 
 performed. I[k] is the array of permutations.
***********************************************************/
t_err	R[],N;
{
LOCAL	t_err	I[MAXSIZE];
	t_err	i,j;
 	t_bool	next;

	for (i = 0; i < N; i++) I[i] = ERROR;
        /* if one line = NE result = NE */
        for (i = 0; i < N; i++) {
            next = TRUE;
            for (j = 0; j < N; j++) {
                if (R[i * N + j] != NE) next = FALSE;
                }
            if (next) return(NE);
            }                
        if (equal(EQ,R,0,N,I)) return(EQ);
        if (equal(LT,R,0,N,I)) return(LT);
        if (equal(GT,R,0,N,I)) return(GT);
        return(NE);
	}

LOCAL	t_err relation_potential_type(s1,s2)
/**********************************************************
 Return the  relation between two potential types ie :
 	s1 =  s2
 	s1 > (GT) s2 s1 include s2
 	s1 < (LT) s2 s1 is included in s2
 	s1 != s2
***********************************************************/
p_string	s1,s2;
{
	if ((s1 == NIL) || (s2 == NIL)) return(ERROR);

	if ((s1[0] == s2[0]) && (s1[1] == s2[1])) return(EQ);

	if ((s2[0] == '*') && (s2[1] == '_')) return(GT);
	if ((s1[0] == '*') && (s1[1] == '_')) return(LT);

	if ((s2[0] == 'X') && (s1[0] != 'h')) return(GT);
	if ((s1[0] == 'X') && (s2[0] != 'h')) return(LT);

	if ((s2[0] == 'x') && (s1[0] != 'h'))
	   if (s1[1] == s2[1]) return(GT);
	   
	if ((s1[0] == 'x') && (s2[0] != 'h'))
	   if (s1[1] == s2[1]) return(LT);

	return(NE);     	
	}

DEFINE	t_err siut_relation_string_signature(hight,dim1,s1,dim2,s2)
/**********************************************************
 Return the  relation between two strings s1,s2 ie :
 	s1 =  s2
 	s1 >  s2 s1 include s2
 	s1 <  s2 s1 is included in s2
 	s1 != s2
 WARNING algorithm not simple !!!!!
 The algorithm works even if the signature doesn't have the
 same dimension.
***********************************************************/
t_integer	dim1,dim2,hight;
p_string	s1,s2;
{
	t_integer	dim,b1,b2,B = ((hight == 0) ? VALENCE : VALENCE-1);
	p_string	ss1,ss2;
	t_err		r;
	t_err		R[VALENCE*VALENCE];

	dim = ((dim1 < dim2) ? dim1 : dim2);
        for (b1 = 0; b1 < VALENCE; b1++) for (b2 = 0; b2 < VALENCE; b2++) R[b1*B+b2] = ERROR;
	if ((s1 == NIL) || (s2 == NIL)) return(ERROR);
	if (hight > dim) return(ERROR);
	r = relation_potential_type(s1,s2);
	if (r == LT) if (dim1 > dim2) r = NE;
	if (r == GT) if (dim1 < dim2) r = NE;
	if (r != EQ) return(r);

/* r == EQ */
        if ((s1[0] == '*') && (s1[1] == '_')) return(r);
	if (hight == dim) {
	   if (dim1 == dim2) return(EQ);
	   if (dim1 <  dim2) return(LT);
	   if (dim1 >  dim2) return(GT);
	   }
          
	for (b1 = 0; b1 < B; b1++) {
	    ss1 = built_atom_signature(dim1,hight+1,b1,s1);
	    for (b2 = 0; b2 < B; b2++) {
                ss2 = built_atom_signature(dim2,hight+1,b2,s2);
                R[b1*B+b2] = siut_relation_string_signature(hight+1,dim1,ss1,dim2,ss2);

                }
            }
        r = siut_minor_major(R,B);
	return(r);     	
	}

DEFINE	t_err siut_relation_atom_signature(as1,as2)
/**********************************************************
 Return the  relation between two atom signatures ie :
 	as1 =  as2
 	as1 <  as2 as1 include as2 (checked first)
 	as1 >  as2 as1 is included in as2 (checked second)
 	as1 != as2
 WARNING the relation is not reflexive 
 ie a GT b not impose b LT a. 
***********************************************************/
atom_signature_p	as1,as2;
{
	t_err	relation;

	if ((as1 == NIL) || (as2 == NIL)) return(ERROR);
        if (as1->dis != as2->dis) return(NE);
	relation = siut_relation_string_signature(0,as1->dim,as1->s,as2->dim,as2->s); 
        if (relation == GT) {
           relation = siut_relation_string_signature(0,as2->dim,as2->s,as1->dim,as1->s);
           relation = ((relation == LT) ? GT : NE);
	   } 
	return(relation); 	
	}

DEFINE	t_err siut_set_signature(signature,i,value,as)
/**********************************************************
 no comment
***********************************************************/
signature_p		signature;
t_integer		i;
t_real			value;
atom_signature_p	as;
{
	if ((signature == NIL) || (i == ERROR)) return(ERROR);
        signature[i].as = as;
	signature[i].value = value;
	return(OK);
	}

DEFINE	t_err siut_inq_signature(signature,i,value,as)
/**********************************************************
 no comment
***********************************************************/
signature_t		signature;
t_integer		i;
t_real			*value;
atom_signature_p	*as;
{
	if ((signature == NIL) || (i == ERROR)) return(ERROR);
	if ((value == NIL) || (as == NIL)) return(ERROR);
	*value = signature[i].value;
        *as = signature[i].as;
	return(OK);
	}

DEFINE	t_real siut_average_signature(s)
/**********************************************************
 no comment
***********************************************************/
signature_p		s;
{
	t_integer	i = 0;
	t_real		value = 0,V = 0;
	t_pointer	as;

	if (s == NIL) return((t_real)0);
	if (s == (signature_p)ERROR) return(INFINI);
        LOOP {
	  siut_inq_signature(s,i,&value,&as);
	  if (as == NIL) break;
	  V += fabs(value);
          i++;
	  }
	if (i == 0) return((t_real)0);
	return((t_real)(V/(t_real)i));
	}

DEFINE	t_err siut_card_signature(s)
/**********************************************************
 Return the number of element of the signature.
 WARNING the absolute value is returned
***********************************************************/
signature_t	s;
{
	t_integer	n = 0,i = 0;

	if (s == NIL) return(0);
	for (i = 0; s[i].as != NIL; i++)
            n += ((s[i].value > 0) ? s[i].value : -s[i].value);
	return(n);
        }

LOCAL	t_bool pt_atom_signature(as,pt,d)
/**********************************************************
 Return TRUE if one of the first d+1 potential types of as 
 is equal to pt.
***********************************************************/
atom_signature_p	as;
t_name			pt;
t_integer		d;
{
	t_integer	i;

	if ((as == NIL) || (pt == NIL)) return(FALSE);
        if (as->s == NIL) return(FALSE);
	for (i = 0; i < d+1; i += 2) { 
	    if (as->s[i] == '\0') break;
	    if (as->s[i+1] == '\0') break;
	    if ((as->s[i] == pt[0]) && (as->s[i+1] == pt[1])) return(TRUE);
	    }
	return(FALSE);
	}

DEFINE t_err	siut_pt_signature(s,pt,d)
/**********************************************************
 Count the number of potential type equal to pt in s, d is
 the size of the signature (0,alpha,...).
***********************************************************/
signature_t	s;
t_name		pt;
t_integer	d;
{
	t_integer	i = 0;
	t_real		n = 0;

	if ((s == NIL) || (pt == NIL)) return(ERROR);
	for (i = 0; s[i].as != NIL; i++)
            if (pt_atom_signature(s[i].as,pt,d)) n += s[i].value;
	return(ROUND(n));
	}



/* SIGNATURE--------------------------------------------------------------------*/

DEFINE	t_real	siut_max_value_signature(s)
/**********************************************************
 Max value for signature s.
***********************************************************/
signature_p	s;
{
	t_integer		i = 0;
	t_real			v, maxv = -INFINI;
	atom_signature_p	as;

	if (s == NIL) return((t_real)ERROR);
	LOOP {
	  siut_inq_signature(s,i,&v,&as);
	  if (as == NIL) break;
	  if (v > maxv) maxv = v; i++;
	  }
	return(maxv);
	}	

DEFINE	t_err	siut_reset_value_signature(s)
/**********************************************************
 Reset the field value.
***********************************************************/
signature_p	s;
{
	t_integer		i = 0;
	t_real			v;
	atom_signature_p	as;

	if (s == NIL) return(ERROR);
	LOOP {
	  siut_inq_signature(s,i,&v,&as);
	  if (as == NIL) break;
          v = 0; siut_set_signature(s,i,v,as); i++;
	  }
	return(OK);
	}	

DEFINE	t_err	siut_print_signature(s)
/**********************************************************
 For debug.
***********************************************************/
signature_p	s;
{
	t_integer		i = 0;
	t_real			value;
	atom_signature_p	as;

	LOOP {
	  siut_inq_signature(s,i,&value,&as);
	  printf("value = %f ",value);
	  if (as == NIL) break;
	  printf("%s\n",as->s); i++;
	  }
	printf("\n");
	return(OK);
	}	

DEFINE	t_err	siut_print_signature_file(FDOUT,s)
/**********************************************************
 For debug.
***********************************************************/
FILE		*FDOUT;
signature_p	s;
{
	t_integer		i = 0;
	t_real			value;
	atom_signature_p	as;

	if (s == NIL) return(ERROR);
	LOOP {
	  siut_inq_signature(s,i,&value,&as);
	  if (as) fprintf(FDOUT,"%.1f ",value);
	  else { /* as == NIL */
	     if (sisi_inq_dimension_min() != sisi_inq_dimension_max() &&
		     sisi_inq_dimension(NIL) != sisi_inq_dimension_max()) break;
		 else fprintf(FDOUT,"%.1f ",value);
		 }
	  if (as == NIL) break;
	  fprintf(FDOUT,"%s\n",as->s); i++;
	  }
	fprintf(FDOUT,"\n");
	return(OK);
	}	
	

