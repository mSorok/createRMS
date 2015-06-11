/***********************************************************
  File of data structure allocation for the signature 
  Jean-loup Faulon Aug. 91 Pennstate
***********************************************************/
#include <general.h>
#include <eps.h>
#include "signature.h"

EXTERN	t_pointer hea_alloc();

DEFINE	atom_signature_p sial_create_atom_signature(HEAP,dim,dis,s,min,max)
/**********************************************************
 Create an atomic signature.
***********************************************************/
t_pointer	HEAP; 
t_integer	dim,dis;
p_string	s;
t_real		min,max;
{
	atom_signature_p	signature = NIL;

       	if (HEAP  == ((t_pointer)ERROR)) return((atom_signature_p)ERROR);
        signature = (atom_signature_p)hea_alloc(HEAP,sizeof(atom_signature_t));
        signature->dim = dim; signature->dis = dis;
 	signature->s = (p_string)hea_alloc(HEAP,2 * strlen(s)); 
        strcpy(signature->s,s);
        signature->min = min; signature->max = max;
	return(signature);
        }

DEFINE	t_err sial_size_signature(signature)
/**********************************************************
 Return the size of signature
***********************************************************/
signature_t	signature;
{
	t_integer	i;

	if (signature == NIL) return(0);
	for (i = 0; signature[i].as != NIL; i++);
	return(i);
        }

DEFINE	signature_p sial_create_signature(HEAP,size)
/**********************************************************
 Create a signature. Initialize to NIL.
 No more than MAXSIG elements
***********************************************************/
t_pointer	HEAP; 
t_integer	size;
{
	signature_p	signature;
	t_integer	i,s;

       	if (HEAP  == ((t_pointer)ERROR)) return((signature_p)ERROR);
       	if (HEAP  == NIL) return((signature_p)ERROR);
        s = (size+1) * sizeof(struct_signature_t);
        if (s < MAXHEAP/2) {
           signature = (signature_p)hea_alloc(HEAP,s);
	       }
	   else {
	       signature = (signature_p)malloc((size+1)*sizeof(struct_signature_t));
	       }
        for (i = 0; i < size+1; i++) siut_set_signature(signature,i,ZERO,NIL);
	    return(signature);
        }

LOCAL	signature_p copy_signature(HEAP,sig,SIG)
/**********************************************************
 sig = SIG. 
 WARNING the string are now copied !!! (Nov 1996)
***********************************************************/
t_pointer	HEAP;
signature_t	sig,SIG;
{
	t_integer		i;
	atom_signature_p	as;
	t_real			value;
	t_string		s;

        for (i = 0; SIG[i].as != NIL; i++) {
	    siut_inq_signature(SIG,i,&value,&as);
	    if (sig[i].as == NIL) {
	       s = (t_string)hea_alloc(HEAP,strlen(as->s)*sizeof(t_string)); strcpy(s,as->s);
               as = sial_create_atom_signature(HEAP,as->dim,as->dis,s,as->min,as->max);
	       }
	    siut_set_signature(sig,i,value,as);
	    } 
	return(sig);
	}

DEFINE	signature_p sial_copy_signature(HEAP,sig,SIG)
/**********************************************************
 sig = SIG.
 WARNING the string are now copied !!! (Nov 1996)
***********************************************************/
t_pointer	HEAP;
signature_p	sig,SIG;
{
	t_integer	SIZE, size;

       	if (SIG  == NIL) return((signature_p)ERROR);
       	if (SIG  == (signature_p)ERROR) return((signature_p)ERROR);
        if ((HEAP == NIL) && (sig == NIL)) return((signature_p)ERROR);
	SIZE = sial_size_signature(SIG);
	if (sig == NIL) size = SIZE; else size = sial_size_signature(sig);
	if (SIZE != size) return((signature_p)ERROR);
	if (sig == NIL) sig = sial_create_signature(HEAP,size);
	return(copy_signature(HEAP,sig,SIG));
	}
