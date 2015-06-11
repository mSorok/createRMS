/**********************************************                                                   
 Data structure for the signature
 program. 
 Jean-loup Faulon May 91 Albuquerque N.M.
 Modified PSU Jan 92 and May 92
**********************************************/
#ifndef	GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include <time.h>
#include <c.h>


/* CONSTANTS                                 */
#define    MAXATOM    100000
#define    VALENCE    4         /* for atoms */
#define    FUNCTIONALITY  4    /* for monomers */
#define    VDW    1.54
#define    MINRING    3
#define    MAXRING    20
#define    MAXBONDINGSITE  2000
#define    IBCREATE  1
#define    IBDELETE  0

#define    MICRO    1
#define    MESO    10
#define    MACRO    100
#define    MICROANDMESO  1000
#define    MICROTOMESO  110
#define    MESOTOMICRO  101


/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
   for dimension 5 one need :
   4 * ( 3xx6 + 3xx5 + ... + 3 + 1) + 1 = 4373 potential types
                                       -> 8746 characters
   4 * (        3xx5 + ... + 3 + 1) + 1 = 1457 pairs of ()
                                       -> 2914 characters
                                 total  = 11660 
          = MAXHEAP/2
         CHANGE OF MAXDIM => CHANGE OF MAXHEAP
ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/ 
 
#define			MAXHEAP    1048576 // PC: Increased the HEAP 48000


/* TYPES                                     */
typedef  struct  ATOM {
    t_name		name;        /* topology */
    t_integer	ID;
	t_pointer	group; /* NO RECURSION IN DECLARATIONS */
    t_pointer	SB;
	struct ATOM  *up,*down;			/* scaling */     
    t_coord		geom;				/* geometry */
    t_name		potential_type;    /* chemistry */
    t_real		charge;
    t_integer	degre;
    t_bufstring	comment;     /* work */
    t_bool		marq;
	} atom_t, *atom_p;

typedef  struct   SET_ATOM {
    atom_p    atom;   
    struct SET_ATOM  *succ;
    } set_atom_t, *set_atom_p;

typedef struct  BOND {
    atom_p    atom1,atom2; /* topology */
    t_coord    geom1,geom2; /* geometry */ 
                t_bool    order;       /* chemistry */ 
    t_bufstring  comment;     /* work */
    t_bool    marq,visit;     
    } bond_t, *bond_p;

typedef struct  SET_BOND {
    bond_p    bond;
    struct SET_BOND  *succ;
    } set_bond_t, *set_bond_p;

typedef struct  GROUP {
    t_name    name;  /* topology */
    t_integer  ID;
    t_pointer  molecule;
    t_integer  degre;
    set_atom_p  SA;
	set_bond_p  SB;
    t_bufstring  comment; /* work */
    t_bool    marq; 
	t_pointer  work;   
    } group_t, *group_p;

typedef struct  SET_GROUP {
    group_p				group;
    struct SET_GROUP	*succ;
    } set_group_t, *set_group_p;

typedef  struct  MOLECULE {
    t_pointer		HEAP; /* topologie */
    t_name			name;
    t_integer		ID;
	t_integer		size;
    set_group_p		SG;
	t_pointer       signature,lmg,lib; /* chemistry */
    t_bufstring		comment; /* work */
    t_bool			marq;
    } molecule_t, *molecule_p;

typedef struct  SET_MOLECULE {
    molecule_p    molecule;
    struct SET_MOLECULE  *succ;
    } set_molecule_t, *set_molecule_p;

#endif
