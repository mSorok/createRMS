/**********************************************                                                    
  Data structure for the signature
        manipulation .
 Jean-loup Faulon Aug. 91 Pennstate.
 Modify May 94.
**********************************************/


/* CONSTANTS                                 */
#define    MAXSIZE				MAXHEAP/2
#define  	MAXSIGSTRING		10000

/* TYPES                                     */

typedef  struct    ATOM_SIGNATURE {
      t_integer    dim,dis;
      t_string    s;
      t_real      min,max;
    } atom_signature_t, *atom_signature_p;

typedef struct    SIGNATURE {
      t_real      value;
            atom_signature_p  as;
      } struct_signature_t, signature_t[], *signature_p;

typedef  struct     SET_SIGNATURE {
      signature_p    signature;   
      struct SET_SIGNATURE  *succ;
    } set_signature_t, *set_signature_p;

typedef  struct     INFO_MOLECULAR_GROUP {
      t_bufstring  name;
      t_real    mass,bp,energy;
    } info_mg_t, *info_mg_p;  

typedef  struct     MOLECULAR_GROUP {
      t_pointer  group;
      signature_p  signature;
      t_real    min,max;
      t_pointer  work;  /* working field */
    } molecular_group_t, *molecular_group_p;   

typedef  struct     INFO_INTERGROUP_BOND {
      t_bufstring  name;
      t_integer  nature; /* New for kinetics purposes  */
      t_real    p1,p2,p3,p4,p5,p6,p7;
    } info_ib_t, *info_ib_p;  

typedef  struct     INTERGROUP_BOND {
      t_pointer  bond;
      signature_p  signature;
      t_real    min,max;
      t_pointer  work;
    } intergroup_bond_t, *intergroup_bond_p; 
 
typedef  struct     MASS_DATA {
      t_bufstring  name;
      t_real    min,max;
      t_pointer  work;
    } mass_data_t, *mass_data_p;  


typedef struct    LIST_MG {
      t_real      value;
            molecular_group_p  mg;
      } struct_lmg_t, list_mg_t[], *list_mg_p;    

typedef struct    LIST_IB {
      t_real      value;
      intergroup_bond_p  ib;
      } struct_lib_t, list_ib_t[], *list_ib_p;

typedef struct    LIST_MA {
      t_real      value;
      mass_data_p    ma;
      } struct_lma_t, list_ma_t[], *list_ma_p;

typedef struct    LIST_BP {
      t_real      mass,bp;
      } struct_lbp_t, list_bp_t[], *list_bp_p;

typedef  struct     ELEM_REACTION { 
      intergroup_bond_p  ibr,ibp; 
    } elem_reaction_t, *elem_reaction_p;  /* elementary reaction */

typedef struct    LIST_ER {
      elem_reaction_p  er;
      } struct_ler_t, list_er_t[], *list_er_p;
  
typedef  struct     SET_LMG {
      list_mg_p    lmg;   
      struct SET_LMG    *succ;
    } set_lmg_t, *set_lmg_p;

typedef  struct     SET_LIB {
      list_ib_p    lib; 
      struct SET_LIB    *succ;
    } *set_lib_p, set_lib_t;
