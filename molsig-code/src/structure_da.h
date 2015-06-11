/*
 *  structure_da.h
 *  
 *
 *  Created by Davide Fichera on 03/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef STRUSTURE_DA_H
#define STRUSTURE_DA_H



typedef struct ARG_EXIST {
	bond_p	bond;
	atom_p	atom;
} arg_exist_t, *arg_exist_p;

t_err	dast_all_bond(molecule_p molecule, group_p group, atom_p atom,f_traitement _function, t_bool marq,t_pointer arg);
bond_p	dast_inq_bond_directed(atom_p atom1,atom_p atom2);
t_err	dast_all_atom(molecule_p molecule, group_p group,f_traitement _function,t_bool marq,t_pointer arg);
t_err	dast_set_marq(molecule_p molecule,group_p group,atom_p atom,t_bool marq);
bond_p	dast_inq_bond(atom_p atom1,atom_p atom2);


#endif



