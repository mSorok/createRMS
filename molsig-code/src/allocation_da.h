/*
 *  allocation_da.h
 *  
 *
 *  Created by Davide Fichera on 03/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ALLOCATION_DA_H
#define ALLOCATION_DA_H

#include <string.h>
#include "hea.h"

molecule_p	daal_create_molecule(t_pointer TAS, t_name name, int ID, int size, set_group_p SG, t_pointer signature,t_pointer lmg,t_pointer lib, t_bufstring comment, 	t_bool marq);
bond_p	daal_create_bond( t_pointer TAS, atom_p a1,atom_p a2,bond_p b1,bond_p b2,t_coord geom1,t_coord geom2, t_bool order, t_bufstring comment, t_bool marq);
atom_p	daal_create_atom(t_pointer	TAS, t_name	name, int	ID, group_p group, set_bond_p SB, atom_p up,atom_p down, t_coord geom, t_name	potential_type, t_real charge, int degre, t_bufstring	comment, t_bool	marq); 
group_p	daal_create_group(t_pointer	TAS,  t_name	name, int	ID, molecule_p	molecule, int	degre, set_atom_p	SA, set_bond_p	SB, t_bufstring	comment, t_bool		marq, t_pointer	work); 

#endif
