/******************************************************************************
 Write MOL format
 ******************************************************************************/
#include <stdio.h>
#include <general.h>
#include <eps.h>

LOCAL	t_integer	BONDID = 0;
LOCAL	atom_p *	ATOM;

LOCAL	t_err	write_bond3000(b,FD)
/******************************************************************************
 Write bond 
 ******************************************************************************/
bond_p		b;
FILE		*FD;
{
	t_integer	order,n1,n2;       
	
	if (b->atom1 == NIL) return(OK);
	if (b->atom2 == NIL) return(OK);
	order = b->order; n1 = b->atom1->ID+1; n2 = b->atom2->ID+1; BONDID++;
	fprintf(FD,"M  V30 %d %d %d %d \n",BONDID,order,n1,n2);
	return(OK);
}

LOCAL	t_err	write_bond(b,FD)
/******************************************************************************
 Write bond 
 ******************************************************************************/
bond_p		b;
FILE		*FD;
{
	t_integer	order;       
	
	if (b->atom1 == NIL) return(OK);
	if (b->atom2 == NIL) return(OK);
	order = b->order;
	/*
	 order = ((b->order > 3) ? 1 : b->order);
	 if (b->atom1->potential_type[strlen(b->atom1->potential_type)-1] == 'p' && 
	 b->atom2->potential_type[strlen(b->atom2->potential_type)-1] == 'p')   
	 order=4; 
	 */      
	//DF add these variables in order to make standard print
	
	char Aux1[4];
	char Aux2[4];
	if(b->atom1->ID < 9)
		sprintf(Aux1,"  %1d",b->atom1->ID+1);
	else{ 
		if(b->atom1->ID < 99)
			sprintf(Aux1," %2d",b->atom1->ID+1);
		else
			sprintf(Aux1,"%3d",b->atom1->ID+1);
	}
	if(b->atom2->ID < 9)
		sprintf(Aux2,"  %1d",b->atom2->ID+1);
	else{ 
		if(b->atom2->ID < 99)
			sprintf(Aux2," %2d",b->atom2->ID+1);
		else
			sprintf(Aux2,"%3d",b->atom2->ID+1);
	}
	Aux1[3] = '\0'; Aux2[3] = '\0';
	//End: DF add these variables in order to make standard print
	//DF replaced the commented one with this.
	//	fprintf(FD,"%s%s  %d  0     0  0\n",Aux1,Aux2,order);
	//DF as above but with fixed order.
	fprintf(FD,"%s%s%3d  0     0  0\n",Aux1,Aux2,order);
	//  fprintf(FD,"%d %d %d 0 0 0\n",b->atom1->ID+1,b->atom2->ID+1,order);
	return(OK);
}

LOCAL	t_err	mol_atom(a,FD)
/******************************************************************************
 Write atom
 ******************************************************************************/
atom_p		a;
FILE		*FD;
{
	t_integer	a_pt = ERROR;
	t_name		element;
	
	if (a->name[0] == '\0') return(OK);
	ATOM[a->ID] = a;
	return(OK);
}

LOCAL	t_err	write_atom3000(FD)
/******************************************************************************
 Write atom
 ******************************************************************************/
FILE		*FD;
{
	t_integer	a_pt = ERROR,i;
	t_name		element;
	
	for (i = 0; ATOM[i]; i++) {
		atom_p	a = ATOM[i];
		strcpy(element,a->potential_type);
		t_integer j,l = strlen(element); element[4] = '\0';
		// remove charge from element
		for (j=0; j < l; j++) if ((element[j] == '+') || (element[j] == '-')) break ; 
		element[j] = '\0';
		fprintf(FD,"M  V30 %d %s %.4f %.4f %.4f 0\n",i+1,element,a->geom.x,a->geom.y,a->geom.z);

	}
	return(OK);
}

LOCAL	t_err	write_atom(FD)
/******************************************************************************
 Write atom
 ******************************************************************************/
FILE		*FD;
{
	t_integer	a_pt = ERROR,i;
	t_name		element;
	
	for (i = 0; ATOM[i]; i++) {
		atom_p	a = ATOM[i];
		strcpy(element,a->potential_type);
		t_integer j,l = strlen(element); element[4] = '\0';
		// remove charge from element
		for (j=0; j < l; j++) if ((element[j] == '+') || (element[j] == '-')) break ; 
		while (j < 4) element[j++] = ' ';
		fprintf(FD,"%10.4f%10.4f%10.4f",a->geom.x,a->geom.y,a->geom.z);
		fprintf(FD," %s",element); 
		fprintf(FD,"0  %1d  0  0  0  0  0  0  0  0  0  0\n",(int)a->charge);
	}
	return(OK);
}


DEFINE	t_err	moou_write3000(file_name,SM,extension)
/******************************************************************************
 Write  MOL file.mol 3000 version
 Output : - file name modified
 Input  : - MOLECULE
 ******************************************************************************/
t_name		file_name;
set_molecule_p	SM;
t_bool		extension;
{
	EXTERN	t_err		geut_translation_atom();
	FILE		*FDMOL;
	set_molecule_p  S;
	t_name		mol_name;
	molecule_p	molecule;
	t_integer	na = 0,nb = 0,i;
	
	
	if (SM == NIL) return(ERROR);
	if (SM->molecule == NIL) return(ERROR);
	
	strcpy(mol_name,file_name);  
	if (extension) strcat(mol_name,"_out.mol"); 
	else           strcat(mol_name,".mol");   
	if ((FDMOL = fopen(mol_name,"w")) == NULL) return(ERROR);      
	S = SM;
	while (S) { 
		na += S->molecule->size; 
		nb += dast_number_link(S->molecule,NIL,NIL);
		S = S->succ; if (S == NIL) break;
		if (S->molecule == NIL) break;
	}
	ATOM = (atom_p *)malloc((na+2)*sizeof(atom_p));
	for (i = 0; i <= na; i++) ATOM[i] = NIL;
	
	fprintf(FDMOL,"\n\n\n  0  0  0     0  0            999 V3000\n");  
	fprintf(FDMOL,"M  V30 BEGIN CTAB\n");
	fprintf(FDMOL,"M  V30 COUNTS %d %d 0 0 0\n",na,nb);

	fprintf(FDMOL,"M  V30 BEGIN ATOM\n");
	S = SM; 
	while (S) {
		molecule = S->molecule;
		if (molecule == NIL) break;
		dast_all_atom(molecule,NIL,mol_atom,FALSE,FDMOL);
		write_atom3000(FDMOL);
		S = S->succ; if (S == NIL) break;
	}
	fprintf(FDMOL,"M  V30 END ATOM\n");
	
	fprintf(FDMOL,"M  V30 BEGIN BOND\n");
	S = SM; 
	while (S) {
		molecule = S->molecule;
		if (molecule == NIL) break;
		dast_all_bond(molecule,NIL,NIL,write_bond3000,FALSE,FDMOL);
		S = S->succ; if (S == NIL) break;
	}
	fprintf(FDMOL,"M  V30 END BOND\n");
	fprintf(FDMOL,"M  V30 END CTAB\n");
	fprintf(FDMOL,"M  END\n");
	fclose(FDMOL);
	free(ATOM);
	return(OK);
}	



DEFINE	t_err	moou_write(file_name,SM,extension)
/******************************************************************************
 Write  MOL file.mol
 Output : - file name modified
 Input  : - MOLECULE
 ******************************************************************************/
t_name		file_name;
set_molecule_p	SM;
t_bool		extension;
{
	EXTERN	t_err		geut_translation_atom();
	FILE		*FDMOL;
	set_molecule_p  S;
	t_name		mol_name;
	molecule_p	molecule;
	t_integer	na = 0,nb = 0,i;
	
	
	if (SM == NIL) return(ERROR);
	if (SM->molecule == NIL) return(ERROR);
	
	strcpy(mol_name,file_name);  
	if (extension) strcat(mol_name,"_out.mol"); 
	else           strcat(mol_name,".mol");   
	S = SM;
	while (S) { 
		na += S->molecule->size; 
		nb += dast_number_link(S->molecule,NIL,NIL);
		S = S->succ; if (S == NIL) break;
		if (S->molecule == NIL) break;
	}
	if (na > 999) return(moou_write3000(file_name,SM,extension));
	
	if ((FDMOL = fopen(mol_name,"w")) == NULL) return(ERROR);      
	ATOM = (atom_p *)malloc((na+2)*sizeof(atom_p));
	for (i = 0; i <= na; i++) ATOM[i] = NIL;
    
	//DF 
	char Aux1[4];
	char Aux2[4];
	if(na < 9)
		sprintf(Aux1,"  %d",na);
	else{ 
		if(na < 99)
			sprintf(Aux1," %d",na);
		else{
			sprintf(Aux1,"%d",na);
			// printf("LLL%s %dLLL",Aux1,na);
		}
	}
	
	// printf("LLL%s %dLLL",Aux1,na);
	
	if(nb < 9)
		sprintf(Aux2,"  %d",nb);
	else{ 
		if(nb < 99)
			sprintf(Aux2," %d",nb);
		else
			sprintf(Aux2,"%d",nb);
	}
	
	//	printf("LLL%s %dLLL",Aux1,na);
	
	//	printf("FLAG %d %d %s %s ",na,nb,Aux1,Aux2);
	fprintf(FDMOL,"\n\n\n%s%s  0  0  1  0  0  0  0  0999 V2000\n",Aux1,Aux2);    
	//	fprintf(FDMOL,"%d %d\n",na,nb);   
	S = SM;
	
	S = SM; 
	while (S) {
		molecule = S->molecule;
		if (molecule == NIL) break;
		dast_all_atom(molecule,NIL,mol_atom,FALSE,FDMOL);
		write_atom(FDMOL);
		S = S->succ; if (S == NIL) break;
	}
	S = SM; 
	while (S) {
		molecule = S->molecule;
		if (molecule == NIL) break;
		dast_all_bond(molecule,NIL,NIL,write_bond,FALSE,FDMOL);
		S = S->succ; if (S == NIL) break;
	}
	fprintf(FDMOL,"M  END\n");
	fclose(FDMOL);
	free(ATOM);
	return(OK);
}	


