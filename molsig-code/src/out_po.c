/******************************************************************************
  From signature tp polygraf
  generation of data structures of signature cf signature.h
  Jean-loup Faulon PSU May 93.
  For modification cf release of biosym and polygtaf.
  For BIOSYM cf  REFERENCE GUIDE DISCOVER  2.7.0/3-31-91
  For POLYGRAF cf Programmer's Manuel for Polygraf, Biograf, NMRgraf and Polaris
               Version 3.1, September 14,1992 page I-22.
  Modify to allow large molecule (i.e. more than 99,999 atoms).
******************************************************************************/
#include <stdio.h>
#include <general.h>
#include <eps.h>

EXTERN  t_bool		MATCH;

typedef	struct	BOX {
	t_coord	min,max;
	} t_box, *p_box;

typedef	struct	ARG	{
        t_bufstring	s;
	t_integer	i;
	atom_p		atom;
	} arg_t, *arg_p;
LOCAL	t_integer	MAXLEN = ERROR;

LOCAL	t_err	po_bond(b,arg)
/******************************************************************************
 Write bond connected to atom
 Output : - string modified
 Input  : - string to be written, bond b , atom, increment 
 ******************************************************************************/
bond_p		b;
arg_p		arg;
{
EXTERN	t_real		log10();
	atom_p		a;
        t_name		id;
	t_integer	order;      

        a = ((b->atom1 == arg->atom) ? b->atom2 : b->atom1);
        if (a == NIL) return(OK);
	if (a->name[0] == '\0') return(OK);
	order = b->order; if (order == 4) order = 1;
	if (arg->s[0] == 'C') sprintf(id,"%d",a->ID+1);
	else sprintf(id,"%d",order);
	strcpy(&arg->s[arg->i+MAXLEN-strlen(id)],id); 
        arg->i += MAXLEN;
	return(OK);
	}

DEFINE	t_err	poou_po_atom_bond(a,FD)
/******************************************************************************
 Write CONECT Field
 Output : - file name modified
 Input  : - atom a ,file FD
 ******************************************************************************/
atom_p		a;
FILE		*FD;
{
	t_bufstring	s;
	t_name		id;
	t_integer	i,n = 0;
	arg_t		arg;

	if (a == NIL) return(OK);
	if (MAXLEN == ERROR) {
	   molecule_p	m = (molecule_p)((group_p)a->group)->molecule;
	   if (m->size > 99999) MAXLEN = ROUND(log10((t_real)m->size)) + 1;
	   else MAXLEN = 6;
	   }

	if (a->name[0] == '\0') return(OK);
	arg.atom = a; arg.i = 0;
	for (i = 0; i < MAXSTRING; i++) arg.s[i] = ' ';
	sprintf(id,"%d",a->ID+1);
	strcpy(&arg.s[arg.i],"CONECT");
	arg.i  = MAXLEN; strcpy(&arg.s[arg.i+MAXLEN-strlen(id)],id); 
        arg.i += MAXLEN; dast_all_bond(NIL,NIL,arg.atom,po_bond,FALSE,&arg);
	for (i = 0; i < arg.i; i++) if (arg.s[i] == '\0') arg.s[i] = ' '; strcpy(&arg.s[arg.i],"\0"); 
        fprintf(FD,"%s\n",arg.s);
	arg.atom = a; arg.i = 0;
	for (i = 0; i < MAXSTRING; i++) arg.s[i] = ' ';
	sprintf(id,"%d",a->ID+1);
	strcpy(&arg.s[arg.i],"ORDER");
	arg.i  = MAXLEN; strcpy(&arg.s[arg.i+MAXLEN-strlen(id)],id); 
        arg.i += MAXLEN; dast_all_bond(NIL,NIL,arg.atom,po_bond,FALSE,&arg);
	for (i = 0; i < arg.i; i++) if (arg.s[i] == '\0') arg.s[i] = ' '; strcpy(&arg.s[arg.i],"\0"); 
        fprintf(FD,"%s\n",arg.s);

	return(OK);
	}


DEFINE	t_err	poou_po_atom(a,FD)
/******************************************************************************
 Write atom
 Output : - file name modified
 Input  : - atom a ,file FD
 ******************************************************************************/
atom_p		a;
FILE		*FD;
{
	t_bufstring	s;
	t_name		a_name,a_id,a_x,a_y,a_z,g_name,g_id,degre,a_pt,a_charge;
	t_integer	i,n = 0,size;
	group_p		g;
	molecule_p	m;
	t_coord		coord,angle;

	for (i = 0; i < MAXSTRING; i++) s[i] = ' ';
	sprintf(a_name,"%s",a->name);
	dast_coord_pbc(a,&coord);
        sprintf(a_x,"%.5f",coord.x);
        sprintf(a_y,"%.5f",coord.y);
        sprintf(a_z,"%.5f",coord.z);
	g = (group_p)a->group; m = (molecule_p)g->molecule;
        sprintf(g_name,"%s",g->name);
	for (i = 0; g_name[i] != '\0'; i++) if (g_name[i] == '/') n = i;
	if (g_name[n] == '/') for (i = 0; i < MAXNAME -n - 2; i++) g_name[i] = g_name[i+1+n];
	for (i = 0; g_name[i] != '\0'; i++) if (g_name[i] == '_') break;

        g_name[i] = '\0'; sprintf(g_id,"%d",g->ID); /* strcpy(g_id,"1"); */
	a_pt[0] = '\0'; pout_potential_type(a_pt,a->potential_type);        
        sprintf(a_charge,"%.5f",a->charge);
        sprintf(a_id,"%d",a->ID+1);
        sprintf(degre,"%d",a->degre);

	if (m->size > 99999) {
	   fprintf(FD,"HETATM 	%s 	%s 	%s 	%s 	%s 	%s 	%s 	%s 	%s 	0 	%s\n",
		                a_id,a_name,g_name,g_id,a_x,a_y,a_z,a_pt,degre,a_charge);
	   return(OK);
	   }

	strcpy(&s[0],"HETATM");
        strcpy(&s[12-strlen(a_id)],a_id);
        strcpy(&s[13],a_name); s[18] = ' ';
        strcpy(&s[19],g_name);
        strcpy(&s[28-strlen(g_id)],g_id);
        strcpy(&s[40-strlen(a_x)],a_x);  
        strcpy(&s[50-strlen(a_y)],a_y);
        strcpy(&s[60-strlen(a_z)],a_z); 
	strcpy(&s[61],a_pt); 
	strcpy(&s[68],degre); 
        strcpy(&s[70],"0");
        strcpy(&s[73],a_charge); 
	for (i = 0; i < 81; i++) if (s[i] == '\0') s[i] = ' '; strcpy(&s[81],"\0");
        fprintf(FD,"%s\n",s);
	return(OK);
	}

DEFINE	t_err	poou_molecule(file_name,po,molecule)
/******************************************************************************
 Write  POLYGRAF files .bgf
 Output : - file names modified
 Input  : - MOLECULE
******************************************************************************/
t_name		file_name,po;
molecule_p	molecule;
{
EXTERN	t_err		dast_set_name_atom(),dast_set_ID_atom();
	FILE		*FDPO;
	t_bufstring	po_name,s;
	set_molecule_p  S;
	t_integer	na = 0, i;
	arg_t		arg;
	t_name		a_x,a_y,a_z;
	t_coord 	coord,angle;


        if (molecule == NIL) return(ERROR);
	/* i = 0; dast_all_atom(molecule,NIL,dast_set_ID_atom,FALSE,&i); */
        if (MATCH == FALSE) dast_all_atom(molecule,NIL,dast_set_name_atom,FALSE,NIL);
	strcpy(po_name,file_name); strcat(po_name,po); 
        if ((FDPO = fopen(po_name,"w")) == NULL) return(ERROR);    
        na += molecule->size;

	/* first string */
	fprintf(FDPO,"BIOGRF  160\n");
	fprintf(FDPO,"DESCRP TRANSLATOR\n");
	fprintf(FDPO,"REMARK This file was created by the TRANSLATOR program\n");
	fprintf(FDPO,"FORCEFIELD DREIDING\n");
dast_set_coord_in_pbc(FALSE);

	if (dast_inq_pbc(&coord,&angle) == OK) {
/*
dast_set_coord_in_pbc(TRUE);
*/
/*
printf("WARNING: PBC has been removed for .bgf files\n");
*/
	  fprintf(FDPO,"PERIOD 111     MOLXTL T\n");
	  fprintf(FDPO,"AXES   ZYX\n");
	  fprintf(FDPO,"SGNAME P1       C1(1)       1    1    1\n");
	  for (i = 0; i < MAXSTRING; i++) s[i] = ' ';
	  strcpy(&s[0],"CRYSTX");
          sprintf(a_x,"%.5f",coord.x);
          sprintf(a_y,"%.5f",coord.y);
          sprintf(a_z,"%.5f",coord.z);
          strcpy(&s[17-strlen(a_x)],a_x);  
          strcpy(&s[28-strlen(a_y)],a_y);
          strcpy(&s[39-strlen(a_z)],a_z);
          sprintf(a_x,"%.5f",angle.x);
          sprintf(a_y,"%.5f",angle.y);
          sprintf(a_z,"%.5f",angle.z);
          strcpy(&s[50-strlen(a_x)],a_x);  
          strcpy(&s[61-strlen(a_y)],a_y);
          strcpy(&s[72-strlen(a_z)],a_z);
	  for (i = 0; i < 73; i++) if (s[i] == '\0') s[i] = ' '; strcpy(&s[73],"\0");
	  fprintf(FDPO,"%s\n",s);
	  fprintf(FDPO,"CELLS     0    0    0    0    0    0\n");

	  }

	fprintf(FDPO,"FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)\n"); 
	dast_all_atom(molecule,NIL,poou_po_atom,FALSE,FDPO);
        fprintf(FDPO,"FORMAT CONECT (a6,12i6)\n");
	dast_all_atom(molecule,NIL,poou_po_atom_bond,FALSE,FDPO);
	fprintf(FDPO,"END\n"); 
	fclose(FDPO);
        return(OK);
	}


LOCAL	t_err	write(file_name,SM,extension)
/******************************************************************************
 Write  POLYGRAF files .bgf (one file for each solution)
 Output : - file names modified
 Input  : - MOLECULE
******************************************************************************/
t_name		file_name;
set_molecule_p	SM;
t_bool		extension;
{
	set_molecule_p  S;
	t_name		number,po;
	molecule_p	molecule;
	t_integer	n = 0,na = 0;

        if (SM == NIL) return(ERROR);
        if (SM->molecule == NIL) return(ERROR);
        S = SM; n = 0;	
        while (S) { 
          molecule = S->molecule;  if (molecule == NIL) break;
          n++; number[0] = '\0'; sprintf(number,"%d",n);
          strcpy(po,"");
	  if (extension) { strcpy(po,"_out"); strcat(po,number);}
	  strcat(po,".bgf");
	  poou_molecule(file_name,po,molecule);
          S = S->succ; if (S == NIL) break;
          }
        return(OK);
	}

DEFINE	t_err	poou_write(file_name,SM,extension)
/******************************************************************************
 Write  POLYGRAF file.bgf
 Output : - file name modified
 Input  : - MOLECULE
******************************************************************************/
t_name		file_name;
set_molecule_p	SM;
t_bool		extension;
{
EXTERN	t_err		dast_set_name_atom(),dast_set_ID_atom();
	FILE		*FDPO;
	set_molecule_p  S;
	t_bufstring	s;
	t_name		po_name,a_x,a_y,a_z;
	molecule_p	molecule;
	t_coord 	t,coord,angle;
	t_integer	n = 0,n2 = 0,na = 0, i;


        if (SM == NIL) return(ERROR);
        if (SM->molecule == NIL) return(ERROR);
	if (xsio_output() > 1) return(write(file_name,SM,extension));

	strcpy(po_name,file_name);  
	if (extension) strcat(po_name,"_out.bgf"); 
	else				 strcat(po_name,".bgf"); 
      
        if ((FDPO = fopen(po_name,"w")) == NULL) return(ERROR);      
        S = SM;
        while (S) { 
          na += S->molecule->size; n++; 
          S = S->succ; if (S == NIL) break;
          if (S->molecule == NIL) break;
          }
        while (n2*n2 < n) n2++;

	/* First strings */
	fprintf(FDPO,"BIOGRF  160\n");
	fprintf(FDPO,"DESCRP TRANSLATOR\n");
	fprintf(FDPO,"REMARK This file was created by the TRANSLATOR program\n");
	fprintf(FDPO,"FORCEFIELD DREIDING\n");
	if (dast_inq_pbc(&coord,&angle) == OK) {
	  fprintf(FDPO,"PERIOD 111     MOLXTL T\n");
	  fprintf(FDPO,"AXES   ZYX\n");
	  fprintf(FDPO,"SGNAME P1       C1(1)       1    1    1\n");
	  for (i = 0; i < MAXSTRING; i++) s[i] = ' ';
	  strcpy(&s[0],"CRYSTX");
          sprintf(a_x,"%.5f",coord.x);
          sprintf(a_y,"%.5f",coord.y);
          sprintf(a_z,"%.5f",coord.z);
          strcpy(&s[17-strlen(a_x)],a_x);  
          strcpy(&s[28-strlen(a_y)],a_y);
          strcpy(&s[39-strlen(a_z)],a_z);
          sprintf(a_x,"%.5f",angle.x);
          sprintf(a_y,"%.5f",angle.y);
          sprintf(a_z,"%.5f",angle.z);
          strcpy(&s[50-strlen(a_x)],a_x);  
          strcpy(&s[61-strlen(a_y)],a_y);
          strcpy(&s[72-strlen(a_z)],a_z);
	  for (i = 0; i < 73; i++) if (s[i] == '\0') s[i] = ' '; strcpy(&s[73],"\0");
	  fprintf(FDPO,"%s\n",s);
	  fprintf(FDPO,"CELLS     0    0    0    0    0    0\n");
	  }
	fprintf(FDPO,"FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)\n"); 
	t.x = 0; t.y = 0; t.z = 0; S = SM; n = 0;
        while (S) {
          n++; molecule = S->molecule;
          if (molecule == NIL) break;
	  /* i = 0; dast_all_atom(molecule,NIL,dast_set_ID_atom,FALSE,&i); */
          if (MATCH == FALSE) dast_all_atom(molecule,NIL,dast_set_name_atom,FALSE,NIL);
	  dast_all_atom(molecule,NIL,poou_po_atom,FALSE,FDPO);
          fprintf(FDPO,"FORMAT CONECT (a6,12i6)\n");
	  dast_all_atom(molecule,NIL,poou_po_atom_bond,FALSE,FDPO);
          S = S->succ; if (S == NIL) break;
          if (n > n2) n = 0;
          }

	fprintf(FDPO,"END\n"); 
	fclose(FDPO);
        return(OK);
	}	


