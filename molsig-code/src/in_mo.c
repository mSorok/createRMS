#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "general.h"
#include "DFgeneral.h"
#include "DFManageMemory.h"
#include "allocation_da.h"
#include "structure_da.h"

DEFINE t_err	moin_read3000(file_name,file_type,operation,molecule,nci,MyGraph)
/******************************************************************************
 Read  .mol file
 Input :  - file name.
 Output : - MOLECULE
 WARNING :  Clorine atom and bond are removed.
 There is only one group = file_name 
 ******************************************************************************/
t_name		file_name,file_type,operation;
molecule_p	molecule;
t_bool		nci;
struct GRAPH *MyGraph;
{
#define MAXCHARLINE     1024 /* MAX NUMBER OF CHAR PER LINE */
	FILE 		*FD;
	t_name		mol_name;
	int	na,nb,deg = 0,i,chirality,charge;
	t_coord	    ZERO = {0,0,0};
	atom_p		*V;
	group_p		group;
	t_bufstring	name,dump,stereo,mass;
	char        line[MAXCHARLINE];
	char* ok_fgets;
	char completeName[MaxChaine];
	
	
	strcpy(mol_name,file_name); strcat(mol_name,".mol");
	if ((FD = fopen(mol_name,"r")) == NULL) return(ERROR); 
	
	
	//BEGIN READ COUNTS LINE
	/* move to the 4th line */		
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 1st line */
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 2nd line */
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 3rd line */
	ok_fgets =  fgets(line,MAXCHARLINE,FD); // skip the 4th line 
	ok_fgets =  fgets(line,MAXCHARLINE,FD); // skip the 5th line 
	
	ok_fgets =  fgets(line,MAXCHARLINE,FD); // skip the 6th line 
	char * ptr;//DF
	ptr = strtok (line+14," ");
	if(ptr!= NULL)	
		sscanf(ptr,"%d\n",&na);
	ptr = strtok(NULL," ");
	if(ptr!= NULL)	
		sscanf(ptr,"%d\n",&nb);
	//END READ COUNTS LINE
	
	ok_fgets = fgets(line,MAXCHARLINE,FD);//drop next line.
	
	group = (group_p)daal_create_group(molecule->HEAP,file_name,1,\
									   molecule,0,NIL,NIL,file_name,FALSE,NIL); 
	V = (atom_p *)calloc(na,sizeof(atom_p));
	
	//We know the size of the molecule, we can allocated the space needed for its 
	//desription as a Graph.
	if (MyGraph!=NIL)
		dfmm_CreateGRAPHMemory(na, nb, MyGraph);
	
	
	//BEGIN READ ALL ATOMS
	int *connectionT_map = calloc(na,sizeof(int));
	int index;
	
	for (i = 0; i < na; i++) {
		fgets(line,MAXCHARLINE,FD);
		float 	x,y,z;
		char * ptr;//DF
		char dump2[12];
		
		ptr = strtok (line,"M V");
		
		if(ptr!= NULL)	
			sscanf(ptr,"%d",&index);
		
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d",&index);
		connectionT_map[i] = index;
		
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%s",&name);
		else
			strcpy(name,"X");
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&x);
		else x= 0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&y);
		else y= 0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&z);
		else z= 0;
		
		charge = 0;
		strcpy(stereo,"0");
		strcpy(mass,"0");
		
		ptr = strtok (NULL," ");
		while(ptr != NULL){
			if(strncmp(ptr,"CHG=",4) == 0){
				sscanf(ptr+4,"%d",&charge);
			}
			if(strncmp(ptr,"CFG=",4) == 0){
				sscanf(ptr+4,"%s",stereo);
			}
			if(strncmp(ptr,"MASS=",5) == 0){
				sscanf(ptr+5,"%s",mass);
			}
			ptr = strtok (NULL," ");
		}
		
		//BEGIN GRAPH PART		
		
		//Atom by atom we store in its identity all the informations about it and we initialize 
		//variable to be used later. 
		if (MyGraph!=NIL) {
			dfmm_MakeIdentity(MyGraph[0].atom_identity,i,x,y,z,name,atoi(mass),charge,atoi(stereo));  //Change name later according to charge
			MyGraph[0].atom_identity[i].Name = Name2Number_atom(name);
		}
		
		if(atoi(mass)){
			sprintf(completeName,"%d%s",atoi(mass)+Name2Number_atom(name),name);//In V3000 we do not need +Name2Number_atom(name)
		}
		else
			sprintf(completeName,"%s",name);//In V3000 we do not need +Name2Number_atom(name)
		if(atoi(mass)){
			char aux[6];
			sprintf(aux,"%d",mass+Name2Number_atom(name));
			strncat(name,",i",3);
			strncat(name,aux,3);
		}
		
		if (charge) {
			if (charge == 1) strcat(completeName,"+++");
			if (charge == 2) strcat(completeName,"++");
			if (charge == 3) strcat(completeName,"+");
			//if (charge == 4) strcat(name,"");
			if (charge == 5) strcat(completeName,"-");
			if (charge == 6) strcat(completeName,"--");
			if (charge == 7) strcat(completeName,"---");
		}
		
		if (MyGraph!=NIL)
			strncpy(MyGraph[0].atom_identity[i].charName,completeName,12);
		//END GRAPH PART		
		
		if (mass[0] != '0') {
			if (mass[0] == '-') { strcat(name,",i"); strcat(name,mass); }
			else { strcat(name,",i+"); strcat(name,mass); }
		}
		if (charge) {
			if (charge == 1) strcat(name,"+3");
			if (charge == 2) strcat(name,"+2");
			if (charge == 3) strcat(name,"+1");
			if (charge == 4) strcat(name,"::");
			if (charge == 5) strcat(name,"-1");
			if (charge == 6) strcat(name,"-2");
			if (charge == 7) strcat(name,"-3");
		}
		V[i] = (atom_p)daal_create_atom(molecule->HEAP,name,i,group,NIL,\
										NIL,NIL,ZERO,name,(t_real)charge,deg,"",FALSE);
		strcpy(V[i]->comment,stereo);
	}
	//END READ ALL ATOMS
	int Max = 0;
	for(i=0;i<na;i++)
		if(connectionT_map[i] > Max)
			Max = connectionT_map[i];
	int* InvMap = calloc(Max + 1,sizeof(int));
	for(i=0;i<na;i++)
		InvMap[connectionT_map[i]] = i;
	
	molecule->size = na;
	
	ok_fgets =  fgets(line,MAXCHARLINE,FD); // skip two lines: 
	ok_fgets =  fgets(line,MAXCHARLINE,FD); //  
	/* read bonds */
	
	for (i = 0; i < nb; i++) {
		fgets(line,MAXCHARLINE,FD);
		int	n1,n2;
		int order = 0;
		strcpy(stereo,"0");
		char* ptr;
		
		ptr = strtok (line,"M V");
		if(ptr!= NULL)	
			sscanf(ptr,"%d",&index);
		
		ptr = strtok(NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d", &n1); //Tag to drop 
		ptr = strtok(NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d", &order);//atom 1  
		
		ptr = strtok(NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d", &n1);
		ptr = strtok(NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d", &n2); /* atom 2 */
		n1 = InvMap[n1];
		n2 = InvMap[n2];
		
		//as for the atoms the identity of the graph is initialized.
		if (MyGraph!=NIL)
			dfmm_AddBond(n1,n2,MyGraph,order,atoi(stereo));
		
		
		if (dast_inq_bond(V[n1],V[n2]) == NIL) { // JLF was n1-1 now is n1
			daal_create_bond(molecule->HEAP,V[n1],V[n2],NIL,NIL,ZERO,ZERO,order,"",FALSE);
			V[n1]->degre += order; V[n2]->degre += order;
		}
	}
	free(connectionT_map);
	free(InvMap);
	
	fclose(FD);	
	
	if (MyGraph!=NIL)
		dfmm_MakeNeigh(MyGraph);
	
	
	return(OK);
}


DEFINE	t_err	moin_read(file_name,file_type,operation,molecule,nci,MyGraph)
/******************************************************************************
 Read  .mol file
 Input :  - file name.
 Output : - MOLECULE
 WARNING :  Clorine atom and bond are remove.
 There is only one group = file_name 
 ******************************************************************************/
t_name		file_name,file_type,operation;
molecule_p	molecule;
t_bool		nci;
struct GRAPH *MyGraph;
{
#define MAXCHARLINE     1024 /* MAX NUMBER OF CHAR PER LINE */
	FILE 		*FD;
	t_name		mol_name;
	int	na,nb,deg = 0,i,chirality,charge;
	t_coord	    ZERO = {0,0,0};
	atom_p		*V;
	group_p		group;
	t_bufstring	name,dump,stereo,mass;
	char        line[MAXCHARLINE];
	char* ok_fgets;
	char completeName[MaxChaine];
	
	strcpy(mol_name,file_name); strcat(mol_name,".mol");
	if ((FD = fopen(mol_name,"r")) == NULL) return(ERROR); 
	
	/* move to the 4th line */		
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 1st line */
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 2nd line */
	ok_fgets = fgets(line,MAXCHARLINE,FD); /* skip the 3rd line */
	
	/* get the number of atoms and bonds and chirality */
	ok_fgets = fgets(dump,4,FD);
	sscanf(dump,"%d", &na); /* number atoms in can be in 3 spaces */
	ok_fgets = fgets(dump,4,FD); 
	sscanf(dump,"%d", &nb); /* number of bond can be in 3 spaces */
	ok_fgets = fgets(dump,7,FD); /* ignore obsolete information */
	ok_fgets = fgets(dump,4,FD);
	sscanf(dump,"%d", &chirality); /* stereochemistry */
	ok_fgets = fgets(line,MAXCHARLINE,FD);
	if (strstr(line,"V3000")) {
		fclose(FD);
		return(moin_read3000(file_name,file_type,operation,molecule,nci,MyGraph));
	}
	/*	printf("na = %d nb = %d chiral = %d\n",na,nb,chirality); */
	
	/* create atoms */
	group = (group_p)daal_create_group(molecule->HEAP,file_name,1,\
									   molecule,0,NIL,NIL,file_name,FALSE,NIL); 
	V = (atom_p *)calloc(na,sizeof(atom_p));
	
	
	//We know the size of the molecule, we can allocated the space needed for its 
	//desription as a Graph.
	if (MyGraph!=NIL)
		dfmm_CreateGRAPHMemory(na, nb, MyGraph);
	
	
	/* read all atoms */
	deg = 0; 
	
	for (i = 0; i < na; i++) {
		fgets(line,MAXCHARLINE,FD);
		float 	x,y,z;
		char * ptr;//DF
		ptr = strtok (line," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&x);
		else 
			x=0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&y);
		else y =0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%f",&z);
		else z = 0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%s",name);
		else
			strcpy(name,"X");
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%s",mass);
		else 
			strcpy(mass,"0");
		ptr = strtok (NULL," ");
		if(ptr!= NULL)	
			sscanf(ptr,"%d",&charge);
		else
			charge = 0;
		ptr = strtok (NULL," ");
		if(ptr!= NULL)		
			sscanf(ptr,"%s",stereo);
		else{strcpy(stereo,"0");}
		
		//Atom by atom we store in its identity all the informations about it and we initialize 
		//variable to be used later. 
		if (MyGraph!=NIL) {
			dfmm_MakeIdentity(MyGraph[0].atom_identity,i,x,y,z,name,atoi(mass),charge,atoi(stereo));  //Change name later according to charge
			MyGraph[0].atom_identity[i].Name = Name2Number_atom(name);
		}
		
		if(atoi(mass)){
			sprintf(completeName,"%d%s",atoi(mass)+Name2Number_atom(name),name);//In V3000 we do not need +Name2Number_atom(name)
		}
		else
			sprintf(completeName,"%s",name);//In V3000 we do not need +Name2Number_atom(name)
		if(atoi(mass)){
			char aux[6];
			sprintf(aux,"%d",mass+Name2Number_atom(name));
			strncat(name,",i",3);
			strncat(name,aux,3);
		}
		if (charge) {
			if (charge == 1) strcat(completeName,"+++");
			if (charge == 2) strcat(completeName,"++");
			if (charge == 3) strcat(completeName,"+");
			//if (charge == 4) strcat(name,"");
			if (charge == 5) strcat(completeName,"-");
			if (charge == 6) strcat(completeName,"--");
			if (charge == 7) strcat(completeName,"---");
		}		 
		if (mass[0] != '0') {
			if (mass[0] == '-') { strcat(name,",i"); strcat(name,mass); }
			else { strcat(name,",i+"); strcat(name,mass); }
		}
		if (charge) {
			if (charge == 1) strcat(name,"+3");
			if (charge == 2) strcat(name,"+2");
			if (charge == 3) strcat(name,"+1");
			if (charge == 4) strcat(name,"::");
			if (charge == 5) strcat(name,"-1");
			if (charge == 6) strcat(name,"-2");
			if (charge == 7) strcat(name,"-3");
		}
		V[i] = (atom_p)daal_create_atom(molecule->HEAP,name,i,group,NIL,\
										NIL,NIL,ZERO,name,(t_real)charge,deg,"",FALSE);
		if (MyGraph!=NIL)
			strncpy(MyGraph[0].atom_identity[i].charName,completeName,12);
		strcpy(V[i]->comment,stereo);
	}
	
	molecule->size = na;
	
	/* read bonds */
	for (i = 0; i < nb; i++) {
		int	n1,n2,order;
		ok_fgets = fgets(dump,4,FD);
		sscanf(dump,"%d", &n1); /* atom 1 */
		ok_fgets = fgets(dump,4,FD);
		sscanf(dump,"%d", &n2); /* atom 2 */
		ok_fgets = fgets(dump,4,FD);
		sscanf(dump,"%d", &order); /* bond order */
		ok_fgets = fgets(dump,4,FD);
		sscanf(dump,"%s",&stereo); /* stereo bond */

		//as for the atoms the identity of the graph is initialized.
		if (MyGraph!=NIL)
			dfmm_AddBond(n1,n2,MyGraph,order,atoi(stereo));
		

		if (dast_inq_bond(V[n1-1],V[n2-1]) == NIL) {
			daal_create_bond(molecule->HEAP,V[n1-1],V[n2-1],NIL,NIL,ZERO,ZERO,order,"",FALSE);
			V[n1-1]->degre += order; V[n2-1]->degre += order;
		}
		
		ok_fgets = fgets(line,MAXCHARLINE,FD); /* omit the rest of the line */
		
	}
	fclose(FD);
	
	if (MyGraph!=NIL)
		dfmm_MakeNeigh(MyGraph);
	
	return(OK);
}



