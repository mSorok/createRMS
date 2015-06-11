/******************************************************************************
  From polygraf to signature and reverse module utility
  generation of data structures of signature cf signature.h
  Jean-loup Faulon PSU May 93.
  For modification cf release of biosym and polygtaf.
  For BIOSYM cf  REFERENCE GUIDE DISCOVER  2.7.0/3-31-91
  For MOL cf Programmer's Manuel for Polygraf, Biograf, NMRgraf and Polaris
               Version 3.1, September 14,1992 page I-22.

******************************************************************************/
#include <stdio.h>
#include <general.h>


DEFINE	t_err	mout_potential_type(pt_mol,pt_signature)
/******************************************************************************
 Correspondance for the potential types.
 The potential type of signature are thoose of biosym
 At -> Cl when outputed
******************************************************************************/
t_name		pt_mol,pt_signature;
{
	if ((strcmp(pt_mol,"") == OK) && (strcmp(pt_signature,"") == OK))
            return(ERROR);

	/* SIGNATURE -> MOL */
        if (strcmp(pt_mol,"") == OK) {
	if      (strcmp(pt_signature,"h_") == OK) strcpy(pt_mol,"H");

        else if (strcmp(pt_signature,"c_") == OK) strcpy(pt_mol,"C");
        else if (strcmp(pt_signature,"cp") == OK) strcpy(pt_mol,"C");
        else if (strcmp(pt_signature,"c=") == OK) strcpy(pt_mol,"C");
        else if (strcmp(pt_signature,"ct") == OK) strcpy(pt_mol,"C");

        else if (strcmp(pt_signature,"o_") == OK) strcpy(pt_mol,"O"); 
        else if (strcmp(pt_signature,"o'") == OK) strcpy(pt_mol,"O");
        else if (strcmp(pt_signature,"o=") == OK) strcpy(pt_mol,"O");
 
	else if (strcmp(pt_signature,"n5") == OK) strcpy(pt_mol,"N");
        else if (strcmp(pt_signature,"n_") == OK) strcpy(pt_mol,"N");
        else if (strcmp(pt_signature,"n2") == OK) strcpy(pt_mol,"N");
        else if (strcmp(pt_signature,"np") == OK) strcpy(pt_mol,"N");
        else if (strcmp(pt_signature,"nt") == OK) strcpy(pt_mol,"N");
  
	else if (strcmp(pt_signature,"s6") == OK) strcpy(pt_mol,"S");
	else if (strcmp(pt_signature,"s4") == OK) strcpy(pt_mol,"S");
        else if (strcmp(pt_signature,"s_") == OK) strcpy(pt_mol,"S");    
	else if (strcmp(pt_signature,"s=") == OK) strcpy(pt_mol,"S");

        else if (strcmp(pt_signature,"si") == OK) strcpy(pt_mol,"Si");

        else if (strcmp(pt_signature,"f_") == OK) strcpy(pt_mol,"F");
        else if (strcmp(pt_signature,"cl") == OK) strcpy(pt_mol,"Cl");
        else if (strcmp(pt_signature,"br") == OK) strcpy(pt_mol,"Br");
        else if (strcmp(pt_signature,"i_") == OK) strcpy(pt_mol,"I");
        else if (strcmp(pt_signature,"at") == OK) strcpy(pt_mol,"AT");
        else if (strcmp(pt_signature,"p5") == OK) strcpy(pt_mol,"P");
	else if (strcmp(pt_signature,"p_") == OK) strcpy(pt_mol,"P");
	else if (strcmp(pt_signature,"na") == OK) strcpy(pt_mol,"Na");

        /* super atom must have different colors */
	else if (pt_signature[0] == 'Z') {
             if (pt_signature[2] == '1')          strcpy(pt_mol,"C");
             else if (pt_signature[2] == '2')     {
                  if (pt_signature[1] == '1')     strcpy(pt_mol,"O");
	          else 				  strcpy(pt_mol,"O");
	          }
             else if (pt_signature[2] == '3')     strcpy(pt_mol,"N");
             else if (pt_signature[2] == '4')     strcpy(pt_mol,"S");
             else                                 strcpy(pt_mol,"C");
	     } /* Z */
        else {  printf("WARNING unknown potential type : %s\n",pt_signature); 
                return(ERROR);
             }
        return(OK);
        }

	/* MOL -> SIGNATURE */
        else {
	pt_mol[3] = '\0';
	if (pt_mol[2] == '_') pt_mol[2] = '\0';
	if      (strcmp(pt_mol,"H") == OK) strcpy(pt_signature,"h_");
        else if (strcmp(pt_mol,"C") == OK) strcpy(pt_signature,"c_");
        else if (strcmp(pt_mol,"O") == OK)  strcpy(pt_signature,"o_"); 
	else if (strcmp(pt_mol,"N") == OK)  strcpy(pt_signature,"n_");	
        else if (strcmp(pt_mol,"S") == OK)  strcpy(pt_signature,"s_");    
        else if (strcmp(pt_mol,"Si") == OK)  strcpy(pt_signature,"si");
        else if (strcmp(pt_mol,"F") == OK) strcpy(pt_signature,"f_");
        else if (strcmp(pt_mol,"Cl") == OK) strcpy(pt_signature,"cl");
        else if (strcmp(pt_mol,"Br") == OK) strcpy(pt_signature,"br");
        else if (strcmp(pt_mol,"I_") == OK) strcpy(pt_signature,"i_");
        else if (strcmp(pt_mol,"At") == OK) strcpy(pt_signature,"at");
        else if (strcmp(pt_mol,"P") == OK) strcpy(pt_signature,"p_");
	else if (strcmp(pt_mol,"Na") == OK) strcpy(pt_signature,"na");
        else {  printf("WARNING unknown potential type : %s\n",pt_mol); 
	        strcpy(pt_signature,"?_");
                return(ERROR);
             }
        return(OK);
        }
	}

