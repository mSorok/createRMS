/******************************************************************************
  Writting scan
  Jean-Loup Faulon 08/29/00 modified 09/09/2009
******************************************************************************/
#include <stdio.h>
#include <general.h>
#include <eps.h>

DEFINE	t_err	scou_write_scan(file_name,SM,extension)
/******************************************************************************
 Write scan files
******************************************************************************/
t_name		file_name;
set_molecule_p	SM;
t_bool		extension;
{
EXTERN	t_err		geut_translation_atom();
EXTERN	t_pointer	sisc_signature();
	FILE		*FD;
	set_molecule_p  S;
	t_pointer	s;
	t_name		crtc_name,dim,dim_min,dim_max;
	molecule_p	molecule;

	if (SM == NIL) return(ERROR);
	if (SM->molecule == NIL) return(ERROR);

	strcpy(crtc_name,file_name);  
	if (extension) strcat(crtc_name,"_out.scan"); 
	else	       strcat(crtc_name,".scan"); 
	sprintf(dim,"%d",sisi_inq_dimension(NIL));
	if (sisi_inq_dimension_min() != sisi_inq_dimension_max()) {
	     sprintf(dim_min,"%d",sisi_inq_dimension_min());
	     sprintf(dim_max,"%d",sisi_inq_dimension_max());
	     strcat(crtc_name,dim_min);
		 strcat(crtc_name,"-");
		 strcat(crtc_name,dim_max);
		 }
	else strcat(crtc_name,dim);
	if ((FD = fopen(crtc_name,"w")) == NULL) return(ERROR);

	S = SM; 
	while (S) {
	  t_integer	h;
	  molecule = S->molecule;
	  if ((molecule == NIL) || (molecule == (molecule_p)ERROR)) break;
	  s = NIL;
	  for (h = sisi_inq_dimension_min(); h <= sisi_inq_dimension_max(); h++) {
	      sisi_set_dimension(h);
	      s = (t_pointer)sisc_signature(molecule->HEAP,molecule,NIL,NIL,NIL);
	      siut_print_signature_file(FD,s);
		  }
	  S = S->succ;
	  }
	fclose(FD);
	return(OK);
	}	

