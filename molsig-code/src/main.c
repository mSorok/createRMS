/***********************************************************
 MAIN PROGRAM SIGNATURE TRANSLATOR
 Jean-loup Faulon
 Major updates: May 1996, July 2002, July 2003, Feb 2006
 Sept 2009
 ***********************************************************/
#include	<stdio.h>
#include	<stdlib.h>
#include	<general.h>
#include	<eps.h>
#include 	<time.h>
#include "signature.h"
#include "DFgeneral.h"
#include "PCStereoSig.h"
EXTERN t_pointer heat_creer_hea(), daal_create_molecule();
DEFINE FILE *FDOUT = NIL;
DEFINE t_bool MATCH = FALSE;
DEFINE t_bool ORIENTED = FALSE;
LOCAL t_real MINGR = 1.0, MAXGR = 3.0;
LOCAL t_bool SETMINGR = FALSE, SETMAXGR = FALSE;
LOCAL t_integer MAXSOLUTION = MAXINTEGER;
LOCAL t_integer MULTIPLE = 2;
LOCAL t_bufstring SCANTYPE;
// PC
LOCAL t_bufstring USERSCAN;
LOCAL t_bool STEREO_MODE = FALSE;
LOCAL t_bool FAST_MODE = FALSE;
t_bool OPT_NEWCIP = TRUE;
//

DEFINE t_err xsio_output()
/**********************************************************************
 no comment.
 **********************************************************************/
{
	return (0);
}

DEFINE p_string xsio_scan_type()
/**********************************************************************
 no comment.
 **********************************************************************/
{
	return (&SCANTYPE[0]);
}

DEFINE t_err xsma_set_maxsolution(value)
	/**************************************************************
	 no comment.
	 **************************************************************/
	t_err value; {
	MAXSOLUTION = value;
	return (OK);
}

DEFINE t_err xsma_inq_maxsolution()
/**************************************************************
 no comment.
 **************************************************************/
{
	return (MAXSOLUTION);
}

DEFINE t_err xsio_file_read(file_name, file_read)
	/******************************************************************************
	 no comment
	 ******************************************************************************/
	t_name file_name;f_traitement *file_read; /* function that reads the molecular file */
{
	return (ERROR);
}

DEFINE t_err xsio_file_write(file_name, SM)
	/******************************************************************************
	 no comment
	 ******************************************************************************/
	t_name file_name;set_molecule_p SM; {
	return (ERROR);
}


LOCAL void print_usage() {
	printf("Usage: sscan <file.mol> <output-type>\n");
	printf("\t<file> mol file\n");
	printf("\t<output-type> scanx, sscanx, sigx, ssigx\n");
	printf("\tx between 0,10000; x omitted: canonicalize structure (scan,sscan)\n");
}

LOCAL t_err write_scan(file_name, SM, extension)
	/******************************************************************************
	 Write scan files
	 ******************************************************************************/
	t_name file_name;set_molecule_p SM;t_bool extension; {
	EXTERN t_err geut_translation_atom();
	EXTERN t_pointer sisc_signature();
	FILE *FD;
	set_molecule_p S;
	t_pointer s;
	t_name crtc_name, dim, dim_min, dim_max;
	molecule_p molecule;

	if (SM == NIL)
		return (ERROR);
	if (SM->molecule == NIL)
		return (ERROR);

	strcpy(crtc_name, file_name);
	if (extension)
		strcat(crtc_name, "_out");
	strcat(crtc_name, ".");
	strcat(crtc_name, USERSCAN);

	sprintf(dim, "%d", sisi_inq_dimension(NIL));
	if (sisi_inq_dimension_min() != sisi_inq_dimension_max()) {
		sprintf(dim_min, "%d", sisi_inq_dimension_min());
		sprintf(dim_max, "%d", sisi_inq_dimension_max());
		strcat(crtc_name, dim_min);
		strcat(crtc_name, "-");
		strcat(crtc_name, dim_max);
	} else {
		if (atoi(dim) != MAXINTEGER) // PC: don't add height if it's canonical
			strcat(crtc_name, dim);
	}
	if ((FD = fopen(crtc_name, "w")) == NULL) return (ERROR);

	S = SM;
	while (S) {
		t_integer h;
		molecule = S->molecule;
		if ((molecule == NIL) || (molecule == (molecule_p) ERROR))
			break;
		s = NIL;
		for (h = sisi_inq_dimension_min(); h <= sisi_inq_dimension_max(); h++) {
			sisi_set_dimension(h);
			s = (t_pointer) sisc_signature(molecule->HEAP, molecule, NIL, NIL,
					NIL);
			siut_print_signature_file(FD, s);
		}
		S = S->succ;
	}
	fclose(FD);
	return (OK);
}

LOCAL void set_diameter(int diameter)
/******************************************************************************
 Set current diameter
 ******************************************************************************/
{
	sisi_set_dimension(diameter);
	sisi_set_dimension_min(diameter);
	sisi_set_dimension_max(diameter);
}

main(argc, argv)
	/**********************************************************************
	 no comment
	 **********************************************************************/
	long argc;char *argv[]; {
	EXTERN set_molecule_p isin_isomer();
	EXTERN t_err biin_read(), pcin_read();
	EXTERN t_real dast_string_real(), dast_molecular_weight();
	t_pointeur TAS = NIL, SIGNATURE = NIL, LG = NIL, LB;
	molecule_p molecule, groups, bonds;
	set_molecule_t SM;
	t_bufstring w, file, file_name, file_atom, operation;
	t_bufstring input, output, ff, option;
	t_integer dim = 1, dim_min = 1, dim_max = 1, col = MAXINTEGER;
	t_integer t1, t2, i, n = 0, err = ERROR, r = 0;
	t_real bp = 0;

	// PC: Start options section
	int options;
	t_bool OPT_GLOBAL_STEREOSIG = TRUE;
	t_bool OPT_PRINT_WEIGHTS = FALSE;
	t_bool OPT_PRINT_ATOMNUMBER = FALSE;
	t_bool OPT_CLOCK = FALSE;
	int OPT_RANDOMIZER = 0;
	int OPT_MOL_VERSION = 2000;
	VERBOSE = FALSE;
	while ((options = getopt(argc, argv, "lmntvwopk:r:")) != -1) {
		switch (options) {
		case 'k': // Mol version ([2000],3000)
			OPT_MOL_VERSION = atoi(optarg);
			break;
		case 'l': // Compute local stereochemistry
			OPT_GLOBAL_STEREOSIG = FALSE;
			break;
		case 'm':
			OPT_PRINT_ATOMNUMBER = TRUE;
			break;
		case 'o': // Oriented
			ORIENTED = TRUE;
			break;
		case 'r': // Call randomizer
			OPT_RANDOMIZER = atoi(optarg);
			break;
		case 't': // Timer
			OPT_CLOCK = TRUE;
			break;
		case 'v':
			VERBOSE = TRUE;
			break;
		case 'p':
			PCINFO = TRUE;
			break;
		case 'w': // Compute local stereochemistry
			OPT_PRINT_WEIGHTS = TRUE;
			break;
		case 'n': // Use old cip computation
			OPT_NEWCIP = FALSE;
			break;
		}
	}
	argc -= optind - 1; // A quick way to skip the options
	argv += optind - 1;
	argv[0] = argv[1 - optind];
	// PC: Ends options sections


	if (VERBOSE) {
		printf("\n--------SCAN PROGRAM Version 1.00-----------\n");
		printf("Jean-loup Faulon - Sept 2009\n");
		printf("email : jfaulon@gmail.com\n");
		printf("http  : //jfaulon.wikispaces.com\n");
		printf("\nLimitations of the program\n");
		printf("-----input:\n");
		printf("     .mol (chemical) \n");
		printf("-----output only:\n");
		printf("     .scan .sig .cip signature (o) \n");
		printf("-----no more than 30,000 atoms\n");
	}
	file_atom[0] = '\0';
	if (argc == 3) {
		strcpy(output, argv[argc - 1]);
		strcpy(file_name, argv[argc - 2]);
	} else {
		print_usage();
		return(0);
	}

	char* ext = strstr(file_name,".mol");
	if (ext == NIL) {
		print_usage();
		return(0);  // PC Print Usage
	}
	ext[0] = '\0';

	// PC : predefined values ([s]scan alone means canonical)
	if ((strcmp(output, "sscan") == OK) || (strcmp(output, "scan") == OK)
			|| (strcmp(output, "cip") == OK) || (strcmp(output, "fscan") == OK)) {
		dim = MAXINTEGER;
		dim_min = MAXINTEGER;
		dim_max = MAXINTEGER;
	} else if ((strcmp(output,"sig") == OK)||(strcmp(output,"ssig") == OK)) {
		dim = 0;
		dim_min = 0;
		dim_max = 0;
	}
	//PC END

	for (i = 0; output[i] != '\0'; i++)
		if (('0' <= output[i]) && (output[i] <= '9')) {
			t_bufstring nb1, nb2;
			t_integer v1, v2, j;
			strcpy(nb1, &output[i]);
			for (j = 0; nb1[j] != '\0'; j++)
				if (nb1[j] < '0' || nb1[j] > '9')
					break;
			if (nb1[j] != '\0')
				strcpy(nb2, &nb1[j + 1]);
			else
				nb2[0] = '\0';
			v1 = atoi(nb1);
			if (nb2[0] != '\0')
				v2 = atoi(nb2);
			else
				v2 = v1;
			output[i] = '\0';
			//PC
			if ( (strncmp(output, "scan", 4) == OK) 	|| (strncmp(output, "sscan", 5) == OK)
					|| (strncmp(output, "cip", 3) == OK) || (strncmp(output, "fscan", 5) == OK)) {
				if (v1>MAXINTEGER) v1 = MAXINTEGER;
				if (v2>MAXINTEGER) 	v2 = MAXINTEGER;
				dim = v1;
				dim_min = v1;
				dim_max = v2;
			} else if ((strncmp(output, "sig",3) == OK)||(strncmp(output,"ssig",4) == OK)) {
				if (v1>12) v1 = 12;
				if (v2>12) v2 = 12;
				dim = v1;
				dim_min = v1;
				dim_max = v2;
			}
			//PC END
			if (VERBOSE)
				printf("signature dimension : %d\n", dim);
			break;
		}

	if ((dim>1) && ((strncmp(output, "sscan", 5) == OK)||(strncmp(output, "fscan", 5) == OK)
			|| (strncmp(output, "ssig", 4) == OK)))  {
				STEREO_MODE = TRUE;
	}
	if (strncmp(output, "fscan", 5) == OK)  {
		FAST_MODE = TRUE;
	}

	strcpy(SCANTYPE, output);
	sisi_set_dimension(dim);
	sisi_set_dimension_min(dim_min);
	sisi_set_dimension_max(dim_max);
	sisi_set_color(col);
	strcpy(operation, "-e1");
	siut_set_directed(FALSE);
	siut_set_multigraph(TRUE);

	if (VERBOSE)
		printf(
				"\n----------------------READ DATA STRUCTURES------------------\n");
	if ((TAS = (t_pointer) heat_creer_hea(MAXHEAP)) == NIL)
		return (ERROR);
	molecule = (molecule_p) daal_create_molecule(TAS, "", 0, 0, NIL, NIL, NIL,NIL, "", 0);

	//PC: Replaced original moin_read with DF's version that creates data structure Graph
	struct GRAPH Graph;
	if (STEREO_MODE)
		err = moin_read(file_name, (molecule_p) "_group",
				operation, molecule, FALSE, &Graph);
	 else
			err = moin_read(file_name, (molecule_p) "_group",
					operation, molecule, FALSE, NIL);
	if (err == ERROR) {
		print_usage();
		return(0);
	}
	//PC END

	SM.molecule = molecule;
	SM.succ = NIL;

	if (OPT_RANDOMIZER > 0) {
		for (i = 0; i < OPT_RANDOMIZER; i++)
			dara_random_integer(OPT_RANDOMIZER);
		dara_random_molecule(molecule);
	}


	strcpy(USERSCAN, SCANTYPE);
	////////////////////////////////////////////////////////
	// PC: BEGIN STEREO
	// PC: Store the scan type for later use
	if (STEREO_MODE) {

		int size = Graph.size;
		int NumbEdges = Graph.NumbEdges;
		signature_p s = NIL;

		if (FAST_MODE) {
			// Switch to DAG (fast) mode
			strcpy(SCANTYPE, "fscancip");
		} else {
			// Switch to CIP mode
			strcpy(SCANTYPE, "cip");
		}
		// Initial diameter
		int init_dim = dim;

		// Search for and classify stereocenters and stereobonds
		dffg_StereoCenters(Graph);
		dffg_StereoBonds(Graph);

		// Flag initial stereo candidates
		int *Parity, *OldParity, *BondParity, *OldBondParity;
//		pcss_InitParities
		pcss_InitParities(&Parity,&BondParity,&OldParity,&OldBondParity,Graph);

		// Start stereo assignment
		t_bool NEWPARITIES_ATOM = FALSE;
		t_bool NEWPARITIES_BOND = FALSE;
		NEWPARITIES_ATOM  = pcss_ParitiesChange(Graph.size,OldParity,Parity);
		if (dim>2) {
			NEWPARITIES_BOND  = pcss_ParitiesChange(Graph.NumbEdges,OldBondParity,BondParity);
		}

		while (NEWPARITIES_ATOM || NEWPARITIES_BOND) { // LOOP UNTIL PARITY IS FIXED
			// We go for a next round if either a stereocenter or a stereobond was updated in this iteration
			if (NEWPARITIES_ATOM) {
				// Set even diameter for computing stereocenters parity
				set_diameter(pcss_DiameterAtomSignature(init_dim));
				// Apply CIP rules to determine stereocenters parity
				pcss_StereocenterParities(Parity,OldParity,Graph,&SM);
				// True if some stereocenters were updated while others were not assigned yet
				NEWPARITIES_ATOM  = pcss_ParitiesChange(Graph.size,OldParity,Parity) && pcss_AtomParitiesChange(Graph.size, Parity,Graph);
				pcss_UpdateAtomParity(Parity,OldParity,Graph);
			}

			// If we go for a next round, we compute the bond parities if not all stereobond candidates have been already assigned
			if ((dim>2)&&(NEWPARITIES_BOND)) {
				// Set odd diameter for computing stereobonds parity
				set_diameter(pcss_DiameterBondSignature(init_dim));
				// Apply CIP rules to determine stereobonds parity
				 pcss_StereobondParities(BondParity,OldBondParity,Graph,&SM);
				// True if some stereobonds were updated while others were not assigned yet
				NEWPARITIES_BOND  = pcss_ParitiesChange(Graph.NumbEdges,OldBondParity,BondParity) && pcss_BondParitiesChange(Graph.NumbEdges, BondParity,Graph);
				pcss_UpdateBondParity(BondParity,OldBondParity,Graph);
			}

		}// Do loop until parity is fixed
		pcss_FreeParities(Parity,BondParity,OldParity,OldBondParity);

		// Return to initial mode
		strcpy(SCANTYPE, USERSCAN);
		set_diameter(init_dim);

		// PC: END STEREO
		////////////////////////////////////////////////////////
	}

	if (VERBOSE)
		printf(
				"\n---------------------WRITE DATA STRUCTURES------------------\n");

	if ((OPT_CLOCK) || (VERBOSE))
		t1 = clock();
	if (dim == MAXINTEGER) {
		strcat(file_name,"."); strcat(file_name,output);
		sisc_canonize(file_name,molecule);
		return(0);
	}
	if ((strncmp(output, "scan", 4) == OK)|| (strncmp(output, "sscan", 5) == OK)
		|| (strncmp(output, "fscan", 5) == OK) || (strncmp(output, "sig", 3) == OK)
		|| (strncmp(output, "ssig", 4) == OK) || (strcmp(output, "cip") == OK))
			write_scan(file_name, &SM, FALSE);

	if ((OPT_CLOCK) || (VERBOSE))
		t2 = clock();
	if ((OPT_CLOCK) || (VERBOSE)) {
		printf("%s vertices %d edges %d cpu time %.4f s\n", file_name,
				molecule->size, dast_number_link(molecule, NIL, NIL),
				(float) (t2 - t1) / (float) 1000000L);
	}
	if (STEREO_MODE)
		dfmm_DeleteGraphMemory(Graph);

}

