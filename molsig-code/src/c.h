/****************************************************************************                                                    
  General Data structure for C language (in french and english).
        
     Jean-loup Faulon May 91 Albuquerque N.M.
*****************************************************************************/

#ifndef C_H
#define C_H


/***************** MACRO instructions **************************************/
  /* majuscule */
# define    MAJ(X)      ( 'a' <= (X) && (X) <= 'z' ? (X) -'a'+'A' : (X))
        /* minuscule */
# define    MINUS(X)      ( 'A' <= (X) && (X) <= 'Z' ? (X) -'A'+'a' : (X))
            
  /* Functions MIN and MAX */
# define    MAX(X,Y)    ( (X) < (Y) ? (Y) : (X) )
# define    MIN(X,Y)    ( (X) < (Y) ? (X) : (Y) )
# define    ABS( n )    ( (n) > 0 ? (n) : -(n) )


/****************** CONSTANTS AND TYPES     ********************************/
  
  /* Unsigned Integer Types  */
typedef unsigned char   us_char;
typedef unsigned short  us_short;
typedef unsigned long   us_int;
typedef unsigned char	t_digit;

      /* The bool type must be a long since sometimes
                           pointer are assinged to it = crapy programing */
typedef long			t_bool;     /*  C  or  PASCAL !! PSU !!*/
typedef unsigned long   t_bool_pf;  /* only FORTRAN  */
           
                                              /**********/
# define VRAI_PF             0xffffffff   /* true fortran or pascal */
# define TRUE_PF             0xffffffff
# define FAUX_PF             (long)0           /* false fortran or pascal */
# define FALSE_PF            (long)0
# define VRAI                1    /* no null !!!! */
# define TRUE                1
# define FAUX                0    /* null         */
# define FALSE         0

/*********************************************/
typedef  char                *t_pointeur;
typedef  char                *t_pointer;
                                              /**********/
# define NIL                  0L
/******************************/
typedef  long                                   t_err;
                                              /**********/
# define OK                   NIL    /* can be also a pointer     */
                                     /* be careful when non zero  */
/*********************************************/
typedef  long                                   t_compteur;
typedef  long                                   t_count;
                                              /**********/
# define ERREUR               -1L     /* ERROR if negativ */
# define ERROR                -1L     
/*********************************************/
typedef  char                          *t_chaine,*p_chaine;
typedef  char                          *t_string,*p_string;
                                              /********************/
# define FIN                  '\0'
# define END                  '\0'
# define MAXCHAINE            256
# define MAXSTRING            256
# define MAXINTEGER           1000000000
typedef  char                          t_bufchaine[MAXCHAINE];
typedef  char                          t_bufstring[MAXCHAINE];
                                              /**********/
/*********************************************/
# define MAXNOM                256
# define MAXNAME               256
               /* 32 car on apollo plus FIN */

typedef   char                            t_nom[MAXNOM],*p_nom;
typedef   char                            t_name[MAXNOM],*p_name;
          
/*********************************************/
typedef     struct{
                unsigned char       jour[2] ;
                unsigned char       mois[2] ;
                unsigned char       an[2] ;
                                                }t_date,*p_date ;
/*********************************************/
typedef  t_err                                (*f_traitement)();
                                              /**********/



/*********************************************/
typedef long  t_integer;  
typedef  double  t_real;
typedef  struct  COORD {
    t_real  x,y,z;
    } t_coord, *p_coord;


# define BOUCLER                0
# define FINBOUCLE             -1
# define ERR_TRAITEMENT        -2
/*********************************************/

# define  SERVICE          extern
# define  UPWARD           extern
# define  FORWARD          extern
# define  LOCAL            static
# define  STATIC           static
# define  DEFINE       
# define  EXTERN           extern
# define  LOOP             for(;;)
# define  ENDLOOP          break

# define  FAULT            -1               /* ERROR call UNIX system  */
# define  OPEN        02000 
# define  CREAT       00400  
          /* flag of open */
# define  READ        0   
# define  WRITE       1
# define  UPDATE      2
          /* mode of open */
# define  FROMBOF     0
# define  FROMCUR     1
# define  FROMEOF     2
          /* for lseek */
# define  EXIST       0
          /* for acces function */


/*********************************************/
/* FOR DEBUG  AND VERBOSE                    */
EXTERN    t_bool      DEBUG;         /* declared in error.c */
EXTERN    t_bool      VERBOSE;            /* global to the progam must be declared in error.c */

#endif
