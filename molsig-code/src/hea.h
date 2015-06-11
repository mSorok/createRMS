# define TAS_TAIL_BLOC_DEFAUT                  4088
                 /* 4K - 8octets pris par le systeme */

# define TAS_WORD_ALIGN(T)                   ((T) &= 0xfffffffe)
           /* alignement sur le mot de 16 bit */

typedef  unsigned long                        t_taille;


typedef  struct CEL {
                      struct CEL   *suc ;
                                              } t_cellule ,*p_cellule;
                        /* les cellule libre std sont chainees et
                           gerees comme une pile*/

typedef  struct     {
                      p_cellule     pile;
                      t_taille      taille;
                             /* taille des cellules */
                                              } t_std_cel ,*p_std_cel;

typedef  struct BLOC {
                      struct BLOC   *suc ;
                              /* les bloc sont chainees et gerees comme
                                 une pile */
                      t_pointeur    deblib;
                      t_pointeur    finlib;
                                              } t_bloc ,*p_bloc;


typedef  struct      {
                      p_bloc         pile;
                      t_taille       taille_bloc;
                      t_taille       nb_std;
                      t_std_cel      *std;
                                              } t_hea ,*p_hea;




