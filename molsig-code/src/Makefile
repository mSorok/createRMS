##########################################################
#
# makefile for module (load and run) Compilation
#
#  Jean-Loup Faulon
#  iSSB U. Evry Genopole
#  09/2009
#
##########################################################

GDB_SUPPORT = 
ifdef GDB_SUPPORT
CFLAGS = -w -g -m32 -I.
else
CFLAGS = -w -O3 -m32 -I.
endif

LINKFLAGS = -m32

cc = /usr/bin/gcc

#--------------------------------------------------------------------
#  Original scan headers and objects, except for files that
#  are modified to run with the new executables.
#--------------------------------------------------------------------

include_h =  c.h  eps.h  general.h  

# data modules #######################################################
da_data_o = allocation_da.o random_da.o structure_da.o \
structure_nb.o structure_tp.o structure_ut.o \
twin_da.o


he_hea_o = hea.o structure_he.o
he_hea_h = hea.h  

si_signature_o = allocation_si.o \
signature.o cansig.o candag.o utility_si.o signature_mo.o
si_signature_h = signature.h

# io modules #########################################################
mo_mol_o       = in_mo.o out_mo.o utility_mo.o 

# stereo modules ########################################################
df_geometry_o = DFConst.o DFManageMemory.o DFFindGeometry.o
pc_stereo_o = PCStereoSig.o PCcip.o

# main module ########################################################
x_main_o = error.o vector.o main.o 


#####################################################################
#
#  EXECUTABLE:
#     sscan  clean         
#####################################################################

all:    sscan 

sscan:    \
	$(include_h) $(he_hea_h) $(si_signature_h) \
	$(da_data_o) $(he_hea_o) $(si_signature_o) \
	$(mo_mol_o) $(df_geometry_o) $(pc_stereo_o) \
	$(x_main_o) $(invar_o)
	$(cc) -o ../bin/sscan $(LINKFLAGS) \
        $(da_data_o) $(he_hea_o) $(si_signature_o) $(mo_mol_o)  \
        $(df_geometry_o) $(pc_stereo_o) \
	$(x_main_o) -lm 
clean:
	rm -f *.o ../bin/sscan; 
#####################################################################
#
#  objects
#
#####################################################################


#da_data
allocation_da.o : allocation_da.c $(include_h)
	$(cc) $(CFLAGS) -c allocation_da.c
random_da.o : random_da.c $(include_h)
	$(cc) $(CFLAGS) -c random_da.c
structure_da.o : structure_da.c $(include_h)
	$(cc) $(CFLAGS) -c structure_da.c
structure_nb.o : structure_nb.c $(include_h)
	$(cc) $(CFLAGS) -c structure_nb.c
structure_tp.o : structure_tp.c $(include_h)
	$(cc) $(CFLAGS) -c structure_tp.c
structure_ut.o : structure_ut.c $(include_h)
	$(cc) $(CFLAGS) -c structure_ut.c
twin_da.o : twin_da.c $(include_h)
	$(cc) $(CFLAGS) -c twin_da.c

#he_hea
hea.o : hea.c $(include_h) $(he_hea_h)
	$(cc) $(CFLAGS) -c hea.c
structure_he.o : structure_he.c $(include_h) $(he_hea_h)
	$(cc) $(CFLAGS) -c structure_he.c

#si_signature
allocation_si.o : allocation_si.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c allocation_si.c
equation.o : equation.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c equation.c
in_out_si.o : in_out_si.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c in_out_si.c
signature.o : signature.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c signature.c
utility_si.o : utility_si.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c utility_si.c
cansig.o : cansig.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c cansig.c
candag.o : candag.c $(include_h) $(si_signature_h)
	$(cc) $(CFLAGS) -c candag.c
signature_mo.o : signature_mo.c $(include_h)
	$(cc) $(CFLAGS) -c signature_mo.c

# main
main.o : main.c $(include_h)
	$(cc) $(CFLAGS) -c main.c
error.o : error.c $(include_h)
	$(cc) $(CFLAGS) -c error.c
vector.o : vector.c $(include_h)
	$(cc) $(CFLAGS) -c vector.c

# io
in_mo.o : in_mo.c $(include_h)
	$(cc) $(CFLAGS) -c in_mo.c
out_mo.o : out_mo.c $(include_h)
	$(cc) $(CFLAGS) -c out_mo.c
utility_mo.o : utility_mo.c $(include_h)
	$(cc) $(CFLAGS) -c utility_mo.c

# stereo
DFConst.o : DFConst.c $(include_h)
	$(cc) $(CFLAGS) -c DFConst.c
DFManageMemory.o : DFManageMemory.c $(include_h)
	$(cc) $(CFLAGS) -c DFManageMemory.c
DFFindGeometry.o : DFFindGeometry.c $(include_h)
	$(cc) $(CFLAGS) -c DFFindGeometry.c
PCStereoSig.o : PCStereoSig.c $(include_h)
	$(cc) $(CFLAGS) -c PCStereoSig.c
PCcip.o : PCcip.c $(include_h)
	$(cc) $(CFLAGS) -c PCcip.c



