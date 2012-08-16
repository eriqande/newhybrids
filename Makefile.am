ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = newhybs newhybsng

newhybs_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) \
	-I$(top_srcdir)/eca-shared/ecalibs \
	-I$(top_srcdir)/src  \
	-I$(top_srcdir)/eca-shared/ranlib/src  \
	-I$(top_srcdir)/eca-shared/gfmcmc


newhybsng_CFLAGS = -I$(top_srcdir)/eca-shared/ecalibs \
	-I$(top_srcdir)/src  \
	-I$(top_srcdir)/eca-shared/ranlib/src  \
	-I$(top_srcdir)/eca-shared/gfmcmc \
	-DCOMPILE_NEW_HYB_WITH_NO_GUI


newhybs_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm

newhybsng_LDFLAGS = -lm

newhybs_SOURCES = $(top_srcdir)/src/GLUT_for_NewHybrids.c  \
	$(top_srcdir)/src/DataInput.c  \
	$(top_srcdir)/src/DataOutput.c  \
	$(top_srcdir)/src/NewHybrids.c  \
	$(top_srcdir)/src/RunWOG_func.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_MemAlloc.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_Opt3.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_ReadData.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_utilities.c  \
	$(top_srcdir)/eca-shared/ecalibs/MathStatRand.c  \
	$(top_srcdir)/eca-shared/ecalibs/MCTypesEtc.c  \
	$(top_srcdir)/eca-shared/gfmcmc/GFMCMC.c \
	$(top_srcdir)/eca-shared/gfmcmc/GFMCMC_DrawingUtilities.c  \
	$(top_srcdir)/eca-shared/ranlib/linpack/linpack.c  \
	$(top_srcdir)/eca-shared/ranlib/src/com.c  \
	$(top_srcdir)/eca-shared/ranlib/src/ranlib.c 


newhybsng_SOURCES = $(top_srcdir)/src/GLUT_for_NewHybrids.c  \
	$(top_srcdir)/src/DataInput.c  \
	$(top_srcdir)/src/DataOutput.c  \
	$(top_srcdir)/src/NewHybrids.c  \
	$(top_srcdir)/src/RunWOG_func.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_MemAlloc.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_Opt3.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_ReadData.c  \
	$(top_srcdir)/eca-shared/ecalibs/ECA_utilities.c  \
	$(top_srcdir)/eca-shared/ecalibs/MathStatRand.c  \
	$(top_srcdir)/eca-shared/ecalibs/MCTypesEtc.c  \
	$(top_srcdir)/eca-shared/ranlib/linpack/linpack.c  \
	$(top_srcdir)/eca-shared/ranlib/src/com.c  \
	$(top_srcdir)/eca-shared/ranlib/src/ranlib.c 

