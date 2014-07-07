ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = eca-shared

bin_PROGRAMS = newhybs newhybsng

newhybs_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
newhybsng_CFLAGS = -Ieca-shared/include	-DCOMPILE_NEW_HYB_WITH_NO_GUI


newhybs_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
newhybsng_LDFLAGS = -lm

newhybs_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a eca-shared/libecaopt3.a
newhybsng_LDADD = eca-shared/libecashared.a  eca-shared/libecaopt3.a


newhybs_SOURCES = src/GLUT_for_NewHybrids.c  \
	src/DataInput.c  \
	src/DataOutput.c  \
	src/NewHybrids.c  \
	src/RunWOG_func.c


newhybsng_SOURCES = src/GLUT_for_NewHybrids.c  \
	src/DataInput.c  \
	src/DataOutput.c  \
	src/NewHybrids.c  \
	src/RunWOG_func.c

