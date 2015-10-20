

gcc src/GLUT_for_NewHybrids.c  \
	src/DataInput.c  \
	src/DataOutput.c  \
	src/NewHybrids.c  \
	src/RunWOG_func.c \
	eca-shared/src/ECA_Opt3.c \
	eca-shared/src/MathStatRand.c \
	eca-shared/src/MCTypesEtc.c \
	eca-shared/src/com.c \
	eca-shared/src/linpack.c \
	eca-shared/src/ranlib.c \
	eca-shared/src/ECA_ReadData.c \
	eca-shared/src/ECA_MemAlloc.c \
	eca-shared/src/ECA_utilities.c \
	-I eca-shared/include \
	-I src \
	-DCOMPILE_NEW_HYB_WITH_NO_GUI