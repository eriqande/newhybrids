


gcc -o simdata_nh src/SimNewHybData.c \
  src/DataInput.c \
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
	-I src