#*************************************************************************
#*                                                                       *
#*     Makefile for IviaCLR (missing value imputation via nonsmooth      * 
#*     clusterwise linear regression, last modified 01.03.2017).         *
#*                                                                       *
#*************************************************************************
# 
FF = gfortran -O2 
#-fbounds-check -Wall
OBJ = parameters.o initiviaclr.o clrfunctionmod.o clrsubpro.o lmbmsubpro.o lmbm.o imput.o iviaclr.o

all:	iviaclr

parameters.o: parameters.f03
	$(FF) -c parameters.f03

initiviaclr.o: r_precision.mod initiviaclr.f03
	$(FF) -c initiviaclr.f03

clrfunctionmod.o: r_precision.mod param.mod initclr.mod clrfunctionmod.f03
	$(FF) -c clrfunctionmod.f03

clrsubpro.o: r_precision.mod param.mod initclr.mod clrsubpro.f03
	$(FF) -c clrsubpro.f03

lmbmsubpro.o: r_precision.mod param.mod lmbmsubpro.f03
	$(FF) -c lmbmsubpro.f03

lmbm.o: r_precision.mod param.mod initclr.mod initlmbm.mod exe_time.mod clrobjfun.mod lmbmsubpro.mod lmbm.f03
	$(FF) -c lmbm.f03

imput.o: r_precision.mod param.mod initimp.mod initlmbm.mod clrobjfun.mod lmbmsubpro.mod lmbm_mod.mod imput.f03
	$(FF) -c imput.f03

iviaclr.o: r_precision.mod param.mod exe_time.mod initimp.mod initclr.mod initlmbm.mod clrsubpro.mod lmbmsubpro.mod lmbm_mod.mod firstimput.mod prediction.mod iviaclr.f03
	$(FF) -c iviaclr.f03

iviaclr: $(OBJ)
	$(FF) -o iviaclr $(OBJ)

clean:	
	rm iviaclr $(OBJ) *.mod
	echo Clean done	