#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Check_XS
CALC    = Calc_XS
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)
XEMCDIR = ./SRC

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include -I$(XEMCDIR)
RTLIBS := -L${ROOTSYS}/lib -lCore -lgfortran # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore ${x2}/mylib/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore /usr/lib64/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore  /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
all: $(PROGRAM) ${CALC}

GFORTRAN := gfortran
#GFORTRAN := f77
#FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

SRC       = ./SRC/

my_objs = $(SRC)F1F209.o \
  $(SRC)hyd_elastic.o $(SRC)deut_elastic.o $(SRC)trit_elastic.o \
  $(SRC)he3_elastic.o $(SRC)he4_elastic.o  $(SRC)c12_elastic.o \
  $(SRC)smott.o $(SRC)interp.o $(SRC)formfact.o


%.o: %.f
	$(GFORTRAN) $(FFLAGS) -c $< -o $@

$(PROGRAM): $(PROGRAM).o ${my_objs}
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) ${my_objs}

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}


$(CALC): $(CALC).o ${my_objs}
	g++ -o $(CALC) $(CALC).o $(CFLAGS) $(LIBS) $(RTLIBS) ${my_objs}

$(CALC).o: $(CALC).C
	g++ $(CFLAGS) -c -o $(CALC).o $(CALC).C ${INCLUDE}

clean:
	rm -f *.o ./SRC/*.o	
