FC = gfortran -ffixed-line-length-none
FFLAGS = -O
CC = gcc-11 -fopenmp
CFLAGS = ${FFLAGS} 

SACHOME =/Users/qb/codes/sac

SUBS = sacio.o fft.o Complex.o radiats.o radpmt.o butstsub.o grid3d.o nrutil.o jacobi.o eigsrt.o

SRC_INV = mpsmp 

all: $(SRC_INV) fault 
	rm *.o

$(SRC_INV): %:%.o $(SUBS)
	$(LINK.f) -o $@ $@.o $(SUBS) $(FGCCLIB) -fopenmp -L${SACHOME}/lib -lsac -lsacio -lm

fault: fault.c
	$(LINK.c) -o fault fault.c -lm

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

fmap: fmap.o
	$(LINK.f) -o $@ $@.o $(FGCCLIB)
clean:
	rm -f *.o
