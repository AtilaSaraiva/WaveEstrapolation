FC=gfortran
FFLAGS=-c
LDFLAGS= -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -O3 -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe  -fopenmp -fopt-info-optimized-omp -pthread -fbounds-check
SOURCES = $(wildcard *.f95)
OBJECTS = $(subst .f95,.o,$(SOURCES))
EXECUTABLE=exe

.PHONY: clean help exec

all: $(SOURCES) $(EXECUTABLE) exec

$(EXECUTABLE): $(SOURCES)
	$(FC) $(SOURCES) $(LDFLAGS) -o $@

%.o : %.f95
	$(FC) $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

help:
	@echo "Valid targets:"
	@echo "  main.exe"
	@echo "  main.o"
	@echo "  sub1.o"
	@echo "  sub2.o"
	@echo "  clean: removes .o and .exe files"


exec:
	/usr/bin/time ./$(EXECUTABLE) vel=campo.rsf wav=wavelet.rsf snaps=snaps.rsf >data.rsf sx=1 sz=0 gxbeg=0.5 gzbeg=0 jgx=0.01 nr=100
