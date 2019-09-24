.PHONY: clean run

run: exe
	#valgrind --tool=callgrind ./exe vel=campo.rsf wav=wavelet.rsf >snaps.rsf
	./exe vel=campo.rsf wav=wavelet.rsf >snaps.rsf
	#nimage perc=98 n1=200 n2=200 <snap.ad

exe: wave.f95
	gfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -fopenacc -fopt-info-optimized-omp -pthread  -O3 -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe
	#pgfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -fast -Minfo=opt -ta=multicore -Minfo=accel -Mfree -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe

clean:
	rm exe *.mod
