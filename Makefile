.PHONY: clean run

run: exe
	#valgrind --tool=callgrind ./exe vel=campo.rsf wav=wavelet.rsf >snaps.rsf
	/usr/bin/time ./exe vel=campo.rsf wav=wavelet.rsf snaps=snaps.rsf >data.rsf sx=1 sz=0 gxbeg=0.5 gzbeg=0 jgx=0.01 nr=100
	sfgrey gainpanel=a <snaps.rsf | sfpen
	#nimage perc=98 n1=200 n2=200 <snap.ad

exe: wave.f95
	gfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -O3 -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe  -fopenmp -fopt-info-optimized-omp -pthread -fbounds-check
	#pgfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -fast -Minfo=opt -ta=multicore -Minfo=accel -Mfree -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe

clean:
	rm exe *.mod
