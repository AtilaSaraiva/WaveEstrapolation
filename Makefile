.PHONY: clean run

run: exe
	./exe vel=campo.rsf wav=wavelet.rsf >snaps.rsf
	#nimage perc=98 n1=200 n2=200 <snap.ad

exe: wave.f95
	gfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -fopenacc -lrsf -lm -ltirpc -lfftw3f -lfftw3 -o exe

clean:
	rm exe *.mod
	rm snap.ad
