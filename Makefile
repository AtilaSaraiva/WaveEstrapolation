.PHONY: clean run

run: exe
	./exe vel=mod1.rsf <pulso2.rsf >snaps.rsf
	#nimage perc=98 n1=200 n2=200 <snap.ad

exe: wave.f95
	gfortran wave.f95 -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsff90 -lrsf -lm -ltirpc -lfftw3f -o exe

clean:
	rm exe *.mod
	rm snap.ad
