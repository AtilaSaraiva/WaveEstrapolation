#!/bin/zsh
.PHONY: clean run

run: exe
	./exe clip=3000 <campovel.rsf >teste.rsf
	sfgrey <teste.rsf | sfpen
	#nimage perc=98 n1=200 n2=200 <snap.ad

exe: clip.f95
	g++ clip.cpp -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsf++ -lrsf -lm -ltirpc -lfftw3f -o exe

clean:
	rm exe *.mod
	rm snap.ad
