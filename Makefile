.PHONY: clean run

run: exe
	./exe
	nimage perc=98 n1=200 n2=200 <snap.ad

exe: wave.f95
	gfortran -o exe wave.f95 -fbounds-check

clean:
	rm exe *.mod
	rm snap.ad
