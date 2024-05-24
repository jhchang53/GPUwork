
GEOM =  Geometry.o TRcoord.o CG.o CGcyl.o CGrpp.o CGtrc.o UVfile.o Box.o Vox.o

CXXFLAGS = -g -fopenmp

# tally
batchtle: batchtle.o $(GEOM) Tally.o Xor128.o
	g++ -o batchtle -g $^ -lm -fopenmp
rantle: rantle.o $(GEOM) Tally.o Xor128.o
	g++ -o rantle -g $^ -lm
tle: tle.o $(GEOM) Tally.o Xor128.o
	g++ -o tle -g $^ -lm

#  tracking
scat: scat.o $(GEOM) Xor128.o
	g++ -o scat -g $^ -lm
collide: collide.o $(GEOM) Xor128.o
	g++ -o collide -g $^ -lm
geomall: geomall.o $(GEOM) Xor128.o
	g++ -o geomall -g $^ -lm
geomcg: geomcg.o $(GEOM) Xor128.o
	g++ -o geomcg -g $^ -lm
#  coordinate translation
tran: tran.o TRcoord.o
	g++ -o tran -g $^ -lm
tranvec: tranvec.o TRcoord.o
	g++ -o tranvec -g $^ -lm

# CG combination tests
cgsera: cgsera.o  CG.o CGcyl.o CGrpp.o CGtrc.o CGset.o Xor128.o
	g++ -o cgsera -g $^ -lm
cgcon: cgcon.o  CG.o CGcyl.o CGrpp.o CGtrc.o CGset.o Xor128.o
	g++ -o cgcon -g $^ -lm
cgset: cgset.o CG.o CGcyl.o CGrpp.o CGtrc.o CGset.o Xor128.o
	g++ -o cgset -g $^ -lm
#  CG test
trc: trc.o CGtrc.o CG.o Xor128.o
	g++ -o trc -g $^ -lm
cyl: cyl.o CGcyl.o CG.o Xor128.o
	g++ -o cyl -g $^ -lm
rpp: rpp.o CGrpp.o CG.o Xor128.o
	g++ -o rpp -g $^ -lm

# UV test

vox: vox.o UVfile.o Box.o Vox.o
	g++ -o vox -g $^ -lm
box: box.o Box.o
	g++ -o box -g $^ -lm
clear:
	rm -f *.o

