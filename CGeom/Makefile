
GEOM =  Geometry.o TRcoord.o CG.o CGcyl.o CGrpp.o CGtrc.o UVfile.o Box.o Vox.o

CXXFLAGS = -g -std=c++11 -fopenmp

# howfar: howfar.o CG.o CGcyl.o CGrpp.o CGtrc.o ../Random/Xor128.o
#	g++ -o howfar -g $^ -lm
batchtle: batchtle.o $(GEOM) Scatter.o Tally.o Xor128.o
	g++ -o batchtle -g $^ -lm -fopenmp
rantle: rantle.o $(GEOM) Scatter.o Tally.o Xor128.o
	g++ -o rantle -g $^ -lm
tle: tle.o $(GEOM) Scatter.o Tally.o Xor128.o
	g++ -o tle -g $^ -lm
scat: scat.o $(GEOM) Scatter.o Xor128.o
	g++ -o scat -g $^ -lm
collide: collide.o $(GEOM) Xor128.o
	g++ -o collide -g $^ -lm
geomall: geomall.o $(GEOM)
	g++ -o geomall -g $^ -lm
geomcg: geomcg.o $(GEOM)
	g++ -o geomcg -g $^ -lm
#  coordinate translation
tran: tran.o TRcoord.o
	g++ -o tran -g $^ -lm
tranvec: tranvec.o TRcoord.o
	g++ -o tranvec -g $^ -lm
#  CG test

cgsera: cgsera.o  CG.o CGcyl.o CGrpp.o CGtrc.o  CGset.o  Xor128.o
	g++ -o cgsera -g $^ -lm
cgcon: cgcon.o  CG.o CGcyl.o CGrpp.o CGtrc.o  CGset.o  Xor128.o
	g++ -o cgcon -g $^ -lm
cgset: cgset.o CG.o CGcyl.o CGrpp.o CGtrc.o CGset.o Xor128.o
	g++ -o cgset -g $^ -lm
trc: trc.o CGtrc.o CG.o Xor128.o
	g++ -o trc -g $^ -lm
cyl: cyl.o CGcyl.o CG.o Xor128.o
	g++ -o cyl -g $^ -lm
rpp: rpp.o CGrpp.o CG.o Xor128.o
	g++ -o rpp -g $^ -lm
# UV tests
uv3d: uv3d.o UVfile.o Vox.o
	g++ -o uv3d -g $^ -lm
box: box.o  Box.o
	g++ -o box -g $^ -lm
vox: vox.o UVfile.o Box.o Vox.o
	$(CXX) -o vox -g $^ -lm
#  scattering
compton: compton.o  Scatter.o Xor128.o
	g++ -o compton -g $^ -lm
inel: inel.o Scatter.o Xor128.o
	g++ -o inel -g $^ -lm
elas: elas.o Scatter.o Xor128.o
	g++ -o elas -g $^ -lm
# random number
random: random.o Xor128.o
	g++ -o random -g $^
clear:
	rm -f *.o
