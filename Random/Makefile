CXXFLAGS = -g

isotro: isotro.o Xor128.o
	g++ -o isotro -g $^
ranxor1: ranxor1.o
	g++ -o ranxor -g $^

rand1: rand1.o
	g++ -o rand1 -g $^
rand2: rand2.o
	g++ -o rand2 -g $^
rand3: rand3.o
	g++ -o rand3 -g $^
rand4: rand4.o
	g++ -o rand4 -g $^

rand: random.o
	g++ -o rand -g $^
clear:
	rm -f *.o
