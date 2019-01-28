all:
	g++ kmeans.cpp
	./a.out < in.txt

clean:
	rm *.out

