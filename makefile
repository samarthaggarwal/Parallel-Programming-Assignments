all:
	g++ kmeans.cpp
	./a.out < in1.txt

clean:
	rm *.out

randomInput:
	g++ randomInput.cpp
	./a.out > in1.txt

