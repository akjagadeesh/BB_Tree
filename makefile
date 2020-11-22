all: main

main: BubbleBucket.h BB_Tree.h
	g++ -std=c++14 main.cpp -o main

clean:
	rm main

run: main
	./main