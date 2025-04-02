all: main.cpp
	g++ -pthread main_other.cpp
	./a.out


clean: 
	rm a.out 
