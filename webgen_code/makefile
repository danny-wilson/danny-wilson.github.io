webgen : main.o
	/usr/bin/g++ -w -o3 -o webgen main.o -L/usr/local/lib  -lstdc++ -static -lm

main.o : main.cpp
	/usr/bin/g++ -w -o3 -c -o main.o main.cpp

clean :
	rm main.o
