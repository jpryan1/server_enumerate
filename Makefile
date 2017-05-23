AFLAGS = 
LFLAGS = -std=c++11 -lpthread 
CFLAGS = -std=c++11 -c -o
FAST_FLAG =

recompile:
	rm build/src/* ;
	make compile ;

configure:
	@cd src/nauty && \
	echo "Configuring nauty..." && \
	./configure &> ../../docs/nauty/config.log  && \
	echo "Making nauty library..." && \
	make &>../../docs/nauty/make.log nauty.a && \
	echo "Cleaning up..." && \
	cd ../.. && \
	cp src/nauty/nauty.a build/lib/ && \
	rm -rf src/nauty/ 


build/src/main.o: src/main.cpp src/Configuration.h src/Bank.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/main.o src/main.cpp 

build/src/process.o: src/process.cpp src/Configuration.h src/Bank.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/process.o src/process.cpp 


build/src/Config.o: src/Configuration.cpp src/Configuration.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/Config.o src/Configuration.cpp

build/src/canonize.o: src/canonize.cpp src/Configuration.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/canonize.o src/canonize.cpp

build/src/project.o: src/project.cpp src/Configuration.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/project.o src/project.cpp

build/src/walk.o: src/walk.cpp src/Configuration.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/walk.o src/walk.cpp

build/src/dimension.o: src/dimension.cpp src/Configuration.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/dimension.o src/dimension.cpp

build/src/Bank.o: src/Bank.cpp src/Configuration.h src/Bank.h
	g++ $(FAST_FLAG) $(CFLAGS) build/src/Bank.o src/Bank.cpp

build/src/Timer.o: src/Timer.cpp src/Timer.h 
	g++ $(FAST_FLAG) $(CFLAGS) build/src/Timer.o src/Timer.cpp

build/src/naugroup.o: src/naugroup.c
	gcc $(FAST_FLAG) -c -o build/src/naugroup.o src/naugroup.c




compile: build/src/main.o build/src/Config.o build/src/walk.o build/src/dimension.o build/src/canonize.o build/src/project.o build/src/Bank.o build/src/Timer.o build/src/naugroup.o
	g++ $(LFLAGS) $(AFLAGS) $(FAST_FLAG) -o build/enumerate_clusters build/src/main.o build/src/Config.o build/src/walk.o build/src/dimension.o build/src/canonize.o build/src/project.o build/src/Bank.o build/src/Timer.o build/src/naugroup.o build/lib/nauty.a


compile_process: build/src/process.o build/src/Config.o build/src/walk.o build/src/dimension.o build/src/canonize.o build/src/project.o build/src/Bank.o build/src/Timer.o
	g++ $(LFLAGS) $(AFLAGS) -O2 -o build/process build/src/process.o build/src/Config.o build/src/walk.o build/src/dimension.o build/src/canonize.o build/src/project.o build/src/Bank.o build/src/Timer.o build/lib/nauty.a

run:
	./build/enumerate_clusters r

process:
	./build/process

debug:
	./build/enumerate_clusters d

fast_compile:
	make compile FAST_FLAG=-O2

fast_compile_process:
	make compile_process FAST_FLAG=-O2


fast_recompile:
	make recompile FAST_FLAG=-O2
