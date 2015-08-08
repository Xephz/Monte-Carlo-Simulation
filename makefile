CC = gcc
OBJ = GEM.o coreMCfunc.o particleIO.o analysis.o
OUT = GEM
FLAGS = -Wall -lm

build: $(OBJ)
	$(CC) $(FLAGS) -o $(OUT) $(OBJ)

%.o: %.c
	$(CC) $(FLAGS) -c $<

clean:
	rm *.o
	rm $(OUT)

clearData:
	rm *.xyz
	rm *.GRAPH
	rm *.DIST

rebuild:
	clean build