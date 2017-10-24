GCC=mpicc
EXEC=filtru
NAME=331AA_Barbu_Florina_tema3.zip

all: build

build: main.c
	$(GCC) $^ -o $(EXEC)

clean:
	rm $(EXEC) 

run:
	mpirun -np 12 ./filtru topologie.in imagini.in statistica.out

pack:
	zip $(NAME) main.c README Makefile 

