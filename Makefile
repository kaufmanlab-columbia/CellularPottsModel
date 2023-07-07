CC=g++
CFLAGS=-I.
DEPS = potts_analysis_.h potts_energy_.h potts_flip_.h potts_print_.h potts_spawn_.h
OBJ = potts.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -lm

potts: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

.PHONY: clean all

clean:
	rm -f *.o *~ core *~ 
  
all:
	make potts
