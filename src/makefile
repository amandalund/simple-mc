OBJECTS = main.o initialize.o utils.o io.o transport.o eigenvalue.o
HEADERS = simple_transport_header.h
CFLAGS = -Wall -g
CC = gcc

default: simple_transport

simple_transport: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ -lm

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) simple_transport
