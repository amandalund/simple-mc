# User options

COMPILER = gnu

# Program and source code

PROGRAM = simple-transport

HEADERS = simple_transport_header.h

SOURCE = \
main.c \
initialize.c \
utils.c \
io.c \
transport.c \
eigenvalue.c

OBJECTS = $(SOURCE:.c=.o)

# Set flags

CFLAGS = -Wall -g
LDFLAGS = -lm

ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# Targets to build

$(PROGRAM): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(PROGRAM)
