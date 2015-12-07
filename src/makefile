# User options

COMPILER = gnu

DEBUG = true

# Program and source code

PROGRAM = transport

HEADERS = header.h

SOURCE = \
main.c \
initialize.c \
utils.c \
io.c \
transport.c \
tally.c \
eigenvalue.c

OBJECTS = $(SOURCE:.c=.o)

# Set flags

CFLAGS = -Wall
LDFLAGS = -lm

ifeq ($(COMPILER),gnu)
  CC = gcc
endif

ifeq ($(DEBUG),true)
  CFLAGS += -g
endif

# Targets to build

$(PROGRAM): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(PROGRAM)
