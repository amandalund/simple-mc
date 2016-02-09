# User options

COMPILER = gnu
OPTIMIZE = yes
DEBUG    = no

# Program and source code

PROGRAM = transport

HEADERS = \
header.h

SOURCE = \
main.c \
initialize.c \
prng.c \
utils.c \
io.c \
transport.c \
tally.c \
eigenvalue.c

OBJECTS = $(SOURCE:.c=.o)

# Set flags

CFLAGS = -Wall
LDFLAGS = -lm

ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

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
