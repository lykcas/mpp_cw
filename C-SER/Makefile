MF=	Makefile

CC=	gcc
CFLAGS=	-O3

LFLAGS= $(CFLAGS)

EXE=	percolate_2

INC= \
	percolate.h

SRC= \
	percolate.c \
	uni.c \

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(OBJ) core
