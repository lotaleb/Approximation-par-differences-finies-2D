CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pthread
LDFLAGS = -lm

SRC_MAIN = main.c
SRC_FUNCTIONS = Functions.c
OBJ_MAIN = $(SRC_MAIN:.c=.o)
OBJ_FUNCTIONS = $(SRC_FUNCTIONS:.c=.o)
EXECUTABLE = my_program

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_MAIN) $(OBJ_FUNCTIONS)
	$(CC) $(LDFLAGS) $^ -o $@

$(OBJ_MAIN): $(SRC_MAIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_FUNCTIONS): $(SRC_FUNCTIONS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_MAIN) $(OBJ_FUNCTIONS) $(EXECUTABLE)
