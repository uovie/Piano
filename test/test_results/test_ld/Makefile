CC      := g++
INCS    := -I ../include
LIBS    := -L ../lib

TARGETS := Piano
OBJ     := test_ld.o process.o model.o ld.o
DEPS    := 

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c $< -o $@ $(INCS) $(LIBS)

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@ $(INCS) $(LIBS)
