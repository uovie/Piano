CC      := g++
INCS    := -I ../include
LIBS    := -L ../lib

TARGETS := Piano
OBJ     := main.o process.o model.o thermostat.o ld.o at.o nhc.o pimd.o
DEPS    := 

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c $< -o $@ $(INCS) $(LIBS)

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@ $(INCS) $(LIBS)
