CC      := g++
INCS    := -I ../include \
           -I ../include/thermostat
LIBS    := -L ../lib \
           -L ../lib/thermostat

TARGETS := Piano
OBJ     := test_nhc.o process.o nhc.o
DEPS    := 

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c $< -o $@ $(INCS) $(LIBS)

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@ $(INCS) $(LIBS)
