CC      := g++
INCS    := -I ../include \
           -I ../include/thermostat
LIBS    := -L ../lib \
           -L ../lib/thermostat

TARGETS := Piano
OBJ     := test_nhc.o process_print.o process_read.o nose_hoover_chain.o
DEPS    := 

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c $< -o $@ $(INCS) $(LIBS)

$(TARGETS): $(OBJ)
	$(CC)  $^ -o $@ $(INCS) $(LIBS)