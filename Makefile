CXX = g++
CXXFLAGS = -O2
# CXXFLAGS += -fopenmp

# Tenta usar pkg-config; se falhar, usa caminhos manuais
YAML_CFLAGS := $(shell pkg-config --cflags yaml-cpp 2>/dev/null || echo "-I/usr/local/include")
YAML_LIBS   := $(shell pkg-config --libs yaml-cpp 2>/dev/null || echo "-L/usr/local/lib -lyaml-cpp")

CXXFLAGS += $(YAML_CFLAGS)
LDFLAGS  += $(YAML_LIBS)

TARGET = main.x
OBJS = CG.o thomas.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

CG.o: CG.cpp CG.h algebra.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

thomas.o: thomas.cpp thomas.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

main.o: main.cpp manufaturadas.cpp CG.h matrizes.h algebra.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Regra genérica para outros .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET) *.csv *.mp4 *.png

.PHONY: all clean
