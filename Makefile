CXXFLAGS := -std=c++11 -Wall -Wno-deprecated-declarations -Wno-unused-command-line-argument
LDFLAGS := /usr/local/lib
LDLIBS := -lglfw -framework Cocoa -framework OpenGL -framework IOKit
INCLUDE := /usr/local/include
SRC := $(wildcard *.cpp)
OBJECTS  := $(SRC:%.cpp=%.o)
TARGET := hw1

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -L$(LDFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -L$(LDFLAGS) $(LDLIBS) -o $(TARGET) $^

.PHONY: clean

clean:
	-@rm -rvf $(OBJECTS)
	-@rm -rvf $(TARGET)
