CXXFLAGS := -std=c++1z -Wall -Wno-deprecated-declarations -Wno-unused-command-line-argument
LDFLAGS := /usr/local/lib
LDLIBS := -lglfw -framework Cocoa -framework OpenGL -framework IOKit
INCLUDE := /usr/local/include
SRC := $(wildcard *.cpp)
OBJECTS  := $(SRC:%.cpp=%.o)
DEPENDS := $(SRC:%.cpp=%.d)
TARGET := hw1

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -I$(INCLUDE) -L$(LDFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -L$(LDFLAGS) $(LDLIBS) -o $(TARGET) $^

-include $(DEPENDS)

.PHONY: clean

clean:
	-@rm -rvf $(OBJECTS)
	-@rm -rvf $(DEPENDS)
	-@rm -rvf $(TARGET)
