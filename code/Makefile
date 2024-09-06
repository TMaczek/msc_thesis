# Define compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -ggdb3

# Source files
SRCS = main.cpp System.cpp Particle.cpp Cell.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = main_exec

# Default rule
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean rule
clean:
	$(RM) $(OBJS) $(TARGET)
