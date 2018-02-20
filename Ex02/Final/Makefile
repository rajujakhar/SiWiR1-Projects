CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wshadow -Werror -O3 -pipe -fopenmp -DNDEBUG


TARGET = rbgs
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

$(TARGET).o: $(TARGET).cpp Makefile
	$(CXX) -c $(CXXFLAGS) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)
