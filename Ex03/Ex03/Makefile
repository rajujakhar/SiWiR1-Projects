CXX = mpic++
#CXXFLAGS = -std=c++11 -Wall -pg -O3 -Werror -Wshadow -Wextra -pipe -DNDEBUG
CXXFLAGS = -std=c++11 -Wall -pg -Wshadow -Werror -O3 -pipe -DNDEBUG

TARGET = cg
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

$(TARGET).o: $(TARGET).cpp Makefile
	$(CXX) -c $(CXXFLAGS) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)
