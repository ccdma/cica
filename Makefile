CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS += single.cpp
SRCS += symbolic.cpp
SRCS += csymbolic.cpp

COMMIT := $(shell git rev-parse --short HEAD | tr -d '\n'; if [ `git status -s -uno | wc -l` -ne 0 ]; then echo "_unstaged"; fi)

$(SRCS:.cpp=.out): clean
	$(CC) $(CFLAGS) $(@:.out=.cpp) -o $@ -DCOMMIT_ID=\"$(COMMIT)\"

clean:
	rm -f *.out 

clean-all: clean
	rm -f *.csv

.PHONY: clean clean-all