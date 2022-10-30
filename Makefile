CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS += single.cpp
SRCS += symbolic.cpp
SRCS += csymbolic.cpp
SRCS += cdma.cpp
SRCS += orth.cpp

COMMIT := $(shell git rev-parse --short HEAD | tr -d '\n'; if [ `git status -s -uno | wc -l` -ne 0 ]; then echo "_unstaged"; fi)

$(SRCS:.cpp=.out): clean commit
	$(CC) $(CFLAGS) $(@:.out=.cpp) -o $@ -DCOMMIT_ID=\"$(COMMIT)\"

clean:
	rm -f *.out 

clean-all: clean
	rm -f *.csv

commit:
	git add -A .
	git commit -m "modified: $(shell git diff master --name-only)"

.PHONY: clean clean-all