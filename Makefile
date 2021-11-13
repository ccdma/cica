CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS := batch.cpp
SRCS += single.cpp
SRCS += symbolic.cpp

COMMIT := $(shell git rev-parse --short HEAD | tr -d '\n'; if [ `git status -s -uno | wc -l` -ne 0 ]; then echo "_unstaged"; fi)

$(SRCS:.cpp=):
	$(CC) $(CFLAGS) $@.cpp -o $@.out -DCOMMIT_ID=\"$(COMMIT)\"

clean:
	rm -f *.out *.csv

RSYNC := rsync -avc --exclude '*.out' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

send:
	$(RSYNC)

send-force:
	$(RSYNC) --delete

.PHONY: clean send send-force