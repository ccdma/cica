CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
COMMIT_ID := $(sh git rev-parse --short HEAD)
UNSTAGED := $(sh if [ `git status --porcelain --untracked-files=no` -n 1 ]; then echo true;)
SRCS := batch.cpp
SRCS += single.cpp
SRCS += symbolic.cpp

$(SRCS:.cpp=):
	$(CC) $(CFLAGS) $@.cpp -o $@.out

clean:
	rm -f *.out *.csv

RSYNC := rsync -avc --exclude '*.out' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

send:
	$(RSYNC)

send-force:
	$(RSYNC) --delete

.PHONY: clean send send-force