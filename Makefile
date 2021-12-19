CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS := batch.cpp
SRCS += single.cpp
SRCS += symbolic.cpp
SRCS += csymbolic.cpp

COMMIT := $(shell git rev-parse --short HEAD | tr -d '\n'; if [ `git status -s -uno | wc -l` -ne 0 ]; then echo "_unstaged"; fi)

$(SRCS:.cpp=.out): clean
	$(CC) $(CFLAGS) $(@:.out=.cpp) -o $@ -DCOMMIT_ID=\"$(COMMIT)\"

adi:
	$(CC) $(CFLAGS) $@.cpp -liio -o $@.out -DCOMMIT_ID=\"$(COMMIT)\"

clean:
	rm -f *.out 

clean-all: clean
	rm -f *.csv

RSYNC := rsync -avc --exclude '*.out' --exclude './*.csv' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

send:
	$(RSYNC)

send-force:
	$(RSYNC) --delete

push:
	git add -A .
	git commit -m "modified: $(shell git diff master --name-only)"
	make send
	git push

.PHONY: clean clean-all send send-force