CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS := batch.cpp
SRCS += single.cpp
SRCS += symbolic.cpp

$(SRCS:.cpp=):
	$(CC) $(CFLAGS) $@.cpp -o $@.out

clean:
	rm -f *.out *.csv

RSYNC := rsync -avc --exclude '.git' --exclude '*.out' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

send:
	$(RSYNC)

send-force:
	$(RSYNC) --delete

.PHONY: clean send send-force