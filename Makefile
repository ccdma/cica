CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11
SRCS := batch.cpp
SRCS += single.cpp
SRCS += symbolic.cpp

$(SRCS:.cpp=):
	$(CC) $(CFLAGS) $@.cpp -o $@.out

send:
	rsync -avc --exclude '.git' --exclude '*.out' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

clean:
	rm -f *.out *.csv

.PHONY: clean send