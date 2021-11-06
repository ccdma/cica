CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -std=c++11

batch: batch.cpp
	$(CC) $(CFLAGS) $< -o $@

ica: ica.cpp
	$(CC) $(CFLAGS) $< -o $@

tssrun:
	tssrun -A t=144:c=72:m=10G batch

send:
	rsync -avc --exclude '.git' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

clean:
	rm -f *.out batch ica

.PHONY: clean send tssrun