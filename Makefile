CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11

single: single.cpp
	$(CC) $(CFLAGS) $< -o $@

batch: batch.cpp
	$(CC) $(CFLAGS) $< -o $@

# 相対パスで指定することを忘れない！
tssrun:
	tssrun -A t=72:c=72 ./batch

send:
	rsync -avc --exclude '.git' --exclude 'batch' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

clean:
	rm -f *.out batch single

.PHONY: clean send tssrun