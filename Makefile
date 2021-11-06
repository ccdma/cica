CC := g++
CFLAGS := -I ./include/ -fopenmp -O3 -mtune=native -march=native -std=c++11

batch: batch.cpp
	$(CC) $(CFLAGS) $< -o $@

ica: ica.cpp
	$(CC) $(CFLAGS) $< -o $@

# 相対パスで指定することを忘れない！
tssrun:
	tssrun -A t=72:c=72 ./batch

send:
	rsync -avc --exclude '.git' --exclude 'batch' ./ b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/cica

clean:
	rm -f *.out batch ica

.PHONY: clean send tssrun