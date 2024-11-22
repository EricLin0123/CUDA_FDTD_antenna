gauss: gauss.c
	gcc -shared -o gaussian.so -fPIC gauss.c -lm