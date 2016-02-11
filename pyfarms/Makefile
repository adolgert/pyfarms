libherdsize.so: herd_size_factor.c
	gcc

sizes: herd_size_factor.c herd_size_factor.h herd_size_main.c
	gcc -o sizes herd_size_factor.c herd_size_main.c -lgsl


runone.h5: /work/myers/run_HPAI_oneunit.out
	python naadsmout.py --singlefarm --input /work/myers/run_HPAI_oneunit.out -v
	python naadsmout.py --singlefarm --input /work/myers/run_HPAI_oneunit_shortinfectious.out -v
