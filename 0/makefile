all: nbody.gcc nbody.clang nbody.icc nbody.icx

nbody.gcc: nbody.c
	gcc -march=native -Wall -Wextra -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

nbody.gccp: nbody.c
	gcc -march=native -O2 $< -o $@ -lm -fopenmp

nbody.clang: nbody.c
	clang -Wall -Wextra -march=native -Ofast $< -o $@ -lm -fopenmp

nbody.clangp: nbody.c
	clang -O2 -ffp-model=precise $< -o $@ -lm -fopenmp

nbody.clangp_f64: nbody_f64.c
	clang -O2 -ffp-model=precise $< -o $@ -lm -fopenmp

nbody.icc: nbody.c
	icc -xhost -Ofast -qopt-report-file=$@.optrpt $< -o $@ -qmkl -qopenmp

nbody.iccp: nbody.c
	icc -xhost -O2 -fp-model=precise -qopt-report-file=$@.optrpt $< -o $@ -qmkl -qopenmp

nbody.icx: nbody.c
	icx -xhost -Ofast -qopt-report-file=$@.optrpt $< -o $@ -qmkl -qopenmp

prof: nbody.c
	gcc -pg -march=znver3 -mtune=znver3 -finline-functions -mavx -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

verif: clean nbody.clang
	rm -f out.dat
	taskset -c 3 ./nbody.clang -n=10000 -d -s=1 -w=0 -c

verif_f64: clean nbody.clangp_f64
	rm -f out.dat
	taskset -c 3 ./nbody.clangp_f64 -n=10000 -d -s=1 -w=0 -c

verifKNL: clean nbody.icc
	rm -f out.dat
	taskset -c 3 ./nbody.icc -n=10000 -d -s=1 -w=0 -c

bench : clean nbody.gcc nbody.clang
	mkdir -p benchs
	taskset -c 3 ./nbody.gcc -s=23 -b=benchs/b_gcc.dat
	taskset -c 3 ./nbody.clang -s=23 -b=benchs/b_clang.dat

benchKNL : nbody.gcc nbody.icc nbody.icx
	mkdir -p benchs
	taskset -c 3 ./nbody.gcc -s=23 -b=benchs/b_gcc.dat
	# taskset -c 3 ./nbody.clang -s=23 -b=benchs/b_clang.dat
	taskset -c 3 ./nbody.icc -s=23 -b=benchs/b_icc.dat
	taskset -c 3 ./nbody.icx -s=23 -b=benchs/b_icx.dat

clean:
	rm -Rf *~ nbody.gcc nbody.gccp nbody.icc nbody.iccp nbody.icx nbody.clang nbody.clangp nbody.clangp_f64 prof gmon.out *.optrpt

