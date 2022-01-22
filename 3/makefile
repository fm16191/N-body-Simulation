all: nbody.gcc nbody.clang nbody.icc nbody.icx

nbody.gcc: nbody.c
	gcc -Wall -Wextra -g3 -march=native -mavx -Ofast -finline-functions -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp
	# gcc -Wall -Wextra -g3 -march=znver3 -mavx -Ofast -finline-functions -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

nbody.clang: nbody.c
	clang -march=native -Wall -Wextra -g3 -Wpadded -march=native -mtune=native -mavx -Ofast -finline-functions -funroll-loops $< -o $@ -lm -fopenmp
	
nbody.aocc: nbody.c
	~/aocc/bin/clang -Wall -Wextra -g3 -march=native -mtune=native -mavx -Ofast -finline-functions -funroll-loops $< -o $@ -lm -fopenmp

nbody.icc: nbody.c
	icc -xhost -fast -fp-model=fast -qopt-report-file=$@.optrpt $< -o $@ -qmkl -qopenmp

nbody.icx: nbody.c
	icx -xhost -fast -fp-model=fast -qopt-report-file=$@.optrpt $< -o $@ -qmkl -qopenmp

prof: nbody.c
	gcc -pg -march=znver3 -mtune=znver3 -finline-functions -mavx -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

verif_last: clean nbody.clang
	rm -f out.dat
	taskset -c 3 ./nbody.clang -n=10000 -d -s=1 -w=0 -c
	python3 precision_checker.py ../2/out.dat out.dat

verif: clean nbody.clang
	rm -f out.dat
	taskset -c 3 ./nbody.clang -n=10000 -d -s=1 -w=0 -c
	python3 precision_checker.py

verifKNL: clean nbody.icc
	rm -f out.dat
	taskset -c 3 ./nbody.icc -n=10000 -d -s=1 -w=0 -c
	python3 precision_checker.py

verif_all: clean nbody.gcc nbody.clang nbody.aocc
	mkdir -p precisions
	touch precisions/avg_p.dat
	taskset -c 3 ./nbody.gcc -n=10000 -s=1 -w=0 -d=precisions/p_gcc.dat
	taskset -c 3 ./nbody.clang -n=10000 -s=1 -w=0 -d=precisions/p_clang.dat
	taskset -c 3 ./nbody.aocc -n=10000 -s=1 -w=0 -d=precisions/p_aocc.dat
	python3 precision_checker.py precisions/

verif_allKNL: nbody.gcc nbody.icc nbody.icx
	mkdir -p precisions
	touch precisions/avg_p.dat
	taskset -c 3 ./nbody.gcc -n=10000 -s=1 -w=0 -d=precisions/p_gcc.dat
	# taskset -c 3 ./nbody.clang -n=10000 -s=1 -w=0 -d=precisions/p_clang.dat
	taskset -c 3 ./nbody.icc -n=10000 -s=1 -w=0 -d=precisions/p_icc.dat
	taskset -c 3 ./nbody.icx -n=10000 -s=1 -w=0 -d=precisions/p_icccx.dat
	python3 precision_checker.py precisions/

bench : clean nbody.gcc nbody.clang nbody.aocc
	mkdir -p benchs
	taskset -c 3 ./nbody.gcc -s=23 -b=benchs/b_gcc.dat
	taskset -c 3 ./nbody.clang -s=23 -b=benchs/b_clang.dat
	taskset -c 3 ./nbody.aocc -s=23 -b=benchs/b_amd.dat

benchKNL : nbody.gcc nbody.icc nbody.icx
	mkdir -p benchs
	taskset -c 3 ./nbody.gcc -s=23 -b=benchs/b_gcc.dat
	# taskset -c 3 ./nbody.clang -s=23 -b=benchs/b_clang.dat
	taskset -c 3 ./nbody.icc -s=23 -b=benchs/b_icc.dat
	taskset -c 3 ./nbody.icx -s=23 -b=benchs/b_icx.dat

clean:
	rm -Rf *~ nbody.gcc nbody.clang nbody.aocc nbody.icc prof gmon.out *.optrpt

run: clean nbody.gcc
	time taskset -c 3 ./nbody.gcc

maqao:
	# rm -rf maqao*
	scp nbody.c makefile fob1:nbody3D/3
	# ssh fob1 rm -rf /users/user2235/nbody3D/3/maqao*
	# ssh fob1 ssh knl05 "module load 2.15.0"
	# ssh fob1 ssh knl05 maqao oneview -R1 --binary="./nbody.clang" --run_command="./nbody.clang 50000"

get_makao:
	rm -rf maqao_last
	scp -r fob1:nbody3D/3/maqao* maqao_last
	firefox maqao_last/RESULTS/*/index.html
	# rm -rf /users/user2235/nbody3D/3/maqao* && maqao oneview -R1 --executable="nbody.gcc" --run_command="nbody.gcc 50000"
	# rm -rf /users/user2235/nbody3D/3/maqao* && maqao oneview -R1 --executable="nbody.clang" --run_command="nbody.clang 50000"
	# rm -rf /users/user2235/nbody3D/3/maqao* && maqao oneview -R1 --executable="nbody.icc" --run_command="<executable> 50000"