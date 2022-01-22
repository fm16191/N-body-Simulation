//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s_SOA {
  f32 *x, *y, *z;
  f32 *vx, *vy, *vz;
} particle_t_SOA;

#define NB_VALUES_PACKED 4

// Horizontal addition of __m128
static inline f32 _hadd_ps(const __m128 a) {
  // return a[0] + a[1] + a[2] + a[3];
  __m128 shuf = _mm_movehdup_ps(a); // [3,2,1,0] -> [2,2,0,0]
  __m128 sums = _mm_add_ps(a, shuf); // [3+2,2+2,1+0,0+0]
  shuf = _mm_movehl_ps(shuf, sums); // [3a,2a,1a,0a][3b,2b,1b,0b] -> [1b,0b,1a,0a] [1+0,0+0,0,0]
  sums = _mm_add_ss(sums, shuf); // [3+2,2+2,1+0,0+0] + [1+0,0+0,0,0] -> [3+2+1+0, 2+2+0+0, 1+0+0, 0+0+0]
  return _mm_cvtss_f32(sums); // return [0] : 3+2+1+0
}

// Random particle initialization
void init_SOA(particle_t_SOA* p, u64 n) {
  for (u64 i = 0; i < n; i++) {
    //  
    u64 r1 = (u64)rand();
    u64 r2 = (u64)rand();
    f32 sign = (r1 > r2) ? 1 : -1;

    //
    p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
    p->y[i] = (f32)rand() / (f32)RAND_MAX;
    p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

    //
    p->vx[i] = (f32)rand() / (f32)RAND_MAX;
    p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
    p->vz[i] = (f32)rand() / (f32)RAND_MAX;
  }
}


// Malloc for particle_t_SOA
void malloc_particle_t_SOA(particle_t_SOA* p, u64 n){
  p->x = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
  p->y = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
  p->z = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
  p->vx = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
  p->vy = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
  p->vz = aligned_alloc(sizeof(f32)*2*NB_VALUES_PACKED, sizeof(f32) * n / NB_VALUES_PACKED * NB_VALUES_PACKED);
}

void free_particle_t_SOA(particle_t_SOA* p){
  free(p->x);
  free(p->y);
  free(p->z);
  free(p->vx);
  free(p->vy);
  free(p->vz);
}

void print_debug(particle_t_SOA* p, u64 n, char* DEBUG_filename, u64 step) {
  FILE* fp;
  if (step == 0) fp = fopen(DEBUG_filename, "w");
  else fp = fopen(DEBUG_filename, "a");
  fprintf(fp, "===%llu===\n", step);
  for (u64 i=0; i < n; i++) {
    fprintf(fp, "|%llu|%e_%e_%e\n", i, p->x[i], p->y[i], p->z[i]);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

//
void move_particles_SOA(particle_t_SOA* p, const f32 dt, u64 n) {
  //
  const f32 softening = 1e-20;

  __m128 fx, fy, fz; 
  __m128 pxi, pyi, pzi; 
  __m128 pvxi, pvyi, pvzi; 
  __m128 dx, dy, dz; 

  __m128 d_2, soften, vdt;

  soften = _mm_set1_ps(softening);
  vdt = _mm_set1_ps(dt);

  //
  for (u64 i = 0; i < n; i ++) {
    //
    fx = _mm_setzero_ps();
    fy = _mm_setzero_ps();
    fz = _mm_setzero_ps();

    pxi = _mm_set1_ps(p->x[i]);
    pyi = _mm_set1_ps(p->y[i]);
    pzi = _mm_set1_ps(p->z[i]);

    for (u64 j = 0; j < n; j += NB_VALUES_PACKED) {
      //Newton's law

      //// dx, dy, dz calculation. 3 Floating operations (3 loop total)
      dx = _mm_sub_ps(_mm_loadu_ps(p->x+j), pxi);
      dy = _mm_sub_ps(_mm_loadu_ps(p->y+j), pyi);
      dz = _mm_sub_ps(_mm_loadu_ps(p->z+j), pzi);

      //// d_2 Calculation : 6 Floating operations (9 loop total)
      d_2 = _mm_mul_ps(dx, dx) + _mm_mul_ps(dy, dy) + _mm_mul_ps(dz, dz) + soften;

      //// d3_inv_sqrt_d_2 calculation. 5 Floating operations (14 loop total)
      // sqrt_d_2 = _mm_sqrt_ps(d_2);
      // inv_sqrt_d_2 = _mm_div_ps(_mm_set1_ps(1.0), sqrt_d_2);
      d_2 = _mm_rsqrt_ps(d_2);
      d_2 = _mm_mul_ps(_mm_mul_ps(d_2, d_2), d_2);

      //// Net force calculation. 6 Floating operations (20 loop total)
      // temp = _mm_mul_ps(dx, d3_inv_sqrt_d_2);
      // fx = _mm_add_ps(temp, fx);
      // temp = _mm_mul_ps(dy, d3_inv_sqrt_d_2);
      // fy = _mm_add_ps(temp, fy);
      // temp = _mm_mul_ps(dz, d3_inv_sqrt_d_2);
      // fz = _mm_add_ps(temp, fz);
      fx = _mm_fmadd_ps(dx, d_2, fx);
      fy = _mm_fmadd_ps(dy, d_2, fy);
      fz = _mm_fmadd_ps(dz, d_2, fz);
    }

    //// Velocity calculation. 6 Floating operationgs (6 loop total)
    // temp = _mm_mul_ps(vdt, _mm_set1_ps(_hadd_ps(fx)));
    // pvxi = _mm_add_ps(temp, _mm_set1_ps(p->vx[i]));
    // temp = _mm_mul_ps(vdt, _mm_set1_ps(_hadd_ps(fy)));
    // pvyi = _mm_add_ps(temp, _mm_set1_ps(p->vy[i]));
    // temp = _mm_mul_ps(vdt, _mm_set1_ps(_hadd_ps(fz)));
    // pvzi = _mm_add_ps(temp, _mm_set1_ps(p->vz[i]));

    p->vx[i] += dt * _hadd_ps(fx);
    p->vy[i] += dt * _hadd_ps(fy);
    p->vz[i] += dt * _hadd_ps(fz);
    // pvxi = _mm_fmadd_ps(vdt, _mm_set1_ps(_hadd_ps(fx)), _mm_set1_ps(p->vx[i]));
    // pvyi = _mm_fmadd_ps(vdt, _mm_set1_ps(_hadd_ps(fy)), _mm_set1_ps(p->vy[i]));
    // pvzi = _mm_fmadd_ps(vdt, _mm_set1_ps(_hadd_ps(fz)), _mm_set1_ps(p->vz[i]));
    // p->vx[i] = _mm_cvtss_f32(pvxi);
    // p->vy[i] = _mm_cvtss_f32(pvyi);
    // p->vz[i] = _mm_cvtss_f32(pvzi);
  }

  for (u64 i = 0; i < n; i += NB_VALUES_PACKED) {

    // Loading velocities
    pvxi = _mm_loadu_ps(p->vx + i);
    pvyi = _mm_loadu_ps(p->vy + i);
    pvzi = _mm_loadu_ps(p->vz + i);

    //// Positions calculation. 6 Floating operations (12 loop total)
    // temp = _mm_mul_ps(vdt, pvxi);
    // pxi = _mm_add_ps(temp, _mm_loadu_ps(p->x + i));
    // temp = _mm_mul_ps(vdt, pvyi);
    // pyi = _mm_add_ps(temp, _mm_loadu_ps(p->y + i));
    // temp = _mm_mul_ps(vdt, pvzi);
    // pzi = _mm_add_ps(temp, _mm_loadu_ps(p->z + i));
    pxi = _mm_fmadd_ps(vdt, pvxi, _mm_loadu_ps(p->x + i));
    pyi = _mm_fmadd_ps(vdt, pvyi, _mm_loadu_ps(p->y + i));
    pzi = _mm_fmadd_ps(vdt, pvzi, _mm_loadu_ps(p->z + i));

    _mm_storeu_ps(p->x + i, pxi);
    _mm_storeu_ps(p->y + i, pyi);
    _mm_storeu_ps(p->z + i, pzi);
  }
}

int usage(char** argv, int argi){
  if (argv[argi][1] != 'h') fprintf(stderr, "Invalid argument : %s.\n\n", argv[argi]);
  fprintf(stderr, "USAGE\n\n");
  fprintf(stderr, "\t-h\t\tPrints out this message.\n");
  fprintf(stderr, "\t-c\t\tPrints configuration.\n");
  fprintf(stderr, "\t-b[=bench.file]\tActivate BENCHMARK mode. Prints average performance to the given file. By default, -b=bench_[prog_name].dat\n");
  fprintf(stderr, "\t-d[=debug.file]\tActivate DEBUG mode. Prints positions to the given file. By default, -d=out.dat\n");
  fprintf(stderr, "\t-n=<nbodies> \tSpecify the number of nbodies. By default, -n=16384\n");
  fprintf(stderr, "\t-s=<steps> \tSpecify the number of steps. By default, -s=10\n");
  fprintf(stderr, "\t-w=<warmups> \tSpecify the number of warmups. By default, -w=3\n");
  fprintf(stderr, "\t-t=<deltaTime> \tSpecify the time delta. By default, -t=0.01\n");
  fprintf(stderr, "\n");
  return 1;
}


//
int main(int argc, char** argv) {
  srand(0);
  
  //
  u64 n = 16384; // = (argc > 1) ? atoll(argv[1]) : 16384;
  int CONFIG = 0;

  u64 steps = 10;
  f32 dt = 0.01;

  // DEBUG mode
  int DEBUG = 0;
  char* DEBUG_filename;
  DEBUG_filename = malloc(sizeof(char) * 8);
  strcpy(DEBUG_filename, "out.dat");

  // BENCH mode
  int BENCH = 0;
  int strl = strlen(argv[0]);
  char* BENCH_filename;
  BENCH_filename = malloc(sizeof(char) * (9 + strl)); // 6 + 4 + strl - 2
  char* test = malloc(sizeof(char) * strl);
  strcpy(test, argv[0]);
  // printf("=%s=\n", test);
  strcpy(BENCH_filename, "bench_");
  strncpy(BENCH_filename+6, test+2, strl - 2);
  strcat(BENCH_filename+6+strl-2, ".dat");
  free(test);

  //
  f64 rate = 0.0, drate = 0.0;

  // Steps to skip for warm up
  u64 warmup = 3;

  if (argc > 1) {
    char* arg;
    int arglen;
    for (int argi = 1; argi < argc; argi++) {
      if (argv[argi][0] == '-') {
        arglen = strlen(argv[argi]);
        if (arglen < 2) return usage(argv, argi);
        if (argv[argi][1] == 'd') {
          DEBUG = 1;
          if (arglen > 3 && argv[argi][2] == '=') {
            DEBUG_filename = malloc(sizeof(char) * (arglen - 2));
            if (!DEBUG_filename) return fprintf(stderr, "Cannot allocate memory for the debug file.");
            strncpy(DEBUG_filename, argv[argi] + 3, arglen - 2);
          }
        }
        else if (argv[argi][1] == 'b') {
          BENCH = 1;
          if (arglen > 3 && argv[argi][2] == '=') {
            BENCH_filename = malloc(sizeof(char) * (arglen - 2));
            if (!BENCH_filename) return fprintf(stderr, "Cannot allocate memory for the bench file.");
            strncpy(BENCH_filename, argv[argi] + 3, arglen - 2);
          }
        }
        else if (argv[argi][1] == 'h') return usage(argv, argi);
        else if (argv[argi][1] == 'c') CONFIG = 1;
        else {
          if (arglen < 3 || argv[argi][2] != '=') return usage(argv, argi);
          arg = malloc(sizeof(char) * (arglen - 2));
          strncpy(arg, argv[argi] + 3, arglen - 2);
          if (atof(arg) < 0) return usage(argv, argi);
          if (argv[argi][1] == 'n') n = atoi(arg);
          else if (argv[argi][1] == 'w') warmup = atoi(arg);
          else if (argv[argi][1] == 's') steps = atoi(arg);
          else if (argv[argi][1] == 't') dt = atof(arg);
          else return usage(argv, argi);
          free(arg);
        }
      }
      else fprintf(stderr, "Ignoring unrecognized argument : %s\n", argv[argi]);
    }
  }

  // Initialisation
  particle_t_SOA p;
  malloc_particle_t_SOA(&p, n);
  init_SOA(&p,n);

  const u64 s = 6 * n * sizeof(f32);

  fprintf(stderr, "N-body 3D Problem version 4.0 - using SSE optimized\n\n");
  if (CONFIG) {
    fprintf(stderr, "CONFIGURATION :\n");
    fprintf(stderr, " - DEBUG mode :\t %s\n", ((DEBUG == 1) ? DEBUG_filename : "0"));
    fprintf(stderr, " - BENCH mode :\t %s\n", ((BENCH == 1) ? BENCH_filename : "0"));
    fprintf(stderr, " - nbodies :\t %llu\n", n);
    fprintf(stderr, " - steps :\t %llu\n", steps);
    fprintf(stderr, " - warmups :\t %llu\n", warmup);
    fprintf(stderr, " - dt :\t\t %f\n\n", dt);
  }

  fprintf(stderr, "\n\033[1mTotal memory size:\033[0m %llu B, %.2lf KiB, %.2lf MiB\n\n", s, (f32)s / 1024.0f, (f32)s / 1048576.0f);
  fprintf(stderr, "\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s");

  //
  for (u64 i = 0; i < steps; i++) {
    // Measure travel computation
    const f64 start = omp_get_wtime();
    move_particles_SOA(&p, dt, n);
    const f64 end = omp_get_wtime();

    if (DEBUG) print_debug(&p, n, DEBUG_filename, i);

    //Number of interactions/iterations
    const f32 h1 = (f32)(n) * (f32)(n);

    //GFLOPS
    const f32 h2 = (20 * h1 + 12 * (f32)n) * 1e-9;

    if (i >= warmup) {
      rate += h2 / (end - start);
      drate += (h2 * h2) / ((end - start) * (end - start));
    }

    //
    fprintf(stderr, "%5llu %10.3e %10.3e %8.1f %s\n",
      i,
      (end - start),
      h1 / (end - start),
      h2 / (end - start),
      (i < warmup) ? "*" : "");
          
    if (BENCH){
      FILE *fp;
      if (i<warmup) continue;
      if(i == warmup) fp = fopen(BENCH_filename, "w"); 
      else fp = fopen(BENCH_filename, "a");
      fprintf(fp, "%llu %e %e %f\n",
      i-warmup,
      (end - start),
      h1 / (end - start),
      h2 / (end - start));
      fclose(fp);
    }
}

  //
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

  fprintf(stderr, "-----------------------------------------------------\n");
  fprintf(stderr, "\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
    "Average performance:", "", rate, drate);
  fprintf(stderr, "-----------------------------------------------------\n");

  //
  free_particle_t_SOA(&p);
  free(DEBUG_filename);
  free(BENCH_filename);

  //
  return 0;
}