//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #define DEBUG 1
int DEBUG = 0;

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

    f64 x, y, z;
    f64 vx, vy, vz;

} particle_t;


//
void init(particle_t* p, u64 n) {
    for (u64 i = 0; i < n; i++) {
        //
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        f64 sign = (r1 > r2) ? 1 : -1;

        //
        p[i].x = sign * (f64)rand() / (f64)RAND_MAX;
        p[i].y = (f64)rand() / (f64)RAND_MAX;
        p[i].z = sign * (f64)rand() / (f64)RAND_MAX;

        //
        p[i].vx = (f64)rand() / (f64)RAND_MAX;
        p[i].vy = sign * (f64)rand() / (f64)RAND_MAX;
        p[i].vz = (f64)rand() / (f64)RAND_MAX;
    }
}

//
void move_particles(particle_t* p, const f64 dt, u64 n) {
    //
    const f64 softening = 1e-20;

    //
    for (u64 i = 0; i < n; i++) {
        //
        f64 fx = 0.0;
        f64 fy = 0.0;
        f64 fz = 0.0;

        //23 floating-point operations
        for (u64 j = 0; j < n; j++) {
            //Newton's law
            const f64 dx = p[j].x - p[i].x; //1
            const f64 dy = p[j].y - p[i].y; //2
            const f64 dz = p[j].z - p[i].z; //3
            const f64 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9
            const f64 d_3_over_2 = pow(d_2, 3.0 / 2.0); //11

            //Net force
            fx += dx / d_3_over_2; //13
            fy += dy / d_3_over_2; //15
            fz += dz / d_3_over_2; //17
        }

        //
        p[i].vx += dt * fx; //1
        p[i].vy += dt * fy; //2
        p[i].vz += dt * fz; //3
    }

    //3 floating-point operations
    for (u64 i = 0; i < n; i++) {
        p[i].x += dt * p[i].vx; // 5
        p[i].y += dt * p[i].vy; // 7
        p[i].z += dt * p[i].vz; // 9

        if (DEBUG) printf("|%llu|%.64lf_%.64lf_%.64lf\n", i, p[i].x, p[i].y, p[i].z);
    }
}

//
int main(int argc, char** argv) {
    srand(0);

    //
    const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
    DEBUG = (argc > 2) ? atoll(argv[2]) : 0;
    const u64 steps = 10;
    const f64 dt = 0.01;

    //
    f64 rate = 0.0, drate = 0.0;

    //Steps to skip for warm up
    const u64 warmup = 3;

    //
    particle_t* p = malloc(sizeof(particle_t) * n);

    //
    init(p, n);

    const u64 s = sizeof(particle_t) * n;

    // printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
    fprintf(stderr,
        "\n\033[1mTotal memory size:\033[0m %llu B, %.2lf KiB, %.2lf MiB\n\n",
        s, (f64)s / 1024.0f, (f64)s / 1048576.0f);

    //
    // printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
    fprintf(stderr, "\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s");

    //
    for (u64 i = 0; i < steps; i++) {
        //Measure
        const f64 start = omp_get_wtime();

        if (DEBUG) printf("\n===%llu===\n", i);
        move_particles(p, dt, n);

        const f64 end = omp_get_wtime();

        //Number of interactions/iterations
        const f64 h1 = (f64)(n) * (f64)(n);

        //GFLOPS
        const f64 h2 = (17 * h1 + 9.0 * (f64)n) * 1e-9;

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
    }

    //
    rate /= (f64)(steps - warmup);
    drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

    fprintf(stderr, "-----------------------------------------------------\n");
    fprintf(stderr, "\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
        "Average performance:", "", rate, drate);
    fprintf(stderr, "-----------------------------------------------------\n");

    //
    free(p);

    //
    return 0;
}
