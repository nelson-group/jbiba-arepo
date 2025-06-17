#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <complex.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

#include "turbulence.h"


#ifdef TURBULENCE
#define INDEX(i, j, k, n) ((i) * (n) * (n) + (j) * (n) + (k))


double turb_random_normal(double mean, double variance) {
    // Generate a standard normal random number using Box-Muller transform
    double u1 = rand() / (double)RAND_MAX;
    double u2 = rand() / (double)RAND_MAX;
    return mean + sqrt(variance) * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}


void turb_OU_check_update() {
    if (All.Time > turb_time) {
        if (ThisTask == 0) {
            mpi_printf("TURBULENCE: Task %d said its time to update the OU phases.\n", ThisTask);
            while (turb_time < All.Time) {
                OU_update();
                turb_time = turb_time + turb_dt;
            }
        }
    }
    mpi_printf("TURBULENCE: Waiting for task 0 to update the OU phases or not.\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(turb_OUphases, 2*turb_nonzero * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&turb_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    mpi_printf("TURBULENCE: Task %d broadcasted the updated phases.\n", ThisTask);
    MPI_Barrier(MPI_COMM_WORLD);
}


void OU_update() {
    double wiener[3] = {0, 0, 0};
    double wiener_decomp[3] = {0, 0, 0};
    for (int i = 0; i < 2*turb_nonzero; i++) {
        for (int l = 0; l < 3; l++) {
            wiener[l] = turb_random_normal(0, 1);
        }
        double k_mag_sq = 0;
        double k_dot = 0;
        for (int l = 0; l < 3; l++) {
            k_mag_sq = k_mag_sq + turb_nonzero_k_modes[i/2][l]*turb_nonzero_k_modes[i/2][l];
            k_dot = k_dot + turb_nonzero_k_modes[i/2][l]*wiener[l];
        }
        // apply spectral decomposition
        for (int l = 0; l < 3; l++) {
            // habe ich das plus da falsch gemacht?
            wiener_decomp[l] = All.turb_zeta*wiener[l] + (1 - 2*All.turb_zeta)*turb_nonzero_k_modes[i/2][l]*k_dot/(k_mag_sq);
        }

        for (int l = 0; l < 3; l++) {
            turb_OUphases[i][l] = exp(-turb_dt/turb_T)*turb_OUphases[i][l] + turb_OUvar * sqrt(1 - exp(-2*turb_dt/turb_T)) * wiener_decomp[l];
        }
    }
}


void get_driving_vec(double pos[3], double acc3d[3]) {
    double ampl, cosx, cosy, cosz, sinx, siny, sinz, real, imag;
    for (int i = 0; i < turb_nonzero; i++) {
        ampl = 2*sqrt(3)/sqrt(1 - 2*All.turb_zeta + 3*All.turb_zeta*All.turb_zeta) * turb_nonzero_amplitudes[i];
        sinx = sin(turb_nonzero_k_modes[i][0]*pos[0]);
        siny = sin(turb_nonzero_k_modes[i][1]*pos[1]);
        sinz = sin(turb_nonzero_k_modes[i][2]*pos[2]);
        cosx = cos(turb_nonzero_k_modes[i][0]*pos[0]);
        cosy = cos(turb_nonzero_k_modes[i][1]*pos[1]);
        cosz = cos(turb_nonzero_k_modes[i][2]*pos[2]);

        real = (cosx*cosy - sinx*siny) * cosz - (sinx*cosy + cosx*siny) * sinz;
        imag = cosx * (cosy*sinz + siny*cosz) + sinx * (cosy*cosz - siny*sinz);

        // real part of the inverse fourier transform, the above is just exp(i*k*x) = cos(k_x *x + ...) + i sin(k_x * x + ...)
        // = real +  i * imag. Then Re(F^-1(f) = Re((a + bi)*(real+i*imag)) = a*real - b*imag)
        acc3d[0] += ampl * (turb_OUphases[2*i][0]*real - turb_OUphases[2*i + 1][0]*imag);
        acc3d[1] += ampl * (turb_OUphases[2*i][1]*real - turb_OUphases[2*i + 1][1]*imag);
        acc3d[2] += ampl * (turb_OUphases[2*i][2]*real - turb_OUphases[2*i + 1][2]*imag);
    }
}


void turb_apply_driving() {
    TIMER_START(CPU_TURBULENCE);
    // printf("WE ARE HERE Task %d\n", ThisTask);
    int count = 0;
    for(int idx = 0; idx < NumGas; idx++) {
        int i = idx; // TimeBinsHydro.ActiveParticleList[idx];
        if (i < 0) {
            continue;
        }
        double acc3d[3] = {0, 0, 0};
        get_driving_vec(P[i].Pos, acc3d);
        for (int l = 0; l < 3; l++) {
            P[i].GravAccel[l] = acc3d[l];
        }
        count++;
    }
    int global_count = 0;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    mpi_printf("TURBULENCE: Updated acceleration of %d particles\n", global_count);
    TIMER_STOP(CPU_TURBULENCE);
}
#endif