/*
 * main_fast.cpp — Standalone Sinkhorn reflector, no MKL, NK=381
 *
 * Compile:  cd BenchmarkCode && g++ -O2 -o main_fast main_fast.cpp -lm
 * Run:      ./main_fast
 *
 * Uses the small (381-point) QMC cloud as the sole problem grid.
 * Skips the UseSmall warm-start (NK == NK_small, so starts from f=g=0).
 * Skips the FinalGrid regular-grid computation.
 * Outputs to output_fast/  (f.txt, g.txt, f_id.txt, g_id.txt, R.txt, Ref.txt).
 */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

using namespace std;

/* ------------------------------------------------------------------
 * MKL replacement functions — trivial loops over the six calls used
 * ------------------------------------------------------------------ */
inline void vdMul(int n, const double* a, const double* b, double* c) {
    for (int i = 0; i < n; i++) c[i] = a[i] * b[i];
}
inline void vdLn(int n, const double* a, double* c) {
    for (int i = 0; i < n; i++) c[i] = log(a[i]);
}
inline void vdAdd(int n, const double* a, const double* b, double* c) {
    for (int i = 0; i < n; i++) c[i] = a[i] + b[i];
}
inline void vdExp(int n, const double* a, double* c) {
    for (int i = 0; i < n; i++) c[i] = exp(a[i]);
}
inline int cblas_idamax(int n, const double* a, int /*inc*/) {
    int idx = 0;
    double mx = fabs(a[0]);
    for (int i = 1; i < n; i++) {
        if (fabs(a[i]) > mx) { mx = fabs(a[i]); idx = i; }
    }
    return idx;
}
inline void cblas_dscal(int n, double s, double* a, int /*inc*/) {
    for (int i = 0; i < n; i++) a[i] *= s;
}

/* ------------------------------------------------------------------
 * Small cloud header — defines NK_small=381,
 *   x_small[NK_small][3]  (upper hemisphere, z > 0)
 *   y_small[NK_small][3]  (same data; z MUST be negated before use)
 * ------------------------------------------------------------------ */
#include "SmallGrid/3D_MonteCarlo_Pointcloud_small.h"

/* ------------------------------------------------------------------
 * Problem dimensions
 * ------------------------------------------------------------------ */
const int dim = 3;
const int NK  = NK_small;   /* 381 */

int getk(int Nk) { return (int)sqrt((double)Nk); }

/* ------------------------------------------------------------------
 * Density functions   (SquareToCircle benchmark)
 * P: uniform on stereographic north-projection square [-0.5,0.5]^2
 * Q: uniform on stereographic south-projection disk  |.| <= 0.5
 * ------------------------------------------------------------------ */
double P(double v[]) {
    double p1 = v[0] / (1.0 + v[2]);
    double p2 = v[1] / (1.0 + v[2]);
    return (-0.5 < p1 && p1 < 0.5 && -0.5 < p2 && p2 < 0.5) ? 1.0 : 0.0;
}
double Q(double v[]) {
    double p1 = v[0] / (1.0 - v[2]);
    double p2 = v[1] / (1.0 - v[2]);
    return (p1 * p1 + p2 * p2 <= 0.25) ? 1.0 : 0.0;
}

/* Wang (2004) reflector cost: c(x,y) = -log(1 - x·y) */
double Cost_Func(double a[], double b[]) {
    double t = 1.0;
    for (int d = 0; d < dim; d++) t -= a[d] * b[d];
    return -log(t);
}

/* ------------------------------------------------------------------
 * Working arrays
 * x = x_small  (upper hemisphere, z > 0)
 * y = x_small with z negated  (lower hemisphere, z < 0)
 * ------------------------------------------------------------------ */
double x[NK][dim];
double y[NK][dim];

double p[NK], q[NK], logp[NK], logq[NK];

double F[NK], G[NK];
double f[NK], g[NK];
double F_id[NK], G_id[NK];
double f_id[NK], g_id[NK];

double fc[NK], gc[NK];
double R[NK];
double Ref[NK][dim];

double tempvecNK_1[NK], tempvecNK_2[NK], tempvecNK_3[NK];

/* threshold: skip summation entries whose contribution is < 1e-6/NK */
double treshold = log(1.e-6 / NK);

int IndJ_begin[NK], IndJ_end[NK];
int IndI_begin[NK], IndI_end[NK];
int IndId_begin[NK], IndId_end[NK];

double maxdif = -1.0;
double Reg_Param;
string outputname = "output_fast";

const int    multiplier    = 8;
const double cap_treshold  = 1.e-5;
const int    cap_iteration = 16;

/* ------------------------------------------------------------------
 * fillall: initialise all working arrays
 * ------------------------------------------------------------------ */
void fillall()
{
    /* Source = x_small (upper hemisphere).
     * Target = x_small with z negated → lower hemisphere antipodal points. */
    for (int i = 0; i < NK; i++) {
        x[i][0] = x_small[i][0];
        x[i][1] = x_small[i][1];
        x[i][2] = x_small[i][2];
        y[i][0] = x_small[i][0];
        y[i][1] = x_small[i][1];
        y[i][2] = -x_small[i][2];
    }

    double sum = 0.0;
    for (int i = 0; i < NK; i++) { p[i] = P(x[i]); sum += p[i]; }
    for (int i = 0; i < NK; i++) {
        p[i] /= sum;
        logp[i] = (p[i] > 0.0) ? log(p[i]) : -1.0e300;
    }
    printf("Source support: %d / %d points\n",
           (int)(sum + 0.5), NK);

    sum = 0.0;
    for (int i = 0; i < NK; i++) { q[i] = Q(y[i]); sum += q[i]; }
    for (int i = 0; i < NK; i++) {
        q[i] /= sum;
        logq[i] = (q[i] > 0.0) ? log(q[i]) : -1.0e300;
    }
    printf("Target support: %d / %d points\n",
           (int)(sum + 0.5), NK);

    fill_n(F,    NK, 1.0);  fill_n(G,    NK, 1.0);
    fill_n(f,    NK, 0.0);  fill_n(g,    NK, 0.0);
    fill_n(fc,   NK, 0.0);  fill_n(gc,   NK, 0.0);
    fill_n(F_id, NK, 1.0);  fill_n(G_id, NK, 1.0);
    fill_n(f_id, NK, 0.0);  fill_n(g_id, NK, 0.0);

    fill_n(IndJ_begin, NK, 0);   fill_n(IndJ_end, NK, NK);
    fill_n(IndI_begin, NK, 0);   fill_n(IndI_end, NK, NK);
    fill_n(IndId_begin, NK, 0);  fill_n(IndId_end, NK, NK);
}

/* ------------------------------------------------------------------
 * One Sinkhorn iteration (main transport: G update then F update)
 * Matches C++ Sinkhorn_axb.
 * ------------------------------------------------------------------ */
void Sinkhorn_axb(int k)
{
    /* G update: G[j] = 1 / sum_i( exp(-k*(C-f[i]-g[j])) * F[i]*p[i] ) */
    vdMul(NK, F, p, tempvecNK_1);   /* tempvecNK_1 = F*p */
    for (int j = 0; j < NK; j++) {
        if (q[j] == 0.0) continue;
        double sum = 0.0;
        int ind_begin = 0, ind_end = 0;
        bool first = true;
        for (int i = IndI_begin[j]; i < IndI_end[j]; i++) {
            double temp = -k * (Cost_Func(x[i], y[j]) - f[i] - g[j]);
            if (logp[i] + temp < treshold) continue;
            if (first) { ind_begin = i; first = false; }
            ind_end = i;
            sum += exp(temp) * tempvecNK_1[i];
        }
        IndI_begin[j] = ind_begin;
        IndI_end[j]   = ind_end + 1;
        G[j] = (sum > 0.0) ? 1.0 / sum : 1.0;
    }

    /* F update: F[i] = 1 / sum_j( exp(-k*(C-f[i]-g[j])) * G[j]*q[j] ) */
    vdMul(NK, G, q, tempvecNK_1);   /* tempvecNK_1 = G*q */
    for (int i = 0; i < NK; i++) {
        if (p[i] == 0.0) continue;
        double sum = 0.0;
        int ind_begin = 0, ind_end = 0;
        bool first = true;
        for (int j = IndJ_begin[i]; j < IndJ_end[i]; j++) {
            double temp = -k * (Cost_Func(x[i], y[j]) - f[i] - g[j]);
            if (temp + logq[j] < treshold) continue;
            if (first) { ind_begin = j; first = false; }
            ind_end = j;
            sum += exp(temp) * tempvecNK_1[j];
        }
        IndJ_begin[i] = ind_begin;
        IndJ_end[i]   = ind_end + 1;
        F[i] = (sum > 0.0) ? 1.0 / sum : 1.0;
    }

    vdLn(NK, F, tempvecNK_3);
    int maxind = cblas_idamax(NK, tempvecNK_3, 1);
    maxdif = fabs(tempvecNK_3[maxind] / k);
}

/* Absorb F,G into f,g and reset to 1. */
void absorbtion(int k)
{
    vdLn(NK, F, tempvecNK_1);
    cblas_dscal(NK, 1.0 / (2 * k), tempvecNK_1, 1);
    vdAdd(NK, f, tempvecNK_1, tempvecNK_2);
    swap(f, tempvecNK_2);

    vdLn(NK, G, tempvecNK_1);
    cblas_dscal(NK, 1.0 / (2 * k), tempvecNK_1, 1);
    vdAdd(NK, g, tempvecNK_1, tempvecNK_2);
    swap(g, tempvecNK_2);

    fill_n(F, NK, 1.0);
    fill_n(G, NK, 1.0);
}

/* ------------------------------------------------------------------
 * Identity Sinkhorn steps (damped, for source marginal)
 * Matches C++ Sinkhorn_identity_F_axb.
 * ------------------------------------------------------------------ */
void Sinkhorn_identity_F_axb(int k)
{
    maxdif = -1.0;
    vdMul(NK, F_id, p, tempvecNK_1);   /* tempvecNK_1 = F_id * p */
    for (int i = 0; i < NK; i++) {
        if (p[i] == 0.0) continue;
        double sum = 0.0;
        /* No threshold contraction for identity steps — iterate all j to
         * match Python's full logsumexp (threshold is too aggressive at
         * k=152 for NK=381, excluding many valid pairs).             */
        for (int j = 0; j < NK; j++) {
            if (p[j] == 0.0) continue;
            double temp = -k * (Cost_Func(x[i], y[j]) - f_id[i] - f_id[j]);
            sum += exp(temp) * tempvecNK_1[j];
        }
        F_id[i] = (sum > 0.0) ? 1.0 / sum : 1.0;
    }
    vdLn(NK, F_id, tempvecNK_3);
    int maxind = cblas_idamax(NK, tempvecNK_3, 1);
    maxdif = fabs(tempvecNK_3[maxind] / k);
}

void absorbtion_f_id(int k)
{
    vdLn(NK, F_id, tempvecNK_1);
    cblas_dscal(NK, 1.0 / (2 * k), tempvecNK_1, 1);
    vdAdd(NK, f_id, tempvecNK_1, tempvecNK_2);
    swap(f_id, tempvecNK_2);
    fill_n(F_id, NK, 1.0);
}

/* Identity Sinkhorn step for target marginal. */
void Sinkhorn_identity_G_axb(int k)
{
    maxdif = -1.0;
    vdMul(NK, G_id, q, tempvecNK_1);   /* tempvecNK_1 = G_id * q */
    for (int j = 0; j < NK; j++) {
        if (q[j] == 0.0) continue;
        double sum = 0.0;
        /* No threshold contraction — iterate all i to match Python. */
        for (int i = 0; i < NK; i++) {
            if (q[i] == 0.0) continue;
            double temp = -k * (Cost_Func(x[i], y[j]) - g_id[i] - g_id[j]);
            sum += exp(temp) * tempvecNK_1[i];
        }
        G_id[j] = (sum > 0.0) ? 1.0 / sum : 1.0;
    }
    vdLn(NK, G_id, tempvecNK_3);
    int maxind = cblas_idamax(NK, tempvecNK_3, 1);
    maxdif = fabs(tempvecNK_3[maxind] / k);
}

void absorbtion_g_id(int k)
{
    vdLn(NK, G_id, tempvecNK_1);
    cblas_dscal(NK, 1.0 / (2 * k), tempvecNK_1, 1);
    vdAdd(NK, g_id, tempvecNK_1, tempvecNK_2);
    swap(g_id, tempvecNK_2);
    fill_n(G_id, NK, 1.0);
}

/* ------------------------------------------------------------------
 * Helper: print a 1-D double array to a file
 * ------------------------------------------------------------------ */
void print_vec(FILE* fp, double v[], int n)
{
    fprintf(fp, "%d\n", n);
    for (int i = 0; i < n; i++)
        fprintf(fp, "%.*e\n", DECIMAL_DIG, v[i]);
    fprintf(fp, "\n");
}

/* Helper: print a 2-D double array (n rows, dim cols) */
void print_mat(FILE* fp, double v[][dim], int n)
{
    fprintf(fp, "%d\n", n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++)
            fprintf(fp, "%.*e ", DECIMAL_DIG, v[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

/* ------------------------------------------------------------------
 * Main algorithm: Sinkhorn divergence (no warm-start, NK=381)
 * Matches do_sinkhorn_subtracted_axb but skips UseSmall.
 * ------------------------------------------------------------------ */
void do_sinkhorn_subtracted_axb(int k)
{
    printf("NK=%d, Reg_Param=%d\n", NK, k);

    clock_t t;
    int i;

    /* --- Multi-scale loop: k_small → k_final ---
     * When NK=381, k_small = 8*getk(381) = 8*19 = 152 = k_final,
     * so this loop body never executes. */
    int regvariable = multiplier * getk(NK_small);   /* 152 */
    i = 0;
    while (regvariable < k) {
        t = clock();
        Sinkhorn_axb(regvariable);
        absorbtion(regvariable);
        t = clock() - t;
        printf("  iter %d, k=%d, %.3fs, maxdif=%e\n",
               ++i, regvariable, (float)t / CLOCKS_PER_SEC, maxdif);
        regvariable += (int)pow((double)k, 1.0 / 3.0);
    }

    /* --- Final loop at k_final ---
     * maxdif is -1 initially; set sentinel so we enter the loop. */
    printf("Final Sinkhorn at k=%d:\n", k);
    maxdif = cap_treshold + 1.0;
    i = 0;
    while (maxdif > cap_treshold) {
        t = clock();
        Sinkhorn_axb(k);
        absorbtion(k);
        t = clock() - t;
        printf("  iter %d, %.3fs, maxdif=%e\n",
               ++i, (float)t / CLOCKS_PER_SEC, maxdif);
        if (i >= cap_iteration) break;
    }
    printf("Final: %d iterations, last change=%e\n", i, maxdif);

    /* --- Identity F loop --- */
    printf("\nIdentity Sinkhorn for source (f_id):\n");
    regvariable = 1;
    i = 0;
    while (regvariable < k) {
        t = clock();
        Sinkhorn_identity_F_axb(regvariable);
        absorbtion_f_id(regvariable);
        t = clock() - t;
        printf("  iter %d, k=%d, %.3fs, maxdif=%e\n",
               ++i, regvariable, (float)t / CLOCKS_PER_SEC, maxdif);
        regvariable += (int)sqrt((double)k);
    }
    maxdif = cap_treshold + 1.0;
    i = 0;
    while (maxdif > cap_treshold) {
        Sinkhorn_identity_F_axb(k);
        absorbtion_f_id(k);
        i++;
        if (i >= cap_iteration) break;
    }
    printf("Identity F: %d final iterations, last change=%e\n", i, maxdif);

    /* Reset IndId bounds for G loop. */
    fill_n(IndId_begin, NK, 0);
    fill_n(IndId_end,   NK, NK);

    /* --- Identity G loop --- */
    printf("\nIdentity Sinkhorn for target (g_id):\n");
    regvariable = 1;
    i = 0;
    while (regvariable < k) {
        t = clock();
        Sinkhorn_identity_G_axb(regvariable);
        absorbtion_g_id(regvariable);
        t = clock() - t;
        printf("  iter %d, k=%d, %.3fs, maxdif=%e\n",
               ++i, regvariable, (float)t / CLOCKS_PER_SEC, maxdif);
        regvariable += (int)sqrt((double)k);
    }
    maxdif = cap_treshold + 1.0;
    i = 0;
    while (maxdif > cap_treshold) {
        Sinkhorn_identity_G_axb(k);
        absorbtion_g_id(k);
        i++;
        if (i >= cap_iteration) break;
    }
    printf("Identity G: %d final iterations, last change=%e\n", i, maxdif);
}

/* ------------------------------------------------------------------
 * Normalise, subtract identity terms, build reflector surface R/Ref.
 * Matches sink_to_Reflector_subtracted_axb.
 * ------------------------------------------------------------------ */
void sink_to_Reflector(int k)
{
    /* Total transport cost (informational) */
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < NK; i++) if (p[i] > 0) sum1 += p[i] * (f[i] - f_id[i]);
    for (int j = 0; j < NK; j++) if (q[j] > 0) sum2 += q[j] * (g[j] - g_id[j]);
    printf("\nTotal cost: %e + %e = %e\n", sum1, sum2, sum1 + sum2);

    /* Shift f_id so max=0, compensate in g_id */
    double mx = -1e300;
    for (int i = 0; i < NK; i++) if (f_id[i] > mx) mx = f_id[i];
    for (int i = 0; i < NK; i++) { f_id[i] -= mx; g_id[i] += mx; }

    /* Shift f so max=0, compensate in g */
    mx = -1e300;
    for (int i = 0; i < NK; i++) if (f[i] > mx) mx = f[i];
    for (int i = 0; i < NK; i++) { f[i] -= mx; g[i] += mx; }

    /* Subtract identity */
    for (int i = 0; i < NK; i++) { f[i] -= f_id[i]; g[i] -= g_id[i]; }

    /* Build reflector: R = exp(f),  Ref = 2 * x * R */
    vdExp(NK, f, R);
    for (int i = 0; i < NK; i++)
        for (int j = 0; j < dim; j++)
            Ref[i][j] = 2.0 * x[i][j] * R[i];
}

/* ------------------------------------------------------------------
 * Write results to output_fast/
 * ------------------------------------------------------------------ */
void Printall()
{
    FILE* fp;
    string path;

    auto open_or_die = [&](const char* name) -> FILE* {
        path = outputname + "/" + name;
        FILE* h = fopen(path.c_str(), "w");
        if (!h) { perror(path.c_str()); exit(1); }
        return h;
    };

    fp = open_or_die("f.txt");       print_vec(fp, f,    NK); fclose(fp);
    fp = open_or_die("g.txt");       print_vec(fp, g,    NK); fclose(fp);
    fp = open_or_die("f_id.txt");    print_vec(fp, f_id, NK); fclose(fp);
    fp = open_or_die("g_id.txt");    print_vec(fp, g_id, NK); fclose(fp);
    fp = open_or_die("R.txt");       print_vec(fp, R,    NK); fclose(fp);
    fp = open_or_die("Ref.txt");     print_mat(fp, Ref,  NK); fclose(fp);

    printf("\nResults written to %s/\n", outputname.c_str());
}

/* ------------------------------------------------------------------
 * main
 * ------------------------------------------------------------------ */
int main()
{
    clock_t total_time = clock();

    /* Create output directory (POSIX mkdir; -p silences "already exists") */
    system("mkdir -p output_fast");

    fillall();

    Reg_Param = multiplier * getk(NK);   /* 8 * 19 = 152 */
    printf("Regularisation parameter: %g\n", Reg_Param);

    do_sinkhorn_subtracted_axb((int)Reg_Param);
    sink_to_Reflector((int)Reg_Param);
    Printall();

    total_time = clock() - total_time;
    printf("Total time: %.3f s\n", (float)total_time / CLOCKS_PER_SEC);

    return 0;
}
