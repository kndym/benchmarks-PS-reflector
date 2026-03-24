/*
 * main_compare.cpp — Sinkhorn reflector comparison harness
 *
 * Generates NK_COMPARE Halton QMC points (same algorithm as Python gen_cloud),
 * runs the full Sinkhorn divergence pipeline (no warm-start, f=g=0 init),
 * and writes results to output_cpp_NK{N}/
 *
 * Compile for NK=300:
 *   g++ -O2 -DNK_COMPARE=300  -o main_compare_300  main_compare.cpp -lm
 * Compile for NK=1600:
 *   g++ -O2 -DNK_COMPARE=1600 -o main_compare_1600 main_compare.cpp -lm
 *
 * Run:
 *   ./main_compare_300   -> output_cpp_NK300/
 *   ./main_compare_1600  -> output_cpp_NK1600/
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
using namespace std;

#ifndef NK_COMPARE
#  define NK_COMPARE 300
#endif

/* ------------------------------------------------------------------
 * MKL replacements
 * ------------------------------------------------------------------ */
inline void vdMul(int n, const double* a, const double* b, double* c) {
    for (int i=0;i<n;i++) c[i]=a[i]*b[i];
}
inline void vdLn(int n, const double* a, double* c) {
    for (int i=0;i<n;i++) c[i]=log(a[i]);
}
inline void vdAdd(int n, const double* a, const double* b, double* c) {
    for (int i=0;i<n;i++) c[i]=a[i]+b[i];
}
inline void vdExp(int n, const double* a, double* c) {
    for (int i=0;i<n;i++) c[i]=exp(a[i]);
}
inline int cblas_idamax(int n, const double* a, int) {
    int idx=0; double mx=fabs(a[0]);
    for (int i=1;i<n;i++) if (fabs(a[i])>mx){mx=fabs(a[i]);idx=i;}
    return idx;
}
inline void cblas_dscal(int n, double s, double* a, int) {
    for (int i=0;i<n;i++) a[i]*=s;
}

/* ------------------------------------------------------------------
 * Halton QMC — matches Python _halton() and gen_cloud()
 * ------------------------------------------------------------------ */
static double halton(int index, int base) {
    double result=0.0, f=1.0;
    int i=index;
    while (i>0) { f/=base; result+=f*(i%base); i/=base; }
    return result;
}

static void sphere_pt(double X, double Y, double pt[3], int upper) {
    double N2=X*X+Y*Y, d=1.0+N2;
    pt[0]=2.0*X/d;
    pt[1]=2.0*Y/d;
    pt[2]=(upper?1.0:-1.0)*(1.0-N2)/d;
}

static void gen_cloud(int n, double pts[][3], int upper, int skip) {
    const double half=0.6;
    int idx=skip, cnt=0;
    while (cnt<n) {
        double X=(halton(idx,2)-0.5)*2.0*half;
        double Y=(halton(idx,3)-0.5)*2.0*half;
        sphere_pt(X,Y,pts[cnt],upper);
        cnt++; idx++;
    }
}

/* ------------------------------------------------------------------
 * Problem size
 * ------------------------------------------------------------------ */
const int dim = 3;
const int NK  = NK_COMPARE;

int getk(int Nk) { return (int)sqrt((double)Nk); }

/* ------------------------------------------------------------------
 * Density functions — SquareToCircle (same as main_fast.cpp)
 * ------------------------------------------------------------------ */
double P(double v[]) {
    double p1=v[0]/(1.0+v[2]), p2=v[1]/(1.0+v[2]);
    return (-0.5<p1&&p1<0.5&&-0.5<p2&&p2<0.5)?1.0:0.0;
}
double Q(double v[]) {
    double p1=v[0]/(1.0-v[2]), p2=v[1]/(1.0-v[2]);
    return (p1*p1+p2*p2<=0.25)?1.0:0.0;
}
double Cost_Func(double a[], double b[]) {
    double t=1.0;
    for (int d=0;d<dim;d++) t-=a[d]*b[d];
    return -log(t);
}

/* ------------------------------------------------------------------
 * Working arrays
 * ------------------------------------------------------------------ */
double x[NK][dim], y[NK][dim];
double p[NK], q[NK], logp[NK], logq[NK];
double F[NK], G[NK], f[NK], g[NK];
double F_id[NK], G_id[NK], f_id[NK], g_id[NK];
double R[NK], Ref[NK][dim];
double tempvecNK_1[NK], tempvecNK_2[NK], tempvecNK_3[NK];

double treshold;
double maxdif=-1.0;
string outputname;

const int    multiplier    = 8;
const double cap_treshold  = 1.e-5;
const int    cap_iteration = 16;

/* ------------------------------------------------------------------
 * fillall — generate clouds + evaluate densities
 * ------------------------------------------------------------------ */
void fillall() {
    gen_cloud(NK, x, 1, 0);   /* source: upper hemisphere, skip=0 */
    gen_cloud(NK, y, 0, 0);   /* target: lower hemisphere, skip=0 */

    double sum=0.0;
    for (int i=0;i<NK;i++) { p[i]=P(x[i]); sum+=p[i]; }
    for (int i=0;i<NK;i++) {
        p[i]/=sum;
        logp[i]=(p[i]>0.0)?log(p[i]):-1.0e300;
    }
    printf("Source support: %d / %d\n",(int)(sum+0.5),NK);

    sum=0.0;
    for (int i=0;i<NK;i++) { q[i]=Q(y[i]); sum+=q[i]; }
    for (int i=0;i<NK;i++) {
        q[i]/=sum;
        logq[i]=(q[i]>0.0)?log(q[i]):-1.0e300;
    }
    printf("Target support: %d / %d\n",(int)(sum+0.5),NK);

    fill_n(F,NK,1.0);  fill_n(G,NK,1.0);
    fill_n(f,NK,0.0);  fill_n(g,NK,0.0);
    fill_n(F_id,NK,1.0); fill_n(G_id,NK,1.0);
    fill_n(f_id,NK,0.0); fill_n(g_id,NK,0.0);

    treshold = log(1.e-6/(double)NK);
}

/* ------------------------------------------------------------------
 * Sinkhorn step (main OT)
 * ------------------------------------------------------------------ */
void Sinkhorn_axb(int k) {
    vdMul(NK,F,p,tempvecNK_1);
    for (int j=0;j<NK;j++) {
        if (q[j]==0.0) continue;
        double sum=0.0;
        for (int i=0;i<NK;i++) {
            double temp=-k*(Cost_Func(x[i],y[j])-f[i]-g[j]);
            if (logp[i]+temp<treshold) continue;
            sum+=exp(temp)*tempvecNK_1[i];
        }
        G[j]=(sum>0.0)?1.0/sum:1.0;
    }
    vdMul(NK,G,q,tempvecNK_1);
    for (int i=0;i<NK;i++) {
        if (p[i]==0.0) continue;
        double sum=0.0;
        for (int j=0;j<NK;j++) {
            double temp=-k*(Cost_Func(x[i],y[j])-f[i]-g[j]);
            if (temp+logq[j]<treshold) continue;
            sum+=exp(temp)*tempvecNK_1[j];
        }
        F[i]=(sum>0.0)?1.0/sum:1.0;
    }
    vdLn(NK,F,tempvecNK_3);
    int maxind=cblas_idamax(NK,tempvecNK_3,1);
    maxdif=fabs(tempvecNK_3[maxind]/k);
}

void absorbtion(int k) {
    vdLn(NK,F,tempvecNK_1); cblas_dscal(NK,1.0/(2*k),tempvecNK_1,1);
    vdAdd(NK,f,tempvecNK_1,tempvecNK_2); swap(f,tempvecNK_2);
    vdLn(NK,G,tempvecNK_1); cblas_dscal(NK,1.0/(2*k),tempvecNK_1,1);
    vdAdd(NK,g,tempvecNK_1,tempvecNK_2); swap(g,tempvecNK_2);
    fill_n(F,NK,1.0); fill_n(G,NK,1.0);
}

/* ------------------------------------------------------------------
 * Identity Sinkhorn step (source marginal)
 * No threshold — full sum over supported j, matching Python logsumexp.
 * maxdif computed only over supported i (p[i]>0), matching Python fix.
 * ------------------------------------------------------------------ */
void Sinkhorn_identity_F_axb(int k) {
    maxdif=-1.0;
    vdMul(NK,F_id,p,tempvecNK_1);
    for (int i=0;i<NK;i++) {
        if (p[i]==0.0) continue;
        double sum=0.0;
        for (int j=0;j<NK;j++) {
            if (p[j]==0.0) continue;  /* same support on both sides */
            double temp=-k*(Cost_Func(x[i],y[j])-f_id[i]-f_id[j]);
            sum+=exp(temp)*tempvecNK_1[j];
        }
        F_id[i]=(sum>0.0)?1.0/sum:1.0;
    }
    /* maxdif only over supported points (F_id[unsupported]=1 → log=0) */
    vdLn(NK,F_id,tempvecNK_3);
    maxdif=0.0;
    for (int i=0;i<NK;i++) {
        if (p[i]==0.0) continue;
        double v=fabs(tempvecNK_3[i]/k);
        if (v>maxdif) maxdif=v;
    }
}

void absorbtion_f_id(int k) {
    vdLn(NK,F_id,tempvecNK_1); cblas_dscal(NK,1.0/(2*k),tempvecNK_1,1);
    vdAdd(NK,f_id,tempvecNK_1,tempvecNK_2); swap(f_id,tempvecNK_2);
    fill_n(F_id,NK,1.0);
}

/* Identity Sinkhorn step (target marginal) */
void Sinkhorn_identity_G_axb(int k) {
    maxdif=-1.0;
    vdMul(NK,G_id,q,tempvecNK_1);
    for (int j=0;j<NK;j++) {
        if (q[j]==0.0) continue;
        double sum=0.0;
        for (int i=0;i<NK;i++) {
            if (q[i]==0.0) continue;
            double temp=-k*(Cost_Func(x[i],y[j])-g_id[i]-g_id[j]);
            sum+=exp(temp)*tempvecNK_1[i];
        }
        G_id[j]=(sum>0.0)?1.0/sum:1.0;
    }
    vdLn(NK,G_id,tempvecNK_3);
    maxdif=0.0;
    for (int j=0;j<NK;j++) {
        if (q[j]==0.0) continue;
        double v=fabs(tempvecNK_3[j]/k);
        if (v>maxdif) maxdif=v;
    }
}

void absorbtion_g_id(int k) {
    vdLn(NK,G_id,tempvecNK_1); cblas_dscal(NK,1.0/(2*k),tempvecNK_1,1);
    vdAdd(NK,g_id,tempvecNK_1,tempvecNK_2); swap(g_id,tempvecNK_2);
    fill_n(G_id,NK,1.0);
}

/* ------------------------------------------------------------------
 * Main Sinkhorn divergence loop
 * ------------------------------------------------------------------ */
void do_sinkhorn(int k) {
    printf("NK=%d, k_final=%d\n", NK, k);
    clock_t t;
    int i;

    /* Multi-scale: k_small → k_final */
    int regvar = multiplier * getk(NK);  /* same k_final, so loop won't run if NK_small==NK */
    /* For general NK we start from k=1 multi-scale like generate_results.py */
    /* Actually: k_small = 8*sqrt(NK), which equals k_final here.
       So multi-scale body never runs — same as original main_fast.cpp. */
    i=0;
    int k_step = (int)round(pow((double)k, 1.0/3.0));
    int k_ms   = multiplier * getk(NK);  /* start = k_final (no multi-scale) */
    while (k_ms < k) {
        Sinkhorn_axb(k_ms); absorbtion(k_ms);
        printf("  multiscale iter %d, k=%d, maxdif=%e\n", ++i, k_ms, maxdif);
        k_ms += k_step;
    }

    /* Final loop */
    printf("Final Sinkhorn at k=%d:\n", k);
    maxdif = cap_treshold + 1.0; i=0;
    while (maxdif > cap_treshold) {
        t=clock(); Sinkhorn_axb(k); absorbtion(k); t=clock()-t;
        printf("  iter %d, %.3fs, maxdif=%e\n", ++i, (float)t/CLOCKS_PER_SEC, maxdif);
        if (i >= cap_iteration) break;
    }
    printf("Final: %d iters, last maxdif=%e\n", i, maxdif);

    /* Identity F */
    printf("\nIdentity F:\n");
    int id_step = (int)sqrt((double)k);
    regvar=1; i=0;
    while (regvar < k) {
        Sinkhorn_identity_F_axb(regvar); absorbtion_f_id(regvar);
        printf("  iter %d, k=%d, maxdif=%e\n", ++i, regvar, maxdif);
        regvar += id_step;
    }
    maxdif = cap_treshold + 1.0; i=0;
    while (maxdif > cap_treshold) {
        Sinkhorn_identity_F_axb(k); absorbtion_f_id(k); i++;
        if (i >= cap_iteration) break;
    }
    printf("Identity F: %d final iters, last maxdif=%e\n", i, maxdif);

    /* Identity G */
    printf("\nIdentity G:\n");
    regvar=1; i=0;
    while (regvar < k) {
        Sinkhorn_identity_G_axb(regvar); absorbtion_g_id(regvar);
        printf("  iter %d, k=%d, maxdif=%e\n", ++i, regvar, maxdif);
        regvar += id_step;
    }
    maxdif = cap_treshold + 1.0; i=0;
    while (maxdif > cap_treshold) {
        Sinkhorn_identity_G_axb(k); absorbtion_g_id(k); i++;
        if (i >= cap_iteration) break;
    }
    printf("Identity G: %d final iters, last maxdif=%e\n", i, maxdif);
}

/* ------------------------------------------------------------------
 * Normalise and build reflector (matching main_fast.cpp sink_to_Reflector)
 * ------------------------------------------------------------------ */
void sink_to_Reflector(int k) {
    double sum1=0.0, sum2=0.0;
    for (int i=0;i<NK;i++) if (p[i]>0) sum1+=p[i]*(f[i]-f_id[i]);
    for (int j=0;j<NK;j++) if (q[j]>0) sum2+=q[j]*(g[j]-g_id[j]);
    printf("\nTotal cost: %e\n", sum1+sum2);

    double mx=-1e300;
    for (int i=0;i<NK;i++) if (f_id[i]>mx) mx=f_id[i];
    for (int i=0;i<NK;i++) { f_id[i]-=mx; g_id[i]+=mx; }

    mx=-1e300;
    for (int i=0;i<NK;i++) if (f[i]>mx) mx=f[i];
    for (int i=0;i<NK;i++) { f[i]-=mx; g[i]+=mx; }

    for (int i=0;i<NK;i++) { f[i]-=f_id[i]; g[i]-=g_id[i]; }

    vdExp(NK,f,R);
    for (int i=0;i<NK;i++)
        for (int j=0;j<dim;j++)
            Ref[i][j]=2.0*x[i][j]*R[i];

    double rmin=R[0], rmax=R[0], rmean=0.0;
    for (int i=0;i<NK;i++) if (p[i]>0) {
        if (R[i]<rmin) rmin=R[i];
        if (R[i]>rmax) rmax=R[i];
        rmean+=R[i];
    }
    int ns=0; for (int i=0;i<NK;i++) if(p[i]>0) ns++;
    printf("R (supported): min=%.4f, max=%.4f, mean=%.4f\n",
           rmin, rmax, rmean/ns);
}

/* ------------------------------------------------------------------
 * Output
 * ------------------------------------------------------------------ */
void print_vec(FILE* fp, double v[], int n) {
    fprintf(fp,"%d\n",n);
    for (int i=0;i<n;i++) fprintf(fp,"%.17e\n",v[i]);
}
void print_mat(FILE* fp, double v[][dim], int n) {
    fprintf(fp,"%d\n",n);
    for (int i=0;i<n;i++) {
        for (int j=0;j<dim;j++) fprintf(fp,"%.17e ",v[i][j]);
        fprintf(fp,"\n");
    }
}

void Printall() {
    FILE* fp;
    auto W=[&](const char* name, auto fn){
        string path=outputname+"/"+name;
        fp=fopen(path.c_str(),"w");
        if (!fp){perror(path.c_str());exit(1);}
        fn(fp); fclose(fp);
    };
    W("f.txt",    [](FILE*fp){print_vec(fp,f,   NK);});
    W("g.txt",    [](FILE*fp){print_vec(fp,g,   NK);});
    W("f_id.txt", [](FILE*fp){print_vec(fp,f_id,NK);});
    W("g_id.txt", [](FILE*fp){print_vec(fp,g_id,NK);});
    W("R.txt",    [](FILE*fp){print_vec(fp,R,   NK);});
    W("Ref.txt",  [](FILE*fp){print_mat(fp,Ref, NK);});
    printf("Results written to %s/\n", outputname.c_str());
}

/* ------------------------------------------------------------------
 * main
 * ------------------------------------------------------------------ */
int main() {
    clock_t t0=clock();

    char buf[64];
    snprintf(buf, sizeof(buf), "output_cpp_NK%d", NK);
    outputname = buf;

    char cmd[128];
    snprintf(cmd, sizeof(cmd), "mkdir -p %s", outputname.c_str());
    system(cmd);

    fillall();

    int k_final = multiplier * getk(NK);
    printf("k_final = %d\n", k_final);

    do_sinkhorn(k_final);
    sink_to_Reflector(k_final);
    Printall();

    printf("Total time: %.3fs\n", (float)(clock()-t0)/CLOCKS_PER_SEC);
    return 0;
}
