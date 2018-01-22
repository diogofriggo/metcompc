#include <stdio.h>
#include <math.h>

#define NDIM 3

#if NDIM == 3

#define VSet(v, sx, sy, sz) \
    (v).x = sx,             \
    (v).y = sy,             \
    (v).z = sz

#define VCopy(v1, v2)   \
    (v1).x = (v2).x,    \
    (v1).y = (v2).y,    \
    (v1).z = (v2).z

#define VScale(v, s)    \
    (v).x *= s,         \
    (v).y *= s,         \
    (v).z *= s

#define VAdd(v1, v2, v3)        \
    (v1).x = (v2).x + (v3).x,   \
    (v1).y = (v2).y + (v3).y,   \
    (v1).z = (v2).z + (v3).z

#define VSub(v1, v2, v3)        \
    (v1).x = (v2).x - (v3).x,   \
    (v1).y = (v2).y - (v3).y,   \
    (v1).z = (v2).z - (v3).z

#define VSAdd(v1, v2, s3, v3)       \
    (v1).x = (v2).x + (s3)*(v3).x,  \
    (v1).y = (v2).y + (s3)*(v3).y,  \
    (v1).y = (v2).z + (s3)*(v3).z

#define VSSAdd(v1, s2, v2, s3, v3)      \
    (v1).x = (s2)*(v2).x + (s3)*(v3).x, \
    (v1).y = (s2)*(v2).y + (s3)*(v3).y, \
    (v1).y = (s2)*(v2).z + (s3)*(v3).z

#define VMul(v1, v2, v3) (v1).x =   \
    (v2).x * (v3).x,                \
    (v1).y = (v2).y * (v3).y,       \
    (v1).z = (v2).z * (v3).z

#define VDiv(v1, v2, v3)        \
    (v1).x = (v2).x / (v3).x,   \
    (v1).y = (v2).y / (v3).y,   \
    (v1).z = (v2).z / (v3).z

#define VDot(v1, v2) ((v1).x*(v2).x + (v1).y*(v2).y + (v1).z*(v2).z)

#define VCross(v1, v2, v3)                      \
    (v1).x = (v2).y * (v3).z - (v2).z * (v3).y, \
    (v1).y = (v2).z * (v3).x - (v2).x * (v3).z, \
    (v1).z = (v2).x * (v3).y - (v2).y * (v3).x

#define VWDot(v1, v2, v3)       \
    (v1).x * (v2).x * (v3).x +  \
    (v1).y * (v2).y * (v3).y +  \
    (v1).z * (v2).z * (v3).z

#define MVMul(v1, m, v2)                                    \
    (v1).x = (m)[0]*(v2).x + (m)[3]*(v2).x + (m)[6]*(v2).x, \
    (v1).y = (m)[1]*(v2).y + (m)[4]*(v2).y + (m)[7]*(v2).y, \
    (v1).z = (m)[2]*(v2).z + (m)[5]*(v2).z + (m)[8]*(v2).z

#define MVMulT(v1, m, v2)                                   \
    (v1).x = (m)[0]*(v2).x + (m)[1]*(v2).y + (m)[2]*(v2).z  \
    (v1).y = (m)[3]*(v2).x + (m)[4]*(v2).y + (m)[5]*(v2).z  \
    (v1).z = (m)[6]*(v2).x + (m)[7]*(v2).y + (m)[8]*(v2).z

#define VSetAll(v, s) VSet(v, s, s, s)

#define VAddCon(v1, v2, s) \
    (v1).x = (v2).x + (s), \
    (v1).y = (v2).y + (s), \
    (v1).z = (v2).z + (s)

#define VProd(v) \
    ((v).x * (v).y * (v).z)

#define VGe(v1, v2) ((v1).x >= (v2).x && (v1).y >= (v2).y && (v1).z >= (v2).z)

#define VLt(v1, v2) ((v1).x < (v2).x && (v1).y < (v2).y && (v1).z < (v2).z)

#define VLinear(p, s) (((p).z*(s).y + (p).y) * (s).x + (p).x)

#define VCSum(v) ((v).x + (v).y + (v).z)

#define VComp(v, k) *((k == 0) ? &(v).x : ((k == 1) ? &(v).y : &(v).z))

#define VToLin(a, n, v) \
a[(n) + 0] = (v).x, \
a[(n) + 1] = (v).y, \
a[(n) + 2] = (v).z

#define VFromLin(v, a, n) VSet(v, a[(n) + 0], a[(n) + 1], a[(n) + 2])

#elif NDIM == 2

#define VCross(v1, v2) (v1).x*(v2).y  - (v1).y*(v2).x

#endif

#define VZero(v) VSetAll(v, 0)

#define VLenSq(v) VDot(v, v)

#define VWLenSq(v1, v2) VWDot(v1, v2, v2)

#define VLen(v) sqrt(VDot(v, v))

#define VVAdd(v1, v2) VAdd(v1, v1, v2)

#define VVSub(v1, v2) VSub(v1, v1, v2)

#define VVSAdd(v1, s2, v2) VSAdd(v1, v1, s2, v2)

#define VInterp(v1, s2, v2, v3) VVSAdd(v1, s2, v2, 1. - (s2), v3)

//

#define MAT(a, n, i, j) (a)[(i) + n * (j)]

#define M3(a, i, j) MAT(a, 3, i, j)

#define M6(a, i, j) MAT(a, 6, i, j)

#define DO(m, n) for(m = 0; m < n; m++)

//

#define Sqr(x) ((x) * (x))
//#define Cube(x) ((x) * (x) * (x))

typedef double real;

typedef struct {
    real x, y, z;
} VecR;

typedef struct {
    real x, y, z;
} VecI;

typedef struct {
    real u[9];
} RMat;

typedef struct {
    RMat rMatT;
    VecR r, rv, omega, omegah, bV, cV, hV;
    real inertiaM[9], fV[6], gV[6], xV[6], yV[6], mass, s, sv, svh, sa, torq;
} Link;

typedef struct {
    Link *L;
    VecR ra, wa;
    int nLink;
} Poly;

typedef struct {
    VecR f, r;
} Site;

typedef struct {
    real val, sum, sum2;
} Prop;

void ComputeLinkCoordsVels();
void BuildLinkRotmatT(RMat *rMat, real dihedA, real bondA);
void BuildLinkPhimatT(real *phiT, int k);
void MulMat(real uk[9], real ukm1[9], real u[9], real a);
void MulMatVec(real vsp[6], real phiT[36], real vs[6], real a);
void BuildLinkInertiaMats();
void ComputeLinkForces();
void BuildLinkXYVecs(int k);
    
Poly P;
Site *site;
real bondAng, bondLen, kinEnVal, totEnVal, twistAng, uCon;
int chainLen, helixPeriod, nDof, nSite;

//Mol *mol
VecR region, vSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum;
int moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit;

int main() {
    bondLen = 1.3;
    chainLen = 80;
    helixPeriod = 8; //6
    uCon = 1.; //5.;
    deltaT = 0.001; //0.004;
    //region 24. 24. 24.;
    stepAvg = 2000;
    stepEquil = 10000;
    stepLimit = 100000;
    temperature = 2.;
    //stepAdjustTemp = 2000;
    //stepReduceTemp = 4000;
    //stepSnap = 10000;
    //tempFinal = 0.001;
    //tempInit = 4.;
    //tempReduceFac = 0.97;
    
    return 0;
}

void ComputeLinkCoordsVels()
{
    RMat rMat;
    VecR bEx, bVp, hVp;
    real phiT[36], vs[6], vsp[6];
    int k;
    
    VSet (bVp, 0., 0., bondLen);
    VSet (hVp, 0., 0., 1.);
    for (k = 0; k < P.nLink; k++) {
        if(k > 0) {
            MVMul(P.L[k].hV, P.L[k-1].rMatT.u, hVp);
            BuildLinkRotmatT(&rMat, P.L[k].s, bondAng);
            MulMat(P.L[k].rMatT.u, P.L[k-1].rMatT.u, rMat.u, 3);
        }
        MVMul(P.L[k].bV, P.L[k].rMatT.u, bVp);
    }
    for (k = 0; k < P.nLink; k++) {
        VToLin(vs, 0, P.L[k].omega);
        VToLin(vs, 3, P.L[k].rv);
        BuildLinkPhimatT(phiT, k);
        MulMatVec(vsp, phiT, vs, 6);
        if(k < P.nLink - 1){
            VFromLin(P.L[k+1].omega, vsp, 0);
            VVSAdd(P.L[k+1].omega, P.L[k+1].sv, P.L[k+1].hV);
        }
        VFromLin(P.L[k+1].rv, vsp, 3);
    }
    for (k = 0; k < P.nLink; k++)
        VAdd(P.L[k+1].r, P.L[k].r, P.L[k].bV);
    for (k = 0; k < P.nLink+1; k++)
        site[k+1].r = P.L[k].r;
    VSet(bEx, 0., -sin(bondAng), -cos(bondAng));
    VScale(bEx, bondLen);
    MVMul(site[0].r, P.L[0].rMatT.u, bEx);
    VVAdd(site[0].r, site[1].r);
}

void BuildLinkRotmatT(RMat *rMat, real dihedA, real bondA){
    real cb, cd, sb, sd;
    cb = cos(bondA);
    sb = sin(bondA);
    cd = cos(dihedA);
    sd = sin(dihedA);
    M3(rMat->u, 0, 0) = cd;
    M3(rMat->u, 1, 0) = sd;
    M3(rMat->u, 2, 0) = 0.;
    
    M3(rMat->u, 0, 1) = sb * cb;
    M3(rMat->u, 1, 1) = cd * cb;
    M3(rMat->u, 2, 1) = sb;
    
    M3(rMat->u, 0, 2) = sd * sb;
    M3(rMat->u, 1, 2) = cd * sb;
    M3(rMat->u, 2, 2) = cb;
}

void BuildLinkPhimatT(real *phiT, int k){
    int i, j;
    DO(i, 6){
        DO(j, 6) M6(phiT, i, j) = (i == j) ? 1. : 0.;
    }
    M6(phiT, 3, 1) = P.L[k].bV.z;
    M6(phiT, 3, 2) = -P.L[k].bV.y;
    M6(phiT, 4, 0) = -P.L[k].bV.z;
    M6(phiT, 4, 2) = P.L[k].bV.x;
    M6(phiT, 5, 0) = P.L[k].bV.y;
    M6(phiT, 5, 1) = -P.L[k].bV.x;
}

void MulMat(real uk[9], real ukm1[9], real u[9], real a){

}

void MulMatVec(real vsp[6], real phiT[36], real vs[6], real a){
    
}

void BuildLinkInertiaMats() {
    VecR d;
    real dd, iBall, inertiaK;
    int k;
    
    inertiaK = 0.1;
    for(k = 0; k < P.nLink; k++) {
        if(k > 0) {
            P.L[k].mass = 1.;
            VSub(P.L[k].cV, site[k+2].r, site[k+1].r);
        }
        else{
            P.L[k].mass = 3.;
            VAdd(P.L[k].cV, site[2].r, site[1].r);
            VVAdd(P.L[k].cV, site[0].r);
            VScale(P.L[k].cV, 1./3.);
            VVSub(P.L[k].cV, site[1].r);
        }
        iBall = inertiaK * P.L[k].mass;
        VSub(d, site[k+2].r, site[k+1].r);
        dd = VLenSq(d);
        M3(P.L[k].inertiaM, 0, 0) = dd - Sqr(d.x) + iBall;
        M3(P.L[k].inertiaM, 1, 1) = dd - Sqr(d.y) + iBall;
        M3(P.L[k].inertiaM, 2, 2) = dd - Sqr(d.z) + iBall;
        M3(P.L[k].inertiaM, 0, 1) = -d.x * d.y;
        M3(P.L[k].inertiaM, 0, 2) = -d.x * d.z;
        M3(P.L[k].inertiaM, 1, 2) = -d.y * d.z;
        if(k == 0){
            VSub(d, site[0].r, site[1].r);
            M3(P.L[k].inertiaM, 0, 0) += dd - Sqr(d.x);
            M3(P.L[k].inertiaM, 0, 1) -= -d.x * d.y;
            //... (similarly for other matrix elements) ...
        }
        M3(P.L[k].inertiaM, 1, 0) = M3(P.L[k].inertiaM, 0, 1);
        M3(P.L[k].inertiaM, 2, 0) = M3(P.L[k].inertiaM, 0, 2);
        M3(P.L[k].inertiaM, 2, 1) = M3(P.L[k].inertiaM, 1, 2);
    }
}

void ComputeLinkForces(){
    VecR d, fc, tq, tq1;
    real ang;
    int k;
    
    for(k = 1; k < P.nLink; k++) {
        ang = P.L[k].s - twistAng;
        P.L[k].torq = -uCon * sin(ang);
        uSum -= uCon * cos(ang);
    }
    for(k = 0; k < P.nLink; k++) {
        fc = site[k+2].f;
        VCross(tq, P.L[k].bV, fc);
        if(k==0){
            VVAdd(fc, site[1].f);
            VVAdd(fc, site[0].f);
            VSub(d, site[0].r, site[1].r);
            VCross(tq1, d, site[0].f);
            VVAdd(tq, tq1);
        }
        VToLin(P.L[k].fV, 0, tq);
        VToLin(P.L[k].fV, 3, fc);
    }
}

void BuildLinkXYVecs(int k){
    VecR dv, w, w1, w2;
    int i;
    if(k > 0){
        VCross(w, P.L[k-1].omega, P.L[k].hV);
        VScale(w, P.L[k].sv);
        VToLin(P.L[k].xV, 0, w);
        VSub(dv, P.L[k].rv, P.L[k-1].rv);
        VCross(w, P.L[k-1].omega, dv);
        VToLin(P.L[k].xV, 3, w);
    }
    else{
        DO(i,6) P.L[k].xV[i] = 0.;
    }
    MVMul(w, P.L[k].inertiaM, P.L[k].omega);
    VCross(w1, P.L[k].omega, w);
    VToLin(P.L[k].yV, 0, w1);
    VCross(w, P.L[k].omega, P.L[k].cV);
    VCross(w2, P.L[k].omega, w);
    VScale(w2, P.L[k].mass);
    VToLin(P.L[k].yV, 3, w2);
}

void ComputeLinkAccels(){
    real as[6], asp[6], h[3];
    real mMat[36], phi[36], phiT[36], pMat[36], tMat1[36], tMat2[36];
    real z[6], zp[6], zt[6], dk, e;
    int i, j, k;
    DO(i, 6) zp[i] = 0;
    for(k = P.nLink - 1; k >= 0; k--){
        BuildLinkPhimatT(phiT, k);
        DO(i, 6){
            DO(j, 6) M6(phi, i, j) = M6(phiT, j, i);
        }
        BuildLinkXYVecs(k);
        BuildLinkMMat((k == P.nLink - 1) ? pMat : mMat, k);
        if(k < P.nLink - 1){
            DO(i, 6){
                DO(j, 6) M6(phi, i, j) = M6(phiT, j, i);
            }
        }
    }
}



























































