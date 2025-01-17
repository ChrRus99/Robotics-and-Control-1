/* This file was automatically generated by CasADi 3.6.5.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) f_g_vect_mex_ ## ID
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_from_mex CASADI_PREFIX(from_mex)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_to_mex CASADI_PREFIX(to_mex)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

#ifdef MATLAB_MEX_FILE
casadi_real* casadi_from_mex(const mxArray* p, casadi_real* y, const casadi_int* sp, casadi_real* w) {
  casadi_int nrow, ncol, is_sparse, c, k, p_nrow, p_ncol;
  const casadi_int *colind, *row;
  mwIndex *Jc, *Ir;
  const double* p_data;
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError",
      "\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
  p_nrow = mxGetM(p);
  p_ncol = mxGetN(p);
  is_sparse = mxIsSparse(p);
  Jc = 0;
  Ir = 0;
  if (is_sparse) {
    Jc = mxGetJc(p);
    Ir = mxGetIr(p);
  }
  p_data = (const double*)mxGetData(p);
  if (p_nrow==1 && p_ncol==1) {
    casadi_int nnz;
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    nnz = sp[ncol];
    casadi_fill(y, nnz, v);
  } else {
    casadi_int tr = 0;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError",
                                 "\"from_mex\" failed: Dimension mismatch. "
                                 "Expected %d-by-%d, got %d-by-%d instead.",
                                 nrow, ncol, p_nrow, p_ncol);
    }
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

#endif

#define casadi_to_double(x) ((double) x)

#ifdef MATLAB_MEX_FILE
mxArray* casadi_to_mex(const casadi_int* sp, const casadi_real* x) {
  casadi_int nrow, ncol, c, k;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int nnz;
#endif
  const casadi_int *colind, *row;
  mxArray *p;
  double *d;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int i;
  mwIndex *j;
#endif /* CASADI_MEX_NO_SPARSE */
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  nnz = sp[ncol];
  if (nnz!=nrow*ncol) {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = casadi_to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    d = (double*)mxGetData(p);
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = casadi_to_double(*x++);
      }
    }
  }
  return p;
}

#endif

#ifndef CASADI_PRINTF
#ifdef MATLAB_MEX_FILE
  #define CASADI_PRINTF mexPrintf
#else
  #define CASADI_PRINTF printf
#endif
#endif

static const casadi_int casadi_s0[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
static const casadi_int casadi_s1[13] = {1, 6, 0, 0, 1, 2, 3, 4, 4, 0, 0, 0, 0};

/* f_g_vect:(i0[6])->(o0[1x6,4nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[0]? arg[0][1] : 0;
  a1=cos(a0);
  a2=arg[0]? arg[0][2] : 0;
  a3=cos(a2);
  a4=arg[0]? arg[0][3] : 0;
  a5=cos(a4);
  a6=-5.8999748394999996e-01;
  a7=arg[0]? arg[0][4] : 0;
  a8=sin(a7);
  a9=(a6*a8);
  a10=-1.3729310000000000e-01;
  a11=6.1232339957367660e-17;
  a8=(a11*a8);
  a8=(a10*a8);
  a9=(a9+a8);
  a8=(a5*a9);
  a12=2.8991448428249997e+00;
  a13=sin(a4);
  a14=(a12*a13);
  a8=(a8+a14);
  a13=(a11*a13);
  a14=5.8999748394999996e-01;
  a15=cos(a7);
  a16=(a14*a15);
  a17=1.3729310000000000e-01;
  a15=(a11*a15);
  a15=(a17*a15);
  a16=(a16+a15);
  a16=(a16+a17);
  a15=(a13*a16);
  a8=(a8-a15);
  a15=1.0049854919999998e+01;
  a8=(a8+a15);
  a15=(a3*a8);
  a18=sin(a2);
  a19=sin(a4);
  a20=(a19*a9);
  a21=-2.8991448428249997e+00;
  a22=cos(a4);
  a23=(a21*a22);
  a20=(a20+a23);
  a22=(a11*a22);
  a23=(a22*a16);
  a20=(a20+a23);
  a23=(a18*a20);
  a15=(a15-a23);
  a23=8.4680422749999991e+01;
  a24=-5.7230000000000003e-01;
  a25=(a24*a3);
  a25=(a23*a25);
  a15=(a15+a25);
  a25=4.7326892899999997e+01;
  a15=(a15+a25);
  a25=-1.2804562518299997e+02;
  a15=(a15+a25);
  a1=(a1*a15);
  a15=sin(a0);
  a25=cos(a2);
  a26=(a25*a20);
  a27=sin(a2);
  a28=(a27*a8);
  a26=(a26+a28);
  a28=(a24*a27);
  a28=(a23*a28);
  a26=(a26+a28);
  a15=(a15*a26);
  a1=(a1-a15);
  if (res[0]!=0) res[0][0]=a1;
  a1=cos(a2);
  a15=cos(a0);
  a26=(a15*a8);
  a28=(a23*a15);
  a28=(a24*a28);
  a26=(a26+a28);
  a1=(a1*a26);
  a26=cos(a2);
  a0=sin(a0);
  a28=(a0*a20);
  a26=(a26*a28);
  a28=sin(a2);
  a20=(a15*a20);
  a28=(a28*a20);
  a26=(a26+a28);
  a1=(a1-a26);
  a2=sin(a2);
  a8=(a0*a8);
  a23=(a23*a0);
  a24=(a24*a23);
  a8=(a8+a24);
  a2=(a2*a8);
  a1=(a1-a2);
  if (res[0]!=0) res[0][1]=a1;
  a1=cos(a4);
  a25=(a15*a25);
  a18=(a0*a18);
  a25=(a25-a18);
  a18=(a25*a9);
  a1=(a1*a18);
  a18=sin(a4);
  a0=(a0*a3);
  a15=(a15*a27);
  a0=(a0+a15);
  a9=(a0*a9);
  a18=(a18*a9);
  a1=(a1-a18);
  a18=cos(a4);
  a12=(a12*a0);
  a9=(a0*a16);
  a9=(a11*a9);
  a12=(a12-a9);
  a18=(a18*a12);
  a1=(a1+a18);
  a4=sin(a4);
  a21=(a21*a25);
  a16=(a25*a16);
  a16=(a11*a16);
  a21=(a21+a16);
  a4=(a4*a21);
  a1=(a1-a4);
  if (res[0]!=0) res[0][2]=a1;
  a1=cos(a7);
  a5=(a0*a5);
  a19=(a25*a19);
  a5=(a5+a19);
  a6=(a6*a5);
  a10=(a10*a5);
  a10=(a11*a10);
  a6=(a6+a10);
  a1=(a1*a6);
  a7=sin(a7);
  a25=(a25*a22);
  a0=(a0*a13);
  a25=(a25-a0);
  a25=(a25+a11);
  a14=(a14*a25);
  a17=(a17*a25);
  a11=(a11*a17);
  a14=(a14+a11);
  a7=(a7*a14);
  a1=(a1-a7);
  if (res[0]!=0) res[0][3]=a1;
  return 0;
}

CASADI_SYMBOL_EXPORT int f_g_vect(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int f_g_vect_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int f_g_vect_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_g_vect_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int f_g_vect_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_g_vect_release(int mem) {
}

CASADI_SYMBOL_EXPORT void f_g_vect_incref(void) {
}

CASADI_SYMBOL_EXPORT void f_g_vect_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int f_g_vect_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int f_g_vect_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real f_g_vect_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_g_vect_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_g_vect_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_g_vect_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_g_vect_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int f_g_vect_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

CASADI_SYMBOL_EXPORT int f_g_vect_work_bytes(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1*sizeof(const casadi_real*);
  if (sz_res) *sz_res = 1*sizeof(casadi_real*);
  if (sz_iw) *sz_iw = 0*sizeof(casadi_int);
  if (sz_w) *sz_w = 0*sizeof(casadi_real);
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_f_g_vect(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  int mem;
  casadi_real w[39];
  casadi_int *iw = 0;
  const casadi_real* arg[1] = {0};
  casadi_real* res[1] = {0};
  if (argc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_g_vect\" failed. Too many input arguments (%d, max 1)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_g_vect\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+10);
  --resc;
  res[0] = w+6;
  f_g_vect_incref();
  mem = f_g_vect_checkout();
  i = f_g_vect(arg, res, iw, w+10, mem);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_g_vect\" failed.");
  f_g_vect_release(mem);
  f_g_vect_decref();
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif

casadi_int main_f_g_vect(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[39];
  const casadi_real* arg[1];
  casadi_real* res[1];
  arg[0] = w+0;
  res[0] = w+6;
  a = w;
  for (j=0; j<6; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = f_g_vect(arg, res, iw, w+10, 0);
  if (flag) return flag;
  r = w+6;
  for (j=0; j<4; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[9];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_f_g_vect(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "f_g_vect")==0) {
    mex_f_g_vect(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'f_g_vect'");
}
#endif
CASADI_SYMBOL_EXPORT int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "f_g_vect")==0) {
    return main_f_g_vect(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'f_g_vect'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
