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
  #define CASADI_PREFIX(ID) f_kin_R_mex_ ## ID
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

static const casadi_int casadi_s0[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const casadi_int casadi_s1[51] = {3, 12, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};

/* f_kin_R:(i0[4])->(o0[3x12]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[0]? arg[0][0] : 0;
  a1=cos(a0);
  if (res[0]!=0) res[0][0]=a1;
  a2=sin(a0);
  if (res[0]!=0) res[0][1]=a2;
  a3=0.;
  if (res[0]!=0) res[0][2]=a3;
  a4=sin(a0);
  a5=(-a4);
  if (res[0]!=0) res[0][3]=a5;
  a0=cos(a0);
  if (res[0]!=0) res[0][4]=a0;
  if (res[0]!=0) res[0][5]=a3;
  if (res[0]!=0) res[0][6]=a3;
  if (res[0]!=0) res[0][7]=a3;
  a5=1.;
  if (res[0]!=0) res[0][8]=a5;
  a6=arg[0]? arg[0][1] : 0;
  a7=cos(a6);
  a8=(a1*a7);
  a9=sin(a6);
  a10=(a4*a9);
  a8=(a8-a10);
  if (res[0]!=0) res[0][9]=a8;
  a7=(a2*a7);
  a9=(a0*a9);
  a7=(a7+a9);
  if (res[0]!=0) res[0][10]=a7;
  if (res[0]!=0) res[0][11]=a3;
  a9=sin(a6);
  a1=(a1*a9);
  a6=cos(a6);
  a4=(a4*a6);
  a1=(a1+a4);
  a4=(-a1);
  if (res[0]!=0) res[0][12]=a4;
  a0=(a0*a6);
  a2=(a2*a9);
  a0=(a0-a2);
  if (res[0]!=0) res[0][13]=a0;
  if (res[0]!=0) res[0][14]=a3;
  if (res[0]!=0) res[0][15]=a3;
  if (res[0]!=0) res[0][16]=a3;
  if (res[0]!=0) res[0][17]=a5;
  if (res[0]!=0) res[0][18]=a8;
  if (res[0]!=0) res[0][19]=a7;
  if (res[0]!=0) res[0][20]=a3;
  if (res[0]!=0) res[0][21]=a4;
  if (res[0]!=0) res[0][22]=a0;
  if (res[0]!=0) res[0][23]=a3;
  if (res[0]!=0) res[0][24]=a3;
  if (res[0]!=0) res[0][25]=a3;
  if (res[0]!=0) res[0][26]=a5;
  a4=arg[0]? arg[0][3] : 0;
  a2=cos(a4);
  a9=(a8*a2);
  a6=sin(a4);
  a10=(a1*a6);
  a9=(a9-a10);
  if (res[0]!=0) res[0][27]=a9;
  a2=(a7*a2);
  a6=(a0*a6);
  a2=(a2+a6);
  if (res[0]!=0) res[0][28]=a2;
  if (res[0]!=0) res[0][29]=a3;
  a2=sin(a4);
  a8=(a8*a2);
  a4=cos(a4);
  a1=(a1*a4);
  a8=(a8+a1);
  a8=(-a8);
  if (res[0]!=0) res[0][30]=a8;
  a0=(a0*a4);
  a7=(a7*a2);
  a0=(a0-a7);
  if (res[0]!=0) res[0][31]=a0;
  if (res[0]!=0) res[0][32]=a3;
  if (res[0]!=0) res[0][33]=a3;
  if (res[0]!=0) res[0][34]=a3;
  if (res[0]!=0) res[0][35]=a5;
  return 0;
}

CASADI_SYMBOL_EXPORT int f_kin_R(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int f_kin_R_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int f_kin_R_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_kin_R_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int f_kin_R_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_kin_R_release(int mem) {
}

CASADI_SYMBOL_EXPORT void f_kin_R_incref(void) {
}

CASADI_SYMBOL_EXPORT void f_kin_R_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int f_kin_R_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int f_kin_R_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real f_kin_R_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_kin_R_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_kin_R_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_kin_R_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_kin_R_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int f_kin_R_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

CASADI_SYMBOL_EXPORT int f_kin_R_work_bytes(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1*sizeof(const casadi_real*);
  if (sz_res) *sz_res = 1*sizeof(casadi_real*);
  if (sz_iw) *sz_iw = 0*sizeof(casadi_int);
  if (sz_w) *sz_w = 0*sizeof(casadi_real);
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_f_kin_R(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  int mem;
  casadi_real w[51];
  casadi_int *iw = 0;
  const casadi_real* arg[1] = {0};
  casadi_real* res[1] = {0};
  if (argc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_kin_R\" failed. Too many input arguments (%d, max 1)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_kin_R\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+40);
  --resc;
  res[0] = w+4;
  f_kin_R_incref();
  mem = f_kin_R_checkout();
  i = f_kin_R(arg, res, iw, w+40, mem);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"f_kin_R\" failed.");
  f_kin_R_release(mem);
  f_kin_R_decref();
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif

casadi_int main_f_kin_R(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[51];
  const casadi_real* arg[1];
  casadi_real* res[1];
  arg[0] = w+0;
  res[0] = w+4;
  a = w;
  for (j=0; j<4; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = f_kin_R(arg, res, iw, w+40, 0);
  if (flag) return flag;
  r = w+4;
  for (j=0; j<36; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[8];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_f_kin_R(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "f_kin_R")==0) {
    mex_f_kin_R(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'f_kin_R'");
}
#endif
CASADI_SYMBOL_EXPORT int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "f_kin_R")==0) {
    return main_f_kin_R(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'f_kin_R'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
