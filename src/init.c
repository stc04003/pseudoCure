#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _pseudoCure_fastDabrowska(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_fastL(SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_fastTau2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_gee(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_pgee(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_pgeeCV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_fastKM(SEXP, SEXP);
extern SEXP _pseudoCure_pseudoKM(SEXP, SEXP);
extern SEXP _pseudoCure_pseudoKM1(SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_pseudoTau(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_pseudoCure_fastDabrowska", (DL_FUNC) &_pseudoCure_fastDabrowska,  7},
    {"_pseudoCure_fastL",         (DL_FUNC) &_pseudoCure_fastL,          4},
    {"_pseudoCure_fastTau2",      (DL_FUNC) &_pseudoCure_fastTau2,       4},
    {"_pseudoCure_gee",           (DL_FUNC) &_pseudoCure_gee,            8},
    {"_pseudoCure_pgee",          (DL_FUNC) &_pseudoCure_pgee,          12},
    {"_pseudoCure_pgeeCV",        (DL_FUNC) &_pseudoCure_pgeeCV,        13},
    {"_pseudoCure_fastKM",        (DL_FUNC) &_pseudoCure_fastKM,       2},
    {"_pseudoCure_pseudoKM",      (DL_FUNC) &_pseudoCure_pseudoKM,       2},
    {"_pseudoCure_pseudoKM1",     (DL_FUNC) &_pseudoCure_pseudoKM1,      3},
    {"_pseudoCure_pseudoTau",     (DL_FUNC) &_pseudoCure_pseudoTau,      4},
    {NULL, NULL, 0}
};

void R_init_pseudoCure(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
