#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void cBase2C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cBaseeco(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cBaseecoX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cBaseecoZ(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cBaseRC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cDPeco(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cDPecoX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cEMeco(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void preBaseX(void *, void *, void *, void *, void *, void *, void *);
extern void preDP(void *, void *, void *, void *, void *, void *, void *);
extern void preDPX(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"cBase2C",   (DL_FUNC) &cBase2C,   22},
    {"cBaseeco",  (DL_FUNC) &cBaseeco,  32},
    {"cBaseecoX", (DL_FUNC) &cBaseecoX, 36},
    {"cBaseecoZ", (DL_FUNC) &cBaseecoZ, 29},
    {"cBaseRC",   (DL_FUNC) &cBaseRC,   23},
    {"cDPeco",    (DL_FUNC) &cDPeco,    36},
    {"cDPecoX",   (DL_FUNC) &cDPecoX,   40},
    {"cEMeco",    (DL_FUNC) &cEMeco,    27},
    {"preBaseX",  (DL_FUNC) &preBaseX,   7},
    {"preDP",     (DL_FUNC) &preDP,      7},
    {"preDPX",    (DL_FUNC) &preDPX,     8},
    {NULL, NULL, 0}
};

void R_init_eco(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
