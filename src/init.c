#include <R_ext/Rdynload.h>
#include "rmargint.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_NativePrimitiveArgType ini_mu_pos_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_cl_pos_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_huber_pos_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType kernel_tukey_pos_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType kernel_cl_lin_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_huber_lin_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType kernel_tukey_lin_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType kernel_cl_alpha_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_huber_alpha_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType kernel_tukey_alpha_multi_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType huber_pos_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static R_NativePrimitiveArgType tukey_pos_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(ini_mu_pos_multi),
    CDEF(kernel_cl_pos_multi),
    CDEF(kernel_huber_pos_multi),
    CDEF(kernel_tukey_pos_multi),
    CDEF(kernel_cl_lin_multi),
    CDEF(kernel_huber_lin_multi),
    CDEF(kernel_tukey_lin_multi),
    CDEF(kernel_cl_alpha_multi),
    CDEF(kernel_huber_alpha_multi),
    CDEF(kernel_tukey_alpha_multi),
    CDEF(huber_pos),
    CDEF(tukey_pos),
    {NULL, NULL, 0}
};


void R_init_rmargint(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}