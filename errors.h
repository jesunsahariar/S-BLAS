/*
    This file is part of HiParTI!.

    HiParTI! is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    HiParTI! is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with HiParTI!.
    If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARTI_ERROR_H
#define PARTI_ERROR_H

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif

/**
 * Check if a value is not zero, print error message and return.
 * @param errcode the value to be checked
 * @param module  the module name of current procedure
 * @param reason  human readable error explanation
 */
#ifndef NDEBUG
#define spt_CheckError(errcode, module, reason) \
    if(unlikely((errcode) != 0)) { \
        spt_ComplainError(module, (errcode), __FILE__, __LINE__, (reason)); \
        return (errcode); \
    }
#else
#define spt_CheckError(errcode, module, reason) \
    if(unlikely((errcode) != 0)) { \
        return (errcode); \
    }
#endif

#ifndef NDEBUG
#define spt_CheckOmpError(errcode, module, reason) \
    if(unlikely((errcode) != 0)) { \
        spt_ComplainError(module, (errcode), __FILE__, __LINE__, (reason)); \
        exit(errcode); \
    }
#else
#define spt_CheckOmpError(errcode, module, reason) \
    if(unlikely((errcode) != 0)) { \
        exit(errcode); \
    }
#endif

/**
 * Check if a condition is true, set the error information as the system error, print error message and return.
 * @param cond   the condition to be checked
 * @param module the module name of current procedure
 */
#define spt_CheckOSError(cond, module) \
    if(unlikely((cond))) { \
        spt_CheckError(errno + SPTERR_OS_ERROR, (module), strerror(errno)); \
    }

/**
 * Check if a condition is true, set the error information as the last CUDA error, print error message and return.
 * @param cond   the condition to be checked
 * @param module the module name of current procedure
 */
#define spt_CheckCudaError(cond, module) \
    if(unlikely((cond))) { \
        cudaError_t _cuda_error = cudaGetLastError(); \
        spt_CheckError(_cuda_error + SPTERR_CUDA_ERROR, (module), cudaGetErrorString(_cuda_error)); \
    }

void spt_ComplainError(const char *module, int errcode, const char *file, unsigned line, const char *reason);

static __thread struct {
    const char *module;
    int errcode;
    const char *file;
    unsigned line;
    char *reason;
} g_last_error = { NULL, 0, NULL, 0, NULL };

/**
 * Set last error information as specified and print the error message.
 * Should not be called directly, use the macro `spt_CheckError`.
 */
void spt_ComplainError(const char *module, int errcode, const char *file, unsigned line, const char *reason) {
    g_last_error.errcode = errcode;
    g_last_error.module = module;
    g_last_error.file = file;
    g_last_error.line = line;
    free(g_last_error.reason);
    if(reason) {
        size_t len = strlen(reason);
        g_last_error.reason = (char *) malloc(len+1);
        if(!g_last_error.reason) {
            abort();
        }
        memcpy(g_last_error.reason, reason, len+1);
    }
    if(g_last_error.reason && g_last_error.reason[0] != '\0') {
        fprintf(stderr, "[%s] error 0x%08x at %s:%u, %s\n",
            g_last_error.module,
            g_last_error.errcode,
            g_last_error.file,
            g_last_error.line,
            g_last_error.reason
        );
    } else {
        fprintf(stderr, "[%s] error 0x%08x at %s:%u\n",
            g_last_error.module,
            g_last_error.errcode,
            g_last_error.file,
            g_last_error.line
        );
    }
}


/**
 * Print out an assertion error and stop the program
 * Should not be called directly, use the macro `sptAssert`.
 */
void spt_Panic(const char *file, unsigned line, const char *expr) {
    fprintf(stderr, "%s:%u: assertion error: \"%s\"\n", file, line, expr);
    abort();
}

#define sptAssert(expr) ((expr) ? (void) true : spt_Panic(__FILE__, __LINE__, #expr))
  
/**
 * Get the last error code and message.
 * @param[out] module store the module name of the last error
 * @param[out] file   store the C source name of the last error
 * @param[out] line   store the line number of the last error
 * @param[out] reason store the human readable error reason
 * @return the error code of the last error
 */
int sptGetLastError(const char **module, const char **file, unsigned *line, const char **reason) {
    if(module) {
        *module = g_last_error.module;
    }
    if(file) {
        *file = g_last_error.file;
    }
    if(line) {
        *line = g_last_error.line;
    }
    if(reason) {
        *reason = g_last_error.reason;
    }
    return g_last_error.errcode;
}


/**
 * Clear the information of the last error.
 * Usually called before a procedure.
 */
void sptClearLastError(void) {
    g_last_error.module = NULL;
    g_last_error.errcode = 0;
    g_last_error.file = NULL;
    g_last_error.line = 0;
    free(g_last_error.reason);
    g_last_error.reason = NULL;
}
  
  
#ifdef __cplusplus
}
#endif

#endif
