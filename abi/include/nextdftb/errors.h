/* 
 * nextDFTB — C ABI: error state
 *
 * Cross-layer error model. Every nextdftb_* function returns an int status:
 *   0            = NEXTDFTB_OK
 *   negative     = error (see enum below)
 *
 * When a function returns non-zero, the last-error state is populated and
 * can be queried via the getters. The state is thread-local (per-thread).
 *
 * Only iso_c_binding-compatible types cross this boundary.
 */
#ifndef NEXTDFTB_ABI_ERRORS_H
#define NEXTDFTB_ABI_ERRORS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum nextdftb_status_e {
    NEXTDFTB_OK                      =   0,
    NEXTDFTB_ERR_INVALID_ARG         =  -1,
    NEXTDFTB_ERR_NOT_INITIALIZED     =  -2,
    NEXTDFTB_ERR_ALREADY_INITIALIZED =  -3,
    NEXTDFTB_ERR_ALLOCATION          =  -4,
    NEXTDFTB_ERR_DIMENSION_MISMATCH  =  -5,
    NEXTDFTB_ERR_NUMERICAL           =  -6,  /* NaN, non-convergence — recoverable */
    NEXTDFTB_ERR_FATAL               =  -7,  /* unrecoverable, clean shutdown */
    NEXTDFTB_ERR_IO                  =  -8,
    NEXTDFTB_ERR_NOT_IMPLEMENTED     = -99
} nextdftb_status_t;

typedef enum nextdftb_severity_e {
    NEXTDFTB_SEV_RECOVERABLE = 0,
    NEXTDFTB_SEV_FATAL       = 1
} nextdftb_severity_t;

/* Getters. Buffer writes are NUL-terminated, truncated to buffer_len. */
int  nextdftb_get_last_error_code(void);
void nextdftb_get_last_error_message (char* buffer, int32_t buffer_len);
void nextdftb_get_last_error_location(char* buffer, int32_t buffer_len);
int  nextdftb_get_last_error_severity(void);
void nextdftb_clear_error(void);

/* Setter — used by layers below the ABI (Fortran kernel, C++ core) to
 * publish an error. Not meant for external callers, but exported for the
 * orchestrator and for C++ exception translators. */
void nextdftb_set_error(int code,
                        int severity,
                        const char* layer,
                        const char* function,
                        const char* message);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* NEXTDFTB_ABI_ERRORS_H */
