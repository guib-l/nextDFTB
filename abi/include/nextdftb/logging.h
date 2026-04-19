/*
 * nextDFTB — C ABI: unified logging sink
 *
 * One logger per language (Fortran, C++, Python) — all route through this
 * single C ABI sink so the output file is unified.
 *
 * Log line format: [TIMESTAMP] [LEVEL] [LAYER:FUNCTION] message
 */
#ifndef NEXTDFTB_ABI_LOGGING_H
#define NEXTDFTB_ABI_LOGGING_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum nextdftb_log_level_e {
    NEXTDFTB_LOG_DEBUG = 0,
    NEXTDFTB_LOG_INFO  = 1,
    NEXTDFTB_LOG_WARN  = 2,
    NEXTDFTB_LOG_ERROR = 3
} nextdftb_log_level_t;

/* Open/close the log file. path==NULL or empty ⇒ stdout only.
 * also_stdout != 0 ⇒ echo every record to stdout as well.
 * Idempotent: calling open twice closes the previous sink first. */
int  nextdftb_log_open (const char* path, int also_stdout);
void nextdftb_log_close(void);

/* Filter threshold. Records below `level` are dropped. */
void nextdftb_log_set_level(int level);
int  nextdftb_log_get_level(void);

/* Emit one record. Thread-safe; the underlying sink serializes writes. */
void nextdftb_log(int level,
                  const char* layer,
                  const char* function,
                  const char* message);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* NEXTDFTB_ABI_LOGGING_H */
