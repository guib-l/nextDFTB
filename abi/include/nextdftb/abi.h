/*
 * nextDFTB — public C ABI
 *
 * This is the ONLY cross-language interface. Everything passed across this
 * boundary is a primitive C type or a raw pointer + explicit dimensions.
 *
 * No C++ STL, no Fortran allocatables / derived types / assumed-shape arrays
 * cross this boundary.
 *
 * Memory ownership: callers own all buffers unless explicitly stated.
 * Shared simulation buffers must be 64-byte aligned (cache line / AVX-512).
 */
#ifndef NEXTDFTB_ABI_H
#define NEXTDFTB_ABI_H

#include <stdint.h>

#include "nextdftb/errors.h"
#include "nextdftb/logging.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ *
 *  Runtime lifecycle
 * ------------------------------------------------------------------ */

/* Initialize the Fortran runtime. Sets OpenMP thread count (<=0 leaves
 * OMP_NUM_THREADS / compiler default untouched). Must be called before any
 * kernel invocation. */
int nextdftb_init(int num_threads);
int nextdftb_finalize(void);
int nextdftb_is_initialized(void);

/* ------------------------------------------------------------------ *
 *  Infrastructure smoke test
 * ------------------------------------------------------------------ */

/* Runs a deterministic computation end-to-end (C++ -> C ABI -> Fortran +
 * OpenMP -> status return). Writes a reference value to *out. */
int nextdftb_test(double* out);

/* ------------------------------------------------------------------ *
 *  Version
 * ------------------------------------------------------------------ */

const char* nextdftb_version_string(void);
int nextdftb_version_major(void);
int nextdftb_version_minor(void);
int nextdftb_version_patch(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* NEXTDFTB_ABI_H */
