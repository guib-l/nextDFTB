#include "nextdftb/abi.h"

// These macros are injected by CMake target_compile_definitions on
// nextdftb_core_cpp (see core_cpp/CMakeLists.txt).
#ifndef NEXTDFTB_VERSION_MAJOR
#  define NEXTDFTB_VERSION_MAJOR 0
#endif
#ifndef NEXTDFTB_VERSION_MINOR
#  define NEXTDFTB_VERSION_MINOR 1
#endif
#ifndef NEXTDFTB_VERSION_PATCH
#  define NEXTDFTB_VERSION_PATCH 0
#endif
#ifndef NEXTDFTB_VERSION_STRING
#  define NEXTDFTB_VERSION_STRING "0.1.0"
#endif

extern "C" {

const char* nextdftb_version_string(void) { return NEXTDFTB_VERSION_STRING; }
int         nextdftb_version_major (void) { return NEXTDFTB_VERSION_MAJOR;  }
int         nextdftb_version_minor (void) { return NEXTDFTB_VERSION_MINOR;  }
int         nextdftb_version_patch (void) { return NEXTDFTB_VERSION_PATCH;  }

} // extern "C"
