#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# AC_PREREQ([2.67])
BUG-REPORT-ADDRESS="maxieds@gmail.com"
AC_INIT([RNABranchIDViz], [1.0], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_HEADERS([rnabranchidviz_config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

#LDFLAGS_LIBPNG="-lpng"
CXXFLAGS_EXTRA="`pkg-config --cflags cairo`"
LDFLAGS_EXTRA="`pkg-config --libs cairo`"

AC_SUBST(CXXFLAGS_EXTRA)
AC_SUBST(LDFLAGS_EXTRA)

AM_INIT_AUTOMAKE([RNABranchIDViz], [1.0])
AC_OUTPUT(Makefile src/Makefile)

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([getcwd memset sqrt strdup strrchr])

AC_OUTPUT
