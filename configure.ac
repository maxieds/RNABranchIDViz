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

CXXFLAGS_FLTK="`fltk-config --use-gl --cxxflags`"
#LDFLAGS_FLTK="`fltk-config --use-gl --ldflags`"
LDFLAGS_FLTK="`find /usr/ -iname '*fltk*.so' | tr "\n" " "`"
LDFLAGS_MESA="`find /usr/lib -iname '*libGL.so.1' | tr "\n" " "`"
LDFLAGS_X11="-lX11"

AC_SUBST(CXXFLAGS_FLTK)
AC_SUBST(LDFLAGS_FLTK)
AC_SUBST(LDFLAGS_MESA)
AC_SUBST(LDFLAGS_X11)

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
