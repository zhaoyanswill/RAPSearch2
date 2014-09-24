#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([RAPSearch2], [2.19])
AC_CONFIG_SRCDIR([Src/main.cpp])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall foreign])

AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create], [], [AC_MSG_ERROR(pthread not found.)])
AC_CHECK_LIB([boost_chrono], [main], [], [AC_MSG_ERROR(boost_chrono not found.)])
AC_CHECK_LIB([boost_serialization], [main], [], [AC_MSG_ERROR(boost_serialization not found.)])
AC_CHECK_LIB([boost_system], [main], [], [AC_MSG_ERROR(boost_system not found.)])
AC_CHECK_LIB([boost_thread], [main], [], [AC_MSG_ERROR(boost_thread not found.)])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h])
AC_CHECK_HEADERS([boost/archive/binary_iarchive.hpp boost/archive/binary_oarchive.hpp boost/bind.hpp boost/chrono/thread_clock.hpp boost/lambda/lambda.hpp boost/lexical_cast.hpp boost/serialization/vector.hpp boost/thread/mutex.hpp boost/thread/thread.hpp])
AC_CHECK_HEADERS([boost/threadpool.hpp])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])

AC_CONFIG_FILES([Makefile
                 Src/Makefile])
AC_OUTPUT
