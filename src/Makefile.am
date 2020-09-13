# "-I$(srcdir)/.." means that I must "include <WFA/gap_affine/etc.h>" (i.e. the WFA dir is needed)
# "-I$(srcdir)/../WFA"  means that "include <gap_affine/etc.h>" works
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../submodules/WFA -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../submodules/WFA/libwfa.la ../biomcmc-lib/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = uvaia 
uvaia_SOURCES = main.c kseq.h
uvaia_LDADD = $(LDADD)
