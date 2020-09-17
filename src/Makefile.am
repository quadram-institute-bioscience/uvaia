# "-I$(srcdir)/.." means that I must "include <WFA/gap_affine/etc.h>" (i.e. the WFA dir is needed)
# "-I$(srcdir)/../WFA"  means that "include <gap_affine/etc.h>" works (obs: top_srcdir is relative to root so dont need the ".."
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/submodules/WFA -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  $(top_builddir)/submodules/WFA/libwfa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = uvaia_align  uvaia_search_u uvaia_search_a
uvaia_align_SOURCES = align.c kseq.h
uvaia_align_LDADD = $(LDADD)

uvaia_search_u_SOURCES = search_u.c kseq.h
uvaia_search_u_LDADD = $(LDADD)

uvaia_search_a_SOURCES = search_a.c utils.c utils.h kseq.h
uvaia_search_a_LDADD = $(LDADD)
