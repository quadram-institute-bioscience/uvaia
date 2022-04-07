# "-I$(srcdir)/.." means that I must "include <WFA/gap_affine/etc.h>" (i.e. the WFA dir is needed)
# "-I$(srcdir)/../WFA"  means that "include <gap_affine/etc.h>" works (obs: top_srcdir is relative to root so dont need the ".."
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/submodules/WFA -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
LDADD = $(top_builddir)/submodules/WFA/libwfa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la $(GTKDEPS_LIBS) @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

bin_PROGRAMS = uvaialign uvaia_search_u uvaia uvaiaclust uvaiaball uvaiann
uvaialign_SOURCES = align.c kseq.h
uvaialign_LDADD = $(LDADD)

uvaia_SOURCES = uvaia.c utils.c utils.h
uvaia_LDADD = $(LDADD)

uvaiaclust_SOURCES = cluster.c utils.c utils.h fastaseq.h fastaseq.c
uvaiaclust_LDADD = $(LDADD)

uvaiaball_SOURCES = ball.c utils.c utils.h fastaseq.h fastaseq.c
uvaiaball_LDADD = $(LDADD)

uvaiann_SOURCES = nearest.c utils.c utils.h fastaseq.h fastaseq.c min_heap.c min_heap.h
uvaiann_LDADD = $(LDADD)
