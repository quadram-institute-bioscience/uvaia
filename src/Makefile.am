# "-I$(srcdir)/.." means that I must "include <WFA/gap_affine/etc.h>" (i.e. the WFA dir is needed)
# "-I$(srcdir)/../WFA"  means that "include <gap_affine/etc.h>" works (obs: top_srcdir is relative to root so dont need the ".."
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/submodules/WFA -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
LDADD = $(top_builddir)/submodules/WFA/libwfa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la $(GTKDEPS_LIBS) @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

bin_PROGRAMS = uvaialign uvaia_search_u uvaia uvaiaclust
uvaialign_SOURCES = align.c kseq.h
uvaialign_LDADD = $(LDADD)

uvaia_search_u_SOURCES = search_u.c kseq.h
uvaia_search_u_LDADD = $(LDADD)

uvaia_SOURCES = uvaia.c utils.c utils.h
uvaia_LDADD = $(LDADD)

uvaiaclust_SOURCES = cluster.c utils.c utils.h fastaseq.h fastaseq.c
uvaiaclust_LDADD = $(LDADD)
