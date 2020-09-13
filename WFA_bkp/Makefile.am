AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ @ZLIB_LIBS@
AM_CFLAGS = @SIMD_FLAGS@ @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @ZLIB_LIBS@ $(AM_LDFLAGS)

common_headers = benchmark/benchmark_edit.h \
benchmark/benchmark_gap_affine.h \
benchmark/benchmark_gap_lineal.h \
benchmark/benchmark_utils.h \
edit/edit_cigar.h \
edit/edit_dp.h \
edit/edit_table.h \
gap_affine/affine_penalties.h \
gap_affine/affine_table.h \
gap_affine/affine_wavefront_align.h \
gap_affine/affine_wavefront_backtrace.h \
gap_affine/affine_wavefront_display.h \
gap_affine/affine_wavefront_extend.h \
gap_affine/affine_wavefront.h \
gap_affine/affine_wavefront_penalties.h \
gap_affine/affine_wavefront_reduction.h \
gap_affine/affine_wavefront_utils.h \
gap_affine/swg.h \
gap_affine/wavefront_stats.h \
gap_lineal/lineal_penalties.h \
gap_lineal/nw.h \
system/mm_allocator.h \
system/profiler_counter.h \
system/profiler_timer.h \
utils/commons.h \
utils/dna_text.h \
utils/string_padded.h \
utils/vector.h 

common_src = 	benchmark/benchmark_edit.c \
benchmark/benchmark_gap_affine.c \
benchmark/benchmark_gap_lineal.c \
benchmark/benchmark_utils.c \
edit/edit_cigar.c \
edit/edit_dp.c \
edit/edit_table.c \
gap_affine/affine_penalties.c \
gap_affine/affine_table.c \
gap_affine/affine_wavefront_align.c \
gap_affine/affine_wavefront_backtrace.c \
gap_affine/affine_wavefront_backtrace.h \
gap_affine/affine_wavefront.c \
gap_affine/affine_wavefront_display.c \
gap_affine/affine_wavefront_extend.c \
gap_affine/affine_wavefront_penalties.c \
gap_affine/affine_wavefront_reduction.c \
gap_affine/affine_wavefront_utils.c \
gap_affine/swg.c \
gap_affine/wavefront_stats.c \
gap_lineal/nw.c \
system/mm_allocator.c \
system/profiler_counter.c \
system/profiler_timer.c \
utils/commons.c \
utils/dna_text.c \
utils/string_padded.c \
utils/vector.c 

noinst_LTLIBRARIES = libwfa.la # convenience library (
libwfa_la_SOURCES = $(common_src) $(common_headers)

otherincludedir = $(includedir)/wfa  ## all headers are installed in this subdir
otherinclude_HEADERS = $(common_headers)
# more info: https://autotools.io/automake/nonrecursive.html
