# ERFA in the tests of myb

Upstream: https://github.com/liberfa/erfa (BSD-3-Clause, see `LICENSE` here; ERFA is derived,
with permission, from the IAU SOFA library). Imported revision:
`1a8044cde5b7763295d472a6443387239127c6c8` (2026-09-03). Files are unmodified.

Used **only by tests** (`tests/erfa/`) as the IAU SOFA reference for calendars and time
scales. Nothing from here is compiled into the library.

Files: `cal2jd.c jd2cal.c dat.c erfadatextra.c dtdb.c dtf2d.c d2dtf.c d2tf.c utctai.c
taiutc.c tttai.c taitt.c tcgtt.c tttcg.c tcbtdb.c tdbtcb.c tttdb.c tdbtt.c` and the headers
`erfa.h erfam.h erfaextra.h erfadatextra.h`.
