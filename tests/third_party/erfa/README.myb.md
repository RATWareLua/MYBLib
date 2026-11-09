# ERFA in the tests of myb

Upstream: https://github.com/liberfa/erfa (BSD-3-Clause, see `LICENSE` here; ERFA is derived,
with permission, from the IAU SOFA library). Imported revision:
`1a8044cde5b7763295d472a6443387239127c6c8` (2026-09-03). Files are unmodified.

Used **only by tests** (`tests/erfa/`) as the IAU SOFA reference for calendars, time scales,
the CIO-based Earth orientation and geodetic coordinates. Nothing from here is compiled into
the library.

Files: `cal2jd.c jd2cal.c dat.c erfadatextra.c dtdb.c dtf2d.c d2dtf.c d2tf.c utctai.c
taiutc.c tttai.c taitt.c tcgtt.c tttcg.c tcbtdb.c tdbtcb.c tttdb.c tdbtt.c`; for frames
`xy06.c s06.c era00.c sp00.c pom00.c c2ixys.c c2ixy.c c2tcio.c fal03.c falp03.c faf03.c
fad03.c faom03.c fame03.c fave03.c fae03.c fama03.c faju03.c fasa03.c faur03.c fane03.c
fapa03.c gmst06.c ir.c rz.c ry.c rx.c rxr.c cr.c cp.c anp.c gd2gc.c gd2gce.c gc2gd.c gc2gde.c eform.c`;
and the headers `erfa.h erfam.h erfaextra.h erfadatextra.h`.
