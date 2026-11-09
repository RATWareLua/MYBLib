# Данные myb

Библиотека не содержит ни одной константы конкретного тела или конвенции: всё приходит
из файлов в формате текстовых ядер NAIF (`myb/io.h`). Каждый файл сам описывает свои
источники и точность; здесь — сводка.

| Файл | Содержание | Источник | Лицензия |
|---|---|---|---|
| `time/time-scales.tk` | дерево шкал TAI, TT, TCG, TDB, TCB, GPS, GST, BDT | резолюции МАС 1991 A4, 2000 B1.9, 2006 B3; IS-GPS-200 | определения (факты) |
| `time/tdb-tt-fb1990.tk` | ряд TDB−TT (787 членов + 5 поправок JPL) | Fairhead & Bretagnon 1990 в редакции SOFA, числа из ERFA `dtdb.c` | BSD-3 (ERFA) |
| `time/utc-tai.tk` | TAI−UTC с 1960 г. (дрейф до 1972), срок действия | USNO tai-utc.dat через ERFA `dat.c`; IERS `Leap_Second.dat` (бюллетень C 72, до 28.06.2027) | BSD-3 (ERFA); данные IERS |
| `constants/codata-2022.tk` | G (CODATA 2022), c | NIST CODATA 2022; определение метра | факты |
| `frames/iau2006-cip.tk` | ряды X, Y, s + XY/2 модели CIP IAU 2006/2000A (2941 член), 14 фундаментальных аргументов, ERA, s′ | IERS Conventions 2010, табл. 5.2a, 5.2b, 5.2d (iers-conventions.obspm.fr; SHA-256 в файле) | конвенции IERS (факты) |
| `frames/eop-tides.tk` | суточные и полусуточные вариации EOP: океанские приливы (71 член) и либрация (10 + 11) | IERS Conventions 2010, табл. 8.2–8.3 и 5.1a/b, таблицы программ interp.f, PMSDNUT2.F, UTLIBR.F (SHA-256 в файле) | конвенции IERS (факты) |
| `grav/iers2010-tides.tk` | приливы поля Земли: числа Лява (табл. 6.3), частотные поправки шага 2 (табл. 6.5a–c, 71 член, ряды Пуассона), постоянный прилив, полюсные приливы, вековой полюс | IERS Conventions 2010, гл. 6 и 7 (редакция 2018; SHA-256 в файле) | конвенции IERS (факты) |
| `frames/naif-inertial.tk` | встроенные инерциальные системы SPICE (B1950, FK4, GALACTIC, ECLIPJ2000, ECLIPB1950) как TK-системы | NAIF Frames Required Reading; Lieske 1977/1979; Explanatory Supplement 1992 | факты |

## Обновление

- **Секунды координации**: при объявлении в бюллетене IERS C добавить строку в
  `time/utc-tai.tk` и сдвинуть `MYB_UTC_VALID_UNTIL_MJD` на новую дату истечения из
  `Leap_Second.dat`. После этой даты преобразования UTC возвращают `MYB_E_RANGE`, а не
  угадывают.
- **Большие наборы** в репозиторий не кладутся: `scripts/fetch-data.cmd` (на Linux —
  `cmake -P cmake/FetchData.cmake`) скачивает их в `download/` (не версионируется).
  Фиксированные выпуски проверяются по SHA-256 из `cmake/FetchData.cmake`: DE440
  (`de440.bsp`, NAIF), `testpo.440` и `header.440` (JPL SSD), ориентация Луны
  `moon_pa_de440_200625.bpc`, `pck00011.tpc` (NAIF), модель поля GOCO06s (ICGEM, CC BY 4.0),
  океанские приливы FES2004 в виде коэффициентов Стокса (IERS). Живые ряды обновляются по
  `fetch-data.cmd refresh`: EOP 20 C04 (`eopc04.1962-now`, IERS EOP Centre), `finals2000A.all`
  (IERS Rapid Service), `earth_latest_high_prec.bpc` (NAIF). Тесты, которым они нужны,
  пропускаются, если файлов нет.

Файлы JPL (`gm_de440.tpc`, `pck00011.tpc`, ...) читаются тем же загрузчиком без
преобразования.
