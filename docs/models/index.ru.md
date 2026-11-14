# Паспорта моделей

[English](index.en.md) | **Русский**

Паспорт — страница, которую читают, прежде чем довериться модели `myb`: что она считает,
откуда она взята, где применима, насколько точна и как это измерено. У всех паспортов одни и
те же девять частей — *Что это*, *Источник*, *Интерфейс*, *Входы и выходы*, *Область
применимости*, *Точность*, *Стоимость*, *Ограничения*, *Данные*, — и каждое число в них
прослеживается до комментария публичного заголовка, до теста (по имени и с эталоном, с которым
он сравнивает) или до плана ([`docs/PLAN.ru.md`](../PLAN.ru.md)). Где ничего не измерялось,
паспорт так и пишет.

В столбце *Проверка* — одна главная проверка; остальные — в паспорте.

## Время, системы отсчёта и файлы

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Шкалы времени и их преобразования, TDB−TT | [time-scales](time-scales.ru.md) | [time.h](../../include/myb/time.h) | против ERFA ≤ 3,2·10⁻¹¹ с; TDB−TT из DE440t — его скорость до 1,4·10⁻¹⁶ (`test_erfa_time`, `test_de440_time`) |
| Ориентация Земли: IAU 2006/2000A, CIO, EOP, таблица нутации, системы от равноденствия | [earth-orientation](earth-orientation.ru.md) | [frames.h](../../include/myb/frames.h) | X, Y против ERFA 1,7·10⁻¹⁷ рад; пункты на Земле против Horizons 0,19 м — разница моделей (0,077 м через PCK ITRF93) (`test_erfa_frames`, `test_earth_frames`) |
| Эфемериды SPK (чтение 15 типов, запись), SP3 | [spk-ephemerides](spk-ephemerides.ru.md) | [kin.h](../../include/myb/kin.h) | `testpo.440` 7,1·10⁻¹⁵ а. е.; положения против CSPICE ≤ 6·10⁻¹⁶ (`test_testpo`, `test_spk`) |
| Двоичные PCK, CK и SCLK (чтение и запись) | [pck-ck-sclk](pck-ck-sclk.ru.md) | [frames.h](../../include/myb/frames.h), [ck.h](../../include/myb/ck.h) | CK против CSPICE: те же ответы, матрица 7·10⁻¹⁶ (типы 1–3), 4,4·10⁻¹⁶ (типы 4–6) (`test_ck`) |
| Повороты, WGCCRE, системы TK, граф систем отсчёта, геодезия, орбитальные системы | [frames-geodesy](frames-geodesy.ru.md) | [frames.h](../../include/myb/frames.h) | граф систем против CSPICE `sxform` 2·10⁻¹⁶; геодезические углы 3,3·10⁻¹⁶ рад (`test_frame_graph`, `test_erfa_frames`) |
| Численные основы и текстовые форматы (корректно округлённые функции, double-double, эллиптические функции, линейная алгебра, числа, CSV, CCSDS) | [numerics-formats](numerics-formats.ru.md) | [cmath.h](../../include/myb/cmath.h), [num.h](../../include/myb/num.h), [io.h](../../include/myb/io.h), [ccsds.h](../../include/myb/ccsds.h) | печать чисел против CPython — ни одного расхождения на 6124 значениях; эллиптические функции против mpmath ≤ 5,6 ulp (`test_format`, `test_elliptic`) |
| События и геометрия: затмения и края тени, пересечение луча с формой | [events-geometry](events-geometry.ru.md) | [events.h](../../include/myb/events.h), [geometry.h](../../include/myb/geometry.h) | контакты затмений против CSPICE `gfoclt` 4,4·10⁻¹⁰ с; луч и эллипсоид против `surfpt` 5·10⁻¹⁶ (`test_events`, `test_geometry`) |

## Два тела и интеграторы

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Уравнение Кеплера, элементы, универсальный пропагатор, матрица перехода двух тел, средние элементы, KS и Зундман | [kepler-regularization](kepler-regularization.ru.md) | [kepler.h](../../include/myb/kepler.h), [regular.h](../../include/myb/regular.h) | пропагатор против MPFR (400 бит) ≤ 4,7·10⁻¹⁵; KS при e = 0,999 — 1,2·10⁻¹² (`test_mpfr_kepler`, `test_regular`) |
| Ламберт, B-плоскость, пристрелка, относительное движение | [transfers-relative-motion](transfers-relative-motion.ru.md) | [lambert.h](../../include/myb/lambert.h), [transfer.h](../../include/myb/transfer.h), [relative.h](../../include/myb/relative.h) | Ламберт пристрелкой 4·10⁻¹² м/с; Клохесси–Уилтшир 3·10⁻¹⁴ (`test_transfer`, `test_relative`) |
| DOP853, RKN12(10) с интерполянтом, Гаусс–Джексон | [classical-integrators](classical-integrators.ru.md) | [integ.h](../../include/myb/integ.h) | RKN12(10): 10 оборотов при e = 0,9 — 4,8·10⁻¹⁵ периода (`test_rkn12`) |
| IAS15 | [ias15](ias15.ru.md) | [integ.h](../../include/myb/integ.h) | 10⁴ оборотов при e = 0,3 — 2,8·10⁻¹¹ периода; закон Брауэра (`test_ias15`) |
| Симплектические композиции (в том числе 7 и 15 стадий), WH/WHFast с корректорами | [symplectic-integrators](symplectic-integrators.ru.md) | [integ.h](../../include/myb/integ.h), [dyn.h](../../include/myb/dyn.h) | коэффициенты бит в бит с MPFR; наблюдаемые порядки 2 … 8 (`test_mpfr_compositions`, `test_whfast`) |

## N тел

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Прямые суммы N тел, граф взаимодействий, потоки; FMM как приближённый уровень | [nbody](nbody.ru.md) | [dyn.h](../../include/myb/dyn.h) | Sum2 ≤ 3,5·10⁻¹⁶ от Σ\|членов\|, биты не зависят от потоков; FMM до округления прямой суммы при p = 16 (`test_nbody`, `test_fmm`) |
| Система тел: уровни, интегрируемое вращение, приливы по полосам, PN-пары, производные | [system-of-bodies](system-of-bodies.ru.md) | [system.h](../../include/myb/system.h) | повторное интегрирование DE440 за 5 лет: Луна 4,9 см, Меркурий 2,1 см (`test_de440_system`) |
| Сближения, столкновения, слияния, поиск соседей | [system-encounters](system-encounters.ru.md) | [system.h](../../include/myb/system.h) | касание двух шаров 10⁻¹² от замкнутой формы; 300 шаров — бильярд до 5,6·10⁻¹³ с (`test_system`) |

## Гравитация и приливы

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Сферические гармоники до степени 2700, градиент | [spherical-harmonics](spherical-harmonics.ru.md) | [grav.h](../../include/myb/grav.h) | против MPFR 256 бит: 7,5·10⁻¹⁶ (степень 120), 4,8·10⁻¹⁶ (степень 2190); градиент против mpmath 2,1·10⁻¹⁵ (`test_mpfr_grav`, `test_gradient`) |
| Многогранник, масконы, эллипсоидальные гармоники | [irregular-bodies-gravity](irregular-bodies-gravity.ru.md) | [grav.h](../../include/myb/grav.h) | многогранник против призмы Nagy et al. 7,6·10⁻¹⁶; масконы 1,2·10⁻¹⁵ (`test_grav`, `test_ellipsoidal`) |
| Фигура–фигура, вращательная динамика, свободный волчок, либрации Луны | [figure-figure-rotation](figure-figure-rotation.ru.md) | [grav.h](../../include/myb/grav.h), [attitude.h](../../include/myb/attitude.h) | свободный волчок 2·10⁻¹⁴ рад; либрации против DE440 0,040″ за год (`test_attitude`, `test_moon_libration`) |
| Приливы твёрдой Земли IERS 2010, океанские приливы FES2004, полюсные приливы, вторичные волны | [tides-iers](tides-iers.ru.md) | [grav.h](../../include/myb/grav.h) | шаг 1 против теоремы сложения 2·10⁻¹⁶; LAGEOS против ILRS 2,05 → 1,57 см (`test_tides`, `test_lageos`) |
| Запаздывающие приливы: постоянное запаздывание, приливы по полосам | [tides-lagged](tides-lagged.ru.md) | [force.h](../../include/myb/force.h), [system.h](../../include/myb/system.h) | удаление Луны 3,79 см в год по формуле Миньяра, десять оборотов интегрирования — до 6 её знаков; Луна DE440 — 19 см за 5 лет (`test_nongrav`, `test_de440_system`) |

## Теория относительности

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| 1PN: Шварцшильд, Лензе–Тирринг, де Ситтер, EIH с β, γ, локальные системы, собственное время | [relativity-1pn](relativity-1pn.ru.md) | [rel.h](../../include/myb/rel.h), [dyn.h](../../include/myb/dyn.h) | Меркурий 42,9804″ за столетие (5·10⁻⁹ от формулы); EIH против независимого лагранжиана 2,7·10⁻¹⁰ (`test_rel`) |
| PN-двойная: 1PN–3,5PN, 4PN с хвостом, спиновые члены, реакция излучения | [pn-binary](pn-binary.ru.md) | [rel.h](../../include/myb/rel.h) | коэффициенты против Blanchet 4·10⁻¹⁶; PSR B1913+16 — Ṗ до 3·10⁻⁵ (`test_pn_higher`, `test_pn_binary`) |

## Негравитационные силы

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Давление излучения, полутень, тени нескольких тел, рефракция, Пойнтинг–Робертсон, альбедо и ИК | [radiation-pressure](radiation-pressure.ru.md) | [force.h](../../include/myb/force.h) | доля диска в тени против замкнутых форм и квадратуры 1,1·10⁻¹⁴; Пойнтинг–Робертсон — 7 знаков усреднённой теории (`test_force`, `test_nongrav`) |
| Свет на поверхностях: законы ориентации, сетки с самозатенением, тепловая отдача и антенна | [radiation-surfaces](radiation-surfaces.ru.md) | [surface.h](../../include/myb/surface.h) | освещённая площадь 9,000000000000 м²; 63 Вт на 241 кг — 8,72·10⁻¹⁰ м/с² (`test_surface`, `test_nongrav`) |
| Ярковский и YORP | [yarkovsky-yorp](yarkovsky-yorp.ru.md) | [force.h](../../include/myb/force.h), [surface.h](../../include/myb/surface.h) | отклик против mpmath 5·10⁻¹³; сезонный дрейф 4·10⁻¹³; YORP меняет знак у зеркального тела (`test_nongrav`, `test_surface`) |
| Кометы и газовыделение | [comets-outgassing](comets-outgassing.ru.md) | [force.h](../../include/myb/force.h) | дрейф от A2 5·10⁻⁷ от формулы; комета 5·10⁻⁶ от усреднённого уравнения Гаусса (`test_nongrav`) |
| Сопротивление, взаимодействие газа с поверхностью, сплошная среда, нагрев при входе и абляция | [aerodynamics](aerodynamics.ru.md) | [force.h](../../include/myb/force.h), [surface.h](../../include/myb/surface.h) | пластины Сентмана по сфере — замкнутая форма до 5·10⁻¹⁵; абляция против замкнутого решения 1,9·10⁻¹⁴ (`test_surface`, `test_entry`) |
| Магнитные поля (IGRF), сила Лоренца, зарядка | [magnetic-field-charging](magnetic-field-charging.ru.md) | [force.h](../../include/myb/force.h) | IGRF-14 против ppigrf 1,4·10⁻¹⁵; OML против mpmath 6·10⁻¹⁶ (`test_charge`, `test_charging`) |
| Контакт и посадка | [contact-landing](contact-landing.ru.md) | [contact.h](../../include/myb/contact.h) | восстановление против замкнутой формы 3,5·10⁻¹³ (`test_contact`) |
| Удары микрометеороидов | [micrometeoroids](micrometeoroids.ru.md) | [meteoroid.h](../../include/myb/meteoroid.h) | число ударов 0,4σ от λT; массы — Колмогоров–Смирнов 1,17 (`test_meteoroid`) |

## Окружение Земли

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| NRLMSISE-00, вход моделей верхней атмосферы, космическая погода; копии пользователя DTM2020, NRLMSIS 2.1, HWM14 (не в репозитории) | [thermosphere-models](thermosphere-models.ru.md) | [myb_nrlmsise00.h](../../ext/nrlmsise00/myb_nrlmsise00.h), [weather.h](../../include/myb/weather.h) | 17 эталонных случаев до 7 напечатанных знаков (`test_nrlmsise00`) |
| JB2008 по статьям | [jb2008](jb2008.ru.md) | [myb_jb2008.h](../../ext/jb2008/myb_jb2008.h) | против программы авторов: T∞ 1,1·10⁻⁶ K, log₁₀ ρ 8,9·10⁻⁴ (СКО) (`test_jb2008`) |
| Тропосфера: модели IERS гл. 9, VMF3/GPT3, градиенты | [troposphere](troposphere.ru.md) | [myb_troposphere.h](../../ext/troposphere/myb_troposphere.h) | GPT3 6·10⁻¹³ и VMF3 1,3·10⁻¹⁵ от программ TU Wien (`test_troposphere`) |
| Ионосфера (IONEX 2-D и 3-D, высшие порядки) и задержка в плазме | [ionosphere-plasma](ionosphere-plasma.ru.md) | [myb_ionosphere.h](../../ext/ionosphere/myb_ionosphere.h), [observe.h](../../include/myb/observe.h) | интегралы вдоль луча против прямых сумм 1,3·10⁻⁹; столб плазмы против mpmath 6·10⁻¹⁶ (`test_ionosphere`, `test_plasma`) |
| SGP4, TLE и их подгонка | [sgp4-tle](sgp4-tle.ru.md) | [myb_sgp4.h](../../ext/sgp4/myb_sgp4.h) | `tcppver.out` 5,0·10⁻⁹ км (половина последнего напечатанного знака); SPK типа 10 против CSPICE 2,7·10⁻¹² км (`test_sgp4`) |

## Наземная станция и измерения

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Смещения пунктов (IERS 2010, гл. 7) | [site-displacement](site-displacement.ru.md) | [site.h](../../include/myb/site.h) | прилив твёрдой Земли 1,1·10⁻⁷ м от DEHANTTIDEINEL; океанская нагрузка 5·10⁻⁷ м от HARDISP (`test_site`) |
| Наблюдаемые и время распространения света | [observables-light-time](observables-light-time.ru.md) | [observe.h](../../include/myb/observe.h), [estim.h](../../include/myb/estim.h) | время света против CSPICE 4·10⁻¹² с; аберрация бит в бит с ERFA (`test_observe`, `test_light_measure`) |

## Приложения

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Определение начальной орбиты: Гиббс, Херрик–Гиббс, Гаусс, Лаплас, двойная r-итерация, Гудинг | [initial-orbit-determination](initial-orbit-determination.ru.md) | [iod.h](../../include/myb/iod.h) | примеры 7-3 и 7-4 Валладо до шести знаков; Гаусс 5·10⁻¹⁵ … 3·10⁻¹⁶ (`test_iod`) |
| Оценивание орбиты: наименьшие квадраты, фильтр, сглаживатель | [estimation](estimation.ru.md) | [estim.h](../../include/myb/estim.h) | LAGEOS против ILRS 1,57 см за неделю; Roadster 306 м за 7 лет, OSIRIS-REx 12,6 м (`test_lageos`, `test_spacecraft`) |
| Сближения и вероятность столкновения (короткие и долгие встречи) | [conjunctions](conjunctions.ru.md) | [conjunction.h](../../include/myb/conjunction.h) | круглая гауссиана 2·10⁻¹⁶; долгая встреча против Монте-Карло 1,3σ (`test_conjunction`) |
| Ограниченная задача трёх тел: точки Лагранжа, периодические орбиты, семейства, многообразия | [cr3bp-periodic](cr3bp-periodic.ru.md) | [cr3bp.h](../../include/myb/cr3bp.h) | каталог JPL (56 орбит): постоянная Якоби 3,6·10⁻¹⁵, неподвижная точка 9·10⁻¹² (`test_cr3bp`) |
| Торы в задаче трёх тел: квазипериодические орбиты, их устойчивость и многообразия | [cr3bp-tori](cr3bp-tori.ru.md) | [cr3bp.h](../../include/myb/cr3bp.h) | окружность возвращается на себя за 20 оборотов до 3,7·10⁻¹¹ (`test_cr3bp`, `test_torus_stability`) |
| Малая тяга: дуги тяги, уравнения Гаусса, Эдельбаум, Q-закон | [low-thrust-qlaw](low-thrust-qlaw.ru.md) | [thrust.h](../../include/myb/thrust.h) | Q-закон 1234,3 м/с против 1232,6 Эдельбаума (`test_thrust`) |
| Оптимальное управление: непрямой метод, НЛП, Sims–Flanagan, цепочки с пролётами | [optimal-control](optimal-control.ru.md) | [optimal.h](../../include/myb/optimal.h) | 486,84 м/с против 486,82 Гомана; Sims–Flanagan +0,11 % от непрямого метода (`test_optimal`, `test_chain`) |

## Пропагация и неопределённость

| Модель | Паспорт | Заголовок | Проверка |
|---|---|---|---|
| Окружение, модель сил, готовые члены, пропагация, уравнения в вариациях, ковариация, бюджет эффектов | [propagation-uncertainty](propagation-uncertainty.ru.md) | [env.h](../../include/myb/env.h), [model.h](../../include/myb/model.h), [orbit.h](../../include/myb/orbit.h), [uncert.h](../../include/myb/uncert.h), [budget.h](../../include/myb/budget.h) | готовые члены против сборки вручную — бит в бит; матрица перехода против замкнутой формы 5,5·10⁻¹⁵ (`test_env`, `test_uncert`) |
