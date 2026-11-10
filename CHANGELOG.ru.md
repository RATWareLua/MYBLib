# История изменений

[English](CHANGELOG.md) | **Русский**

Версии — major.minor.patch.tweak (`MYB_VERSION_*` в `include/myb/core.h`). В пределах старшей
версии API только растёт. Подробности, проверки и измеренные точности каждой фазы —
[`docs/PLAN.ru.md`](docs/PLAN.ru.md), §12.

## 1.2.0.0 — в работе (без тега)

### Добавлено
- `myb/env.h`: окружение модели сил (тела с провайдерами состояния и ориентации, значения на
  момент — один раз, шкалы времени провайдеров по таблице) и готовые члены: точечная масса,
  третье тело, поле тяготения с приливами, давление излучения с тенями, излучение планеты,
  сопротивление, ОТО, эмпирические ускорения.
- `myb/orbit.h`: распространение одного объекта — состояния на запрошенные моменты (не зависят
  от выборки, бит в бит), события с фильтром по направлению, матрица перехода.
- Провайдеры: `myb_spk_body` (`myb/kin.h`); `myb_orientation_function` с ориентациями CIP + EOP,
  PCK, WGCCRE и постоянной (`myb/frames.h`).
- `myb_density_function` и слоистая экспоненциальная атмосфера (`myb/force.h`).
- `myb_gravity_field_gradient`: градиент ускорения поля гармоник; аналитические якобианы
  готовых членов и релятивистских ускорений (`myb_term_*_jacobian`, `myb_rel_*_jacobian`);
  `myb_term_no_jacobian`, `myb_force_model_bind_parameter`.
- `myb/estim.h`: наименьшие квадраты по дуге (`myb_fit_orbit`), расширенный фильтр Калмана
  (`myb_filter`), измерения с производными (положение, дальность, скорость дальности, RA/Dec,
  азимут и высота от станции `myb_station`).
- Вывод: точная печать чисел (`myb_format_double`, `myb_format_fixed`,
  `myb_format_scientific`); писатель SPK `myb_spk_writer` (типы 9, 13, 2, 3); `myb/ccsds.h` —
  чтение и запись CCSDS OEM, OPM и OMM; `myb_sp3_write`; таблицы чисел (`myb_table_*`); в
  `ext/sgp4`: `myb_tle_format`, `myb_tle_from_omm`, `myb_omm_from_tle`.
- Разрывы шага (`myb_force_model_add_break`, `myb_term_radiation_breaks`): `myb_orbit` и
  оцениватели заканчивают шаг на каждом крае тени; `myb_orbit_breaks`.
- `myb/events.h`: готовые события (`myb_event_range`, `myb_event_range_rate`,
  `myb_event_plane`, `myb_event_altitude`, `myb_event_elevation`, `myb_event_eclipse`,
  `myb_event_shape`); `myb_shadow_margin` (`myb/force.h`).
- `myb/geometry.h`: `myb_ray_ellipsoid`, `myb_ellipsoid_nearest`, `myb_ellipsoid_normal`,
  `myb_illumination_angles`, `myb_polyhedron_ray`, `myb_polyhedron_solid_angle`.
- Малые силы (`myb/force.h`): `myb_force_radiation_moving` (Пойнтинг–Робертсон), `myb_nongrav`
  с `myb_nongrav_marsden` и `myb_nongrav_power`, `myb_yarkovsky_drift`, `myb_force_emission`,
  `myb_force_tidal_lag`, `myb_wind_function`; члены `myb_term_nongrav`, `myb_term_tidal_lag`
  (`myb/env.h`).
- `myb/surface.h`: законы ориентации (`myb_attitude_fixed`, `_orbital`, `_sun_nadir`, `_spin`),
  свободномолекулярное обтекание (`myb_force_gas_plates`, `myb_gas_sphere_cd`), свет на сетке
  с самозатенением и его момент (`myb_force_radiation_mesh`), члены
  `myb_term_surface_radiation`, `myb_term_surface_drag`, `myb_term_recoil`.
- `myb_orbit_kick`: импульсные манёвры; `myb_thrust_mass_at`.
- Установка (`cmake --install`), пакет CMake `find_package(myb)`, `myb.pc`, разделяемая
  библиотека по запросу (`MYB_BUILD_SHARED`), `examples/`.
- Оценивание: систематические сдвиги измерений как неизвестные (`myb_fit_options.bias_count`,
  `myb_observation.bias`, `myb_filter_set_bias`, `myb_filter_bias`); сглаживатель
  (`myb_filter_smooth`, `myb_filter_steps`, `myb_filter_step`); измеряемая точка со смещением
  от центра масс (`myb_measure_offset`).
- Силы: `myb_force_nongrav_delayed`; `myb_yarkovsky_response`, `myb_force_yarkovsky`,
  `myb_yarkovsky_season` (сезонный эффект на эксцентричной орбите) и член
  `myb_term_yarkovsky`; `myb_torque_tidal_lag`; `myb_term_field_indirect`;
  `myb_shadow_fraction_many`; `myb_refraction`, `myb_refraction_bending`,
  `myb_shadow_fraction_refracted`; `myb_force_planet_radiation_gradient`; `myb_gas_knudsen`,
  `myb_gas_bridge`; `myb_force_yorp`; якобианы `myb_term_nongrav_jacobian`,
  `myb_term_tidal_lag_jacobian`, `myb_term_empirical_jacobian`,
  `myb_term_planet_radiation_jacobian`, `myb_term_yarkovsky_jacobian` (с
  `myb_force_yarkovsky_gradient` и `myb_force_yarkovsky_season_gradient`),
  `myb_term_surface_radiation_jacobian`, `myb_term_surface_drag_jacobian`,
  `myb_term_recoil_jacobian`.
- `myb_orbit_options.integrator`: DOP853 в `myb_orbit`; `myb_env_clone`.
- Геометрия: `myb_ellipsoid_limb`, `myb_ellipsoid_terminator`, `myb_view_body`,
  `myb_view_intercept`; лучи по сетке идут через дерево ограничивающих объёмов.
- Файлы: `myb_spk_writer_add_raw`; `myb_pck_writer` (бинарный PCK); в `ext/sgp4` —
  `myb_sgp4_spk10_constants`, `myb_sgp4_spk10_packet`, `myb_sgp4_spk10_size`,
  `myb_sgp4_spk10_words`; часы, коды точности и флаги SP3 (`myb_sp3_clock`, `myb_sp3_sigma`,
  `myb_sp3_accuracy`, `myb_sp3_flags`).
- CCSDS: нотация XML (чтение — автоматически, запись — `myb_oem_write_as`,
  `myb_opm_write_as`, `myb_omm_write_as`); ковариации OEM (`myb_oem_covariance`), OPM и OMM;
  сообщение OCM (`myb_ocm_*`).
- Упаковка: версии символов ELF (`cmake/myb.map`); файлы pkg-config библиотек `ext/`.
- Привязки (`bindings/`): C++17 — заголовок `myb/myb.hpp`, цель `myb::cpp`; Python — пакет
  `myb` над разделяемой библиотекой (cffi); Rust — крейты `myb-sys` и `myb`. В каждой весь
  C-интерфейс и классы для обычного пути; результаты совпадают с вызовами C до бита. Пример
  `examples/low_orbit.cpp`.

### Изменено
- `myb_state_function` объявлен в `myb/kin.h` (`myb/observe.h` его включает).
- `myb/force.h` включает `myb/frames.h`.
- Заголовки `ext/` помечают функции `MYB_EXT_API` вместо `MYB_API`.
- `ext/sgp4/myb_sgp4.h` включает `myb/ccsds.h`.
- Новые поля в конце структур (ноль сохраняет прежнее поведение): `myb_object.mass_function`
  и `mass_context`, `myb_radiation_term.moving` и `atmospheres`, `myb_drag_term.wind` и
  `wind_context`, `myb_yarkovsky.conductivity`, `myb_nongrav_term.delay`,
  `myb_sp3_data.clocks` … `flags`, `myb_oem_segment_data.covariance_*`, ковариации в
  `myb_opm` и `myb_omm`.
- Раздел ковариаций OEM теперь читается и проверяется: матрица не из 21 числа —
  `MYB_E_FORMAT` (раньше раздел пропускался).
- Тень от двух и более тел на диске источника в членах давления излучения — объединение
  покрытий вместо произведения долей; с одним телом результат прежний до бита.
- `myb/geometry.h` включает `myb/env.h` и `myb/frames.h`; `myb/estim.h` — `myb/surface.h`.
- `myb_term_drag_jacobian` учитывает градиент ветра (разностями функции ветра); без ветра
  результат прежний до бита.

## 1.1.0.0 — 2026-11-10 (без тега)

### Добавлено
- Компактные двойные до 3PN и 3,5PN со спинами: спин-орбитальные и спин-спиновые члены,
  квадруполи вращения, прецессия спинов (`myb_pn_binary_evaluate`).
- Столб электронов степенных профилей и дисперсионная задержка (`myb/observe.h`).
- `ext/troposphere` (IERS Conventions 2010, гл. 9) и `ext/ionosphere` (IONEX, однослойная
  модель).
- Зарядка сферы в плазме: токи OML, фотоэмиссия, плавающий потенциал (`myb/force.h`).

## 1.0.0.0 — 2026-11-10

Первый выпуск: числовой фундамент (double-double, компенсированные суммы, корректно
округлённые элементарные функции), чтение данных (числа, текстовые ядра NAIF), шкалы времени,
задача двух тел и задача Ламберта, интеграторы (DOP853, RKN12(10), IAS15, Гаусс–Джексон,
симплектические композиции), задача N тел (прямое суммирование, граф взаимодействий, потоки,
семейство Уиздома–Холмана), эфемериды SPK и системы отсчёта (IAU 2006/2000A, EOP, WGCCRE, PCK,
равноденственные системы), ОТО в первом пост-ньютоновском приближении, поля тяготения
протяжённых тел и приливы, негравитационные силы, модели сил с вариационными уравнениями и
ковариациями, ограниченная задача трёх тел, наведение по B-плоскости, бюджет эффектов,
ориентация и вращательная динамика, регуляризация, малая тяга, компактные двойные до 2,5PN,
наблюдаемые, магнитные поля и сила Лоренца; `ext/nrlmsise00`, `ext/sgp4`.
