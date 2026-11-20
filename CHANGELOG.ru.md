# История изменений

[English](CHANGELOG.md) | **Русский**

Версии — major.minor.patch.tweak (`MYB_VERSION_*` в `include/myb/core.h`). В пределах старшей
версии API только растёт. Подробности, проверки и измеренные точности каждой фазы —
[`docs/PLAN.ru.md`](docs/PLAN.ru.md), §12.

## 1.2.2.1 — 2026-11-20

То, что нашёл первый прогон CI, по тегу 1.2.2.0. Интерфейс и результаты — те же, что в 1.2.2.0
(тот же эталон `test_repro`).

### Исправлено
- Публичные заголовки под MSVC 2026: он игнорирует `[[nodiscard]]` после
  `__declspec(dllimport)` (C5240, ошибка под /WX), а C++ допускает стандартный атрибут только в
  начале объявления или после его последнего спецификатора. Каждое объявление теперь
  `MYB_NODISCARD MYB_API …` (их 528); в C23 и C++ атрибут действует у любого компилятора.
- macOS: в репозитории не было `cmake/toolchains/clang-macos.cmake`, который называют пресеты
  `macos-*` (его захватило локальное правило исключения), и они не конфигурировались.
- Скрипты тестов, запускаемые `cmake -P`, задают политики CMake 3.28, минимума проекта: под
  CMake 3.x они шли со старым поведением, и `test_exports` падал на `IN_LIST`.
- `test_go_binding` на Windows с Clang, у которого компоновщик для MinGW по умолчанию GNU ld
  (Clang 20): Go спрашивал компилятор о компоновщике без флагов cgo и отдавал lld сценарий
  компоновщика GNU; компоновщик сборки теперь идёт и в `-extldflags`.
- `test_contact` строит склоны корректно округлёнными `myb_sin`, `myb_cos` и `myb_atan`
  библиотеки вместо libm, у которой последние биты различаются между версиями среды выполнения
  C: брусок, удерживаемый на тысячную ниже угла трения, шёл за ними за границу своей скорости.

## 1.2.2.0 — 2026-11-20

Первый выпуск с тегом после 1.0.0.0: всё, что сделано после 1.1.0.0 (фазы 12–22, привязки).
По пути версия была 1.2.0.0 и отправлялась в `main` без тега.

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
- `myb_budget_linear` (`myb/budget.h`): бюджет эффектов по линейной оценке — одна пропагация
  с вариацией на каждый член вместо пропагации на член; без порога округления.
- `myb/regular.h`: элемент времени в уравнениях KS (`MYB_REGULAR_KS_ELEMENT`,
  `myb_ks_element_from_cartesian`, `myb_ks_element_to_cartesian`, `myb_ks_element_rhs`);
  преобразование Сундмана с любым показателем α ≥ 0.
- XML-нотация сообщений CCSDS: значения в разделах CDATA; объявление типа документа с
  внутренним подмножеством пропускается.
- `myb/system.h`: система тел с моделями (`myb_system`) — над точечными массами `myb_nbody`
  уровни графа: поле тела на точку, квадруполи двух тел, прилив по полосам с запаздыванием,
  искажение тела, члены выше 1PN у пары (со спинами), произвольные действия и свет на сетке
  тела; вращение тел (мантия с жидким ядром), интегрируемое вместе с движением; сближения,
  столкновения, отскок и слияние.
- `myb_ias15_options.groups`, `group_count`: группы координат разных видов, каждая задаёт
  шаг и судит сходимость по себе.
- `myb/iod.h`: начальное определение орбиты — `myb_iod_gibbs`, `myb_iod_herrick_gibbs`,
  `myb_iod_gauss` (по направлениям, с уточнением методом Ньютона).
- `myb/relative.h`: относительное движение — состояние в системе ведущего, линейная карта
  для любой коники (`myb_relative_stm`), Клохесси–Уилтшир (`myb_cw_stm`), двухимпульсная
  встреча, относительные элементы.
- `myb/conjunction.h`: ближайшие сближения двух объектов (`myb_approach_find`,
  `myb_approach_orbits`) и вероятность столкновения (`myb_encounter_plane`,
  `myb_collision_probability`, `_chan`, `_max`).
- `myb/kepler.h`: `myb_mean_elements_average`; `ext/sgp4`: `myb_tle_fit`, `myb_tle_from_state`.
- `myb/cr3bp.h`: многократная стрельба и продолжение семейств периодических орбит
  псевдодлиной дуги (`myb_cr3bp_orbit_from_state`, `_correct`, `_continue`).
- `myb/thrust.h`: закон Q-law малой тяги — `myb_qlaw`, `myb_qlaw_direction`,
  `myb_qlaw_arrived`, член модели `myb_term_qlaw` с порогом эффективности; тот же двигатель
  с массой в состоянии — `myb_qlaw_flight_rhs`, `myb_qlaw_fly`.
- `myb/conjunction.h`: `myb_conjunction_find` — сближения с ковариациями, которые переносят
  матрицы перехода самих орбит, и вероятность столкновения в каждом.
- `myb/optimal.h` (новый): непрямой метод для малой тяги — `myb_indirect_solve` (минимум
  топлива за заданное время, минимум времени), `myb_indirect_sample`; прямой метод
  Sims–Flanagan — `myb_sims_flanagan_mismatch`, `myb_sims_flanagan_solve`; решатель
  нелинейных программ — `myb_nlp_solve`.
- `myb/cr3bp.h`: квазипериодические орбиты как инвариантные окружности потока —
  `myb_cr3bp_torus_from_orbit`, `myb_cr3bp_torus_continue` (семейства постоянной Якоби,
  периода или угла поворота), `myb_cr3bp_torus_state`.
- Файлы и системы отсчёта: SPK типа 15; граф систем отсчёта (`myb_frame_graph_*`) и
  состояние через разные системы (`myb_spk_state_in`); `myb/ck.h` (новый) — часы аппарата
  SCLK (`myb_sclk_*`), чтение CK типов 1–3 (`myb_ck_pointing_at`, `myb_orientation_ck`) и их
  запись (`myb_ck_writer_*`); таблица нутации `myb_cip_tabulate`; разность шкал функцией и
  таблицей (`myb_timesys_set_difference`, `myb_timesys_tabulate`), TDB−TT из эфемериды
  (`myb_spk_time_difference`).
- `myb/system.h`: тела с заданным движением (`myb_system_add_guided`), меняющаяся масса
  (`myb_system_set_mass_rate`), эффект Лензе–Тирринга (`myb_system_add_frame_dragging`).
- `myb/weather.h` (новый): таблица космической погоды (`myb_space_weather_*`) и вход
  моделей верхней атмосферы (`myb_thermosphere_point_at`); в `ext/nrlmsise00` —
  `myb_nrlmsise00_density`, модель как `myb_density_function`.
- `myb/site.h` (новый): смещения пункта по IERS 2010 гл. 7 — `myb_site_displacement`
  (прилив твёрдой Земли, полюсный прилив, нагрузки), океанская нагрузка из таблиц BLQ
  (`myb_ocean_loading_*`), сетка коэффициентов (`myb_site_grid_*`); ядро
  `data/frames/iers2010-site.tk`. Станция, которая движется: `myb_station.motion`.
- `myb/grav.h`: вторичные волны океанского прилива по адмиттансу (`MYB_TIDES_ADMITTANCE`,
  `myb_tides_secondary_waves`), добавление волн к загруженным (`myb_tides_add_waves`).
- `myb/estim.h`: измерения со световым временем — `myb_light_measure`,
  `myb_measure_light_range`, `_range_rate`, `_radec`, `_azel`.
- `ext/troposphere`: `myb_tropo_vmf3`, GPT3 (`myb_gpt3_*`), `myb_tropo_zwd_askne`,
  `myb_tropo_gradient_delay_parts`. `ext/ionosphere`: интегралы вдоль луча
  (`myb_iono_ray_integrals`), плотность слоем над картами и сеткой (`myb_iono_layer_density`,
  `myb_iono_grid_density`), члены высших порядков (`myb_iono_orders_of`,
  `myb_iono_phase_delay`, `myb_iono_code_delay`, `myb_iono_bending`).
- `cmake/FetchData.cmake`: сетка GPT3; большие файлы по `-DLARGE=ON`.
- Фаза 22. `myb/ck.h`: вычисление CK типов 4, 5, 6 и их запись
  (`myb_ck_writer_add_chebyshev`, `myb_ck_writer_add_packets`,
  `myb_ck_writer_add_minisegments`).
- `myb/site.h`: атмосферная нагрузка S1/S2 из файлов — `myb_site_atmosphere_load`,
  `_load_file` (GGFC, станции TU Wien), `myb_site_atmosphere_grid_load`, `_grid_load_file`
  (сетки TU Wien); сетка s1_s2_s3_cm_noib_grid.dat в `cmake/FetchData.cmake`.
- `ext/ionosphere`: трёхмерные карты IONEX — `myb_ionex_dimension`,
  `myb_ionex_density_grid`.
- `ext/jb2008` (новый): JB2008, реализованная по статьям, — `myb_jb2008`, чтение файлов
  индексов SOLFSMY и DTCFILE (`myb_jb2008_solar_*`, `myb_jb2008_storm_*`) и
  `myb_jb2008_density` как `myb_density_function`.
- `myb/iod.h`: `myb_iod_laplace`, `myb_iod_double_r`, `myb_iod_gooding`.
- `myb/integ.h`: композиции из 7 и 15 стадий (`MYB_SYMPLECTIC_6_S7`,
  `MYB_SYMPLECTIC_8_S15`); интерполянт шагов RKN12(10) (`myb_rkn12_options.interpolant`).
- `myb/orbit.h`: `MYB_ORBIT_RKN12` для моделей сил из членов, зависящих только от положения
  (`myb_force_model_set_position_only`, `myb_force_model_position_only` в `myb/model.h`).
- `myb/conjunction.h`: вероятность долгого сближения — `myb_collision_rate`,
  `myb_collision_probability_long`, `myb_conjunction_long`.
- `myb/force.h`: вторичная эмиссия и тонкий слой в заряжении — `myb_secondary_emission`,
  `myb_secondary_yield`, `myb_charging_model`, `myb_charging_current_model`,
  `myb_charging_equilibria`, `myb_debye_length`; вход в атмосферу — `myb_air`,
  `myb_air_heat_ratio`, `myb_newtonian_cp_max`, `myb_heating_sutton_graves`, `myb_ablation`,
  `myb_ablation_rate`, `myb_ablation_rhs`, `myb_ablation_mass`.
- `myb/surface.h`: ньютоновское сплошное обтекание пластин и сеток
  (`myb_force_newtonian_plates`, `myb_force_newtonian_mesh`, `myb_continuum`); свет в данный
  момент с тепловой инерцией (`myb_force_radiation_mesh_lagged`).
- `myb/system.h`: слияния вращающихся тел и тел PN-пар; производные по состоянию и
  параметрам (`myb_system_partials`); поиск соседей в `myb_system_advance` (те же события,
  рассматриваются только пары, чьи пути за шаг могут встретиться).
- `myb/contact.h` (новый): контакт тела с плоскостью, эллипсоидом или сеткой — закон
  Ханта–Кроссли с гладким трением Кулона (`myb_contact_force`, `myb_contact_surface_nearest`,
  `myb_contact_restitution`, `myb_contact_damping`), член `myb_term_contact` и зазор
  `myb_event_contact`; в `myb/system.h` — `myb_system_add_contact`; в `myb/geometry.h` —
  `myb_polyhedron_nearest`.
- `myb/meteoroid.h` (новый): случайные удары микрометеороидов по Grün et al. —
  `myb_meteoroid_flux`, `myb_meteoroid_rate`, `myb_meteoroid_mass`, `myb_meteoroid_next`,
  `myb_meteoroid_propagate`.
- `myb/cr3bp.h`: устойчивость и многообразия торов — `myb_cr3bp_torus_stability`,
  `myb_cr3bp_torus_manifold`; `myb_cr3bp_propagate` с поверхностью сечения.
- `myb/dyn.h`: быстрый метод мультиполей как приближённый уровень графа —
  `myb_nbody_set_approximation`, `myb_nbody_approximation` (по умолчанию выключен).
- `myb/rel.h`: `MYB_PN_4` (4PN с консервативным хвостом), `MYB_PN_4RR` (диссипативный хвост),
  `MYB_PN_SO25`, `MYB_PN_SO35` (спин-орбитальные члены двух следующих порядков),
  `MYB_PN_SORR` (спин-орбитальная радиационная реакция); `myb_system_add_pn_pair` их принимает.
- `myb/optimal.h`: выведение на орбиту со свободной фазой в непрямом методе
  (`myb_indirect_arrival`, `myb_indirect_problem.arrival`, `phase`); поле модели сил в методе
  Симса–Фланагана (`myb_sims_flanagan.field`); цепочки участков с пролётами — `myb_chain`,
  `myb_chain_size`, `myb_chain_guess`, `myb_chain_evaluate`, `myb_chain_solve`;
  `myb_nlp_options.weight`.
- Фаза 21. Массовые формы вызовов (`*_many`, 26 функций: один вызов на ряд, биты одиночного
  вызова в цикле, индекс отказавшего элемента): в `myb/time.h` — `myb_time_from_seconds_many`,
  `myb_time_seconds_many`, `myb_calendar_to_time_many`, `myb_time_to_calendar_many`,
  `myb_timesys_convert_many`, `myb_timesys_from_utc_many`, `myb_timesys_to_utc_many`; в
  `myb/kin.h` — `myb_spk_state_many`, `myb_spk_state_in_many`; в `myb/frames.h` —
  `myb_cip_terrestrial_many`, `myb_cip_xys_many`, `myb_eop_at_many`,
  `myb_equinox_rotation_many`, `myb_equinox_teme_pef_many`, `myb_wgccre_rotation_many`,
  `myb_pck_rotation_many`, `myb_orientation_many`, `myb_frame_graph_rotation_many`,
  `myb_rotation_apply_position_many`, `myb_rotation_apply_state_many`; в `myb/kepler.h` —
  `myb_kepler_propagate_many`, `myb_kepler_stm_many`, `myb_pv_from_elements_many`,
  `myb_elements_from_pv_many`, `myb_pv_from_equinoctial_many`, `myb_equinoctial_from_pv_many`.
- Привязки к Go (`bindings/go`: cgo над статической библиотекой, крупный слой — один вызов,
  одна задача; пакет `myb/sgp4`) и Haskell (`bindings/haskell`: весь C-интерфейс `Myb.C`,
  порождённый из заголовков, верхний слой `Myb` с `ForeignPtr` и `Either`; на Windows — над
  разделяемой библиотекой). Модули pkg-config `myb-static` и `myb_sgp4-static` ставятся вместе с
  библиотекой и называют архивы путями.
- Верхние слои привязок C++, Python, Rust и Haskell для остального интерфейса — N тел и
  системы, оценивание, CR3BP, наведение и оптимальные перелёты, ориентация, наблюдаемые,
  файлы, массовые формы — и те же области в крупном слое Go; в Rust ещё системы отсчёта
  `myb/frames.h` и `OemSegment.covariance_frame`, `covariances`.
- `myb/contact.h`: трение с залипанием — `myb_contact_friction`
  (`MYB_CONTACT_FRICTION_SMOOTH`, `MYB_CONTACT_FRICTION_STICK`), `myb_contact_force_stick`; в
  `myb/system.h` — прогибы его пружин в состоянии: `myb_system_pack_contact`,
  `myb_system_unpack_contact`.
- `myb/integ.h`: масштаб времени шага IAS15 — `myb_ias15_timescale`
  (`MYB_IAS15_TIMESCALE_ACCELERATION`, `MYB_IAS15_TIMESCALE_JERK`); с ним — контакт под IAS15 в
  `myb_orbit`.
- `myb/model.h`: схема разрывов шага — `myb_break_scheme` (`MYB_BREAKS_RETAKE` по умолчанию и
  `MYB_BREAKS_PREDICT`), `myb_force_model_set_break_scheme`, `myb_force_model_break_scheme`.
- `myb/rel.h`: `MYB_PN_SS3` (спин-спиновые члены следующего порядка), `MYB_PN_SORR5`
  (спин-орбитальная радиационная реакция следующего порядка), `MYB_PN_SOTAIL` (диссипативный
  спин-орбитальный хвост на почти круговых орбитах), `MYB_PN_TAILX` (хвосты 4PN на любой
  соприкасающейся конике).
- `myb/optimal.h`: пролёты с импульсом в цепочках — `myb_chain_flyby`,
  `myb_chain.impulse_exhaust`, `myb_chain_result.impulse`.
- `myb/cr3bp.h`: нормальное вращение эллиптических торов — поля `normal_rotation`,
  `normal_eigenvalue`, `normal_frequency`, `normal_residual`, `normal_re`, `normal_im` в
  `myb_cr3bp_torus_spectrum`; торы около вертикальных ляпуновских орбит.
- Фаззинг: 27 обвязок читателей в `tests/fuzz` (`MYB_FUZZ`; libFuzzer с ASan и UBSan в
  пресете `linux-fuzz`, `tests/fuzz/run_fuzzers.cmake`; в остальных сборках — повтор затравок
  и повреждённых вариантов как тесты ctest, метка `fuzz`). Бенчмарки: `bench/` (`MYB_BENCH`,
  программа `myb_bench`, тест `bench_smoke`). Покрытие: `MYB_COVERAGE`, пресеты
  `linux-coverage` и `windows-coverage`, `tests/tools/coverage.cmake`.
- macOS: набор `cmake/toolchains/clang-macos.cmake` (LLVM из Homebrew), пресеты `macos-dev`,
  `macos-release`, `macos-asan`. CI: Linux на ARM64, macOS, фаззинг, покрытие, Go и GHC,
  потребители, собранные MSVC.
- Тесты: `test_msvc_consumer` (заголовки как C17, C++17 и C++20 под `cl`, пример и лицо C++,
  собранные `cl` против `libmyb.dll`), `test_gcc_consumer` (то же gcc и g++ против `libmyb.so`,
  заголовки ещё и как C23 и C++23), `test_rust_binding_msvc`, `test_go_binding`,
  `test_haskell_binding`; на Linux под ASan — `test_python_binding_asan`,
  `test_rust_binding_asan`.
- `docs/models/`: 45 паспортов моделей и указатель, по-английски и по-русски.

### Изменено
- `MYB_TIDES_ALL` включает `MYB_TIDES_ADMITTANCE`: с ядром, где есть таблица адмиттанса, и
  загруженными океанскими волнами результат меняется (вторичные волны добавляются).
- Новые поля в конце структур (ноль сохраняет прежнее поведение): `myb_station.motion`,
  `myb_spk_body.graph`.
- `myb_measure_offset`: производные включают поворот смещения вслед за орбитой (раньше —
  производные измерения в самой точке).
- `myb_ionex_vtec`: значение, отсутствующее в соседней карте, не мешает запросу точно в
  эпоху карты.
- `myb/kin.h` включает `myb/frames.h`; `myb/system.h` — `myb/kin.h`; `myb/estim.h` —
  `myb/site.h`; `myb/weather.h` — `myb/frames.h` и `myb/time.h`.
- `myb_state_function` объявлен в `myb/kin.h` (`myb/observe.h` его включает).
- `myb/force.h` включает `myb/frames.h`.
- `myb/thrust.h` включает `myb/model.h`.
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
- `myb_sundman_rhs` при показателе, отличном от 1, 3/2 и 2, больше не возвращает
  `MYB_E_DOMAIN` (кроме отрицательного и не конечного); для 1, 3/2 и 2 результат прежний до
  бита. В описании `myb/regular.h` исправлено: шагам, равномерным по истинной аномалии,
  отвечает α = 2, а не 3/2.
- Заголовок C++: функции, которые вызывает C, перенесены в `myb::detail` и закреплены
  `static_assert(noexcept(…))`.
- `myb_orbit` с DOP853 принимает `variational` (раньше — `MYB_E_DOMAIN`): вариации
  интегрируются вместе с состоянием.
- Фаза 22. Новые поля в конце структур (ноль сохраняет прежнее поведение):
  `myb_offset_measure.gm`, `myb_system_light.theta`, `spin_axis`, `phases`,
  `myb_surface_drag_term.continuum`, `myb_rkn12_options.interpolant`; `myb_system_merger.spin`.
- `myb_system_merge` принимает тела с интегрируемым вращением и тела PN-пар (раньше —
  `MYB_E_DOMAIN`); `MYB_E_DOMAIN` остаётся для b с интегрируемым вращением при a без него.
- CK: сегменты типов 4, 5, 6 вычисляются (раньше загружались и пропускались).
- `myb_ionex_load` читает трёхмерные карты (раньше — `MYB_E_DOMAIN`); `myb_iono_grid`:
  замкнутые по долготе сетки идут по кругу, пропущенный узел — `MYB_E_RANGE`.
- `myb_ellipsoid_nearest` больше не возвращает `MYB_E_NOCONV` для точек у самой сферы (шаги
  Ньютона падали ниже разрешения); прежние значения те же.
- `myb/myb.h` включает `myb/contact.h` и `myb/meteoroid.h`.
- `myb_wh_start` отказывает (`MYB_E_DOMAIN`) системе с приближением.
- `myb_pn_binary.r0` — также калибровочная длина логарифмов 4PN (обе постоянные r′₀ = r″₀).
- Новые поля в конце структур: `myb_indirect_problem.arrival`, `phase`,
  `myb_indirect_solution.phase`, `myb_indirect_point.costate`, `myb_sims_flanagan.field`,
  `myb_nlp_options.weight`.
- `myb_nlp_solve`: вес штрафа больше не растёт, когда ограничения уже выполнены (результаты
  прежних тестов те же до бита).
- `myb/optimal.h` включает `myb/kin.h` и `myb/time.h`.
- Фаза 21. Новые поля в конце структур (ноль сохраняет прежнее поведение):
  `myb_ias15_options.timescale`, `myb_orbit_options.timescale`, `myb_contact_law.kind`,
  `tangential_stiffness`, `tangential_damping`, `myb_chain.impulse_exhaust`,
  `myb_chain_result.impulse`, поля `normal_*` в `myb_cr3bp_torus_spectrum`.
- `myb_system_partials`: производные дают вариационные уравнения, интегрируемые вместе с
  состоянием (сигнатура та же; два тела против замкнутой формы — 2,3·10⁻¹⁶ вместо 1,8·10⁻¹¹ у
  разностей); система с контактами по-прежнему берётся разностями целых пропагаций.
- FMM (`myb_nbody_set_approximation`) с 2048 тел идёт на пуле `myb_nbody_set_pool` с битами
  вычисления без пула.
- `myb_nlp_solve`: где убыль функции меньше 10⁻¹² от её значений, шаг принимается по наклону
  (приближённые условия Вульфа); задача 71 Хока и Шиттковски сдвигается на 10⁻¹⁴, участок
  Симса–Фланагана берёт 436 вычислений вместо 435, а цепочки, не доходившие до
  стационарности, сходятся.
- `myb_system_add_contact` принимает трение с залипанием (3 координаты на точку после всех
  остальных); `myb_contact_force` и `myb_term_contact` его отвергают (`MYB_E_DOMAIN`).
- IAS15: шаг, чьи проходы расходятся (состояние в узле не конечно), отвергается и повторяется
  вчетверо короче, без вызова ускорения с таким состоянием, при любом масштабе времени.
- `myb_sclk_decode` печатает поля целыми числами (поле длиннее 18 цифр — `MYB_E_RANGE`) и
  принимает пустой текст (NULL) размера 0, чтобы сообщить длину (`MYB_E_OVERFLOW`).
- Эталон `test_repro` — `0xaa258e26488b41eb`; префиксы 22, 20 и 19 — `0xe234728b95d060b4`,
  `0x7fc24de1c1d525be`, `0x4352d2d0bc185d9e` (наклон решателя).
- Описания в заголовках исправлены по сверке паспортов моделей (отказы `myb_wh_start` и
  системы, читатели XML, свет на сетке, диапазоны JB2008 и Мендеса–Павлиса, высшие порядки
  ионосферы, ошибки Ламберта); поведение не менялось.

### Исправлено
- Публичные заголовки компилируются GCC в C23, его режиме по умолчанию с GCC 15: GCC отвергает
  стандартный атрибут после атрибута GNU, а в каждом объявлении
  `__attribute__((visibility("default")))` стоял перед `[[nodiscard]]` (Clang это принимает,
  поэтому первыми ошибку увидели привязки через hsc2hs на Linux). Теперь в C23 видимость
  разделяемой библиотеки пишется как `[[gnu::visibility("default")]]`; `test_gcc_consumer`
  проверяет заголовки и две программы, собранные gcc и g++.
- Фаза 21, найдено при проверке кода и в аудите после неё (у каждого случая регрессионный тест,
  который без исправления падает). Чтение DAF: слово вне файла в релизной сборке даёт NaN вместо
  чтения за образом. CK: окно типа 5 выходило на тик за конец (момент на тике между интервалами,
  сегмент из одного пакета); нулевое смещение пакета типа 4 указывало перед сегментом; бесконечные
  концы интервалов типа 2, число секунд на тик без требования положительности и поворот |w| t с
  переполнением (NaN при `MYB_OK`); начала интервалов типа 3 не на тиках; скорости и кватернионы с
  переполнением норм; рекуррентная схема Кленшоу типа 4, переполнявшаяся на пути к конечному
  значению; границы типа 6, не покрывавшие охват; `myb_ck_pointing_at` отвергает нечисловой ответ
  и не даёт скорости для сегмента без угловых скоростей. SCLK: разбиения с суммой за пределом
  double, таблица на отрицательных тиках, отрицательные поля, которые decode писал, а encode не
  читал. SPK и PCK: состояния с переполнением в метрах или по цепочке, вращения PCK с
  переполнением производной, полюса типа 15 из длинных или коротких векторов, NaN в опорной эпохе
  типа 10 и запросы типа 10 дальше 2³² с от их набора элементов (теперь `MYB_E_RANGE`;
  повреждённый охват интегрировал резонанс глубокого космоса сутками). Текстовые ядра и парсеры:
  длинные мантиссы с насыщенным показателем, TAI−UTC больше полусуток и метки за концом таблицы,
  значения EOP и приливов, доходившие до `MYB_CHECK` или до переполнения, переполнения рядов
  Пуассона, CIP, WGCCRE и ICGEM, кватернион системы TK с очень большими или малыми компонентами,
  таблицы IGRF с равными эпохами или больше 64, космическая погода вне физических пределов.
  Форматы: значения СИ CCSDS с переполнением, ключи ковариации не по порядку, сущности в
  COV_REF_FRAME и в классификации OMM; эпохи SP3 после EOF возвращались неинициализированными,
  спутник, указанный дважды, цикл обрезки перед началом буфера; день 366 невисокосного года в TLE;
  поля IONEX, приводившиеся до проверки, отсутствующий сосед с нулевым весом, из-за которого падал
  верный запрос; GPT3 и сетка S1/S2 читали перед узлами при долготах огромной величины. Элементы
  SGP4 далеко за пределами каталога дают `MYB_E_DOMAIN` (давали NaN). Писатели CK, SPK, PCK, SPK
  типа 10, CCSDS, SP3, TLE и таблиц отвергают то, что отвергают их читатели или что не читается
  обратно: повторённые начала интервалов, имена с пробелами в конце, значения, которые не
  печатаются обратно, эпохи и интервалы вне диапазона.
- Фаза 21. Читатели, найдено фаззингом: `myb_sclk_decode` мог выйти за буфер текста; моменты
  вне диапазона из ядер шкал времени; эпохи WGCCRE; эпохи SPK и PCK дальше 2⁵¹ с и пакеты SPK
  типа 10; непроверенные приведения в читателе EOP; длина имён XML CCSDS считалась до раскрытия
  сущностей; поля TLE, которые читатель брал, а писатель писал без должной строгости
  (`ext/sgp4`).
- `myb_sclk_decode` неверно печатал большие поля с рантаймом C Windows (`%.0f` из msvcrt даёт
  17 значащих цифр).
- `myb_fit_orbit` при точных наблюдениях никогда не сообщала о сходимости (взвешенное СКО
  падает до округления пропагации); теперь она останавливается и при СКО ниже допуска.
- Описание `myb_filter_start`: ковариация, которую она читает, — (6 + np + bias_count)²
  значений.
- Python: писатель SP3 хранил указатель, прочитанный обратно из поля структуры, — тот не
  удерживал массив спутников (запись случайно падала после работы сборщика).
- Rust: с набором MSVC крейты, зависящие от `myb-sys`, не находили библиотеку импорта;
  `build.rs` копирует её в `OUT_DIR` как `myb.lib`.
- Go: тест ключует кеш cgo хешем заголовков (устаревший объект связки читал за концом
  выросшей структуры); после установки другой библиотеки выполняйте `go clean -cache`.

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
