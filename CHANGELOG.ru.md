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
