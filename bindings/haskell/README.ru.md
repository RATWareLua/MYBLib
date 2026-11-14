# myb для Haskell

Лицо [myb](https://github.com/RATWareLua/MYBLib), библиотеки на C о движении объектов в
пространстве, для Haskell. Пакет cabal с двумя слоями:

- `Myb.C` — весь интерфейс C: каждая функция, структура, перечисление и константа
  `include/myb` под именами заголовков C, сгенерированные из них (`tools/Generate.hs`).
  Библиотеки `ext/` — публичные подбиблиотеки `myb:sgp4`, `myb:nrlmsise00`,
  `myb:troposphere`, `myb:ionosphere`, `myb:jb2008`.
- `Myb` — интерфейс в манере языка: ядра и шкалы времени, Кеплер, эфемериды SPK,
  гравитационные поля, окружение и модель сил с готовыми членами, члены, написанные на Haskell,
  орбита с состояниями, событиями, импульсами и матрицей перехода; а также N тел и системы тел,
  оценивание, ограниченная задача трёх тел, наведение и оптимальные перелёты, тяга, начальное
  определение орбиты, относительное движение, сближения, неопределённость, ориентация,
  наблюдаемые, файлы (числа, таблицы, запись SPK и CK, SP3, OEM, OPM, OMM, OCM, бортовые часы),
  формы для массивов; двухстрочные элементы и SGP4 — в подбиблиотеке `myb:tle`. Объект
  библиотеки освобождает его финализатор (`ForeignPtr`); вызов, который может не удаться,
  возвращает `Result a`, то есть `Either Error a` со статусом библиотеки; члены, поставщики,
  события, действия, измерения, моменты сил, невязки и функции нелинейной программы — функции
  Haskell.

На Haskell не делается никакой арифметики: число имеет биты тех же вызовов, сделанных на C.
Тест строит орбиту `examples/low_orbit.c` из `Myb` и из функций C и сравнивает 15 восходящих
узлов, конечное состояние и 776 шагов бит в бит, а в каждой области остального сравнивает
верхний слой с теми же вызовами через `Myb.C`.

Пакету нужны GHC и cabal-install и больше ничего: он зависит только от загрузочных пакетов GHC
(`base` и `bytestring`, а также `containers`, `directory`, `filepath`, `process`,
`transformers` для генератора, теста и примера), а `cabal.project` не называет ни одного
репозитория, так что не нужны ни `cabal update`, ни сеть. Проверено с GHC 9.10.3 и cabal
3.16.1.0.

## Компоновка

Пакет компонует `libmyb` из каталога, указанного через `--extra-lib-dirs`: статическую
библиотеку или библиотеку импорта разделяемой — ту, что компоновщик там найдёт (если в каталоге
обе, он берёт разделяемую; чтобы скомпоновать статически, укажите каталог с одной `libmyb.a`).

- **Linux и другие**: статическая библиотека. Программе во время работы ничего не нужно.
- **Windows**: разделяемая библиотека `libmyb.dll` (её собирает пресет `windows-release`).
  GHC под Windows приносит свой набор инструментов mingw-w64, собранный для UCRT (`-D_UCRT`,
  DLL `api-ms-win-crt-*`), а библиотеку собирает Clang с корнем MSYS2 MINGW64, для
  `msvcrt.dll`. Статический архив в программу GHC не компонуется: ему нужна защита стека
  (`__stack_chk_guard`, `__stack_chk_fail`) более нового mingw-w64 из MSYS2, которой нет в
  инструментах GHC, да и объекты, собранные для одной среды выполнения C, с другой смешивать не
  следует. DLL несёт свою среду выполнения C (она импортирует только `KERNEL32.dll` и
  `msvcrt.dll`), и через интерфейс ничего от среды C не проходит: заголовки используют только
  `stddef.h` и `stdint.h`, никакого `FILE`, а что библиотека выделила, она и освобождает. Чтобы
  запустить программу, добавьте каталог `libmyb.dll` (и `libmyb_sgp4.dll`, ...) в `PATH`.

## Использование

Соберите библиотеку (`cmake --preset linux-release`, `cmake --build --preset linux-release`;
под Windows — пресет `windows-release`), затем из `bindings/haskell`:

    cabal build --extra-lib-dirs=/путь/к/библиотеке
    cabal test  --extra-lib-dirs=/путь/к/библиотеке    # с MYB_DATA_DIR, MYB_TEST_DATA_DIR

GHC и cabal берутся из `PATH`; другой GHC указывается через `-w /путь/к/ghc` (его `hsc2hs`
находится рядом с ним). Программа другого проекта берёт пакет по пути в своём `cabal.project`
(`packages: путь/к/bindings/haskell`) и `build-depends: myb` (и `myb:sgp4`, ... для библиотек
`ext/`, они компонуются как `libmyb_sgp4`, ...). cabal предупреждает, что каталоги заголовков
пакета лежат вне его (`../../include`): пакет собирается из репозитория, как и другие привязки.

Пример `examples/LowOrbit.hs` (исполняемый файл `myb-example`, его собирает и запускает тест):

```haskell
import Control.Monad.IO.Class (liftIO)
import Control.Monad.Trans.Except (ExceptT (..), runExceptT)
import Myb (State (..), Vec3 (..))
import qualified Myb
import System.Exit (exitFailure)

main :: IO ()
main = do
  let gm = 3.986004418e14
      radius = 6378137.0
  done <- runExceptT $ do
    field <- ExceptT (Myb.newGravityField gm radius 2 [1, 0, 0, -1.08263e-3 / sqrt 5, 0, 0] (replicate 6 0))
    env <- ExceptT (Myb.newEnvironment (Myb.timeFromSeconds 0))
    planet <- ExceptT (Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 1, Myb.bodyGm = gm})
    model <- ExceptT Myb.newForceModel
    _ <- ExceptT (Myb.addField model "gravity" env planet field 0 Nothing)
    _ <- ExceptT (Myb.addDrag model "drag" env planet (Myb.ExponentialAtmosphere radius 0 [(400e3, 3e-12, 58e3)])
                    (Myb.Object 500 4 1.3 2.2))
    orbit <- ExceptT (Myb.newOrbit model Myb.defaultOrbitOptions)
    y0 <- ExceptT (pure (Myb.stateFromElements gm (Myb.Elements 6.778e6 0.001 0.9 0.5 1.0 0)))
    ExceptT (Myb.start orbit 0 y0 86400)
    node <- ExceptT (Myb.nextEvent orbit (\_ y -> pure (let Vec3 _ _ z = position y in z)) 86400 Myb.Rising)
    liftIO (print node)                                 -- Just (время, состояние) или Nothing
    ExceptT (Myb.stateAt orbit 86400)
  either (\e -> putStrLn (Myb.errorText e) >> exitFailure) print done
```

(`ExceptT` — из `transformers`, загрузочного пакета GHC; самой библиотеке пакета нужен только
`base`. `ExceptT` лишь несёт `Either` каждого вызова: ничего не бросается.)

## Модули

| модуль | |
|---|---|
| `Myb` | весь верхний слой |
| `Myb.Error` | `Status`, `Error`, `Result`, `statusName`, `errorText` |
| `Myb.Types` | `Vec3`, `State`, `Elements`, `Calendar`, `Time` и их структуры C |
| `Myb.Time` | `KernelPool`, `TimeSystem`: `convertTime`, `fromUtc`, `toUtc` |
| `Myb.Kepler` | `stateFromElements`, `elementsFromState`, `keplerPropagate` (чистые) |
| `Myb.Ephemeris` | `Ephemeris`: `loadEphemerisFile`, `ephemerisState`, `segments` |
| `Myb.Gravity` | `GravityField`: `newGravityField`, `loadIcgemFile` |
| `Myb.Environment` | `Environment`, `Body`: `addOrigin`, `addFromEphemeris`, `addWith` (поставщик на Haskell), `bodyState` |
| `Myb.Force` | `ForceModel`: `addPointMass`, `addThirdBody`, `addField`, `addRadiation`, `addDrag`, `addSchwarzschild`, `addTerm` (член на Haskell), `addRaw`, `acceleration`, `setBreakScheme` (`BreaksRetake`, `BreaksPredict`) |
| `Myb.Orbit` | `Orbit`: `start`, `stateAt`, `statesAt`, `nextEvent` (событие на Haskell), `kick`, `partials`; `OrbitOptions` с `timescale` IAS15 (`TimescaleJerk` для контакта) |
| `Myb.Message` | OEM: запись (с ковариациями) и чтение (заголовок, состояния, ускорения, ковариации, из текста или файла); разбор и запись OPM и OMM (`parseOpm`, `opmText`, `parseOmm`, `ommText`); OCM: чтение (блоки, значения, строки, состояния траектории) и запись (`writeOcm`, `ocmFormatLine`) |
| `Myb.Source` | поставщики: `fromEphemeris` (тело набора SPK, передаётся как его собственный поставщик на C), `fromFunction` (функция Haskell), `fromC`; `orientationFunction`, `orientationC` |
| `Myb.NBody` | `NBody`: `accelerations`, `accelerationsXV`, `setInteraction`, `setRelativity`, `invariants`, `pnInvariants`, `setApproximation` (FMM), `setPool`; `Wh` (Уиздом-Холман): `whStart`, `whStep`, `whIntegrate`, `whState`, `whDense`, `whStatistics` |
| `Myb.System` | `System`: `setOrientation`, `setField`, `setRotation`, `setTide`, `addFieldPair`, `addFigures`, `addTide`, `addDistortion`, `addPnPair`, `addAction` (действие на Haskell), `addLightAction`, `addFrameDragging`, `setMassRate`, `addGuided`, `addContact`, `systemAcceleration`, `systemInertia`, `systemPartials`, упаковка вращений, спинов и пружин трения, которое залипает (`packContact`, `unpackContact`); `pnOrders` (члены пары PN по именам, среди них `PnSs3`, `PnSorr5`, `PnSotail`, `PnTailx`); `SystemIntegrator` (IAS15); `advance`, `merge`, `bounce` |
| `Myb.Estimation` | `Station`, `Measure` (готовые, со временем света, вне центра масс или функция Haskell), `evaluateMeasure`, `Observation`, `fitOrbit`, `Filter` (`filterStart`, `filterUpdate`, `filterPredict`, `filterSmooth`, `filterStep`, сдвиги) |
| `Myb.Cr3bp` | точки Лагранжа, Якоби, `correct`, `monodromy`, `stability`, `manifold`, `richardson`, семейства (`orbitCorrect`, `orbitContinue`), торы (`torusFromOrbit`, `torusContinue`, `torusState`, `torusStability`, `torusManifold`), `propagate` до сечения, системы отсчёта реальных тел, `addCr3bp` (член) |
| `Myb.Transfer` | `bplane`, сопряжённые конические сечения, `shoot` (невязка на Haskell), `lambert` |
| `Myb.Optimal` | `indirectSolve`, `indirectSample`, `nlpSolve` (программа на Haskell), Симс и Фланаган, цепочки участков (`Chain`, пролёты с импульсом, `chainFlyby`, `chainSize`, `chainGuess`, `chainEvaluate`, `chainSolve`) |
| `Myb.Thrust` | дуги тяги (`thrustInit`, `addThrust`), `equinoctialRates`, Эдельбаум, Q-закон (`qlawDirection`, `addQlaw`, `qlawFly`) |
| `Myb.Iod` | `gibbs`, `herrickGibbs`, `gauss`, `laplace`, `doubleR`, `gooding` |
| `Myb.Relative` | относительные состояния, их матрицы перехода, `rendezvous`, относительные орбитальные элементы |
| `Myb.Conjunction` | `approachFind` (два поставщика), `approachOrbits`, плоскость сближения, вероятности, `conjunctionFind`, долгие сближения |
| `Myb.Uncertainty` | `Variational`, `covarianceTransform`, `covarianceRtn`, `Rng`, `sampleStatistics` |
| `Myb.Contact` | ближайшая точка поверхности, `contactForce`, `contactForceStick` (трение, которое залипает: прогибы на входе, их скорости на выходе), восстановление и демпфирование; `Friction` для поля `kind` закона |
| `Myb.Attitude` | кватернионы, углы Эйлера, главные оси, скорости твёрдого тела и тела с ядром, моменты сил, точное свободное вращение, `rotationSplit` (момент на Haskell) |
| `Myb.Observe` | `lightTime`, `roundTrip`, задержка Шапиро, аберрация, отклонение, направления, плазма |
| `Myb.Files` | чтение и печать чисел, таблицы, `SpkWriter` (типы 2, 3, 9, 13, сырые), чтение и запись SP3 |
| `Myb.Ck` | `Sclk`, `Ck` (`ckPointing`, `ckOrientation` как поставщик), `CkWriter` (типы 1–6) |
| `Myb.Frames` | формы для массивов из myb/frames.h над её объектами из `Myb.C` |
| `Myb.Tle` | (подбиблиотека `myb:tle`) двухстрочные элементы, OMM, SGP4, подгонка, SPK типа 10 |
| `Myb.C` | весь интерфейс C; `Myb.C.<Заголовок>` — по модулю на заголовок (`Myb.C.Orbit` — это `myb/orbit.h`) |
| `Myb.C.Ext.<Имя>` | библиотека `ext/` в подбиблиотеке `myb:<имя>` |

Формы для массивов (`timeFromSecondsMany`, `convertTimeMany`, `keplerPropagateMany`,
`ephemerisStateMany`, ...) — функции над списками рядом с одиночными формами; неудача называет
элемент, на котором она случилась ("element 3"). Структуры параметров и результатов остального
интерфейса (орбита CR3BP, результат цепочки, путь света, ...) — записи `Myb.C`, которые
реэкспортируют модули, принимающие или отдающие их; перечисления из нескольких значений,
которые принимает вызов, — типы-суммы (`Interaction`, `WhScheme`, `Symmetry`, ...).

Слои свободно смешиваются: каждый объект отдаёт свой `ForeignPtr` (`forceModelPtr`,
`environmentPtr`, `systemPtr`, ...), а `forceModelKeep`, `environmentKeep`, `systemKeep`
привязывают время жизни того, на что указывает член, тело или уровень, добавленные через `Myb.C`.

### Что остаётся на нижнем слое

Каждая публичная функция dyn.h, system.h, estim.h, cr3bp.h, transfer.h, lambert.h, optimal.h,
thrust.h, iod.h, relative.h, conjunction.h, uncert.h, attitude.h, observe.h, io.h, ccsds.h,
ck.h, записи SPK и SP3 из kin.h, ext/sgp4 и форм для массивов time.h, kin.h, frames.h и
kepler.h достижима из `Myb`: вызывается или передаётся как указатель на функцию, которым она и
задумана (измерения, `myb_cr3bp_term`, `myb_thrust_acceleration`, `myb_term_qlaw`,
`myb_system_light_action`, `myb_orientation_ck`, `myb_spk_body_state`). То, у чего нет верхнего
слоя, берётся через `Myb.C`: заголовки вне этого списка (объекты myb/frames.h — CIP, EOP,
равноденствие, WGCCRE, PCK, граф систем отсчёта, — их формы для массивов принимают их как
`ForeignPtr`; модели поверхности и сетки surface.h; функции сил force.h; член контакта из contact.h, site.h,
weather.h, events.h, integ.h сверх интегратора системы, grav.h и env.h сверх пути выше),
библиотеки ext/, кроме sgp4, и распределители памяти: верхний слой всегда берёт распределитель по
умолчанию.

### Нижний слой

- Функция — внешний импорт под её именем C: `myb_orbit_state :: Ptr Myb_orbit -> CDouble
  -> Ptr CDouble -> IO Myb_status`; её адрес — `p'myb_orbit_state`. Каждый вызов — в IO и
  `safe`, поскольку любой из них может вызвать Haskell обратно.
- Параметр-указатель — `Ptr` его элемента; параметр, записанный как массив массивов
  (`double points[5][3]`), — `Ptr CDouble` на его числа подряд.
- Структура `myb_x` — запись `Myb_x` с полями `myb_x'<поле>` (массив — список), экземпляром
  `Storable`, размер, выравнивание и смещения которого hsc2hs берёт из заголовка, и `zero` —
  значением из одних нулевых битов: `zero {myb_point_mass_term'gm = 3.986004418e14}`. Объект,
  известный только по указателю (`myb_orbit`), — пустой тип.
- Перечисление — newtype над целым типом из заголовка; его константы, константы анонимных
  перечислений и целые `#define` — синонимы образцов под именами C (`MYB_OK`, `MYB_ODM_XML`,
  `MYB_ENV_MAX_OCCULTERS`).
- Тип указателя на функцию `myb_t` — это `Myb_t = FunPtr (...)` с `mk'myb_t` (функция Haskell
  в роли такого указателя; освобождается через `freeHaskellFunPtr`) и `call'myb_t`.
- Внешний импорт не передаёт структуру по значению. Функции, которые принимают или возвращают
  её (`myb_time_from_seconds`, `myb_timesys_convert`, ...), — функции Haskell с тем же именем
  поверх функции C пакета (`generated/cbits`), которая принимает структуру по указателю и
  вызывает библиотеку; ничего там не вычисляется. Обратный вызов, принимающий структуру по
  значению (`myb_state_function`, ...), идёт через трамплин пакета: `mk'myb_state_function`
  возвращает функцию C и контекст, который ей передать.

### Владение и неудачи

Объект уничтожает его финализатор, когда объект становится недостижим; модель сил держит живыми
окружения и поля, на которые указывают её члены (структуры членов — её собственные), окружение —
свои эфемериды и функции Haskell своих поставщиков, орбита — свою модель, так что любой из них
можно отпустить в любом порядке. Финализатор GHC не обязательно выполняется при завершении
программы; тогда память забирает система. У объекта библиотеки рабочая память на один поток:
пользуйтесь им из одного потока за раз.

Система держит данные ей поля, поставщики, действия и точки контакта; орбита, фильтр,
вариационный объект и интегратор — свою модель или систему; станция — своё окружение. Объекты
файлов (таблицы, SP3, OEM, OCM, CK, записывающие) ничего не держат.

Пакет ничего не бросает и не ловит. Функция Haskell, которую вызывает библиотека (член,
поставщик, событие, действие, измерение, момент силы, невязка, нелинейная программа, гауссово
распределение долгого сближения), сообщает о неудаче своим результатом (`Left Range`), и он возвращается как
статус шедшего вызова. Исключение, покинувшее такую функцию, не ловится: оно завершает программу
(обработчик GHC для вызова из C печатает его и выходит), как нарушенный контракт сразу
останавливает библиотеку (`MYB_CHECK`); через кадры C оно никогда не раскручивается.

## Генератор

`tools/Generate.hs` (исполняемый файл пакета `myb-generate`) пишет `generated/` и списки
модулей и подбиблиотек в `myb.cabal` между метками `-- BEGIN myb-generate` и
`-- END myb-generate`. Сгенерированные файлы хранятся в репозитории; чтобы написать их заново
после изменения заголовков, из `bindings/haskell`:

    cabal run myb-generate -- ../..                    # clang в PATH, или --clang <clang>, или MYB_CLANG
    cabal run myb-generate -- --check ../..           # выходит с 1, если файлы устарели

Заголовки проходят через препроцессор Clang (`-E`, freestanding, с его маркерами строк), так
что каждое объявление известно по своему заголовку; небольшой разборщик объявлений C читает их
в структуры, перечисления, указатели на функции и функции, и каждый заголовок становится
модулем hsc2hs. В тексте нет ни размеров, ни смещений, ни путей: он одинаков на любой системе, а
структуры раскладывает hsc2hs по заголовкам при сборке пакета. Объявление, которого читатель не
знает, останавливает его с сообщением; ничего не пропускается молча.

Он написан на Haskell, а не на Python, как `make_cdef.py`: тогда проверке того, что файлы
согласны с заголовками, нужны лишь инструменты, которые пакету нужны и так (GHC, cabal и Clang,
собирающий библиотеку), и ей хватает загрузочных пакетов GHC.

## Тест

`ctest -R test_haskell_binding` (регистрируется, когда найдены GHC и cabal: `MYB_GHC`,
`MYB_CABAL`, их ищут в `PATH`, в `toolchains/` дерева исходников или домашнего каталога и в
`~/.ghcup/bin`) выполняет `tools/test.cmake`: генератор с `--check`, набор тестов, затем
пример. Всё, что пишет cabal — его настройки, хранилище, журналы и сборка, — остаётся в
`build/<пресет>/haskell` (`CABAL_DIR`). `MYB_HASKELL_LINK` (`static` или `shared`, под Windows
по умолчанию) выбирает библиотеку; она копируется одна в `build/<пресет>/haskell/lib`, а то, что
пишет тест (файл SPK, читаемый обратно), попадает в `build/<пресет>/haskell/output`.

Набор тестов (333 проверки) содержит: низкую орбиту примера бит в бит с вызовами C; неудачи как
значения со статусом и строкой разборщика; члены, поставщики и события как функции Haskell и их
неудачи; члены давления излучения и релятивистской поправки; то, что модель держит своё
окружение живым и оба финализируются, когда их отпустили (наблюдается через их финализаторы);
шкалы времени, набор SPK и тело из него, сообщения OEM в обеих нотациях; нижний слой
(константы, структура, переданная и возвращённая по значению, обратный вызов, принимающий её,
`myb:sgp4`); и, вторым процессом, то, что исключение в члене или событии завершает процесс.
Затем, область за областью, верхний слой против тех же вызовов через `Myb.C` (числа — по битам,
структуры — по байтам, тексты и образы файлов — по символам и байтам) и неудачи как `Left` со
статусом: точечные массы, их инварианты, метод быстрых мультиполей, интегратор Уиздома-Холмана;
система, её интегратор, продвижение, действие; подгонка, фильтр, дальность со станции,
измерение на Haskell; точки Лагранжа, исправленная гало-орбита, её монодромия, сечение, участки,
тор с его точками и спектром; Ламберт, Гиббс, Клохесси-Уилтшир, относительное движение,
Эдельбаум, Q-закон, пролёт с импульсом, Симс и Фланаган, нелинейная программа и пристрелка на
Haskell, сближение и его вероятность, сближения двух поставщиков на Haskell, случайные числа,
выборочная статистика; кватернионы, углы Эйлера, свободное вращение, расщепляющий интегратор;
аберрация, направление, время света до тела SPK; числа, таблица, файлы SPK и CK, записанные
байт в байт и прочитанные обратно, SP3, OPM, OCM, OEM с ковариацией, ориентация CK, часы,
переменные пула, SGP4, TLE через OMM; формы для массивов против одиночных форм и элемент, на
котором случилась неудача; схему разрывов и орбиту со шкалой времени рывка, силу трения,
которое залипает, и гладкого, упаковку пружин контакта в системе, новые члены PN.
