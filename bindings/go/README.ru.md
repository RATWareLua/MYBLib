# myb для Go

Лицо [myb](https://github.com/RATWareLua/MYBLib), библиотеки на C о движении тел в космосе,
для Go: пакет `github.com/RATWareLua/MYBLib/bindings/go/myb` — простой cgo над **статической
библиотекой** `libmyb.a`, которую cgo находит через pkg-config (модуль `myb-static`). Модуль ни
от какого другого модуля не зависит.

## Один вызов — одна задача

Переход из Go в C стоит десятки наносекунд, сколько бы работы ни было внутри, поэтому пакет —
только крупный слой библиотеки: каждый вызов — целая задача.

- **Модель сил, собранная из данных**: `ModelSpec` перечисляет тела (их GM, радиусы,
  светимости; состояния из `Ephemeris`; ориентации: постоянная, земная с наблюдёнными ПВЗ,
  двоичный PCK, модель WGCCRE) и члены — `PointMass`, `ThirdBody`, `Field` (`GravityField`),
  `Drag` (экспоненциальная атмосфера), `Radiation` (с тенями), `Schwarzschild`,
  `LenseThirring`, `DeSitter`, `Empirical`, — а также параметры, которые оценивает подгонка;
  `NewModel` превращает её в окружение и модель библиотеки за один вызов.
- **Прогноз** (`Model.Propagate`: состояния во многие моменты, матрицы перехода и
  чувствительности, импульсные манёвры), **поиск событий** (`Model.Events`: узлы, апсиды,
  расстояния, высоты, угол места над станцией, затмения), **подгонка орбиты** (`Model.Fit`:
  положения, состояния, дальность, скорость её изменения, углы со станций, со временем
  распространения света, систематическими сдвигами, априорной информацией) — по одному вызову.
- **Файлы**: сообщения OEM пишутся и читаются в обеих нотациях (`WriteOEM`, `ParseOEM`,
  `ReadOEMFile`), файлы SPK пишутся (`WriteSPK`, `EncodeSPK`) и читаются (`Ephemeris`).
- **Всё «по эпохам» — одним вызовом над срезом**: формы над массивами самой библиотеки C,
  дающие биты одиночных вызовов в цикле: `Ephemeris.States`, `TimeSystem.ConvertMany`,
  `FromUTCMany`, `ToUTCMany`, `TimesFromSeconds`, `CIP.TerrestrialMany`, `CIP.XYSMany`,
  `EOP.AtMany`, `Equinox.RotationMany`, `WGCCRE.RotationMany`, `PCK.RotationMany`,
  `EarthOrientation.Sample` (и прочие ориентации), `ApplyRotations`, `KeplerPropagateMany`,
  `KeplerSTMMany`, `StatesFromElements`, `ElementsFromStates`, … Одиночные вызовы тоже есть.

**Арифметики в Go нет.** Числа копируются между Go и структурами библиотеки, а её функции
вызываются в том порядке, в каком их вызывает программа на C, поэтому у числа те же биты, что у
тех же вызовов на C (компилятор Go на некоторых процессорах может слить `x*y + z` в одну
операцию — пакет не даёт ему такого случая). Тест проверяет это на орбите
`examples/low_orbit.c` против программы на C `tests/bindings/low_orbit_reference.c`: моменты 15
восходящих узлов, состояния в них, состояние через сутки и число шагов совпадают бит в бит.

**Функции Go** могут быть членами модели (`GoTerm`) и событиями (`GoEvent`): библиотека
вызывает их обратно через cgo при каждом вычислении — отдельный и откровенно медленный путь
(переход на каждый вызов, сотни тысяч вызовов на орбиту). Член сообщает о неудаче, возвращая
ошибку; вызов библиотеки останавливается с её статусом (статусом `*myb.Error`, `MYB_E_DOMAIN`
для любой другой ошибки), а возвращаемая ошибка оборачивает ошибку члена. Паника в такой
функции не перехватывается: она завершает программу.

**Владение.** Каждый объект библиотеки принадлежит значению Go с методом `Close`;
`runtime.AddCleanup` освободит его, если `Close` забыли, — только как страховка. Модель держит
то, чем пользуется (эфемериду, поле, систему времени), пока не закрыта сама: закрытие такого
объекта раньше сразу делает его недоступным из Go (`ErrClosed`) и освобождает, когда модель его
отпустит. Значениями можно пользоваться из нескольких горутин; модель обслуживает один вызов за
раз.

**Ошибки** — значения `*myb.Error`: статус библиотеки (`errors.Is(err, myb.ErrRange)`), номер
элемента среза, на котором случилась неудача (`Index`), строка, которую называет разборщик
(`Line`). Ничто не паникует из-за неудачи библиотеки; аргументы, которые библиотека сочла бы
нарушением контракта (момент вне диапазона, срезы разной длины, тело, которого нет в модели),
отвергаются с `MYB_E_DOMAIN` до того, как дойдут до неё.

## Использование

Установите библиотеку (статическую; `-DMYB_BUILD_SHARED=ON` может быть включён тоже) и
укажите на неё pkg-config:

    cmake --preset linux-release && cmake --build --preset linux-release
    cmake --install build/linux-release --prefix /opt/myb
    export PKG_CONFIG_PATH=/opt/myb/lib/pkgconfig       # там myb-static.pc

cgo нужен компилятор C в духе GNU; возьмите тот Clang, которым собрана библиотека. В Linux:

    export CC=clang

В Windows библиотека собрана для цели mingw-w64 со средой выполнения C из sysroot MSYS2
MINGW64; cgo должен компилировать и компоновать для того же:

    set CC=clang --target=x86_64-w64-windows-gnu --sysroot=C:/msys64/mingw64
    set CGO_LDFLAGS=-fuse-ld=lld
    set PKG_CONFIG_PATH=C:/myb/lib/pkgconfig

Модуль не опубликован (библиотека живёт только на GitHub): сошлитесь на рабочую копию в
`go.mod` программы,

    require github.com/RATWareLua/MYBLib/bindings/go v0.0.0
    replace github.com/RATWareLua/MYBLib/bindings/go => /path/to/MYBLib/bindings/go

```go
package main

import (
	"fmt"
	"log"

	"github.com/RATWareLua/MYBLib/bindings/go/myb"
)

func main() {
	const gm, radius, day = 3.986004418e14, 6378137.0, 86400.0
	c20 := -4.841695e-4 // -J2 / sqrt(5), вычисленный там, где биты важны (Go может слить)
	field, err := myb.NewGravityField(gm, radius, 2, []float64{1, 0, 0, c20, 0, 0}, make([]float64, 6))
	if err != nil {
		log.Fatal(err)
	}
	defer field.Close()
	model, err := myb.NewModel(myb.ModelSpec{
		Bodies: []myb.Body{{ID: 1, GM: gm}},
		Terms: []myb.Term{
			myb.Field{Name: "gravity", Body: 0, Field: field, NMax: -1},
			myb.Drag{Name: "drag", Body: 0, Object: myb.Object{Mass: 500, Area: 4, Cd: 2.2},
				EquatorialRadius: radius, Layers: []myb.Layer{{Base: 400e3, Density: 3e-12, Scale: 58e3}}},
		},
	})
	if err != nil {
		log.Fatal(err)
	}
	defer model.Close()
	y0, err := myb.StateFromElements(gm, myb.Elements{P: 6.778e6, E: 0.001, I: 0.9, RAAN: 0.5, ArgP: 1})
	if err != nil {
		log.Fatal(err)
	}
	nodes, err := model.Events(myb.EventSearch{T0: 0, Y0: y0, TEnd: day, Events: []myb.EventSpec{
		myb.PlaneCrossing{Body: myb.Origin, Normal: [3]float64{0, 0, 1}, Direction: myb.Rising},
	}})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(len(nodes.Events), "восходящих узлов; через сутки", nodes.Final)
}
```

## Тест

`test_go_binding` в ctest (регистрируется, когда CMake находит Go — переменную кэша `MYB_GO`,
иначе `go` в `PATH` или в `toolchains/go*/bin` дерева исходников — и pkg-config) собирает пакет
против дерева сборки: файл pkg-config `build/<preset>/bindings/go/pkgconfig/myb-static.pc`,
компилятор и компоновщик сборки передаются cgo (`tests/tools/run_go_test.cmake`). Go держит
свои кэши, рабочие и временные файлы тестов в `build/<preset>/go`; ничего не скачивает
(`GOTOOLCHAIN=local`, `GOPROXY=off`) и не ведёт телеметрии (`GOTELEMETRY=off`). Вручную, из
`bindings/go`:

    go test ./...       # с PKG_CONFIG_PATH, CC (и CGO_LDFLAGS в Windows), как выше

Тест низкой орбиты бит в бит выполняется, когда `MYB_GO_REFERENCE` называет программу
`go_low_orbit_reference` из дерева сборки, иначе пропускается.
