# FetchData.cmake — downloads the large data sets that are not kept in the repository into
# data/download (git-ignored). Files with a SHA-256 are fixed releases and are verified; files
# without one are living series (updated daily or weekly by their producers) and are fetched
# again with -DREFRESH=ON.
#
# Usage:  cmake [-DREFRESH=ON] [-DLARGE=ON] [-DONLY=<name>] -P cmake/FetchData.cmake
#         (Windows: scripts\fetch-data.cmd)
# The files of the second list (a gigabyte of ephemerides, the tests that use them skip
# without them) are fetched only with -DLARGE=ON or when named by -DONLY.
cmake_minimum_required(VERSION 3.28)

get_filename_component(_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(_dir "${_root}/data/download")
file(MAKE_DIRECTORY "${_dir}")

# name | URL | SHA-256 (empty: living series) | purpose
set(_files
  "de440.bsp|https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp|a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5|JPL DE440 planetary ephemeris (NAIF SPK, 1550-2650)"
  "testpo.440|https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/testpo.440|45084a7f406b8e41032e7d99dd8042d4df97a9f752a68c06d1a27af16693429f|JPL acceptance test values of DE440"
  "header.440|https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440|0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359|DE440 constants (AU, GM, ...)"
  "moon_pa_de440_200625.bpc|https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc|60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f|lunar principal-axis orientation of DE440 (binary PCK)"
  "pck00011.tpc|https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc|3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1|WGCCRE rotation models and radii (text PCK)"
  "GOCO06s.gfc|https://icgem.gfz-potsdam.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc|351d9d20b84cd2c0f52ce77146b1e3b774f408200b579ffaf98593cf3d271819|GOCO06s gravity field, degree 300 with trends and annual terms (ICGEM, CC BY 4.0)"
  "fes2004_Cnm-Snm.dat|https://iers-conventions.obspm.fr/content/chapter6/additional_info/tidemodels/fes2004_Cnm-Snm.dat|620dc48fd4d97baa14f532287b19788bd55dc2056e0b478ab0f928b951c194cf|FES2004 ocean tides as Stokes coefficients (IERS Conventions 2010, ch. 6)"
  "gpt3_5.grd|https://vmf.geo.tuwien.ac.at/codes/gpt3_5.grd|3bf35c611e7cd5c42a0650ee14652a8f36430127ec213ed7b7b867bd29245b3d|GPT3 empirical troposphere model, 5 degree grid (TU Wien)"
  "s1_s2_s3_cm_noib_grid.dat|https://vmf.geo.tuwien.ac.at/APL_products/TIDAL/s1_s2_s3_cm_noib_grid.dat|0b54049c98a28ad0116c2ae022dcc444fbce5dd666a486679a6584f4c7ad3940|S1, S2, S3 atmospheric tide loading, 1 degree grid, centre of mass (TU Wien, VMF data server)"
  "igrf14coeffs.txt|https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt|8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c|IGRF-14 geomagnetic field, Gauss coefficients 1900-2025 with secular variation (IAGA, NOAA NCEI)"
  "earth_latest_high_prec.bpc|https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc||ITRF93 Earth orientation from current EOP (binary PCK, updated daily)"
  "eopc04.1962-now|https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now||IERS EOP 20 C04 daily series (updated weekly)"
  "finals2000A.all|https://datacenter.iers.org/data/9/finals2000A.all||IERS Bulletin A with predictions (updated daily)"
)
set(_large
  "de440t.bsp|https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440t.bsp|4f003b4437a8b479aaca262716a779e3683413200adfd6353c4901aab5de800e|JPL DE440 with TT-TDB at the geocentre (152 MB)"
  "sb441-n373s.bsp|https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n373s.bsp|2143113282bfc2b2a0b0b4626125d4f84362339b5a8ae7eea40f4120ca8da10b|the 373 small bodies that perturb DE440 and DE441 (982 MB)"
  "SW-Last5Years.csv|https://celestrak.org/SpaceData/SW-Last5Years.csv||space weather: F10.7, Kp, Ap of the last five years with predictions (updated daily)"
)
foreach(_entry IN LISTS _large)
  string(REPLACE "|" ";" _f "${_entry}")
  list(GET _f 0 _name)
  if(LARGE OR (ONLY AND ONLY STREQUAL _name))
    list(APPEND _files "${_entry}")
  endif()
endforeach()

set(_failed "")
foreach(_entry IN LISTS _files)
  string(REPLACE "|" ";" _f "${_entry}")
  list(GET _f 0 _name)
  list(GET _f 1 _url)
  list(GET _f 2 _hash)
  list(GET _f 3 _what)
  if(ONLY AND NOT ONLY STREQUAL _name)
    continue()
  endif()
  set(_path "${_dir}/${_name}")
  if(EXISTS "${_path}")
    if(_hash)
      file(SHA256 "${_path}" _have)
      if(_have STREQUAL _hash)
        message(STATUS "${_name}: present, SHA-256 verified")
        continue()
      endif()
      message(STATUS "${_name}: SHA-256 mismatch, downloading again")
    elseif(NOT REFRESH)
      message(STATUS "${_name}: present (living series; -DREFRESH=ON to update)")
      continue()
    endif()
  endif()
  message(STATUS "${_name}: downloading ${_what}")
  set(_tmp "${_path}.part")
  if(_hash)
    file(DOWNLOAD "${_url}" "${_tmp}" TLS_VERIFY ON STATUS _status INACTIVITY_TIMEOUT 120
         EXPECTED_HASH SHA256=${_hash})
  else()
    file(DOWNLOAD "${_url}" "${_tmp}" TLS_VERIFY ON STATUS _status INACTIVITY_TIMEOUT 120)
  endif()
  list(GET _status 0 _code)
  if(NOT _code EQUAL 0)
    list(GET _status 1 _text)
    message(WARNING "${_name}: download failed (${_code}: ${_text})")
    file(REMOVE "${_tmp}")
    list(APPEND _failed "${_name}")
    continue()
  endif()
  file(RENAME "${_tmp}" "${_path}")
  file(SIZE "${_path}" _size)
  message(STATUS "${_name}: ${_size} bytes")
endforeach()

if(_failed)
  message(FATAL_ERROR "not downloaded: ${_failed}")
endif()
