/* myb_troposphere.h — delay of radio and optical signals in the Earth's neutral atmosphere as an
 * optional module of myb (target myb_troposphere; see NOTICE).
 *
 * The models of the IERS Conventions (2010), chapter 9 (G. Petit, B. Luzum (eds.), IERS
 * Technical Note 36, with its updates), written from the published equations:
 *   - zenith hydrostatic delay of radio signals (J. Saastamoinen 1972 in the form of J. L. Davis
 *     et al., Radio Sci. 20, 1593, 1985) and Saastamoinen's zenith wet delay;
 *   - zenith delays at optical wavelengths (V. B. Mendes, E. C. Pavlis, Geophys. Res. Lett. 31,
 *     L14602, 2004) and the mapping function FCULa (V. B. Mendes et al., Geophys. Res. Lett. 29,
 *     1414, 2002), for laser ranging;
 *   - mapping functions of radio signals: VMF1 from the coefficients a_h, a_w of the Vienna
 *     service (J. Boehm, B. Werl, H. Schuh, J. Geophys. Res. 111, B02406, 2006) and the
 *     empirical GMF (J. Boehm et al., Geophys. Res. Lett. 33, L07304, 2006), with the height
 *     correction of A. E. Niell (J. Geophys. Res. 101, 3227, 1996);
 *   - the empirical pressure and temperature GPT (J. Boehm, R. Heinkelmann, H. Schuh, J. Geod.
 *     81, 679, 2007) for stations without meteorological data;
 *   - the gradient mapping function of G. Chen, T. A. Herring (J. Geophys. Res. 102, 20489,
 *     1997) for horizontal gradients estimated elsewhere.
 * The slant delay is m_h(e) ZHD + m_w(e) ZWD + gradient term; it lengthens the light time of
 * myb/observe.h by delay / c. Angles in radians, heights in metres, pressures in hPa,
 * temperatures in kelvin unless stated, time as a modified Julian date (days, the seasonal
 * terms of the models need no better than a day). */
#ifndef MYB_TROPOSPHERE_H_
#define MYB_TROPOSPHERE_H_

#include "myb/core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Zenith hydrostatic delay of radio signals (m) for the surface pressure (hPa) at geodetic
 * latitude and ellipsoidal height: 0.0022768 P / (1 - 0.00266 cos 2 phi - 0.28e-6 h). */
MYB_API MYB_NODISCARD myb_status myb_tropo_zhd_saastamoinen(double pressure, double latitude,
                                                            double height, double* zhd);

/* Saastamoinen's zenith wet delay (m) for the surface temperature (K) and water-vapour pressure
 * (hPa): 0.002277 (1255 / T + 0.05) e. Accurate to a few centimetres; estimate it from the
 * data when better is needed. */
MYB_API MYB_NODISCARD myb_status myb_tropo_zwd_saastamoinen(double temperature,
                                                            double vapour_pressure, double* zwd);

/* Mendes-Pavlis zenith hydrostatic and non-hydrostatic delays (m) at optical wavelength
 * (micrometres; 0.3-1.2) for pressure and water-vapour pressure (hPa), geodetic latitude,
 * ellipsoidal height and CO2 content (ppm; the IERS Conventions use 375). */
MYB_API MYB_NODISCARD myb_status myb_tropo_zenith_optical(double latitude, double height,
                                                          double pressure, double vapour_pressure,
                                                          double wavelength, double co2,
                                                          double* zhd, double* zwd);

/* FCULa mapping function of optical signals (both components) for the elevation, geodetic
 * latitude, height above mean sea level and surface temperature (K). */
MYB_API MYB_NODISCARD myb_status myb_tropo_fcula(double elevation, double latitude, double height,
                                                 double temperature, double* mapping);

/* VMF1 hydrostatic and wet mapping functions for the coefficients ah, aw of the site at the
 * epoch, the modified Julian date, geodetic latitude, ellipsoidal height (0 gives the
 * function without the height correction, as for coefficients given at the site height
 * already) and zenith distance. */
MYB_API MYB_NODISCARD myb_status myb_tropo_vmf1(double ah, double aw, double mjd, double latitude,
                                                double height, double zenith_distance,
                                                double* mh, double* mw);

/* GMF hydrostatic and wet mapping functions at the modified Julian date, geodetic latitude,
 * longitude, ellipsoidal height and zenith distance. */
MYB_API MYB_NODISCARD myb_status myb_tropo_gmf(double mjd, double latitude, double longitude,
                                               double height, double zenith_distance, double* mh,
                                               double* mw);

/* GPT pressure (hPa), temperature (degrees Celsius) and geoid undulation (m) at the modified
 * Julian date, geodetic latitude, longitude and ellipsoidal height. */
MYB_API MYB_NODISCARD myb_status myb_tropo_gpt(double mjd, double latitude, double longitude,
                                               double height, double* pressure,
                                               double* temperature, double* undulation);

/* Delay (m) of the horizontal gradients gn (north) and ge (east) in m at elevation and azimuth
 * (from north towards east): m_g(e) (gn cos a + ge sin a), m_g = 1 / (sin e tan e + 0.0032). */
MYB_API MYB_NODISCARD myb_status myb_tropo_gradient_delay(double gn, double ge, double elevation,
                                                          double azimuth, double* delay);

#ifdef __cplusplus
}
#endif

#endif /* MYB_TROPOSPHERE_H_ */
