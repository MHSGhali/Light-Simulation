/* cie_data.h — declarations for the tabulated CIE data in cie_data.c. */
#ifndef LIGHTSIM_CIE_DATA_H
#define LIGHTSIM_CIE_DATA_H

/* CIE 1931 2-degree standard observer, 360-830 nm at 5 nm.
 * ls_cie_ybar IS the luminous efficiency function V(lambda). */
extern const int    ls_cie_count;
extern const double ls_cie_lambda[95];
extern const double ls_cie_xbar[95];
extern const double ls_cie_ybar[95];
extern const double ls_cie_zbar[95];

/* CIE daylight basis, 300-830 nm at 10 nm. */
extern const int    ls_cie_daylight_count;
extern const double ls_cie_daylight_lambda[54];
extern const double ls_cie_s0[54];
extern const double ls_cie_s1[54];
extern const double ls_cie_s2[54];

#endif /* LIGHTSIM_CIE_DATA_H */
