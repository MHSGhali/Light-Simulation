#include "lightsim/units.h"
#include "lightsim/color.h"

ls_real ls_radiometric(const Spectrum *s) {
    return ls_spectrum_integrate(s);
}

ls_real ls_photometric(const Spectrum *s) {
    /* V(lambda) is the CIE ybar colour matching function, by definition. */
    return LS_KM_LM_PER_W * ls_spectrum_integrate_weighted(s, ls_cmf_ybar());
}

ls_real ls_luminous_efficacy_band(const Spectrum *s) {
    ls_real p = ls_radiometric(s);
    return p > 0.0 ? ls_photometric(s) / p : 0.0;
}

ls_real ls_luminous_efficacy_total(const Spectrum *s, ls_real total_radiant_power) {
    return total_radiant_power > 0.0 ? ls_photometric(s) / total_radiant_power : 0.0;
}

ls_real ls_quantity_value(const Spectrum *s, LsUnitSystem sys) {
    return sys == LS_UNITS_PHOTOMETRIC ? ls_photometric(s) : ls_radiometric(s);
}

const char *ls_quantity_unit(LsQuantity q, LsUnitSystem sys) {
    const bool p = (sys == LS_UNITS_PHOTOMETRIC);
    switch (q) {
        case LS_Q_FLUX:       return p ? "lm"      : "W";
        case LS_Q_IRRADIANCE: return p ? "lx"      : "W/m^2";
        case LS_Q_INTENSITY:  return p ? "cd"      : "W/sr";
        case LS_Q_RADIANCE:   return p ? "cd/m^2"  : "W/(m^2 sr)";
        case LS_Q_EXITANCE:   return p ? "lm/m^2"  : "W/m^2";
    }
    return "?";
}

const char *ls_quantity_symbol(LsQuantity q, LsUnitSystem sys) {
    const bool p = (sys == LS_UNITS_PHOTOMETRIC);
    switch (q) {
        case LS_Q_FLUX:       return p ? "Phi_v" : "Phi_e";
        case LS_Q_IRRADIANCE: return p ? "E_v"   : "E_e";
        case LS_Q_INTENSITY:  return p ? "I_v"   : "I_e";
        case LS_Q_RADIANCE:   return p ? "L_v"   : "L_e";
        case LS_Q_EXITANCE:   return p ? "M_v"   : "M_e";
    }
    return "?";
}
