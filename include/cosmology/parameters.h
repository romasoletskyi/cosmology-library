#pragma once

struct PhysParameters {
    double me = 511000; // electron mass in ev
    double delta = 13.6; // hydrogen ionization energy in ev
    double alpha = 1 / 137.036; // fine structure constant

    double evInKelvin = 11605; // translate from K to ev as T -> T / ev_in_kelvin
    double critInEv = 8.616e-20; // proton concentration at a = 1, \omega_b = 1 in ev^3
    double thomsonInEv = 1.708e-15; // Thomson scattering cross-section in ev^-2
    double evInCosmoFrequency = 4.6905e32; // 1 ev corresponds to evInCosmoFrequency measured in 100 km / (Mpc s) units
    double cosmoFrequencyInHz = 3.2404e-18; // 100 km / (Mpc s) in 1 / s units
};

struct CosmoParameters {
    double omegaBaryon = 0.0223; // \omega_b = \Omega_b h^2 - physical baryon density parameter
    double omegaCold = 0.1188; // dark matter
    double omegaLambda = 0.3177; // dark energy
    double omegaPhoton = 2.475e-5; // photon
    double omegaNeutrino = 1.686e-5; // neutrino

    double h = 0.6774; // H0 = 100 h km / (Mpc s), where H0 - Hubble constant
    double cmb = 2.726; // cmb current temperature in kelvin
    double ns = 1; // scalar spectral index

    double cmbScale = 3e+4; // spectrum proportionality coefficient
};
