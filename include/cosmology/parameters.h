#pragma once

struct PhysParameters {

};

struct CosmoParameters {
    float omegaBaryon = 0.0223; // \omega_b = \Omega_b h^2 - physical baryon density parameter
    float omegaCold = 0.1188; // dark matter
    float omegaLambda = 0.3177; // dark energy
    float omegaPhoton = 2.475e-5; // photon
    float omegaNeutrino = 1.686e-5; // neutrino
    float h = 0.6774; // H0 = 100 h km / (Mpc s), where H0 - Hubble constant
};
