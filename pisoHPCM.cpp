#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <cmath>

bool warnings = false;

#include "liquid_sodium.h"

// =======================================================================
//                        TDMA SOLVER
// =======================================================================

std::vector<double> solveTridiagonal(const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d) {

    const int n = b.size();
    std::vector<double> c_star(n), d_star(n), x(n);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = c[i] / m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = d_star[i] - c_star[i] * x[i + 1];

    return x;
}

// =======================================================================
//                        GLOBAL / NUMERICAL SETUP
// =======================================================================

constexpr int    N = 20;
constexpr double L = 1.0;
constexpr double dz = L / N;
constexpr double dt = 1.0e-3;

constexpr int tot_time_steps = 200;

constexpr int tot_outer_x = 20;
constexpr int tot_inner_x = 10;

constexpr double outer_tol_x = 1e-6;
constexpr double inner_tol_x = 1e-6;

constexpr double K = 1.0e-5;               // permeability
constexpr double CF = 0.0;                  // Forchheimer
constexpr double rhie_chow_on_off_x = 1.0;

// =======================================================================
//                        FIELD VARIABLES
// =======================================================================

std::vector<double>
    u_x(N, -0.01),
    u_x_old = u_x,
    T_x_bulk(N, 700.0),
    T_x_bulk_old = T_x_bulk,
    rho_x(N, liquid_sodium::rho(700.0)),
    p_x(N, 10000.0);

std::vector<double> p_storage_x(N + 2, 1000.0);         /// Wick padded pressure vector for R&C correction [Pa]
double* p_padded_x = &p_storage_x[1];           /// Poìnter to work on the wick pressure padded storage with the same indes

std::vector<double> p_prime_x(N, 0.0);    /// Wick correction pressure field [Pa]

std::vector<double> Gamma(N, 0.0);
std::vector<double> Q(N, 0.0);

constexpr double z_evap_start = 0.0;   // [m]
constexpr double z_evap_end = 0.3;     // [m]
constexpr double z_cond_start = 0.7;   // [m]
constexpr double z_cond_end = 1.0;     // [m]

constexpr double Q_tot = 1e8;         // total heat input [W]
const double m_dot = 100.0;               // [kg/s]

// uniform volumetric mass source in evap/cond
const double L_evap = z_evap_end - z_evap_start;
const double L_cond = z_cond_end - z_cond_start;

// BCs
constexpr double u_inlet_x = 0.0;
constexpr double u_outlet_x = 0.0;
constexpr double p_outlet_x = 0.0;

/// The coefficient bXU is needed in momentum predictor loop and pressure correction to estimate the velocities at the faces using the Rhie and Chow correction
std::vector<double>
aXU(N, 0.0),                                                /// Lower tridiagonal coefficient for wick velocity
bXU(N, liquid_sodium::rho(T_x_bulk[0])* dz / dt
    + 2 * liquid_sodium::mu(T_x_bulk[0]) / dz),                    /// Central tridiagonal coefficient for wick velocity
    cXU(N, 0.0),                                                /// Upper tridiagonal coefficient for wick velocity
    dXU(N, 0.0);                                                /// Known vector coefficient for wick velocity

std::ofstream v_out("velocity.dat");
std::ofstream p_out("pressure.dat");
std::ofstream T_out("temperature.dat");

// =======================================================================
//                                MAIN
// =======================================================================

int main() {

    for (int i = 0; i < N; ++i) {

        const double z = i * dz;

        if (z >= z_evap_start && z <= z_evap_end) {
            Q[i] = Q_tot;
            Gamma[i] = m_dot;
        }
        else if (z >= z_cond_start && z <= z_cond_end) {
            Q[i] = -Q_tot;
            Gamma[i] = -m_dot;
        }
    }

    for (int n = 0; n < tot_time_steps; ++n) {

        u_x_old = u_x;
        T_x_bulk_old = T_x_bulk;

        double u_error_x = 1.0;
        int outer_x = 0;

        while (outer_x < tot_outer_x && u_error_x > outer_tol_x) {

            // ===========================================================
            // MOMENTUM PREDICTOR (UPWIND + RC, bXU = a_P completo)
            // ===========================================================
            for (int i = 1; i < N - 1; ++i) {

                // properties
                const double rho_P = liquid_sodium::rho(T_x_bulk[i]);
                const double rho_L = liquid_sodium::rho(T_x_bulk[i - 1]);
                const double rho_R = liquid_sodium::rho(T_x_bulk[i + 1]);
                const double rho_P_old = liquid_sodium::rho(T_x_bulk_old[i]);

                const double mu_P = liquid_sodium::mu(T_x_bulk[i]);
                const double mu_L = liquid_sodium::mu(T_x_bulk[i - 1]);
                const double mu_R = liquid_sodium::mu(T_x_bulk[i + 1]);

                const double D_l = 0.5 * (mu_P + mu_L) / dz;
                const double D_r = 0.5 * (mu_P + mu_R) / dz;

                const double invbX_L = 1.0 / bXU[i - 1] + 1.0 / bXU[i];
                const double invbX_R = 1.0 / bXU[i + 1] + 1.0 / bXU[i];

                // Rhie–Chow corrections for face velocities (4th-order stencil like your "second")
                const double rc_l = -invbX_L / (8.0 * dz) *
                    (p_padded_x[i - 2] - 3.0 * p_padded_x[i - 1] + 3.0 * p_padded_x[i] - p_padded_x[i + 1]);

                const double rc_r = -invbX_R / (8.0 * dz) *
                    (p_padded_x[i - 1] - 3.0 * p_padded_x[i] + 3.0 * p_padded_x[i + 1] - p_padded_x[i + 2]);

                // face velocities (avg + RC)
                const double u_l_face = 0.5 * (u_x[i - 1] + u_x[i]) + rhie_chow_on_off_x * rc_l;
                const double u_r_face = 0.5 * (u_x[i] + u_x[i + 1]) + rhie_chow_on_off_x * rc_r;

                // upwind densities at faces
                const double rho_l = (u_l_face >= 0.0) ? rho_L : rho_P;
                const double rho_r = (u_r_face >= 0.0) ? rho_P : rho_R;

                const double F_l = rho_l * u_l_face; // [kg/(m2 s)]
                const double F_r = rho_r * u_r_face; // [kg/(m2 s)]

                aXU[i] =
                    - std::max(F_l, 0.0)
                    - D_l;
                cXU[i] =
                    - std::max(-F_r, 0.0)
                    - D_r;
                bXU[i] =
                    + std::max(F_r, 0.0) + std::max(-F_l, 0.0)
                    + rho_P * dz / dt
                    + D_l + D_r
                    + mu_P / K * dz
                    + CF * mu_P * dz / sqrt(K) * abs(u_x[i]);
                dXU[i] =
                    -0.5 * (p_x[i + 1] - p_x[i - 1])
                    + rho_P_old * u_x_old[i] * dz / dt;
            }

            // BCs: same structure as your "second"
            const double D_first = liquid_sodium::mu(T_x_bulk[0]) / dz;
            const double D_last = liquid_sodium::mu(T_x_bulk[N - 1]) / dz;

            aXU[0] = 0.0; cXU[0] = 0.0;
            bXU[0] = liquid_sodium::rho(T_x_bulk[0]) * dz / dt + 2.0 * D_first;
            dXU[0] = bXU[0] * u_inlet_x;

            aXU[N - 1] = 0.0; cXU[N - 1] = 0.0;
            bXU[N - 1] = liquid_sodium::rho(T_x_bulk[N - 1]) * dz / dt + 2.0 * D_last;
            dXU[N - 1] = bXU[N - 1] * u_outlet_x;

            u_x = solveTridiagonal(aXU, bXU, cXU, dXU);

            // ===========================================================
            // INNER PISO ITERATIONS
            // ===========================================================
            double p_error_x = 1.0;
            int inner_x = 0;

            while (inner_x < tot_inner_x && p_error_x > inner_tol_x) {

                // -------------------------------------------------------
                // CONTINUITY SATISFACTOR: assemble pressure correction
                // -------------------------------------------------------
                std::vector<double> aXP(N, 0.0), bXP(N, 0.0), cXP(N, 0.0), dXP(N, 0.0);

                for (int i = 1; i < N - 1; ++i) {

                    const double rho_P = liquid_sodium::rho(T_x_bulk[i]);
                    const double rho_L = liquid_sodium::rho(T_x_bulk[i - 1]);
                    const double rho_R = liquid_sodium::rho(T_x_bulk[i + 1]);

                    const double d_l_face = 0.5 * (1.0 / bXU[i - 1] + 1.0 / bXU[i]) / dz;
                    const double d_r_face = 0.5 * (1.0 / bXU[i] + 1.0 / bXU[i + 1]) / dz;

                    // RC corrections consistent with your second block
                    const double rc_l = -d_l_face / 4.0 *
                        (p_padded_x[i - 2] - 3.0 * p_padded_x[i - 1] + 3.0 * p_padded_x[i] - p_padded_x[i + 1]);

                    const double rc_r = -d_r_face / 4.0 *
                        (p_padded_x[i - 1] - 3.0 * p_padded_x[i] + 3.0 * p_padded_x[i + 1] - p_padded_x[i + 2]);

                    const double u_l_star = 0.5 * (u_x[i - 1] + u_x[i]) + rhie_chow_on_off_x * rc_l;
                    const double u_r_star = 0.5 * (u_x[i] + u_x[i + 1]) + rhie_chow_on_off_x * rc_r;

                    const double phi_l = (u_l_star > 0.0) ? rho_L * u_l_star : rho_P * u_l_star;
                    const double phi_r = (u_r_star > 0.0) ? rho_P * u_r_star : rho_R * u_r_star;

                    const double mass_imbalance = (phi_r - phi_l);          // [kg/(m2 s)]
                    const double mass_flux = Gamma[i] * dz;    // [kg/(m2 s)] positive out of wick

                    const double rho_l = 0.5 * (rho_L + rho_P);
                    const double rho_r = 0.5 * (rho_P + rho_R);

                    const double E_l = rho_l * d_l_face; // [s/m]
                    const double E_r = rho_r * d_r_face; // [s/m]

                    aXP[i] = -E_l;
                    cXP[i] = -E_r;
                    bXP[i] = E_l + E_r;
                    dXP[i] = -mass_flux - mass_imbalance;
                }

                // BCs for p' (same as your second)
                aXP[0] = 0.0; bXP[0] = 1.0; cXP[0] = -1.0; dXP[0] = 0.0;
                aXP[N - 1] = 0.0; bXP[N - 1] = 1.0; cXP[N - 1] = 0.0; dXP[N - 1] = 0.0;

                p_prime_x = solveTridiagonal(aXP, bXP, cXP, dXP);

                // -------------------------------------------------------
                // PRESSURE CORRECTOR
                // -------------------------------------------------------
                p_error_x = 0.0;

                for (int i = 0; i < N; ++i) {
                    const double p_prev = p_x[i];
                    p_x[i] += p_prime_x[i];
                    p_storage_x[i + 1] = p_x[i];
                    p_error_x = std::max(p_error_x, std::fabs(p_x[i] - p_prev));
                }

                p_storage_x[0] = p_storage_x[1];
                p_storage_x[N + 1] = p_storage_x[N];

                // -------------------------------------------------------
                // VELOCITY CORRECTOR
                // -------------------------------------------------------
                u_error_x = 0.0;
                for (int i = 1; i < N - 1; ++i) {
                    const double u_prev = u_x[i];
                    u_x[i] -= (p_prime_x[i + 1] - p_prime_x[i - 1]) / (2.0 * dz * bXU[i]);
                    u_error_x = std::max(u_error_x, std::fabs(u_x[i] - u_prev));
                }

                inner_x++;
            }

            outer_x++;
        }

        // ===============================================================
        // TEMPERATURE SOLVER
        // ===============================================================
        std::vector<double>
            aXT(N, 0.0),
            bXT(N, 0.0),
            cXT(N, 0.0),
            dXT(N, 0.0);

        for (int i = 1; i < N - 1; i++) {

            const double rho = liquid_sodium::rho(T_x_bulk[i]);
            const double cp = liquid_sodium::cp(T_x_bulk[i]);
            const double kL = liquid_sodium::k(T_x_bulk[i - 1]);
            const double kR = liquid_sodium::k(T_x_bulk[i + 1]);

            aXT[i] = -kL / dz;
            cXT[i] = -kR / dz;
            bXT[i] = (kL + kR) / dz + rho * cp * dz / dt;
            dXT[i] = rho * cp * dz / dt * T_x_bulk_old[i]
                + Q[i] * dz;
        }

        bXT[0] = bXT[N - 1] = 1.0;
        cXT[0] = -1.0;
        aXT[N - 1] = -1.0;

        T_x_bulk = solveTridiagonal(aXT, bXT, cXT, dXT);

        for (int i = 0; i < N; ++i)
            rho_x[i] = liquid_sodium::rho(T_x_bulk[i]);

        // ===============================================================
        // OUTPUT
        // ===============================================================
        for (int i = 0; i < N; ++i) {

            v_out << u_x[i] << " ";
            p_out << p_x[i] << " ";
            T_out << T_x_bulk[i] << " ";
        }

        v_out << "\n";
        p_out << "\n";
        T_out << "\n";

        v_out.flush();
        p_out.flush();
        T_out.flush();
    }

    v_out.close();
    p_out.close();
    T_out.close();

    return 0;
}
