#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_H_

#include <AMReX_REAL.H>

/* \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumBoris(
    amrex::Real& ux, amrex::Real& uy, amrex::Real& uz,
    const amrex::Real Ex, const amrex::Real Ey, const amrex::Real Ez,
    const amrex::Real Bx, const amrex::Real By, const amrex::Real Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    const amrex::Real econst = 0.5*q*dt/m;

    // First half-push for E
    ux += econst*Ex;
    uy += econst*Ey;
    uz += econst*Ez;
    // Compute temporary gamma factor
    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);
    const amrex::Real inv_gamma = 1./std::sqrt(1. + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Magnetic rotation
    // - Compute temporary variables
    const amrex::Real tx = econst*inv_gamma*Bx;
    const amrex::Real ty = econst*inv_gamma*By;
    const amrex::Real tz = econst*inv_gamma*Bz;
    const amrex::Real tsqi = 2./(1. + tx*tx + ty*ty + tz*tz);
    const amrex::Real sx = tx*tsqi;
    const amrex::Real sy = ty*tsqi;
    const amrex::Real sz = tz*tsqi;
    const amrex::Real ux_p = ux + uy*tz - uz*ty;
    const amrex::Real uy_p = uy + uz*tx - ux*tz;
    const amrex::Real uz_p = uz + ux*ty - uy*tx;
    // - Update momentum
    ux += uy_p*sz - uz_p*sy;
    uy += uz_p*sx - ux_p*sz;
    uz += ux_p*sy - uy_p*sx;
    // Second half-push for E
    ux += econst*Ex;
    uy += econst*Ey;
    uz += econst*Ez;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_H_
