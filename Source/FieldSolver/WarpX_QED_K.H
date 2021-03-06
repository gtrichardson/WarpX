#ifndef WarpX_QED_K_h
#define WarpX_QED_K_h

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_hybrid_QED_push (int j, int k, int l, Array4<Real> const& Ex, Array4<Real>
    const& Ey, Array4<Real> const& Ez, Array4<Real> const& Bx, Array4<Real> const& By,
    Array4<Real const> const& Bz, Real dtsdx, Real dtsdy, Real dtsdy, std:string direction)
{
#if (AMREX_SPACEDIM == 3)

    // Picking out points for stencil to be used in curl function
    const amrex::Real Mpx[3] = {
        
        -2*PhysConst::xi*pow(Physconst::c^2) * (2*(pow(Ex(j+1,k,l),2)+pow(Ey(j+1,k,l),2)+pow(Ez(j+1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*Bx(j+1,k,l)
       -7*(Ex(j+1,k,l)*Bx(j+1,k,l)+Ey(j+1,k,l)*By(j+1,k,l)+Ez(j+1,k,l)*Bz(j+1,k,l))*Ex(j+1,k,l)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j+1,k,l),2)+pow(Ey(j+1,k,l),2)+pow(Ez(j+1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*By(j+1,k,l)
        -7*(Ex(j+1,k,l)*Bx(j+1,k,l)+Ey(j+1,k,l)*By(j+1,k,l)+Ez(j+1,k,l)*Bz(j+1,k,l))*Ey(j+1,k,l)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j+1,k,l),2)+pow(Ey(j+1,k,l),2)+pow(Ez(j+1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j+1,k,l),2)+pow(By(j+1,k,l),2)+pow(Bz(j+1,k,l),2)))*Bz(j+1,k,l)
        -7*(Ex(j+1,k,l)*Bx(j+1,k,l)+Ey(j+1,k,l)*By(j+1,k,l)+Ez(j+1,k,l)*Bz(j+1,k,l))*Ez(j+1,k,l))};
    
    const amrex::Real Mnx[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,l),2)+pow(Ey(j-1,k,l),2)+pow(Ez(j-1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*Bx(j-1,k,l)
        -7*(Ex(j-1,k,l)*Bx(j-1,k,l)+Ey(j-1,k,l)*By(j-1,k,l)+Ez(j-1,k,l)*Bz(j-1,k,l))*Exs(j-1,k,l)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,l),2)+pow(Ey(j-1,k,l),2)+pow(Ez(j-1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*By(j-1,k,l)
        -7*(Ex(j-1,k,l)*Bx(j-1,k,l)+Ey(j-1,k,l)*By(j-1,k,l)+Ez(j-1,k,l)*Bz(j-1,k,l))*Ey(j-1,k,l)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,l),2)+pow(Ey(j-1,k,l),2)+pow(Ez(j-1,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j-1,k,l),2)+pow(By(j-1,k,l),2)+pow(Bz(j-1,k,l),2)))*Bz(j-1,k,l)
        -7*(Ex(j-1,k,l)*Bx(j-1,k,l)+Ey(j-1,k,l)*By(j-1,k,l)+Ez(j-1,k,l)*Bz(j-1,k,l))*Ez(j-1,k,l))};
    
    const amrex::Real Mpy[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,l),2)+pow(Ey(j,k+1,l),2)+pow(Ez(j,k+1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*Bx(j,k+1,l)
        -7*(Ex(j,k+1,l)*Bx(j,k+1,l)+Ey(j,k+1,l)*By(j,k+1,l)+Ez(j,k+1,l)*Bz(j,k+1,l))*Ex(j,k+1,l)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,l),2)+pow(Ey(j,k+1,l),2)+pow(Ez(j,k+1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*By(j,k+1,l)
        -7*(Ex(j,k+1,l)*Bx(j,k+1,l)+Ey(j,k+1,l)*By(j,k+1,l)+Ez(j,k+1,l)*Bz(j,k+1,l))*Ey(j,k+1,l)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,l),2)+pow(Ey(j,k+1,l),2)+pow(Ez(j,k+1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k+1,l),2)+pow(By(j,k+1,l),2)+pow(Bz(j,k+1,l),2)))*Bz(j,k+1,l)
        -7*(Ex(j,k+1,l)*Bx(j,k+1,l)+Ey(j,k+1,l)*By(j,k+1,l)+Ez(j,k+1,l)*Bz(j,k+1,l))*Ez(j,k+1,l))};
    
    const amrex::Real Mny[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,l),2)+pow(Ey(j,k-1,l),2)+pow(Ez(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*Bx(j,k-1,l)
        -7*(Ex(j,k-1,l)*Bx(j,k-1,l)+Ey(j,k-1,l)*By(j,k-1,l)+Ez(j,k-1,l)*Bz(j,k-1,l))*Ex(j,k-1,l)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,l),2)+pow(Ey(j,k-1,l),2)+pow(Ez(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*By(j,k-1,l)
        -7*(Ex(j,k-1,l)*Bx(j,k-1,l)+Ey(j,k-1,l)*By(j,k-1,l)+Ez(j,k-1,l)*Bz(j,k-1,l))*Ey(j,k-1,l)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,l),2)+pow(Ey(j,k-1,l),2)+pow(Ez(j,k-1,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k-1,l),2)+pow(By(j,k-1,l),2)+pow(Bz(j,k-1,l),2)))*Bz(j,k-1,l)
        -7*(Ex(j,k-1,l)*Bx(j,k-1,l)+Ey(j,k-1,l)*By(j,k-1,l)+Ez(j,k-1,l)*Bz(j,k-1,l))*Ez(j,k-1,l))};
    
    const amrex::Real Mpz[3] = {
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l+1),2)+pow(Ey(j,k,l+1),2)+pow(Ez(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*Bx(j,k,l+1)
        -7*(Ex(j,k,l+1)*Bx(j,k,l+1)+Ey(j,k,l+1)*By(j,k,l+1)+Ez(j,k,l+1)*Bz(j,k,l+1))*Ex(j,k,l+1)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l+1),2)+pow(Ey(j,k,l+1),2)+pow(Ez(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*By(j,k,l+1)
        -7*(Ex(j,k,l+1)*Bx(j,k,l+1)+Ey(j,k,l+1)*By(j,k,l+1)+Ez(j,k,l+1)*Bz(j,k,l+1))*Ey(j,k,l+1)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l+1),2)+pow(Ey(j,k,l+1),2)+pow(Ez(j,k,l+1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l+1),2)+pow(By(j,k,l+1),2)+pow(Bz(j,k,l+1),2)))*Bz(j,k,l+1)
        -7*(Ex(j,k,l+1)*Bx(j,k,l+1)+Ey(j,k,l+1)*By(j,k,l+1)+Ez(j,k,l+1)*Bz(j,k,l+1))*Ez(j,k,l+1))};
    
    const amrex::Real Mnz[3] = {
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l-1),2)+pow(Ey(j,k,l-1),2)+pow(Ez(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*Bx(j,k,l-1)
        -7*(Ex(j,k,l-1)*Bx(j,k,l-1)+Ey(j,k,l-1)*By(j,k,l-1)+Ez(j,k,l-1)*Bz(j,k,l-1))*Ex(j,k,l-1)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l-1),2)+pow(Ey(j,k,l-1),2)+pow(Ez(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*By(j,k,l-1)
        -7*(Ex(j,k,l-1)*Bx(j,k,l-1)+Ey(j,k,l-1)*By(j,k,l-1)+Ez(j,k,l-1)*Bz(j,k,l-1))*Ey(j,k,l-1)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k,l-1),2)+pow(Ey(j,k,l-1),2)+pow(Ez(j,k,l-1),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l-1),2)+pow(By(j,k,l-1),2)+pow(Bz(j,k,l-1),2)))*Bz(j,k,l-1)
        -7*(Ex(j,k,l-1)*Bx(j,k,l-1)+Ey(j,k,l-1)*By(j,k,l-1)+Ez(j,k,l-1)*Bz(j,k,l-1))*Ez(j,k,l-1))};

    // Calculating necessary curls
    const amrex::Real VxM[3] = {
        0.5*( dtsdy*(Mpy[2]-Mny[2]) - dtsdz*( Mpz[1]-Mnz[1] ) ),
        0.5*( dtsdz*(Mpz[0]-Mnz[0]) - dtsdx*( Mpx[2]-Mnx[2] ) ),
        0.5*( dtsdx*(Mpx[1]-Mnx[1]) - dtsdy*( Mpy[0]-Mny[0] ) );
    
    const amrex::Real VxE[3] = {
        0.5*(dtsdy*(Ez(j,k+1,l)-Ez(j,k-1.l)) - dtsdz*(Ey(j,k,l+1)-Ey(j,k,l-1))),
        0.5*(dtsdz*(Ex(j,k,l+1)-Ex(j,k,l-1)) - dtsdx*(Ez(j+1,k,l)-Ez(j-1,k,l))),
        0.5*(dtsdx*(Ey(j+1,k,l)-Ey(j-1,k,l)) - dtsdy*(Ex(j,k+1,l)-Ex(j,k-1,l))),};
    
    const amrex::Real VxB[3] = {
        0.5*(dtsdy*(Bz(j,k+1,l)-Bz(j,k-1.l)) - dtsdz*(By(j,k,l+1)-By(j,k,l-1))),
        0.5*(dtsdz*(Bx(j,k,l+1)-Bx(j,k,l-1)) - dtsdx*(Bz(j+1,k,l)-Bz(j-1,k,l))),
        0.5*(dtsdx*(By(j+1,k,l)-By(j-1,k,l)) - dtsdy*(Bx(j,k+1,l)-Bx(j,k-1,l))),};

    // Defining comapct values for QED corrections
    const amrex::Real beta = 4*(PhysConst::xi/PhysConst::ep0)*(pow(Ex(j,k,l),2)+pow(Ey(j,k,l),2)+pow(Ez(j,k,l),2) -
        pow(PhysConst::c,2)*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))+PhysConst::ep0;
    
    const amrex::Real Alpha[3] = {
    
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,l)*VxE[0]+Ey(j,k,l)*VxE[1]+Ez(j,k,l)*VxE[2])*Bx(j,k,l)
            -7*(Ex(j,k,l)*Bx(j,k,l)+Ey(j,k,l)*By(j,k,l)+Ez(j,k,l)*Bz(j,k,l))*VxE[0]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*Ex(j,k,l)) + VxM[0],
        
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,l)*VxE[0]+Ey(j,k,l)*VxE[1]+Ez(j,k,l)*VxE[2])*By(j,k,l)
            -7*(Ex(j,k,l)*Bx(j,k,l)+Ey(j,k,l)*By(j,k,l)+Ez(j,k,l)*Bz(j,k,l))*VxE[1]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*Ey(j,k,l)) + VxM[1]
        
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,l)*VxE[0]+Ey(j,k,l)*VxE[1]+Ez(j,k,l)*VxE[2])*Bz(j,k,l)
            -7*(Ex(j,k,l)*Bx(j,k,l)+Ey(j,k,l)*By(j,k,l)+Ez(j,k,l)*Bz(j,k,l))*VxE[2]
            +4*(Bx(j,k,l)*VxE[0]+By(j,k,l)*VxE[1]+Bz(j,k,l)*VxE[2])*Ez) + VxM[2]};
        
    const amrex::Real Omega[3] = {
    
        alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(Ex(j,k,l)*VxB[0]+Ey(j,k,l)*VxB[1]+Ez(j,k,l)*VxB[2])*Ex(j,k,l)
            +2*(pow(Ex(j,k,l),2)+pow(Ey(j,k,l),2)+pow(Ez(j,k,l),2)-c^2*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[0]
            +7*c^2*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*Bx(j,k,l)),
        
        alpha[1] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(Ex(j,k,l)*VxB[0]+Ey(j,k,l)*VxB[1]+Ez(j,k,l)*VxB[2])*Ey(j,k,l)
            +2*(pow(Ex(j,k,l),2)+pow(Ey(j,k,l),2)+pow(Ez(j,k,l),2)-c^2*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[1]
            +7*c^2*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*By(j,k,l)),
        
        alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) * (4*(Ex(j,k,l)*VxB[0]+Ey(j,k,l)*VxB[1] +
            Ez(j,k,l)*VxB[2])*Ez(j,k,l)+2*(pow(Ex(j,k,l),2)+pow(Ey(j,k,l),2)+pow(Ez(j,k,l),2) - c^2*(pow(Bx(j,k,l),2)+pow(By(j,k,l),2)+pow(Bz(j,k,l),2)))*VxB[2] + 7*c^2*(VxB[0]*Bx(j,k,l)+VxB[1]*By(j,k,l)+VxB[2]*Bz(j,k,l))*Bz(j,k,l))};
    
    // Defining matrix for the algorithm
    const amrex::Real a00 = beta + (PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(Ex(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(Bx(j,k,l),2) );
    
    const amrex::Real a11 = beta +(PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(Ey(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(By(j,k,l),2) );
    
    const amrex::Real a22 = beta + (PhysConst::xi/PhysConst::ep0) *
        ( 8*pow(Ez(j,k,l),2) + 14*pow(PhysConst::c,2)*pow(Bz(j,k,l),2) );
    
    const amrex::Real a01 = (PhysConst::xi/PhysConst::ep0) *
         ( 8*Ex(j,k,l)*Ey(j,k,l) + 14*pow(PhysConst::c,2)*Bx(j,k,l)*By(j,k,l) );
    
    const amrex::Real a02 = (PhysConst::xi/PhysConst::ep0) *
        ( 8*Ex(j,k,l)*Ez(j,k,l) + 14*pow(PhysConst::c,2)*Bx(j,k,l)*Bz(j,k,l) );
    
    const amrex::Real a12 = (PhysConst::xi/PhysConst::ep0) *
        ( 8*Ez(j,k,l)*Ey(j,k,l) + 14*pow(PhysConst::c,2)*Bz(j,k,l)*By(j,k,l) );
    
    const amrex::Real detA = a00*(a11*a22-a12^2)-a01*(a01*a22-a02*a12)+a02*(a01*a12-a02*a11);
    
    // Inverting the matrix
    const amrex::Real invAx[3] = {a22*a11-a12^2, a12*a02-a22*a01, a12*a01-a11*a02};
    const amrex::Real invAy[3] = {a02*a12-a22*a01, a00*a22-a02^2, a01*a02-a12*a00};
    const amrex::Real invAz[3] = {a12*a01-a02*a11, a02*a01-a12*a00, a11*a00-a01^2};
    
    
// Performing push
#if (direction == "x")
    Ex(j,k,l) = Ex(j,k,l) - 0.5*(invAx[0]*Omega[0]+invAx[1]*Omega[1]+invAx[2]*Omega[2])/detA;

#elif (direction == "y")
    Ey(j,k,l) = Ey(j,k,l) - 0.5*(invAy[0]*Omega[0]+invAy[1]*Omega[1]+invAy[2]*Omega[2])/detA;

#else
    Ez(j,k,l) = Ez(j,k,l) - 0.5*(invAz[0]*Omega[0]+invAz[1]*Omega[1]+invAz[2]*Omega[2])/detA;

#endif



#else
        const amrex::Real Mpx[3] = {
        
        -2*PhysConst::xi*Physconst::c^2 * (2*(pow(Ex(j+1,k,0),2)+pow(Ey(j+1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*Bx(j+1,k,0)
        -7*(Ex(j+1,k,0)*Bx(j+1,k,0)+Ey(j+1,k,0)*By(j+1,k,0)+Ez(j+1,k,0)*Bz(j+1,k,0))*Ex(j+1,k,0)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j+1,k,0),2)+pow(Ey(j+1,k,0),2)+pow(Ez(j+1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*By(j+1,k,0)
        -7*(Ex(j+1,k,0)*Bx(j+1,k,0)+Ey(j+1,k,0)*By(j+1,k,0)+Ez(j+1,k,0)*Bz(j+1,k,0))*Ey(j+1,k,0)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j+1,k,0),2)+pow(Ey(j+1,k,0),2)+pow(Ez(j+1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j+1,k,0),2)+pow(By(j+1,k,0),2)+pow(Bz(j+1,k,0),2)))*Bz(j+1,k,0)
        -7*(Ex(j+1,k,0)*Bx(j+1,k,0)+Ey(j+1,k,0)*By(j+1,k,0)+Ez(j+1,k,0)*Bz(j+1,k,0))*Ez(j+1,k,0))};
    
    const amrex::Real Mnx[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,0),2)+pow(Ey(j-1,k,0),2)+pow(Ez(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*Bx(j-1,k,0)
        -7*(Ex(j-1,k,0)*Bx(j-1,k,0)+Ey(j-1,k,0)*By(j-1,k,0)+Ez(j-1,k,0)*Bz(j-1,k,0))*Exs(j-1,k,0)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,0),2)+pow(Ey(j-1,k,0),2)+pow(Ez(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*By(j-1,k,0)
        -7*(Ex(j-1,k,0)*Bx(j-1,k,0)+Ey(j-1,k,0)*By(j-1,k,0)+Ez(j-1,k,0)*Bz(j-1,k,0))*Ey(j-1,k,0)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j-1,k,0),2)+pow(Ey(j-1,k,0),2)+pow(Ez(j-1,k,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j-1,k,0),2)+pow(By(j-1,k,0),2)+pow(Bz(j-1,k,0),2)))*Bz(j-1,k,0)
        -7*(Ex(j-1,k,0)*Bx(j-1,k,0)+Ey(j-1,k,0)*By(j-1,k,0)+Ez(j-1,k,0)*Bz(j-1,k,0))*Ez(j-1,k,0))};
    
     const amrex::Real Mpz[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,0),2)+pow(Ey(j,k+1,0),2)+pow(Ez(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*Bx(j,k+1,0)
        -7*(Ex(j,k+1,0)*Bx(j,k+1,0)+Ey(j,k+1,0)*By(j,k+1,0)+Ez(j,k+1,0)*Bz(j,k+1,0))*Ex(j,k+1,0)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,0),2)+pow(Ey(j,k+1,0),2)+pow(Ez(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*By(j,k+1,0)
        -7*(Ex(j,k+1,0)*Bx(j,k+1,0)+Ey(j,k+1,0)*By(j,k+1,0)+Ez(j,k+1,0)*Bz(j,k+1,0))*Ey(j,k+1,0)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k+1,0),2)+pow(Ey(j,k+1,0),2)+pow(Ez(j,k+1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k+1,0),2)+pow(By(j,k+1,0),2)+pow(Bz(j,k+1,0),2)))*Bz(j,k+1,0)
        -7*(Ex(j,k+1,0)*Bx(j,k+1,0)+Ey(j,k+1,0)*By(j,k+1,0)+Ez(j,k+1,0)*Bz(j,k+1,0))*Ez(j,k+1,0))};
    
    const amrex::Real Mnz[3] = {
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,0),2)+pow(Ey(j,k-1,0),2)+pow(Ez(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*Bx(j,k-1,0)
        -7*(Ex(j,k-1,0)*Bx(j,k-1,0)+Ey(j,k-1,0)*By(j,k-1,0)+Ez(j,k-1,0)*Bz(j,k-1,0))*Ex(j,k-1,0)),
        
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,0),2)+pow(Ey(j,k-1,0),2)+pow(Ez(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*By(j,k-1,0)
        -7*(Ex(j,k-1,0)*Bx(j,k-1,0)+Ey(j,k-1,0)*By(j,k-1,0)+Ez(j,k-1,0)*Bz(j,k-1,0))*Ey(j,k-1,0)),
    
        -2*PhysConst::xi*pow(PhysConst::c,2) * (2*(pow(Ex(j,k-1,0),2)+pow(Ey(j,k-1,0),2)+pow(Ez(j,k-1,0),2)
        -pow(PhysConst::c,2)*(pow(Bx(j,k-1,0),2)+pow(By(j,k-1,0),2)+pow(Bz(j,k-1,0),2)))*Bz(j,k-1,0)
        -7*(Ex(j,k-1,0)*Bx(j,k-1,0)+Ey(j,k-1,0)*By(j,k-1,0)+Ez(j,k-1,0)*Bz(j,k-1,0))*Ez(j,k-1,0))};
    
    const amrex::Real VxM[3] = {
        -0.5*(dtsdz*( Mpz[1]-Mnz[1]) ),
        0.5*( dtsdz*(Mpz[0]-Mnz[0]) - dtsdx*( Mpx[2]-Mnx[2] ) ),
        0.5*( dtsdx*(Mpx[1]-Mnx[1])) );
    
    const amrex::Real VxE[3] = {
        -0.5*(dtsdz*(Ey(j,k+1,0)-Ey(j,k-1,0))),
        0.5*(dtsdz*(Ex(j,k+1,0)-Ex(j,k-1,0)) - dtsdx*(Ez(j+1,k,0)-Ez(j-1,k,0))),
        0.5*(dtsdx*(Ey(j+1,k,0)-Ey(j-1,k,0)))};
    
    const amrex::Real VxB[3] = {
        -0.5*dtsdz*(By(j,k+1,0)-By(j,k-1,0)),
        0.5*(dtsdz*(Bx(j,k+1,0)-Bx(j,k-1,0)) - dtsdx*(Bz(j+1,k,0)-Bz(j-1,k,0))),
        0.5*(dtsdx*(By(j+1,k,0)-By(j-1,k,0)))};

    
    const amrex::Real beta =
        4*(PhysConst::xi/PhysConst::ep0)*(pow(Ex(j,k,0),2)+pow(Ey(j,k,0),2)+pow(Ez(j,k,0),2)-
        pow(PhysConst::c,2)*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))+PhysConst::ep0;
        
    const amrex::Real Alpha[3] = {
    
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,0)*VxE[0]+Ey(j,k,0)*VxE[1]+Ez(j,k,0)*VxE[2])*Bx(j,k,0)
            - 7*(Ex(j,k,0)*Bx(j,k,0)+Ey(j,k,0)*By(j,k,0)+Ez(j,k,0)*Bz(j,k,0))*VxE[0]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*Ex(j,k,0)) + VxM[0],
        
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,0)*VxE[0]+Ey(j,k,0)*VxE[1]+Ez(j,k,0)*VxE[2])*By(j,k,0)
            - 7*(Ex(j,k,0)*Bx(j,k,0)+Ey(j,k,0)*By(j,k,0)+Ez(j,k,0)*Bz(j,k,0))*VxE[1]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*Ey(j,k,0)) + VxM[1]
        
        2*PhysConst::xi*pow(PhysConst::c,2) *
            ((-1)*7*(Ex(j,k,0)*VxE[0]+Ey(j,k,0)*VxE[1]+Ez(j,k,0)*VxE[2])*Bz(j,k,0)
            - 7*(Ex(j,k,0)*Bx(j,k,0)+Ey(j,k,0)*By(j,k,0)+Ez(j,k,0)*Bz(j,k,0))*VxE[2]
            + 4*(Bx(j,k,0)*VxE[0]+By(j,k,0)*VxE[1]+Bz(j,k,0)*VxE[2])*Ez(j,k,0)) + VxM[2]};
        
    const amrex::Real Omega[3] = {
    
        alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(Ex(j,k,0)*VxB[0]+Ey(j,k,0)*VxB[1]+Ez(j,k,0)*VxB[2])*Ex(j,k,0)
            + 2*(pow(Ex(j,k,0),2)+pow(Ey(j,k,0),2)+pow(Ez(j,k,0),2)-c^2*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[0]
            + 7*c^2*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*Bx(j,k,0)),
        
        alpha[1] + 2*PhysConst::xi*pow(PhysConst::c,2) *
            (4*(Ex(j,k,0)*VxB[0]+Ey(j,k,0)*VxB[1]+Ez(j,k,0)*VxB[2])*Ey(j,k,0)
            +2*(pow(Ex(j,k,0),2)+pow(Ey(j,k,0),2)+pow(Ez(j,k,0),2)-c^2*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[1]
            +7*c^2*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*By(j,k,0)),
        
        alpha[0] + 2*PhysConst::xi*pow(PhysConst::c,2)*
            (4*(Ex(j,k,0)*VxB[0]+Ey(j,k,0)*VxB[1]+Ez(j,k,0)*VxB[2])*Ez(j,k,0)
            +2*(pow(Ex(j,k,0),2)+pow(Ey(j,k,0),2)+pow(Ez(j,k,0),2)-c^2*(pow(Bx(j,k,0),2)+pow(By(j,k,0),2)+pow(Bz(j,k,0),2)))*VxB[2]
            +7*c^2*(VxB[0]*Bx(j,k,0)+VxB[1]*By(j,k,0)+VxB[2]*Bz(j,k,0))*Bz(j,k,0))};
    
    const amrex::Real a00 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(Ex(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(Bx(j,k,0),2) );
    const amrex::Real a11 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(Ey(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(By(j,k,0),2) );
    const amrex::Real a22 = beta + (PhysConst::xi/PhysConst::ep0) * ( 8*pow(Ez(j,k,0),2) + 14*pow(PhysConst::c,2)*pow(Bz(j,k,0),2) );
    
    const amrex::Real a01 = (PhysConst::xi/PhysConst::ep0) * ( 8*Ex(j,k,0)*Ey(j,k,0) + 14*pow(PhysConst::c,2)*Bx(j,k,0)*By(j,k,0) );
    const amrex::Real a02 = (PhysConst::xi/PhysConst::ep0) * ( 8*Ex(j,k,0)*Ez(j,k,0) + 14*pow(PhysConst::c,2)*Bx(j,k,0)*Bz(j,k,0) );
    const amrex::Real a12 = (PhysConst::xi/PhysConst::ep0) * ( 8*Ez(j,k,0)*Ey(j,k,0) + 14*pow(PhysConst::c,2)*Bz(j,k,0)*By(j,k,0) );
    
    const amrex::Real detA = a00*(a11*a22-a12^2)-a01*(a01*a22-a02*a12)+a02*(a01*a12-a02*a11);
    
    const amrex::Real invAx[3] = {a22*a11-a12^2, a12*a02-a22*a01, a12*a01-a11*a02};
    const amrex::Real invAy[3] = {a02*a12-a22*a01, a00*a22-a02^2, a01*a02-a12*a00};
    const amrex::Real invAz[3] = {a12*a01-a02*a11, a02*a01-a12*a00, a11*a00-a01^2};
    
#if (direction == "x")
    Ex(j,k,0) = Ex(j,k,0) - 0.5*(invAx[0]*Omega[0]+invAx[1]*Omega[1]+invAx[2]*Omega[2])/detA;

#elif (direction == "y")
    Ey(j,k,0) = Ey(j,k,0) - 0.5*(invAy[0]*Omega[0]+invAy[1]*Omega[1]+invAy[2]*Omega[2])/detA;

#else
    Ez(j,k,0) = Ez(j,k,0) - 0.5*(invAz[0]*Omega[0]+invAz[1]*Omega[1]+invAz[2]*Omega[2])/detA;

#endif

#endif
