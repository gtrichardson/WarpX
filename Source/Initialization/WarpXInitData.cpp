/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 * Axel Huebl, Burlen Loring, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <WarpX.H>
#include <BilinearFilter.H>
#include <NCIGodfreyFilter.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <GpuParser.H>
#include <WarpXUtil.H>

#include <cmath>

namespace QED
{
	// Setting the order of the paraxial approximation
	int order = 5;
	
	// Defining beam paramters
	amrex::Real E0 = 1.e12;
	amrex::Real w0 = 1.;
	amrex::Real psi0 = 0.;
	amrex::Real k = 6283185307.18;
	amrex::Real c = 299792458.;
		
    // Defining Gaussian constants
    amrex::Real Zr =0.5*k*w0*w0;
    amrex::Real eps = w0/Zr;
    /*
    std::string Antonins_Ex(double x, double y, double z){
    
		amrex::Real Ex = 0; // Intializing field
		
		// New coordinate system
		amrex::Real xi = z/w0;
		amrex::Real nu = -y/w0;
			
		// misclanious
		amrex::Real R = x+Zr*Zr/x;
		amrex::Real w = w0*std::sqrt(1+x*x/(Zr*Zr));
		amrex::Real r2 = z*z+y*y;
		amrex::Real rho2 = r2/(w0*w0);
		amrex::Real rho4 = rho2*rho2;
		amrex::Real rho6 = rho2*rho2*rho2;
		amrex::Real rho8 = rho2*rho2*rho2*rho2;
		amrex::Real psip = -k*x;
		amrex::Real psir = k*r2/(2*R);
		amrex::Real psig = std::atan(2*x/(k*w0*w0));
		amrex::Real psi = psi0 + psip - psir + psig;
		amrex::Real E = E0*w0/w*std::exp(-1*r2/(w*w));
		amrex::Real B = E/c;
			
		// Sine and Cosine terms
		amrex::Real S0 = std::sin(psi);
		amrex::Real S1 = (w0/w)*std::sin(psi+psig);
		amrex::Real S2 = std::pow(w0,2)/std::pow(w,2)*std::sin(psi + 2*psig);
		amrex::Real S3 = std::pow(w0,3)/std::pow(w,3)*std::sin(psi + 3*psig);
		amrex::Real S4 = std::pow(w0,4)/std::pow(w,4)*std::sin(psi + 4*psig);
		amrex::Real S5 = std::pow(w0,5)/std::pow(w,5)*std::sin(psi + 5*psig);
		amrex::Real S6 = std::pow(w0,6)/std::pow(w,6)*std::sin(psi + 6*psig);
		  
		amrex::Real C0 = std::sin(psi);
		amrex::Real C1 = (w0/w)*std::cos(psi+psig);
		amrex::Real C2 = std::pow(w0,2)/std::pow(w,2)*std::cos(psi + 2*psig);
		amrex::Real C3 = std::pow(w0,3)/std::pow(w,3)*std::cos(psi + 3*psig);
		amrex::Real C4 = std::pow(w0,4)/std::pow(w,4)*std::cos(psi + 4*psig);
		amrex::Real C5 = std::pow(w0,5)/std::pow(w,5)*std::cos(psi + 5*psig);
		amrex::Real C6 = std::pow(w0,6)/std::pow(w,6)*std::cos(psi + 6*psig);
		amrex::Real C7 = std::pow(w0,7)/std::pow(w,7)*std::cos(psi + 7*psig);
        
	    // Defining fields by gaussian orders
	    amrex::Real Ex_ord0 = 0;
		amrex::Real Ex_ord1 = C1;
		amrex::Real Ex_ord2 = 0;
		amrex::Real Ex_ord3 = -0.5*C2+rho2*C3-0.25*rho2*rho2*C4;
		amrex::Real Ex_ord4 = 0;
		amrex::Real Ex_ord5 = -3./8.*C3-3./8.*rho2*C4+17./16.*rho4*C5-3./8.*rho6*C6+rho8*C7/32.;

		if(order >= 1){Ex = Ex + std::pow(eps,1)*Ex_ord1;}
		if(order >= 3){Ex = Ex + std::pow(eps,3)*Ex_ord3;}
		if(order >= 5){Ex = Ex + std::pow(eps,5)*Ex_ord5;}
		
		Ex = E*xi*Ex;
		return Ex;
		}
	
	amrex::Real Antonins_Ey(double x, double y, double z){
		
		// New coordinate system
		amrex::Real xi = z/w0;
		amrex::Real nu = -y/w0;
			
		// misclanious
		amrex::Real R = x+Zr*Zr/x;
		amrex::Real w = w0*std::sqrt(1+x*x/(Zr*Zr));
		amrex::Real r2 = z*z+y*y;
		amrex::Real rho2 = r2/(w0*w0);
		amrex::Real rho4 = rho2*rho2;
		amrex::Real rho6 = rho2*rho2*rho2;
		amrex::Real rho8 = rho2*rho2*rho2*rho2;
		amrex::Real psip = -k*x;
		amrex::Real psir = k*r2/(2*R);
		amrex::Real psig = std::atan(2*x/(k*w0*w0));
		amrex::Real psi = psi0 + psip - psir + psig;
		amrex::Real E = E0*w0/w*std::exp(-1*r2/(w*w));
		amrex::Real B = E/c;
			
		// Sine and Cosine terms
		amrex::Real S0 = std::sin(psi);
		amrex::Real S1 = (w0/w)*std::sin(psi+psig);
		amrex::Real S2 = std::pow(w0,2)/std::pow(w,2)*std::sin(psi + 2*psig);
		amrex::Real S3 = std::pow(w0,3)/std::pow(w,3)*std::sin(psi + 3*psig);
		amrex::Real S4 = std::pow(w0,4)/std::pow(w,4)*std::sin(psi + 4*psig);
		amrex::Real S5 = std::pow(w0,5)/std::pow(w,5)*std::sin(psi + 5*psig);
		amrex::Real S6 = std::pow(w0,6)/std::pow(w,6)*std::sin(psi + 6*psig);
		  
		amrex::Real C0 = std::sin(psi);
		amrex::Real C1 = (w0/w)*std::cos(psi+psig);
		amrex::Real C2 = std::pow(w0,2)/std::pow(w,2)*std::cos(psi + 2*psig);
		amrex::Real C3 = std::pow(w0,3)/std::pow(w,3)*std::cos(psi + 3*psig);
		amrex::Real C4 = std::pow(w0,4)/std::pow(w,4)*std::cos(psi + 4*psig);
		amrex::Real C5 = std::pow(w0,5)/std::pow(w,5)*std::cos(psi + 5*psig);
		amrex::Real C6 = std::pow(w0,6)/std::pow(w,6)*std::cos(psi + 6*psig);
		amrex::Real C7 = std::pow(w0,7)/std::pow(w,7)*std::cos(psi + 7*psig);
		
		// Intalizing field
		amrex::Real Ey = 0;
		
		// Defining fields by gaussian orders
		amrex::Real Ey_ord0 = 0;
		amrex::Real Ey_ord1 = 0;
		amrex::Real Ey_ord2 = S2;
		amrex::Real Ey_ord3 = 0;
		amrex::Real Ey_ord4 = rho2*S4-0.25*rho2*rho2*S5;
		amrex::Real Ey_ord5 = 0;
		
		if(order >= 2){Ey = Ey + std::pow(eps,2)*Ey_ord2;}
		if(order >= 4){Ey = Ey + std::pow(eps,4)*Ey_ord4;}

		Ey = E*xi*nu*Ey;
		return Ey;
	}
	
	amrex::Real Antonins_Ez(double x, double y, double z){
		// New coordinate system
		amrex::Real xi = z/w0;
		amrex::Real nu = -y/w0;
			
		// misclanious
		amrex::Real R = x+Zr*Zr/x;
		amrex::Real w = w0*std::sqrt(1+x*x/(Zr*Zr));
		amrex::Real r2 = z*z+y*y;
		amrex::Real rho2 = r2/(w0*w0);
		amrex::Real rho4 = rho2*rho2;
		amrex::Real rho6 = rho2*rho2*rho2;
		amrex::Real rho8 = rho2*rho2*rho2*rho2;
		amrex::Real psip = -k*x;
		amrex::Real psir = k*r2/(2*R);
		amrex::Real psig = std::atan(2*x/(k*w0*w0));
		amrex::Real psi = psi0 + psip - psir + psig;
		amrex::Real E = E0*w0/w*std::exp(-1*r2/(w*w));
		amrex::Real B = E/c;
			
		// Sine and Cosine terms
		amrex::Real S0 = std::sin(psi);
		amrex::Real S1 = (w0/w)*std::sin(psi+psig);
		amrex::Real S2 = std::pow(w0,2)/std::pow(w,2)*std::sin(psi + 2*psig);
		amrex::Real S3 = std::pow(w0,3)/std::pow(w,3)*std::sin(psi + 3*psig);
		amrex::Real S4 = std::pow(w0,4)/std::pow(w,4)*std::sin(psi + 4*psig);
		amrex::Real S5 = std::pow(w0,5)/std::pow(w,5)*std::sin(psi + 5*psig);
		amrex::Real S6 = std::pow(w0,6)/std::pow(w,6)*std::sin(psi + 6*psig);
		  
		amrex::Real C0 = std::sin(psi);
		amrex::Real C1 = (w0/w)*std::cos(psi+psig);
		amrex::Real C2 = std::pow(w0,2)/std::pow(w,2)*std::cos(psi + 2*psig);
		amrex::Real C3 = std::pow(w0,3)/std::pow(w,3)*std::cos(psi + 3*psig);
		amrex::Real C4 = std::pow(w0,4)/std::pow(w,4)*std::cos(psi + 4*psig);
		amrex::Real C5 = std::pow(w0,5)/std::pow(w,5)*std::cos(psi + 5*psig);
		amrex::Real C6 = std::pow(w0,6)/std::pow(w,6)*std::cos(psi + 6*psig);
		amrex::Real C7 = std::pow(w0,7)/std::pow(w,7)*std::cos(psi + 7*psig);
		
		// Intalizing field
		amrex::Real Ez = 0;
		
		// Defining fields by gaussian orders
		amrex::Real Ez_ord0 = S0;
		amrex::Real Ez_ord1 = 0;
		amrex::Real Ez_ord2 = xi*xi*S2-0.25*rho4*S3;
		amrex::Real Ez_ord3 = 0;
		amrex::Real Ez_ord4 = .125*S2-0.25*rho2*S3-0.0625*rho2*(rho2-16*xi*xi)*S4-0.125*rho4*(rho2+2*xi*xi)*S5+0.03125*rho8*S6;
		amrex::Real Ez_ord5 = 0;

		
		if(order >= 2){Ez = Ez + std::pow(eps,2)*Ez_ord2;}
		if(order >= 4){Ez = Ez + std::pow(eps,4)*Ez_ord4;}
		
		Ez = E*Ez;
		return Ez;
	}
	
	amrex::Real Antonins_Bx(double x, double y, double z){
		// New coordinate system
		amrex::Real xi = z/w0;
		amrex::Real nu = -y/w0;
			
		// misclanious
		amrex::Real R = x+Zr*Zr/x;
		amrex::Real w = w0*std::sqrt(1+x*x/(Zr*Zr));
		amrex::Real r2 = z*z+y*y;
		amrex::Real rho2 = r2/(w0*w0);
		amrex::Real rho4 = rho2*rho2;
		amrex::Real rho6 = rho2*rho2*rho2;
		amrex::Real rho8 = rho2*rho2*rho2*rho2;
		amrex::Real psip = -k*x;
		amrex::Real psir = k*r2/(2*R);
		amrex::Real psig = std::atan(2*x/(k*w0*w0));
		amrex::Real psi = psi0 + psip - psir + psig;
		amrex::Real E = E0*w0/w*std::exp(-1*r2/(w*w));
		amrex::Real B = E/c;
			
		// Sine and Cosine terms
		amrex::Real S0 = std::sin(psi);
		amrex::Real S1 = (w0/w)*std::sin(psi+psig);
		amrex::Real S2 = std::pow(w0,2)/std::pow(w,2)*std::sin(psi + 2*psig);
		amrex::Real S3 = std::pow(w0,3)/std::pow(w,3)*std::sin(psi + 3*psig);
		amrex::Real S4 = std::pow(w0,4)/std::pow(w,4)*std::sin(psi + 4*psig);
		amrex::Real S5 = std::pow(w0,5)/std::pow(w,5)*std::sin(psi + 5*psig);
		amrex::Real S6 = std::pow(w0,6)/std::pow(w,6)*std::sin(psi + 6*psig);
		  
		amrex::Real C0 = std::sin(psi);
		amrex::Real C1 = (w0/w)*std::cos(psi+psig);
		amrex::Real C2 = std::pow(w0,2)/std::pow(w,2)*std::cos(psi + 2*psig);
		amrex::Real C3 = std::pow(w0,3)/std::pow(w,3)*std::cos(psi + 3*psig);
		amrex::Real C4 = std::pow(w0,4)/std::pow(w,4)*std::cos(psi + 4*psig);
		amrex::Real C5 = std::pow(w0,5)/std::pow(w,5)*std::cos(psi + 5*psig);
		amrex::Real C6 = std::pow(w0,6)/std::pow(w,6)*std::cos(psi + 6*psig);
		amrex::Real C7 = std::pow(w0,7)/std::pow(w,7)*std::cos(psi + 7*psig);
		
		// Intalizing field
		amrex::Real Bx = 0;
		
		// Defining fields by gaussian orders
		amrex::Real Bx_ord0 = 0;
		amrex::Real Bx_ord1 = C1;
		amrex::Real Bx_ord2 = 0;
		amrex::Real Bx_ord3 = 0.5*C2+0.5*rho2*C3-0.25*rho2*rho2*C4;
		amrex::Real Bx_ord4 = 0;
		amrex::Real Bx_ord5 = 3./8.*C3+3./8.*rho2*C4+3./16.*rho4*C5-rho6*C6/4.+rho8*C7/32.;

		if(order >= 1){Bx = Bx + std::pow(eps,1)*Bx_ord1;}
		if(order >= 3){Bx = Bx + std::pow(eps,3)*Bx_ord3;}
		if(order >= 5){Bx = Bx + std::pow(eps,5)*Bx_ord5;}
		
		Bx = nu*B*Bx;
		return Bx;
	}
	
	amrex::Real Antonins_By(double x, double y, double z){
		// New coordinate system
		amrex::Real xi = z/w0;
		amrex::Real nu = -y/w0;
			
		// misclanious
		amrex::Real R = x+Zr*Zr/x;
		amrex::Real w = w0*std::sqrt(1+x*x/(Zr*Zr));
		amrex::Real r2 = z*z+y*y;
		amrex::Real rho2 = r2/(w0*w0);
		amrex::Real rho4 = rho2*rho2;
		amrex::Real rho6 = rho2*rho2*rho2;
		amrex::Real rho8 = rho2*rho2*rho2*rho2;
		amrex::Real psip = -k*x;
		amrex::Real psir = k*r2/(2*R);
		amrex::Real psig = std::atan(2*x/(k*w0*w0));
		amrex::Real psi = psi0 + psip - psir + psig;
		amrex::Real E = E0*w0/w*std::exp(-1*r2/(w*w));
		amrex::Real B = E/c;
			
		// Sine and Cosine terms
		amrex::Real S0 = std::sin(psi);
		amrex::Real S1 = (w0/w)*std::sin(psi+psig);
		amrex::Real S2 = std::pow(w0,2)/std::pow(w,2)*std::sin(psi + 2*psig);
		amrex::Real S3 = std::pow(w0,3)/std::pow(w,3)*std::sin(psi + 3*psig);
		amrex::Real S4 = std::pow(w0,4)/std::pow(w,4)*std::sin(psi + 4*psig);
		amrex::Real S5 = std::pow(w0,5)/std::pow(w,5)*std::sin(psi + 5*psig);
		amrex::Real S6 = std::pow(w0,6)/std::pow(w,6)*std::sin(psi + 6*psig);
		  
		amrex::Real C0 = std::sin(psi);
		amrex::Real C1 = (w0/w)*std::cos(psi+psig);
		amrex::Real C2 = std::pow(w0,2)/std::pow(w,2)*std::cos(psi + 2*psig);
		amrex::Real C3 = std::pow(w0,3)/std::pow(w,3)*std::cos(psi + 3*psig);
		amrex::Real C4 = std::pow(w0,4)/std::pow(w,4)*std::cos(psi + 4*psig);
		amrex::Real C5 = std::pow(w0,5)/std::pow(w,5)*std::cos(psi + 5*psig);
		amrex::Real C6 = std::pow(w0,6)/std::pow(w,6)*std::cos(psi + 6*psig);
		amrex::Real C7 = std::pow(w0,7)/std::pow(w,7)*std::cos(psi + 7*psig);
		
		// Intalizing field
		amrex::Real By = 0;
		
		// Defining fields by gaussian orders
		amrex::Real By_ord0 = S0;
		amrex::Real By_ord1 = 0;
		amrex::Real By_ord2 = 0.5*rho2*S2-0.25*rho2*rho2*S3;
		amrex::Real By_ord3 = 0;
		amrex::Real By_ord4 = -0.125*S2+0.25*rho2*S3+5*0.0625*rho4*S4-0.25*std::pow(rho2,3)*S5+0.03125*rho8*S6;
		amrex::Real By_ord5 = 0;

		if(order >= 0){By = By + std::pow(eps,0)*By_ord0;}
		if(order >= 2){By = By + std::pow(eps,2)*By_ord2;}
		if(order >= 4){By = By + std::pow(eps,4)*By_ord4;}
		
		By = B*By;
		return By;
	}
	
	amrex::Real Antonins_Bz(double x, double y, double z){
		// Intalizing field
		amrex::Real Bz = 0;
		return Bz;
	}
	*/
	std::string Ex = "1";
	std::string Ey = "1";
	std::string Ez = "1";
	std::string Bx = "1";
	std::string By = "1";
	std::string Bz = "1";
}

using namespace amrex;

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    if (restart_chkfile.empty())
    {
        ComputeDt();
        InitFromScratch();
    }
    else
    {
        InitFromCheckpoint();
        if (is_synchronized) {
            ComputeDt();
        }
        PostRestart();
    }

    ComputePMLFactors();

    if (WarpX::use_fdtd_nci_corr) {
        WarpX::InitNCICorrector();
    }

    if (WarpX::use_filter) {
        WarpX::InitFilter();
    }

    BuildBufferMasks();

    InitDiagnostics();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nGrids Summary:\n";
        printGridSummary(std::cout, 0, finestLevel());
    }

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge = new amrex::AmrMeshInSituBridge;
    insitu_bridge->setEnabled(insitu_int > 0 ? 1 : 0);
    insitu_bridge->setConfig(insitu_config);
    insitu_bridge->setPinMesh(insitu_pin_mesh);
    if (insitu_bridge->initialize())
    {
        amrex::ErrorStream()
            << "WarpX::InitData : Failed to initialize the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
    insitu_bridge->setFrequency(1);
#endif

    if (restart_chkfile.empty())
    {
        if (plot_int > 0)
            WritePlotFile();

        if (openpmd_int > 0)
            WriteOpenPMDFile();

        if (check_int > 0)
            WriteCheckPointFile();

        if ((insitu_int > 0) && (insitu_start == 0))
            UpdateInSitu();

        // Write reduced diagnostics before the first iteration.
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->ComputeDiags(-1);
            reduced_diags->WriteToFile(-1);
        }
    }
}

void
WarpX::InitDiagnostics () {
    if (do_back_transformed_diagnostics) {
        const Real* current_lo = geom[0].ProbLo();
        const Real* current_hi = geom[0].ProbHi();
        Real dt_boost = dt[0];
        // Find the positions of the lab-frame box that corresponds to the boosted-frame box at t=0
        Real zmin_lab = current_lo[moving_window_dir]/( (1.+beta_boost)*gamma_boost );
        Real zmax_lab = current_hi[moving_window_dir]/( (1.+beta_boost)*gamma_boost );
        myBFD.reset(new BackTransformedDiagnostic(zmin_lab,
                                               zmax_lab,
                                               moving_window_v, dt_snapshots_lab,
                                               num_snapshots_lab,
                                               dt_slice_snapshots_lab,
                                               num_slice_snapshots_lab,
                                               gamma_boost, t_new[0], dt_boost,
                                               moving_window_dir, geom[0],
                                               slice_realbox,
                                               particle_slice_width_lab));
    }
}

void
WarpX::InitFromScratch ()
{
    const Real time = 0.0;

    AmrCore::InitFromScratch(time);  // This will call MakeNewLevelFromScratch

    mypc->AllocData();
    mypc->InitData();

    // Loop through species and calculate their space-charge field
    bool const reset_fields = false; // Do not erase previous user-specified values on the grid
    ComputeSpaceChargeField(reset_fields);

    InitPML();

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        getLevelMasks(masks);

        // the plus one is to convert from num_cells to num_nodes
        getLevelMasks(gather_masks, n_buffer + 1);
    }
#endif // WARPX_DO_ELECTROSTATIC
}

void
WarpX::InitPML ()
{
    if (do_pml)
    {
        amrex::IntVect do_pml_Lo_corrected = do_pml_Lo;

#ifdef WARPX_DIM_RZ
        do_pml_Lo_corrected[0] = 0; // no PML at r=0, in cylindrical geometry
#endif
        pml[0].reset(new PML(boxArray(0), DistributionMap(0), &Geom(0), nullptr,
                             pml_ncell, pml_delta, 0,
#ifdef WARPX_USE_PSATD
                             dt[0], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                             do_dive_cleaning, do_moving_window,
                             pml_has_particles, do_pml_in_domain,
                             do_pml_Lo_corrected, do_pml_Hi));
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            amrex::IntVect do_pml_Lo_MR = amrex::IntVect::TheUnitVector();
#ifdef WARPX_DIM_RZ
            //In cylindrical geometry, if the edge of the patch is at r=0, do not add PML
            if ((max_level > 0) && (fine_tag_lo[0]==0.)) {
                do_pml_Lo_MR[0] = 0;
            }
#endif
            pml[lev].reset(new PML(boxArray(lev), DistributionMap(lev),
                                   &Geom(lev), &Geom(lev-1),
                                   pml_ncell, pml_delta, refRatio(lev-1)[0],
#ifdef WARPX_USE_PSATD
                                   dt[lev], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                                   do_dive_cleaning, do_moving_window,
                                   pml_has_particles, do_pml_in_domain,
                                   do_pml_Lo_MR, amrex::IntVect::TheUnitVector()));
        }
    }
}

void
WarpX::ComputePMLFactors ()
{
    if (do_pml)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            pml[lev]->ComputePMLFactors(dt[lev]);
        }
    }
}

void
WarpX::InitNCICorrector ()
{
    if (WarpX::use_fdtd_nci_corr)
    {
        for (int lev = 0; lev <= max_level; ++lev)
        {
            const Geometry& gm = Geom(lev);
            const Real* dx = gm.CellSize();
            amrex::Real dz, cdtodz;
            if (AMREX_SPACEDIM == 3){
                dz = dx[2];
            }else{
                dz = dx[1];
            }
            cdtodz = PhysConst::c * dt[lev] / dz;

            // Initialize Godfrey filters
            // Same filter for fields Ex, Ey and Bz
            const bool nodal_gather = (l_lower_order_in_v == 0);
            nci_godfrey_filter_exeybz[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Ex_Ey_Bz, cdtodz, nodal_gather) );
            // Same filter for fields Bx, By and Ez
            nci_godfrey_filter_bxbyez[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Bx_By_Ez, cdtodz, nodal_gather) );
            // Compute Godfrey filters stencils
            nci_godfrey_filter_exeybz[lev]->ComputeStencils();
            nci_godfrey_filter_bxbyez[lev]->ComputeStencils();
        }
    }
}

void
WarpX::InitFilter (){
    if (WarpX::use_filter){
        WarpX::bilinear_filter.npass_each_dir = WarpX::filter_npass_each_dir;
        WarpX::bilinear_filter.ComputeStencils();
    }
}

void
WarpX::PostRestart ()
{
#ifdef WARPX_USE_PSATD
    amrex::Abort("WarpX::PostRestart: TODO for PSATD");
#endif
    mypc->PostRestart();
}


void
WarpX::InitLevelData (int lev, Real time)
{

    ParmParse pp("warpx");

    // default values of E_external_grid and B_external_grid
    // are used to set the E and B field when "constant" or
    // "parser" is not explicitly used in the input.
    pp.query("B_ext_grid_init_style", B_ext_grid_s);
    std::transform(B_ext_grid_s.begin(),
                   B_ext_grid_s.end(),
                   B_ext_grid_s.begin(),
                   ::tolower);

    pp.query("E_ext_grid_init_style", E_ext_grid_s);
    std::transform(E_ext_grid_s.begin(),
                   E_ext_grid_s.end(),
                   E_ext_grid_s.begin(),
                   ::tolower);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (B_ext_grid_s == "constant")
        pp.getarr("B_external_grid", B_external_grid);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (E_ext_grid_s == "constant")
        pp.getarr("E_external_grid", E_external_grid);

    for (int i = 0; i < 3; ++i) {
        current_fp[lev][i]->setVal(0.0);
        if (lev > 0)
           current_cp[lev][i]->setVal(0.0);

        if (B_ext_grid_s == "constant" || B_ext_grid_s == "default") {
           Bfield_fp[lev][i]->setVal(B_external_grid[i]);
           if (lev > 0) {
              Bfield_aux[lev][i]->setVal(B_external_grid[i]);
              Bfield_cp[lev][i]->setVal(B_external_grid[i]);
           }
        }
        if (E_ext_grid_s == "constant" || E_ext_grid_s == "default") {
           Efield_fp[lev][i]->setVal(E_external_grid[i]);
           if (lev > 0) {
              Efield_aux[lev][i]->setVal(E_external_grid[i]);
              Efield_cp[lev][i]->setVal(E_external_grid[i]);
           }
        }
    }

    // if the input string for the B-field is "parse_b_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (B_ext_grid_s == "parse_b_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort("E and B parser for external fields does not work with RZ -- TO DO");
#endif
       std::string Glenn_Bx;
       std::string Glenn_By;
       std::string Glenn_Bz;
       
       Store_parserString(pp, QED::Bx, Glenn_Bx);
       Store_parserString(pp, QED::By, Glenn_By);
       Store_parserString(pp, QED::Bz, Glenn_Bz);

       Bxfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_Bx,{"x","y","z"})));
       Byfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_By,{"x","y","z"})));
       Bzfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_Bz,{"x","y","z"})));

       // Initialize Bfield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Bfield_fp[lev][0].get(),
                                                 Bfield_fp[lev][1].get(),
                                                 Bfield_fp[lev][2].get(),
                                                 Bxfield_parser.get(),
                                                 Byfield_parser.get(),
                                                 Bzfield_parser.get(),
                                                 Bx_nodal_flag, By_nodal_flag,
                                                 Bz_nodal_flag, lev);
       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Bfield_aux[lev][0].get(),
                                                    Bfield_aux[lev][1].get(),
                                                    Bfield_aux[lev][2].get(),
                                                    Bxfield_parser.get(),
                                                    Byfield_parser.get(),
                                                    Bzfield_parser.get(),
                                                    Bx_nodal_flag, By_nodal_flag,
                                                    Bz_nodal_flag, lev);

          InitializeExternalFieldsOnGridUsingParser(Bfield_cp[lev][0].get(),
                                                    Bfield_cp[lev][1].get(),
                                                    Bfield_cp[lev][2].get(),
                                                    Bxfield_parser.get(),
                                                    Byfield_parser.get(),
                                                    Bzfield_parser.get(),
                                                    Bx_nodal_flag, By_nodal_flag,
                                                    Bz_nodal_flag, lev);
       }
    }

    // if the input string for the E-field is "parse_e_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (E_ext_grid_s == "parse_e_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort("E and B parser for external fields does not work with RZ -- TO DO");
#endif

       std::string Glenn_Ex;
       std::string Glenn_Ey;
       std::string Glenn_Ez;
       
       Store_parserString(pp, QED::Ex, Glenn_Ex);
       Store_parserString(pp, QED::Ey, Glenn_Ey);
       Store_parserString(pp, QED::Ez, Glenn_Ez);

       Exfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_Ex,{"x","y","z"})));
       Eyfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_Ey,{"x","y","z"})));
       Ezfield_parser.reset(new ParserWrapper<3>(
                                makeParser(Glenn_Ez,{"x","y","z"})));

       // Initialize Efield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Efield_fp[lev][0].get(),
                                                 Efield_fp[lev][1].get(),
                                                 Efield_fp[lev][2].get(),
                                                 Exfield_parser.get(),
                                                 Eyfield_parser.get(),
                                                 Ezfield_parser.get(),
                                                 Ex_nodal_flag, Ey_nodal_flag,
                                                 Ez_nodal_flag, lev);
       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Efield_aux[lev][0].get(),
                                                    Efield_aux[lev][1].get(),
                                                    Efield_aux[lev][2].get(),
                                                    Exfield_parser.get(),
                                                    Eyfield_parser.get(),
                                                    Ezfield_parser.get(),
                                                    Ex_nodal_flag, Ey_nodal_flag,
                                                    Ez_nodal_flag, lev);

          InitializeExternalFieldsOnGridUsingParser(Efield_cp[lev][0].get(),
                                                    Efield_cp[lev][1].get(),
                                                    Efield_cp[lev][2].get(),
                                                    Exfield_parser.get(),
                                                    Eyfield_parser.get(),
                                                    Ezfield_parser.get(),
                                                    Ex_nodal_flag, Ey_nodal_flag,
                                                    Ez_nodal_flag, lev);
       }
    }

    if (F_fp[lev]) {
        F_fp[lev]->setVal(0.0);
    }

    if (rho_fp[lev]) {
        rho_fp[lev]->setVal(0.0);
    }

    if (F_cp[lev]) {
        F_cp[lev]->setVal(0.0);
    }

    if (rho_cp[lev]) {
        rho_cp[lev]->setVal(0.0);
    }

    if (costs[lev]) {
        costs[lev]->setVal(0.0);
    }
}

#ifdef WARPX_USE_PSATD_HYBRID

void
WarpX::InitLevelDataFFT (int lev, Real time)
{

    Efield_fp_fft[lev][0]->setVal(0.0);
    Efield_fp_fft[lev][1]->setVal(0.0);
    Efield_fp_fft[lev][2]->setVal(0.0);
    Bfield_fp_fft[lev][0]->setVal(0.0);
    Bfield_fp_fft[lev][1]->setVal(0.0);
    Bfield_fp_fft[lev][2]->setVal(0.0);
    current_fp_fft[lev][0]->setVal(0.0);
    current_fp_fft[lev][1]->setVal(0.0);
    current_fp_fft[lev][2]->setVal(0.0);
    rho_fp_fft[lev]->setVal(0.0);

    if (lev > 0)
    {
        Efield_cp_fft[lev][0]->setVal(0.0);
        Efield_cp_fft[lev][1]->setVal(0.0);
        Efield_cp_fft[lev][2]->setVal(0.0);
        Bfield_cp_fft[lev][0]->setVal(0.0);
        Bfield_cp_fft[lev][1]->setVal(0.0);
        Bfield_cp_fft[lev][2]->setVal(0.0);
        current_cp_fft[lev][0]->setVal(0.0);
        current_cp_fft[lev][1]->setVal(0.0);
        current_cp_fft[lev][2]->setVal(0.0);
        rho_cp_fft[lev]->setVal(0.0);
    }

}

#endif

void
WarpX::InitializeExternalFieldsOnGridUsingParser (
       MultiFab *mfx, MultiFab *mfy, MultiFab *mfz,
       ParserWrapper<3> *xfield_parser, ParserWrapper<3> *yfield_parser,
       ParserWrapper<3> *zfield_parser, IntVect x_nodal_flag,
       IntVect y_nodal_flag, IntVect z_nodal_flag,
       const int lev)
{

    const auto dx_lev = geom[lev].CellSizeArray();
    const RealBox& real_box = geom[lev].ProbDomain();
    for ( MFIter mfi(*mfx, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& tbx = mfi.tilebox(x_nodal_flag);
       const Box& tby = mfi.tilebox(y_nodal_flag);
       const Box& tbz = mfi.tilebox(z_nodal_flag);

       auto const& mfxfab = mfx->array(mfi);
       auto const& mfyfab = mfy->array(mfi);
       auto const& mfzfab = mfz->array(mfi);

       auto const& mfx_IndexType = (*mfx).ixType();
       auto const& mfy_IndexType = (*mfy).ixType();
       auto const& mfz_IndexType = (*mfz).ixType();

       // Initialize IntVect based on the index type of multiFab
       // 0 if cell-centered, 1 if node-centered.
       IntVect mfx_type(AMREX_D_DECL(0,0,0));
       IntVect mfy_type(AMREX_D_DECL(0,0,0));
       IntVect mfz_type(AMREX_D_DECL(0,0,0));

       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           mfx_type[idim] = mfx_IndexType.nodeCentered(idim);
           mfy_type[idim] = mfy_IndexType.nodeCentered(idim);
           mfz_type[idim] = mfz_IndexType.nodeCentered(idim);
       }

       amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Shift required in the x-, y-, or z- position
                // depending on the index type of the multifab
                Real fac_x = (1.0 - mfx_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#else
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfx_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the x-component of the field.
                mfxfab(i,j,k) = (*xfield_parser)(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfy_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif (AMREX_SPACEDIM==3)
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfx_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the y-component of the field.
                mfyfab(i,j,k)  = (*yfield_parser)(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfz_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif (AMREX_SPACEDIM==3)
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfz_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the z-component of the field.
                mfzfab(i,j,k) = (*zfield_parser)(x,y,z);
            }
        );
    }

}
