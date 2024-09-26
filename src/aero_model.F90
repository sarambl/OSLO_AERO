module aero_model

  !===============================================================================
  ! Oslo Aerosol Model
  ! Note: SPCAM not supported here
  !===============================================================================

  use shr_kind_mod,             only: r8 => shr_kind_r8
  use spmd_utils,               only: mpicom, mstrid=>masterprocid, masterproc
  use spmd_utils,               only: mpi_logical, mpi_real8, mpi_character, mpi_integer,  mpi_success
  use namelist_utils,           only: find_group_name
  use constituents,             only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,                   only: pcols, pver, pverp
  use phys_control,             only: phys_getopts, cam_physpkg_is
  use cam_abortutils,           only: endrun
  use cam_logfile,              only: iulog
  use perf_mod,                 only: t_startf, t_stopf
  use camsrfexch,               only: cam_in_t, cam_out_t
  use aerodep_flx,              only: aerodep_flx_prescribed
  use physics_types,            only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,           only: physics_buffer_desc, pbuf_get_field, pbuf_get_index, pbuf_set_field
  use physconst,                only: gravit, rair, rhoh2o, pi
  use time_manager,             only: get_nstep
  use cam_history,              only: outfld, fieldname_len, addfld, add_default, horiz_only
  use chem_mods,                only: gas_pcnst, adv_mass
  use mo_tracname,              only: solsym
  use mo_setsox,                only: setsox
  use mo_mass_xforms,           only: vmr2mmr, mmr2vmr, mmr2vmri
  use mo_chem_utls,             only: get_rxt_ndx, get_spc_ndx
  use ref_pres,                 only: top_lev => clim_modal_aero_top_lev
  use wv_saturation,            only: qsat_water
  !
  use oslo_aero_share,          only: nmodes_oslo=>nmodes, nbmodes
  use oslo_aero_share,          only: originalSigma, originalNumberMedianRadius
  use oslo_aero_share,          only: init_interp_constants
  use oslo_aero_share,          only: rTabMin, rTabMax, nk, normnk, rBinEdge, rBinMidpoint
  use oslo_aero_share,          only: volumeToNumber, numberToSurface, nBinsTab, rMinAquousChemistry
  use oslo_aero_share,          only: calculateLognormalCDF, calculatedNdLogR, calculateNumberConcentration
  use oslo_aero_share,          only: calculateEquivalentDensityOfFractalMode
  use oslo_aero_share,          only: numberFractionAvailableAqChem, isTracerInMode
  use oslo_aero_share,          only: chemistryIndex, physicsIndex
  use oslo_aero_share,          only: qqcw_get_field, numberOfProcessModeTracers
  use oslo_aero_share,          only: lifeCycleNumberMedianRadius, rhopart, lifeCycleSigma
  use oslo_aero_share,          only: l_so4_a2, l_bc_n, l_bc_ax
  use oslo_aero_share,          only: MODE_IDX_BC_NUC, MODE_IDX_BC_EXT_AC
  !smb++
  use oslo_aero_share,          only: oslo_aero_share_readnl
  use oslo_aero_condtend,       only: oslo_aero_condtend_readnl
  !smb--
  use oslo_aero_control,        only: oslo_aero_ctl_readnl, use_aerocom
  use oslo_aero_depos,          only: oslo_aero_depos_init, oslo_aero_depos_readnl
  use oslo_aero_depos,          only: oslo_aero_depos_dry, oslo_aero_depos_wet, oslo_aero_wetdep_init
  use oslo_aero_coag,           only: initializeCoagulation, coagtend, clcoag
  use oslo_aero_condtend,       only: N_COND_VAP, COND_VAP_ORG_SV, COND_VAP_ORG_LV, COND_VAP_H2SO4
  use oslo_aero_condtend,       only: initializeCondensation, condtend
  use oslo_aero_seasalt,        only: oslo_aero_seasalt_init, oslo_aero_seasalt_emis, seasalt_active
  use oslo_aero_dust,           only: oslo_aero_dust_init, oslo_aero_dust_emis, dust_active
  use oslo_aero_ocean,          only: oslo_aero_ocean_init, oslo_aero_dms_emis
  use oslo_aero_share,          only: getNumberofTracersInMode, getCloudTracerIndexDirect, getCloudTracerName
  use oslo_aero_share,          only: getCloudTracerName, getTracerIndex, aero_register
  use oslo_aero_sox_cldaero,    only: sox_cldaero_init
  use oslo_aero_microp,         only: oslo_aero_microp_readnl
  use oslo_aero_sw_tables,      only: initopt
  use oslo_aero_aerodry_tables, only: initdry
  use oslo_aero_aerocom_tables, only: initaeropt
  use oslo_aero_logn_tables,    only: initlogn

  use modal_aero_wateruptake, only: modal_strat_sulfate

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch     ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep         ! aerosol dry deposition and sediment
  public :: aero_model_wetdep         ! aerosol wet removal
  public :: aero_model_emissions      ! aerosol emissions
  public :: aero_model_surfarea       ! tropopspheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea ! stratospheric aerosol wet surface area for chemistry

  private :: aero_model_constants

  ! Misc private data
  integer :: pblh_idx= 0
  integer :: ndx_h2so4, ndx_soa_lv, ndx_soa_sv ! for surf_area_dens
  logical :: convproc_do_aer

  ! Namelist variables
  real(r8) :: sol_facti_cloud_borne   = 1._r8
  real(r8) :: sol_factb_interstitial  = 0.1_r8
  real(r8) :: sol_factic_interstitial = 0.4_r8

!=============================================================================
contains
!=============================================================================

  subroutine aero_model_readnl(nlfilename)

    ! read aerosol namelist options

    character(len=*), intent(in) :: nlfilename  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    namelist /aerosol_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial, modal_strat_sulfate
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open(newunit=unitn, file=trim(nlfilename), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if
    call mpi_bcast(sol_facti_cloud_borne, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_facti_cloud_borne")
    call mpi_bcast(sol_factb_interstitial, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_factb_interstitial")
    call mpi_bcast(sol_factic_interstitial, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_factic_interstitial")
    call mpi_bcast(modal_strat_sulfate, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: modal_strat_sulfate")
    !smb++
    call oslo_aero_share_readnl(nlfilename)
    call oslo_aero_condtend_readnl(nlfilename)
    !smb --
    call oslo_aero_ctl_readnl(nlfilename)
    call oslo_aero_microp_readnl(nlfilename)
    call oslo_aero_depos_readnl(nlfilename)

  end subroutine aero_model_readnl

  !=============================================================================
  subroutine aero_model_register()

    call aero_register()

  end subroutine aero_model_register

  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    integer           :: icnst,id
    character(len=20) :: dummy
    logical           :: history_aerosol ! Output MAM or SECT aerosol tendencies
    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    !------------------------------------

    call phys_getopts(history_aerosol_out=history_aerosol, convproc_do_aer_out=convproc_do_aer)

    call aero_model_constants
    call init_interp_constants() ! table initialization constants
    call initopt()               ! table initialization
    call initlogn()              ! table initialization
    if (use_aerocom) then
       call initdry()               ! table initialization
       call initaeropt()            ! table initialization
    end if
    call initializeCondensation()
    call oslo_aero_ocean_init()
    call oslo_aero_depos_init(pbuf2d)
    call oslo_aero_dust_init()
    call oslo_aero_seasalt_init()
    call oslo_aero_wetdep_init()

    dummy = 'RAM1'
    call addfld (dummy,horiz_only, 'A','frac','RAM1')
    if ( history_aerosol ) then
       call add_default (dummy, 1, ' ')
    endif

    dummy = 'airFV'
    call addfld (dummy,horiz_only, 'A','frac','FV')
    if ( history_aerosol ) then
       call add_default (dummy, 1, ' ')
    endif

    ! Get height of boundary layer for boundary layer nucleation
    pblh_idx = pbuf_get_index('pblh')

    call cnst_get_ind ( "H2SO4", ndx_h2so4, abort=.true. )
    ndx_h2so4 = chemistryIndex(ndx_h2so4)
    call cnst_get_ind ( "SOA_LV", ndx_soa_lv,abort=.true.)
    ndx_soa_lv = chemistryIndex(ndx_soa_lv)
    call cnst_get_ind ( "SOA_SV", ndx_soa_sv, abort=.true.)
    ndx_soa_sv = chemistryIndex(ndx_soa_sv)

    do icnst = 1,gas_pcnst
       unit_basename = 'kg'  ! Units 'kg' or '1'

       call addfld( 'GS_'//trim(solsym(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
            trim(solsym(icnst))//' gas chemistry/wet removal (for gas species)')

       call addfld( 'AQ_'//trim(solsym(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
            trim(solsym(icnst))//' aqueous chemistry (for gas species)')

       if(physicsIndex(icnst)<=pcnst) then
          if (getCloudTracerIndexDirect(physicsIndex(icnst)) > 0)then
             call addfld( 'AQ_'//getCloudTracerName(physicsIndex(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
                  trim(solsym(icnst))//' aqueous chemistry (for cloud species)')
          end if
       end if

       if ( history_aerosol ) then
          call add_default( 'GS_'//trim(solsym(icnst)), 1, ' ')
          call add_default( 'AQ_'//trim(solsym(icnst)), 1, ' ')
          if(physicsIndex(icnst)<=pcnst) then
             if(getCloudTracerIndexDirect(physicsIndex(icnst))>0)then
                call add_default( 'AQ_'//getCloudTracerName(physicsIndex(icnst)),1,' ')
             end if
          end if
       endif
    enddo

    call addfld ('NUCLRATE',(/'lev'/), 'A','#/cm3/s','Nucleation rate')
    call addfld ('FORMRATE',(/'lev'/), 'A','#/cm3/s','Formation rate of 12nm particles')
    call addfld ('COAGNUCL',(/'lev'/), 'A', '/s','Coagulation sink for nucleating particles')
    call addfld ('GRH2SO4',(/'lev'/), 'A', 'nm/hour','Growth rate H2SO4')
    call addfld ('GRSOA',(/'lev'/),'A','nm/hour','Growth rate SOA')
    call addfld ('GR',(/'lev'/), 'A', 'nm/hour','Growth rate, H2SO4+SOA')
    call addfld ('NUCLSOA',(/'lev'/),'A','kg/kg','SOA nucleate')
    call addfld ('ORGNUCL',(/'lev'/),'A','kg/kg','Organic gas available for nucleation')

    if(history_aerosol)then
       call add_default ('NUCLRATE', 1, ' ')
       call add_default ('FORMRATE', 1, ' ')
       call add_default ('COAGNUCL', 1, ' ')
       call add_default ('GRH2SO4', 1, ' ')
       call add_default ('GRSOA', 1, ' ')
       call add_default ('GR', 1, ' ')
       call add_default ('NUCLSOA', 1, ' ')
       call add_default ('ORGNUCL', 1, ' ')
    end if

    call addfld( 'XPH_LWC',    (/ 'lev' /), 'A','kg/kg',   'pH value multiplied by lwc')
    call addfld ('AQSO4_H2O2', horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to H2O2')
    call addfld ('AQSO4_O3',   horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to O3')

    if ( history_aerosol ) then
       call add_default ('XPH_LWC', 1, ' ')
       call add_default ('AQSO4_H2O2', 1, ' ')
       call add_default ('AQSO4_O3', 1, ' ')
    endif

  end subroutine aero_model_init

  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    ! args
    type(physics_state),    intent(in)    :: state    ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:) ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in   ! import state
    real(r8),               intent(in)    :: dt       ! time step
    type(cam_out_t),        intent(inout) :: cam_out  ! export state
    type(physics_ptend),    intent(out)   :: ptend    ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

    ! local vars
    real(r8) :: oslo_dgnumwet(pcols, pver, 0:nmodes_oslo)
    real(r8) :: oslo_wetdens(pcols, pver, 0:nmodes_oslo)
    real(r8) :: oslo_dgnumwet_processmodes(pcols, pver, numberOfProcessModeTracers)
    real(r8) :: oslo_wetdens_processmodes(pcols, pver, numberOfProcessModeTracers)

    oslo_wetdens(:,:,:) = 0._r8
    call calcaersize_sub(state%ncol, state%t, state%q(1,1,1), state%pmid, state%pdel, &
         oslo_dgnumwet, oslo_wetdens, oslo_dgnumwet_processmodes, oslo_wetdens_processmodes)

    call oslo_aero_depos_dry(state%lchnk, state%ncol, state%psetcols, &
         state%t, state%pmid, state%pdel, state%pint, state%q, &
         cam_in%landfrac, cam_in%icefrac, cam_in%ocnfrac, cam_in%fv, cam_in%ram1, cam_in%cflx, &
         pbuf, obklen, ustar, dt, &
         oslo_dgnumwet, oslo_wetdens, oslo_dgnumwet_processmodes, oslo_wetdens_processmodes, &
         cam_out, ptend)

  endsubroutine aero_model_drydep

  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    call oslo_aero_depos_wet(state%lchnk, state%ncol, state%psetcols, state%pmid, state%pdel, state%q, state%t, &
         dt, dlf, cam_out, ptend, pbuf)

  endsubroutine aero_model_wetdep

  !=============================================================================
  subroutine aero_model_surfarea(mmr, radmean, relhum, pmid, temp, strato_sad, sulfate, rho, ltrop, &
       dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_trop, reff_trop )

    !-------------------------------------------------------------------------
    ! provides wet tropospheric aerosol surface area info for modal aerosols
    ! called from mo_usrrxt
    !-------------------------------------------------------------------------

    ! arguments
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: rho(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_trop(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)

    ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)
  
    integer :: i,k

    beglev(:ncol)=ltrop(:ncol)+1
    endlev(:ncol)=pver

    ! diameter left out as an argument                                                       
      call surf_area_dens(ncol, mmr, pmid, temp, beglev, endlev, sad_trop, reff_trop, sfc=sfc) 

      do i = 1,ncol                                                                            
         do k = ltrop(i)+1, pver                                                               
            ! djlo : do not use the 0-th mode                                                  
            dm_aer(i,k,:) = 2._r8 * lifeCycleNumberMedianRadius(1:nmodes_oslo) * 1.e-2_r8 ! radius ==> diameter, m ==> cm
         enddo                                                                                 
      enddo

  end subroutine aero_model_surfarea

  !=============================================================================
  subroutine aero_model_strat_surfarea( ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    !-------------------------------------------------------------------------
    ! provides WET stratospheric aerosol surface area info for modal aerosols
    ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
    !-------------------------------------------------------------------------

    ! arguments
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

   ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)   

    reff_strat = 0._r8
    strato_sad = 0._r8

    if (.not. modal_strat_sulfate) return

    beglev(:ncol)=top_lev
    endlev(:ncol)=ltrop(:ncol)
    call surf_area_dens(ncol, mmr, pmid, temp, beglev, endlev, strato_sad, reff_strat)

    return

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, troplev, delt, reaction_rates, &
       tfld, pmid, pdel, mbar, relhum, &
       zm,  qh2o, cwat, cldfr, cldnum, &
       airdens, invariants, del_h2so4_gasprod,  &
       vmr0, vmr, pbuf )

    ! arguments
    integer,  intent(in)    :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in)    :: ncol                   ! number columns in chunk
    integer,  intent(in)    :: lchnk                  ! chunk index
    integer,  intent(in)    :: troplev(pcols)
    real(r8), intent(in)    :: delt                   ! time step size (sec)
    real(r8), intent(in)    :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in)    :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in)    :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in)    :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in)    :: relhum(:,:)            ! relative humidity
    real(r8), intent(in)    :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in)    :: invariants(:,:,:)
    real(r8), intent(in)    :: zm(:,:)
    real(r8), intent(in)    :: qh2o(:,:)
    real(r8), intent(in)    :: cwat(:,:)              ! cloud liquid water content (kg/kg)
    real(r8), intent(in)    :: cldfr(:,:)
    real(r8), intent(in)    :: cldnum(:,:)            ! droplet number concentration (#/kg)
    real(r8), intent(inout) :: del_h2so4_gasprod(:,:) ! [molec/molec/sec]
    real(r8), intent(in)    :: vmr0(:,:,:)            ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)             ! mixing ratios ( vmr )
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars
    integer, parameter :: nmodes_aq_chem = 1
    integer  :: icol,ilev
    integer  :: imode,icnst,itrac
    integer  :: nstep
    real(r8) :: wrk(ncol)
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)   ! cloud-borne aerosol (vmr)
    real(r8) :: del_h2so4_aeruptk(ncol,pver)
    real(r8) :: del_h2so4_aqchem(ncol,pver)
    real(r8) :: mmr_cond_vap_start_of_timestep(pcols,pver,N_COND_VAP)
    real(r8) :: mmr_cond_vap_gasprod(pcols,pver,N_COND_VAP)
    real(r8) :: del_soa_lv_gasprod(ncol,pver)
    real(r8) :: del_soa_sv_gasprod(ncol,pver)
    real(r8) :: dvmrdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: dvmrcwdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: mmr_tend_ncols(ncol, pver, gas_pcnst)
    real(r8) :: mmr_tend_pcols(pcols, pver, gas_pcnst)
    integer  :: cond_vap_idx
    real(r8) :: aqso4(ncol,nmodes_aq_chem)   ! aqueous phase chemistry
    real(r8) :: aqh2so4(ncol,nmodes_aq_chem) ! aqueous phase chemistry
    real(r8) :: aqso4_h2o2(ncol)             ! SO4 aqueous phase chemistry due to H2O2
    real(r8) :: aqso4_o3(ncol)               ! SO4 aqueous phase chemistry due to O3
    real(r8) :: xphlwc(ncol,pver)            ! pH value multiplied by lwc
    real(r8) :: delt_inverse                 ! 1 / timestep
    real(r8), pointer :: pblh(:)
    character(len=32) :: name

    nstep = get_nstep()

    delt_inverse = 1.0_r8 / delt

    ! Get height of boundary layer (needed for boundary layer nucleation)
    call pbuf_get_field(pbuf, pblh_idx, pblh)

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do icnst = 1, gas_pcnst
       wrk(:) = 0._r8
       do ilev = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
       end do
       name = 'GS_'//trim(solsym(icnst))
       call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

    ! Get mass mixing ratios at start of time step
    call vmr2mmr( vmr0, mmr_tend_ncols, mbar, ncol )
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_H2SO4)  = mmr_tend_ncols(1:ncol,:,ndx_h2so4)
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_LV) = mmr_tend_ncols(1:ncol,:,ndx_soa_lv)
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_SV) = mmr_tend_ncols(1:ncol,:,ndx_soa_sv)
    !
    ! Aerosol processes ...
    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - vmr0(1:ncol,:,ndx_h2so4)
    endif

    del_soa_lv_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_soa_lv) - vmr0(1:ncol,:,ndx_soa_lv)
    del_soa_sv_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_soa_sv) - vmr0(1:ncol,:,ndx_soa_sv)

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

    !Save intermediate concentrations
    dvmrdt_sv1 = vmr
    dvmrcwdt_sv1 = vmrcw

    ! aqueous chemistry ...
    call setsox( ncol, lchnk, loffset, delt, pmid, pdel, tfld, mbar, cwat, &
         cldfr, cldnum, airdens, invariants, vmrcw, vmr, xphlwc, &
         aqso4, aqh2so4, aqso4_h2o2, aqso4_o3)

    call outfld( 'AQSO4_H2O2', aqso4_h2o2(:ncol), ncol, lchnk)
    call outfld( 'AQSO4_O3',   aqso4_o3(:ncol),   ncol, lchnk)
    call outfld( 'XPH_LWC',    xphlwc(:ncol,:),   ncol, lchnk )

    ! vmr tendency from aqchem and soa routines
    dvmrdt_sv1 = (vmr - dvmrdt_sv1)/delt
    dvmrcwdt_sv1 = (vmrcw - dvmrcwdt_sv1)/delt

    if (ndx_h2so4 > 0)then
       del_h2so4_aqchem(:ncol,:) = dvmrdt_sv1(:ncol,:,ndx_h2so4)*delt !"production rate" of H2SO4
    else
       del_h2so4_aqchem(:ncol,:) = 0.0_r8
    end if

    do icnst = 1,gas_pcnst
       wrk(:ncol) = 0._r8
       do ilev = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt_sv1(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
       end do
       name = 'AQ_'//trim(solsym(icnst))
       call outfld( name, wrk(:ncol), ncol, lchnk )

       !In oslo aero also write out the tendencies for the
       !cloud borne aerosols...
       itrac = physicsIndex(icnst)
       if (itrac <= pcnst) then
          if (getCloudTracerIndexDirect(itrac) > 0)then
             name = 'AQ_'//trim(getCloudTracerName(itrac))
             wrk(:ncol)=0.0_r8
             do ilev=1,pver
                wrk(:ncol) = wrk(:ncol) + dvmrcwdt_sv1(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
             end do
             call outfld( name, wrk(:ncol), ncol, lchnk )
          end if
       end if
    enddo

    ! condensation
    call vmr2mmr( vmr, mmr_tend_ncols, mbar, ncol )
    do ilev = 1,pver
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_H2SO4) = adv_mass(ndx_h2so4) &
            * (del_h2so4_gasprod(:ncol,ilev)+del_h2so4_aqchem(:ncol,ilev)) / mbar(:ncol,ilev)/delt
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_ORG_LV) = adv_mass(ndx_soa_lv) &
            * del_soa_lv_gasprod(:ncol,ilev) / mbar(:ncol,ilev)/delt
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_ORG_SV) = adv_mass(ndx_soa_sv) &
            * del_soa_sv_gasprod(:ncol,ilev) / mbar(:ncol,ilev)/delt
    end do

    ! This should not happen since there are only production terms for these gases! !
    do cond_vap_idx=1,N_COND_VAP
       where(mmr_cond_vap_gasprod(:ncol,:,cond_vap_idx) < 0.0_r8)
          mmr_cond_vap_gasprod(:ncol,:,cond_vap_idx) = 0.0_r8
       end where
    end do
    mmr_tend_ncols(:ncol,:,ndx_h2so4)  = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_H2SO4)
    mmr_tend_ncols(:ncol,:,ndx_soa_lv) = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_LV)
    mmr_tend_ncols(:ncol,:,ndx_soa_sv) = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_SV)

    ! Rest of microphysics have pcols dimension
    mmr_tend_pcols(:ncol,:,:) = mmr_tend_ncols(:ncol,:,:)

    ! Condensation
    ! Note use of "zm" here. In CAM5.3-implementation "zi" was used..
    ! zm is passed through the generic interface, and it should not change much
    ! to check if "zm" is below boundary layer height instead of zi
    call condtend( lchnk, mmr_tend_pcols, mmr_cond_vap_gasprod,tfld, pmid, &
         pdel, delt, ncol, pblh, zm, qh2o)  ! cka

    ! Coagulation
    ! OS 280415  Concentratiions in cloud water is in vmr space and as a
    ! temporary variable  (vmrcw) Coagulation between aerosol and cloud
    ! droplets moved to after vmrcw is moved into qqcw (in mmr spac)
    call coagtend( mmr_tend_pcols, pmid, pdel, tfld, delt_inverse, ncol, lchnk)

    ! Convert cloud water to mmr again ==> values in buffer
    call vmr2qqcw( vmrcw, mbar, ncol, loffset, pbuf )

    ! Call cloud coagulation routines (all in mass mixing ratios)
    call clcoag( mmr_tend_pcols, pmid, pdel, tfld, cldnum ,cldfr, delt_inverse, ncol, lchnk, loffset, pbuf)

    ! Make sure mmr==> vmr is done correctly
    mmr_tend_ncols(:ncol,:,:) = mmr_tend_pcols(:ncol,:,:)

    ! Go back to volume mixing ratio for chemistry
    call mmr2vmr( mmr_tend_ncols, vmr, mbar, ncol )

  end subroutine aero_model_gasaerexch

  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )

    ! Arguments:
    type(physics_state), intent(in)    :: state   ! Physics state variables
    type(cam_in_t),      intent(inout) :: cam_in  ! import state

    if (dust_active) then
       call oslo_aero_dust_emis( state%lchnk, state%ncol, cam_in%dstflx, cam_in%cflx)
    endif

    if (seasalt_active) then
       call oslo_aero_seasalt_emis(state%ncol, state%lchnk, &
            state%u(:,pver), state%v(:,pver), state%zm(:,pver), &
            cam_in%ocnfrac, cam_in%icefrac, cam_in%sst, cam_in%cflx)
    endif

    ! Pick up correct DMS emissions (replace values from file if requested)
    call oslo_aero_dms_emis(state%ncol, state%lchnk, &
         state%u(:,pver), state%v(:,pver), state%zm(:,pver), &
         cam_in%ocnfrac, cam_in%icefrac, cam_in%sst, cam_in%fdms, cam_in%cflx)

  end subroutine aero_model_emissions

  !=============================================================================
  ! private methods
  !=============================================================================

  subroutine surf_area_dens( ncol, mmr, pmid, temp, beglev, endlev, sad, reff, sfc )
    ! The input, diam, used by MAM, is not used in OSLO_AERO
    use mo_constants     , only : pi

    ! dummy args
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: mmr(:,:,:)
    real(r8), intent(in)  :: pmid(:,:)
    real(r8), intent(in)  :: temp(:,:)
    ! real(r8), intent(in)  :: diam(:,:,:) ! diam, used by MAM, is not used in OSLO_AERO
    integer,  intent(in)  :: beglev(:)
    integer,  intent(in)  :: endlev(:)
    real(r8), intent(out) :: sad(:,:)
    real(r8), intent(out) :: reff(:,:)
    real(r8),optional, intent(out) :: sfc(:,:,:)

    ! local vars
    !HAVE TO GET RID OF THIS MODE 0!! MESSES UP EVERYTHING!!
    real(r8)         :: numberConcentration(pcols,pver,0:nmodes_oslo)
    real(r8), target :: sad_mode(pcols,pver, nmodes_oslo)
    real(r8)         :: vol_mode(pcols,pver, nmodes_oslo)
    real(r8)         :: vol(pcols,pver)                  
    real(r8) :: rho_air(pcols,pver)
    integer :: m
    integer :: i,k

    ! Compute surface aero for each mode.
    ! Total over all modes as the surface area for chemical reactions.

    !Get air density in all layers
    do k=1,pver
       do i=1,ncol
          rho_air(i,k) = pmid(i,k)/(temp(i,k)*rair)
       end do
    end do
    !    
    !Get number concentrations in all layers
    call calculateNumberConcentration(ncol, mmr, rho_air, numberConcentration, isChemistry=.true.)

    !Convert to area using lifecycle-radius
    sad_mode = 0._r8
    sad      = 0._r8
    vol_mode = 0._r8
    vol      = 0._r8
    reff     = 0._r8

    do i=1,ncol ! djlo : added explicit loop over i (because of ltrop(i) dependence)
       do k=beglev(i),endlev(i)
          do m=1,nmodes_oslo
             if ( m .eq. 1 &
             .or. m .eq. 2 &
             .or. m .eq. 4 &
             .or. m .eq. 5 ) then
!               currently only over modes 1,2,4 and 5   
!               might be extended in the future with modes 12 and 14
                sad_mode(i,k,m) = numberConcentration(i,k,m)*numberToSurface(m)*1.e-2_r8 !m2/m3 ==> cm2/cm3

                vol_mode(i,k,m) = numberConcentration(i,k,m) &
                                * 4._r8 / 3._r8 * pi * lifeCycleNumberMedianRadius(m)**3._r8 &
                                * dexp(4.5_r8 * log(lifeCycleSigma(m)) *log(lifeCyclesigma(m))) ! m3/m3 = cm3/cm3  
             endif
          end do

          sad(i,k) = sum(sad_mode(i,k,:))
          vol(i,k) = sum(vol_mode(i,k,:))
          reff(i,k) = 3._r8 * vol(i,k) / sad(i,k) ! djlo : maybe if-test needed when sad=0?
       end do
    end do

    if ( present(sfc) ) then
       sfc(:,:,:) = sad_mode(:,:,:)
    endif

  end subroutine surf_area_dens

  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)

    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: ilev, icnst
    real(r8), pointer :: fldcw(:,:)

    do icnst=1,gas_pcnst
       if( adv_mass(icnst) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, icnst+im)
          if(associated(fldcw)) then
             do ilev = 1,pver
                vmr(:ncol,ilev,icnst) = mbar(:ncol,ilev) * fldcw(:ncol,ilev) / adv_mass(icnst)
             end do
          else
             vmr(:,:,icnst) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr

  !=============================================================================
  subroutine vmr2qqcw( vmr, mbar, ncol, im, pbuf )

    !-----------------------------------------------------------------
    ! Convert from volume to mass mixing ratio
    !-----------------------------------------------------------------

    ! Arguments
    integer , intent(in) :: ncol, im
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! Local variables
    integer :: ilev ! level index
    integer :: icnst ! tracer index
    real(r8), pointer :: fldcw(:,:)

    ! The non-group species
    do icnst = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, icnst+im)
       if( adv_mass(icnst) /= 0._r8 .and. associated(fldcw)) then
          do ilev = 1,pver
             fldcw(:ncol,ilev) = adv_mass(icnst) * vmr(:ncol,ilev,icnst) / mbar(:ncol,ilev)
          end do
       end if
    end do
  end subroutine vmr2qqcw

  !=============================================================================
  subroutine aero_model_constants()

    ! A number of constants used in the emission and size-calculation in CAM-Oslo

    ! local variables
    integer  :: imode,ibin
    real(r8) :: rhob(0:nmodes_oslo) !density of background aerosol in mode
    real(r8) :: rhorbc         !This has to do with fractal dimensions of bc, come back to this!!
    real(r8) :: sumnormnk
    real(r8) :: totalLogDelta
    real(r8) :: logDeltaBin
    real(r8) :: logNextEdge
    !---------------------------------------

    rhob(:)            =-1.0_r8
    volumeToNumber(:)  =-1.0_r8
    numberToSurface(:) =-1.0_r8

    !Prepare modal properties
    do imode =0,nmodes_oslo
       if(getNumberOfTracersInMode(imode) > 0)then

          !Approximate density of mode
          !density of mode is density of first species in mode
          rhob(imode)  = rhopart(getTracerIndex(imode,1,.false.))

          !REPLACE THE EFACT-VARIABLE WITH THIS!!
          volumeToNumber(imode) = 1.0_r8 / &
               ( DEXP ( 4.5_r8 * ( log(originalSigma(imode)) * log(originalSigma(imode)) ) ) &
               *(4.0_r8/3.0_r8)*pi*(originalNumberMedianRadius(imode))**3 )

          numberToSurface(imode) = 4.0_r8*pi*lifeCycleNumberMedianRadius(imode)*lifeCycleNumberMedianRadius(imode)&
               *DEXP(log(lifeCycleSigma(imode))*log(lifeCycleSigma(imode)))
       end if
    end do

    !Find radius in edges and midpoints of bin
    rBinEdge(1) = rTabMin
    totalLogDelta = log(rTabMax/rTabMin)
    logDeltaBin = totalLogDelta / nBinsTab
    do ibin = 2,nBinsTab+1
       logNextEdge = log(rBinEdge(ibin-1)) + logDeltaBin
       rBinEdge(ibin) = DEXP(logNextEdge)
       rBinMidPoint(ibin-1) = sqrt(rBinEdge(ibin)*rBinEdge(ibin-1))
    end do

    !Calculate the fraction of a mode which goes to aquous chemstry
    numberFractionAvailableAqChem(:)=0.0_r8
    do imode = 1,nbmodes
       if(isTracerInMode(imode,l_so4_a2))then
          numberFractionAvailableAqChem(imode) =  1.0_r8 -  &
               calculateLognormalCDF(rMinAquousChemistry,originalNumberMedianRadius(imode), originalSigma(imode))
       end if
    end do

    !Set the density of the fractal mode ==> we get lesser density
    !than the emitted density, so for a given mass emitted, we get
    !more number-concentration!! This is a way of simulating that the
    !aerosols take up more space
    rhorbc = calculateEquivalentDensityOfFractalMode(    &
         rhopart(l_bc_n),                                & !emitted density
         originalNumberMedianRadius(MODE_IDX_BC_NUC),    & !emitted size
         2.5_r8,                                         & !fractal dim
         originalNumberMedianRadius(MODE_IDX_BC_EXT_AC), & !diameter of mode
         originalSigma(MODE_IDX_BC_EXT_AC))                !sigma mode

    rhopart(l_bc_ax) = rhorbc
    !fxm: not the right place for this change of value,
    !but anyway.. this re-calculateion of tracer density
    !influences density of mode used in coagulation
    rhob(MODE_IDX_BC_EXT_AC)=rhorbc

    !Size distribution of the modes!
    !Unclear if this should use the radii assuming growth or not!
    !Mostly used in code where it is sensible to assume some growth has
    !happened, so it is used here
    do imode = 0,nmodes_oslo
       do ibin=1,nBinsTab
          !dN/dlogR (does not sum to one over size range)
          nk(imode,ibin) = calculatedNdLogR(rBinMidPoint(ibin), lifeCycleNumberMedianRadius(imode), lifeCycleSigma(imode))

          !dN (sums to one) over the size range
          normnk(imode,ibin) =logDeltaBin*nk(imode,ibin)
       enddo
    enddo  ! imode

    ! Normalized size distribution must sum to one (accept 2% error)
    do imode=0,nmodes_oslo
       sumNormNk = sum(normnk(imode,:))
       if(abs(sum(normnk(imode,:)) - 1.0_r8) > 2.0e-2_r8)then
          print*, "sum normnk", sum(normnk(imode,:))
          call endrun()
       endif
    enddo

    ! Initialize coagulation including Calculate the coagulation coefficients
    ! Note: Inaccurate density used!
    call initializeCoagulation(rhob, lifeCycleNumberMedianRadius)

  end subroutine aero_model_constants

  !=============================================================================
  subroutine calcaersize_sub( ncol, t, h2ommr, pmid, pdel,wetnumberMedianDiameter,wetrho,   &
       wetNumberMedianDiameter_processmode, wetrho_processmode)

    ! Seland Calculates mean volume size and hygroscopic growth for use in  dry deposition

    use oslo_aero_share

    integer,  intent(in)  :: ncol               ! number of columns
    real(r8), intent(in)  :: t(pcols,pver)      ! layer temperatures (K)
    real(r8), intent(in)  :: h2ommr(pcols,pver) ! layer specific humidity
    real(r8), intent(in)  :: pmid(pcols,pver)   ! layer pressure (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)  ! layer pressure thickness (Pa)
    real(r8), intent(out) :: wetNumberMedianDiameter(pcols,pver,0:nmodes_oslo)
    real(r8), intent(out) :: wetrho(pcols,pver,0:nmodes_oslo) ! wet aerosol density
    real(r8), intent(out) :: wetNumberMedianDiameter_processmode(pcols,pver,numberOfProcessModeTracers)
    real(r8), intent(out) :: wetrho_processmode(pcols,pver,numberOfProcessModeTracers)

    !     local variables
    integer  :: icol,itrac,ilev,imode,irelh,mm, tracerCounter
    real(r8) :: xrh(pcols,pver)
    real(r8) :: relhum(pcols,pver) ! Relative humidity
    real(r8) :: qs(pcols,pver)     ! saturation specific humidity
    real(r8) :: rmeanvol           ! Mean radius with respect to volume
    integer  :: irh1(pcols,pver),irh2(pcols,pver)
    integer  :: t_irh1,t_irh2
    real(r8) :: t_rh1,t_rh2,t_xrh,rr1,rr2
    real(r8) :: volumeFractionAerosol   !with respect to total (aerosol + water)
    real(r8) :: tmp1
    real(r8) :: wetrad_tmp(max_tracers_per_mode)
    real(r8) :: dry_rhopart_tmp(max_tracers_per_mode)
    real(r8) :: mixed_dry_rho


    !Get the tabulated rh in all grid cells
    do ilev=1,pver
       do icol=1,ncol
          call qsat_water(t(icol,ilev), pmid(icol,ilev), tmp1, qs(icol,ilev))
          xrh(icol,ilev) = h2ommr(icol,ilev)/qs(icol,ilev)
          xrh(icol,ilev) = max(xrh(icol,ilev),0.0_r8)
          xrh(icol,ilev) = min(xrh(icol,ilev),1.0_r8)
          relhum(icol,ilev)=xrh(icol,ilev)
          xrh(icol,ilev)=min(xrh(icol,ilev),rhtab(10))
       end do
    end do

    !Find the relh-index in all grid-points
    do irelh=1,SIZE(rhtab) - 1
       do ilev=1,pver
          do icol=1,ncol
             if (xrh(icol,ilev) >= rhtab(irelh) .and. xrh(icol,ilev) <= rhtab(irelh+1)) then
                irh1(icol,ilev)=irelh                !lower index
                irh2(icol,ilev)=irelh+1              !higher index
             end if
          end do
       end do
    end do

    do ilev=1,pver
       do icol=1,ncol

          !Get the indexes out as floating point single numbers
          t_irh1 = irh1(icol,ilev)
          t_irh2 = irh2(icol,ilev)
          t_rh1  = rhtab(t_irh1)
          t_rh2  = rhtab(t_irh2)
          t_xrh  = xrh(icol,ilev)

          do imode = 0, nmodes_oslo
             !Do some weighting to mass mean property
             !weighting by 1.5 is number median ==> volumetric mean
             !http://dust.ess.uci.edu/facts/psd/psd.pdf
             rmeanvol = lifeCycleNumberMedianRadius(imode)*DEXP(1.5_r8*(log(lifeCycleSigma(imode)))**2)
             wetNumberMedianDiameter(icol,ilev,imode) =  0.1e-6_r8 !Initialize to something..
             mixed_dry_rho = 1.e3_r8

             tracerCounter = 0
             do itrac = 1,getNumberOfBackgroundTracersInMode(imode)

                tracerCounter = tracerCounter + 1

                !which tracer is this?
                mm = getTracerIndex(imode,itrac,.false.)

                !radius of lower rh-bin for this tracer
                rr1=rdivr0(t_irh1,mm)

                !radius of upper rh-bin for this tracer
                rr2=rdivr0(t_irh2,mm)

                !linear interpolate dry ==> wet radius for this tracer
                wetrad_tmp(tracerCounter) = (((t_rh2-t_xrh)*rr1+(t_xrh-t_rh1)*rr2)/ &
                     (t_rh2-t_rh1))*rmeanvol

                !mixed density of dry particle
                dry_rhopart_tmp(tracerCounter) = getDryDensity(imode,itrac)

             end do

             !Find the average growth of this mode
             !(still not taking into account how much we have!!)
             if(TracerCounter > 0)then

                !Convert to diameter and take average (note: This is MASS median diameter)
                wetNumberMedianDiameter(icol,ilev,imode) = 2.0_r8 * SUM(wetrad_tmp(1:tracerCounter))/dble(tracerCounter)

                !Take average density
                mixed_dry_rho = SUM(dry_rhopart_tmp(1:tracerCounter))/dble(tracerCounter)

                !At this point the radius is in "mass mean" space
                volumeFractionAerosol = MIN(1.0_r8, ( 2.0_r8*rmeanVol / wetNumberMedianDiameter(icol,ilev,imode) )**3)

                !wet density
                wetrho(icol,ilev,imode) = mixed_dry_rho * volumeFractionAerosol + (1._r8-volumeFractionAerosol)*rhoh2o

                !convert back to number median diameter (wet)
                wetNumberMedianDiameter(icol,ilev,imode) = &
                     wetNumberMedianDiameter(icol,ilev,imode)*DEXP(-1.5_r8*(log(lifeCycleSigma(imode)))**2)
             endif

          end do     !modes

          !Same thing for the process modes
          do itrac = 1,numberOfProcessModeTracers

             mm = tracerInProcessMode(itrac)   !process mode tracer (physics space)

             !weighting by 1.5 is number median ==> volumetric mean
             !http://dust.ess.uci.edu/facts/psd/psd.pdf
             rmeanvol = processModeNumberMedianRadius(itrac)*DEXP(1.5_r8*(log(processModeSigma(itrac)))**2)

             !radius of lower rh-bin for this tracer
             rr1=rdivr0(t_irh1,mm)

             !radius of upper rh-bin for this tracer
             rr2=rdivr0(t_irh2,mm)

             !Note this is MASS median diameter
             wetNumberMedianDiameter_processmode(icol,ilev,itrac) = (((t_rh2-t_xrh)*rr1+(t_xrh-t_rh1)*rr2)/ &
                  (t_rh2-t_rh1))*rmeanvol*2.0_r8

             volumeFractionAerosol = MIN(1.0, (2.0_r8*rmeanVol/wetnumberMedianDiameter_processmode(icol,ilev,itrac))**3)

             wetrho_processmode(icol,ilev,itrac) = volumeFractionAerosol*rhopart(mm) &
                  + (1.0_r8 - volumeFractionAerosol)*rhoh2o

             !convert back to number median diameter (wet)
             wetNumberMedianDiameter_processMode(icol,ilev,itrac) = &
                  wetNumberMedianDiameter_processMode(icol,ilev,itrac)*DEXP(-1.5_r8*(log(processModeSigma(itrac)))**2)
          end do !process modes
       end do !horizontal points
    end do !layers

  end subroutine calcaersize_sub

end module aero_model
