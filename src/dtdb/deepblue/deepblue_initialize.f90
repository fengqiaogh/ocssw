!-----------------------------------------------------------------------
! deepblue_initialize()
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:  none
!
!-----------------------------------------------------------------------
subroutine deepblue_initialize ( config_file, nc4_file, lines, pixels, lat, lon,  &
    year, month, day )

    use landcover

    use calendars, only:            &
        gdatetime,                  &
        doy_from_gregorian,         &
        season_from_doy

    use deepblue_config, only:      &
        viirs_config_type,          &
        load_viirs_config
                            
    use modis_surface, only:        &
        load_terrainflg_tables,     &
        load_seasonal_desert,       &
        load_brdf,                  &
        set_limits,                 &
        set_limits6,                &
        load_hdfLER,                &
        get_LER412,                 &
        get_LER470,                 &
        get_LER650,                 &
        latlon_to_index_ler,        &
        get_geographic_zone,        &
        get_sfc_elev_std,           &
        load_swir_coeffs
  
    use viirs_aerosol_luts, only:   &
        load_viirs_aerosol_luts
  
    use viirs_ler_luts

    implicit none
  include 'db.inc'
  include 'newaottbl90.inc'
  include 'sfc21tbl90.inc'

    character(len=255) ::  config_file, nc4_file
    INTEGER, PARAMETER ::  VLINES=16*203, VPIXELS=3200
    real, dimension(VPIXELS,VLINES):: lat, lon
  
    integer, parameter              ::  R_DBL = kind(1.0d0)
    type (viirs_config_type)        ::  config
    type (gdatetime)                ::  gdt1
 
    real, dimension(:,:), allocatable     ::  minmax, minmax_ler
    integer, dimension(:,:), allocatable  ::  lcld_mask, tmp_mask!,ocld_mask
    integer, dimension(:,:), allocatable  ::  lskip_mask
    integer, dimension(:,:), allocatable  ::  oskip_mask
    integer, dimension(:,:), allocatable  ::  smoke_mask
    integer, dimension(:,:), allocatable  ::  smoke_ae_mask
    integer, dimension(:,:), allocatable  ::  pyrocb_mask
    integer, dimension(:,:), allocatable  ::  high_alt_smoke_mask
    integer, dimension(:,:), allocatable  ::  snow_mask, snow_mask2
    real, dimension(:,:), allocatable     ::  sr650
    real, dimension(:,:), allocatable     ::  sfcstd
    integer, dimension(:,:), allocatable  ::  gzflg
    real, dimension(:,:), allocatable     ::  wv
    real, dimension(:,:), allocatable     ::  oz
    real, dimension(:,:), allocatable     ::  windsp
    real, dimension(:,:), allocatable     ::  winddir
    integer, dimension(:,:), allocatable  ::  lc    ! land cover
    integer, dimension(:,:), allocatable  ::  bathy ! bathymetry
    real, dimension(:,:), allocatable     ::  chl   ! log-10 Chl climatology

    integer, dimension(:,:), allocatable  ::  n_total_pixels
    integer                               ::  cell_resolution
    integer, dimension(2)                 ::  cell_dims
  
    real                                  ::  ndsi
    real                                  ::  ndvi_lower, ndvi_upper
    real                                  ::  dd
    integer                               ::  i,j
    integer                               ::  i1, i2, j1, j2
    integer                               ::  status
    integer                               ::  lines, pixels
    integer                               ::  year, month, day, doy
    integer                               ::  season
    integer, dimension(2)                 ::  dims2
    real, dimension(:,:), allocatable     ::  to_iof
    real, parameter    ::  d2r = 3.14159/180.0   ! convert degrees to radians
    character(len=255) ::  val_fname
    character(len=255) ::  env, out_nc4, dataroot, varroot
    integer            ::  nt, nr

    common  /xday/ doy

    status = -1

    config = load_viirs_config(config_file, status)
    if (status /= 0) then
        print *, "ERROR: Failed to read in VIIRS configuration file: ", status
        stop
    end if

    out_nc4 = trim(config%lut_nc4)
    if (out_nc4 == "") then
        out_nc4 = nc4_file
    endif

    call get_environment_variable ("OCDATAROOT", dataroot)
    env = trim('$OCDATAROOT')
    nt = len_trim(env)
    nr = len_trim(dataroot)
    do
       i = INDEX(out_nc4, env(:nt)) ; IF (i == 0) EXIT
       out_nc4 = out_nc4(:i-1) // dataroot(:nr) // out_nc4(i+nt:)
    end do

    call get_environment_variable ("OCVARROOT", varroot)
    env = trim('$OCVARROOT')
    nt = len_trim(env)
    nr = len_trim(varroot)
    do
       i = INDEX(out_nc4, env(:nt)) ; IF (i == 0) EXIT
       out_nc4 = out_nc4(:i-1) // varroot(:nr) // out_nc4(i+nt:)
    end do

    status = load_viirs_aerosol_luts(out_nc4)
    if (status /= 0) then
        print *, "ERROR: Failed to load VIIRS aerosol LUTS: ", status
        stop
    end if
!    print *, 'done.'

    ! -- LER luts
    call load_viirs_ler_luts(out_nc4, status)
    if (status /= 0) then
        print *, "ERROR: Failed to load VIIRS LER LUTS: ", status
        stop
    end if

    ! -- land cover
    status = load_landcover(out_nc4)
    if (status /= 0) then
        print *, "ERROR: Unable to load land cover input file: ", status
        stop
    end if
        
    ! -- geo zones
    gdt1 = gdatetime(year, month, day, 0, 0, 0, 0, 0)
    doy = doy_from_gregorian(gdt1)

    season = season_from_doy(2005,doy)
    print *, "season: ", season
    status = load_terrainflg_tables(out_nc4,season)
    if (status /= 0) then
        print *, "ERROR: Unable to load geozone table: ", status
        stop
    end if
     
    ! -- seasonal desert data
    status = load_seasonal_desert(out_nc4)
    if (status /= 0) then
        print *, "ERROR: Unable to load seasonal deserts file: ", status
    end if

    ! -- 670nm BRDF data
    status = load_brdf(out_nc4)
    if (status /= 0) then
        print *, "ERROR: Unable to load BRDF base input file: ", status
        stop
    end if

    ! -- surface database and BRDF coefficients
    dims2 = (/pixels, lines/)
    status = set_limits(dims2, lat, lon)
    if (status /= 0) then
        print *, "ERROR: Failure to set limits on surface coefficients table: ", status
        stop
    end if

    status = load_hdfLER(out_nc4, season)
    if (status /= 0) then
        print *, "ERROR: Unable to load surface BRDF coefficients: ", status
        stop
    end if

    ! -- swir vs. vis surface coeffs
    dims2 = (/pixels, lines/)
    status = set_limits6(dims2, lat, lon)
    if (status /= 0) then
        print *, "ERROR: Failure to set limits on 2.2 um surface coefficients table: ", status
        stop
    end if

    status = load_swir_coeffs(out_nc4, season)
    if (status /= 0) then
        print *, "ERROR: Unable to load swir vs. vis surface coeffs file: ", status
    end if

    ! -- vegetated retrieval landcover
    call get_lut_igbp_land_cover(out_nc4, status)
    if (status /= 0) then
        print *, "ERROR: Failed to read in landcover input for vegetated retrieval: ", status
        stop
    end if
  
    call get_lut_211sfc(out_nc4, status)
    if (status /= 0) then
        print *, "ERROR: Failed to read in landcover input for vegetated 2.1um sfc table: ", status
        stop
    end if
  
    allocate(wv(pixels,lines), oz(pixels,lines), windsp(pixels,lines), &
        &        winddir(pixels,lines), lc(pixels,lines), stat=status)
    if (status /= 0) then
        print *, "ERROR: Failed to allocate water vapor, ozone arrays, and wind speed arrays: ", status
        stop
    end if

    return
end
  
!-----------------------------------------------------------------------------------------
! -- deepblue_cleanup
!-----------------------------------------------------------------------------------------
subroutine deepblue_cleanup()

    use landcover

    use modis_surface, only:        &
        unload_brdf

    use viirs_aerosol_luts, only:   &
        unload_viirs_aerosol_luts

    integer                               ::  status
  
    call unload_landcover(status)
    if (status /= 0) then
        print *, "ERROR: Unable to unload landcover data. Continuing: ", status
    end if
  
    call unload_brdf(status)
    if (status /= 0) then
        print *, "ERROR: Unable to unload BRDF input file. Continuing: ", status
    end if
  
    call unload_viirs_aerosol_luts(status)
    if (status /= 0) then
        print *, "ERROR: Unable to unload VIIRS aerosol LUTS. Continuing: ", status
    end if

    return
end

