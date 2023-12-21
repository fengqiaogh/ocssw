  real, dimension(:,:), allocatable               ::  brdf650
  character(len=255), dimension(:), allocatable   ::  aero_sites
  integer, dimension(:), allocatable              ::  aero_zones
  integer, dimension(:), allocatable              ::  aero_types
  integer, dimension(:), allocatable              ::  aero_elev
  real, dimension(:,:), allocatable ::  aero_sr412, aero_sr470, aero_sr650, aero_bgaod

  integer,dimension(3600,1800)  ::  terrain_flag_new
  real, dimension(3600,1800)    ::  terrain_flag, sfc_elev_std
  real, dimension(360,180)    ::  bg_aod! 9 January 2018 JLee

  integer  :: LERstart(2), LERedge(2), dateline, doy
  integer  :: LERstart6(2), LERedge6(2), dateline6
  real, dimension(:,:,:,:), allocatable  :: coefs650_fwd, coefs470_fwd, coefs412_fwd
  real, dimension(:,:,:,:), allocatable  :: coefs650_all, coefs470_all, coefs412_all

  real, dimension(:,:), allocatable  ::  gref412_all, gref412_fwd
  real, dimension(:,:), allocatable  ::  gref470_all, gref470_fwd
  real, dimension(:,:), allocatable  ::  gref650_all, gref650_fwd
  real, dimension(:,:), allocatable  ::  gref865_all

  real, dimension(:,:), allocatable  ::  vgref412_all
  real, dimension(:,:), allocatable  ::  vgref488_all
  real, dimension(:,:), allocatable  ::  vgref670_all

  real, dimension(:,:,:), allocatable  ::  swir_coeffs412, swir_coeffs470
  real, dimension(:,:), allocatable    ::  swir_stderr412, swir_stderr470
  real, dimension(:,:), allocatable    ::  swir_min, swir_max

  common /db/   brdf650, aero_sites, aero_zones, aero_types, aero_elev, &
  aero_sr412, aero_sr470, aero_sr650, aero_bgaod, terrain_flag_new, &
  terrain_flag, sfc_elev_std, bg_aod, LERstart, LERedge, dateline, doy, &
  LERstart6, LERedge6, dateline6, coefs650_fwd, coefs470_fwd, coefs412_fwd, &
  coefs650_all, coefs470_all, coefs412_all, gref412_all, gref412_fwd, &
  gref470_all, gref470_fwd, gref650_all, gref650_fwd, gref865_all, &
  vgref412_all, vgref488_all, vgref670_all, swir_coeffs412, swir_coeffs470, &
  swir_stderr412, swir_stderr470, swir_min, swir_max

