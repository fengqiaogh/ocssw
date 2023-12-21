module deepblue_config
    !
    !    deepblue_config.f95
    !
    !    Description:       Configuration module for Deep Blue aerosol retrievals.
    !
    !    Modules:              None
    !
    !    Author:                Corey Bettenhausen
    !                              Science Systems and Applications, Inc.
    !                              NASA Goddard Space Flight Center
    !                              corey_bettenhausen@ssaihq.com
    !
    !    Last Changed:
    !
    !    History:             2012-08    Original version
    !
    !------------------------------------------------------------------------------
    implicit none

    ! Restrict access to module components.
    private

    public  ::  load_viirs_config
    public  ::  viirs_config_type
  
    type    ::  viirs_config_type
        integer                 ::  year
        integer                 ::  month
        integer                 ::  day
        integer                 ::  hour
        integer                 ::  minute
        integer                 ::  second
        character(len=255)      ::  lut_nc4
        character(len=255)      ::  aerosol_land_file
        character(len=255)      ::  aerosol_dust_file
        character(len=255)      ::  aerosol_ocean_dust_file
        character(len=255)      ::  aerosol_ocean_fine_file
        character(len=255)      ::  aerosol_ocean_mari_file
        character(len=255)      ::  aerosol_ocean_mix_file
        character(len=255)      ::  bathymetry_lut_file
        character(len=255)      ::  chl_lut_file
        character(len=255)      ::  ler_lut_file
        character(len=255)      ::  landcover_file
        character(len=255)      ::  surfpressure_file
        character(len=255)      ::  geozone_file
        character(len=255)      ::  seasonal_deserts_file
        character(len=255)      ::  brdfbase_file
        character(len=255)      ::  modis_surfdb_file
        character(len=255)      ::  viirs_surfdb_file
        character(len=255)      ::  surfcoeffs_file
        character(len=255)      ::  swir_vis_surfcoeffs_file
        character(len=255)      ::  veg_landcover_file
        character(len=255)      ::  veg_sfc21_file
    
        character(len=255)      ::  rayl412_file
        character(len=255)      ::  rayl488_file
        character(len=255)      ::  rayl670_file
        character(len=255)      ::  xcal412_file
        character(len=255)      ::  xcal488_file
        character(len=255)      ::  xcal670_file
    
        character(len=255)      ::  gmtco_file
        character(len=255)      ::  iicmo_file
        character(len=255)      ::  svm01_file
        character(len=255)      ::  svm02_file
        character(len=255)      ::  svm03_file
        character(len=255)      ::  svm04_file
        character(len=255)      ::  svm05_file
        character(len=255)      ::  svm07_file
        character(len=255)      ::  svm08_file
        character(len=255)      ::  svm09_file
        character(len=255)      ::  svm10_file
        character(len=255)      ::  svm11_file
        character(len=255)      ::  svm14_file
        character(len=255)      ::  svm15_file
        character(len=255)      ::  svm16_file
        character(len=255)      ::  gdas_file1
        character(len=255)      ::  gdas_file2
    
        character(len=255)      ::  l1b_m
        character(len=255)      ::  geo_m
    
        character(len=255)      ::  output_l2

    
    end type viirs_config_type

contains

    !
    !    load_viirs_config()
    !
    !        Loads configuration file, config_file, into viirs_config object.
    !
    !-------------------------------------------------------------------
    type(viirs_config_type) function load_viirs_config(config_file, status) result(vcfg)
        implicit none

        character(len=255), intent(in)    :: config_file
        integer, intent(inout)          :: status

        character(255) :: c_line, var, val
        integer        :: eq_index

        !   -- initialize our config variables.
        vcfg%year     = 1900
        vcfg%month    = 1
        vcfg%day      = 1
        vcfg%hour     = 0
        vcfg%minute   = 0
        vcfg%second   = 0
    
        vcfg%lut_nc4                  = ''
        vcfg%aerosol_land_file        = ''
        vcfg%aerosol_dust_file        = ''
        vcfg%aerosol_ocean_dust_file  = ''
        vcfg%aerosol_ocean_fine_file  = ''
        vcfg%aerosol_ocean_mari_file  = ''
        vcfg%aerosol_ocean_mix_file   = ''
        vcfg%bathymetry_lut_file      = ''
        vcfg%chl_lut_file             = ''
        vcfg%ler_lut_file             = ''
        vcfg%landcover_file           = ''
        vcfg%surfpressure_file        = ''
        vcfg%geozone_file             = ''
        vcfg%seasonal_deserts_file    = ''
        vcfg%brdfbase_file            = ''
        vcfg%modis_surfdb_file        = ''
        vcfg%viirs_surfdb_file        = ''
        vcfg%surfcoeffs_file          = ''
        vcfg%swir_vis_surfcoeffs_file = ''
        vcfg%veg_landcover_file       = ''
        vcfg%veg_sfc21_file           = ''
        
        vcfg%gmtco_file               = ''
        vcfg%svm01_file               = ''
        vcfg%svm02_file               = ''
        vcfg%svm03_file               = ''
        vcfg%svm04_file               = ''
        vcfg%svm05_file               = ''
        vcfg%svm07_file               = ''
        vcfg%svm08_file               = ''
        vcfg%svm09_file               = ''
        vcfg%svm10_file               = ''
        vcfg%svm11_file               = ''
        vcfg%svm14_file               = ''
        vcfg%svm15_file               = ''
        vcfg%svm16_file               = ''
        vcfg%iicmo_file               = ''

        vcfg%l1b_m                    = ''
        vcfg%geo_m                    = ''

        vcfg%gdas_file1               = ''
        vcfg%gdas_file2               = ''
    
        vcfg%rayl412_file             = ''
        vcfg%rayl488_file             = ''
        vcfg%rayl670_file             = ''
        vcfg%xcal412_file             = ''
        vcfg%xcal488_file             = ''
        vcfg%xcal670_file             = ''
        vcfg%output_l2                = ''

        open(100, file=trim(config_file), status='OLD', FORM='FORMATTED',  ACTION='READ', iostat=status)
        if (status /= 0) then
            print *, 'ERROR: failed to open configuration file: ', status
            return
        end if

        do
            ! Read a line from the config file.  Die on error.
            read(100, fmt='(A)', iostat=status) c_line
            if (status /= 0) then
                if (status < 0) then
                    status = 0
                    return
                else if (status > 0) then
                    print *, "ERROR: Failed to read configuration file: ", status
                    status = -1
                    return
                end if
            end if

            ! Clean up whitespace.
            c_line =rm_whitespace(c_line)

            ! Detect comment lines or empty lines and skip.
            if (index(c_line, '!') /= 0) then
                cycle
            end if

            if (len_trim(c_line) == 0) then
                cycle
            end if

            ! Begin parsing for variable and value.
            ! If no '=' is detected, we have a funky line, cycle.
            eq_index = index(c_line,'=')
            if (eq_index /= 0) then
                var = c_line(1:eq_index-1)
                val = c_line(eq_index+1:len(c_line))

                select case (trim(var))
                    case ('lut_nc4_db')
                        vcfg%lut_nc4 = trim(val)
                    case ('year')
                        read(val, fmt='(I4)') vcfg%year
                    case ('month')
                        read(val, fmt='(I2)') vcfg%month
                    case ('day')
                        read(val, fmt='(I2)') vcfg%day
                    case ('hour')
                        read(val, fmt='(I2)') vcfg%hour
                    case ('minute')
                        read(val, fmt='(I2)') vcfg%minute
                    case ('second')
                        read(val, fmt='(I2)') vcfg%second
                    case ('ofile')
                        vcfg%output_l2 = trim(val)
                    case ('lut_aero_land_fine')
                        vcfg%aerosol_land_file = trim(val)
                    case ('lut_aero_land_dust')
                        vcfg%aerosol_dust_file = trim(val)
                    case ('lut_aero_ocean_dust')
                        vcfg%aerosol_ocean_dust_file = trim(val)
                    case ('lut_aero_ocean_fine')
                        vcfg%aerosol_ocean_fine_file = trim(val)
                    case ('lut_aero_ocean_mari')
                        vcfg%aerosol_ocean_mari_file = trim(val)
                    case ('lut_aero_ocean_mix')
                        vcfg%aerosol_ocean_mix_file = trim(val)
                    case ('lut_bathymetry')
                        vcfg%bathymetry_lut_file = trim(val)
                    case ('lut_chl')
                        vcfg%chl_lut_file = trim(val)
                    case ('lut_ler_table')
                        vcfg%ler_lut_file = trim(val)
                    case ('lut_landcover')
                        vcfg%landcover_file = trim(val)
                    case ('lut_geozone')
                        vcfg%geozone_file = trim(val)
                    case ('lut_seasonal_deserts')
                        vcfg%seasonal_deserts_file = trim(val)
                    case ('lut_brdf')
                        vcfg%brdfbase_file = trim(val)
                    case ('lut_modis_surfdb')
                        vcfg%modis_surfdb_file = trim(val)
                    case ('lut_viirs_surfdb')
                        vcfg%viirs_surfdb_file = trim(val)
                    case ('lut_surfcoeff')
                        vcfg%surfcoeffs_file = trim(val)
                    case ('lut_swir')
                        vcfg%swir_vis_surfcoeffs_file = trim(val)
                    case ('lut_veg_landcover')
                        vcfg%veg_landcover_file = trim(val)
                    case ('lut_veg_21sfc')
                        vcfg%veg_sfc21_file = trim(val)

                    case ('lut_rayl_412')
                        vcfg%rayl412_file = trim(val)
                    case ('lut_rayl_488')
                        vcfg%rayl488_file = trim(val)
                    case ('lut_rayl_670')
                        vcfg%rayl670_file = trim(val)
                    case ('lut_viirs_xcal_412')
                        vcfg%xcal412_file = trim(val)
                    case ('lut_viirs_xcal_488')
                        vcfg%xcal488_file = trim(val)
                    case ('lut_viirs_xcal_670')
                        vcfg%xcal670_file = trim(val)

                    case ('gmtco')
                        vcfg%gmtco_file = trim(val)
                    case ('iicmo')
                        vcfg%iicmo_file = trim(val)
                    case ('svm01')
                        vcfg%svm01_file = trim(val)
                    case ('svm02')
                        vcfg%svm02_file = trim(val)
                    case ('svm03')
                        vcfg%svm03_file = trim(val)
                    case ('svm04')
                        vcfg%svm04_file = trim(val)
                    case ('svm05')
                        vcfg%svm05_file = trim(val)
                    case ('svm07')
                        vcfg%svm07_file = trim(val)
                    case ('svm08')
                        vcfg%svm08_file = trim(val)
                    case ('svm09')
                        vcfg%svm09_file = trim(val)
                    case ('svm10')
                        vcfg%svm10_file = trim(val)
                    case ('svm11')
                        vcfg%svm11_file = trim(val)
                    case ('svm14')
                        vcfg%svm14_file = trim(val)
                    case ('svm15')
                        vcfg%svm15_file = trim(val)
                    case ('svm16')
                        vcfg%svm16_file = trim(val)
                    case ('gdas1')
                        vcfg%gdas_file1 = trim(val)
                    case ('gdas2')
                        vcfg%gdas_file2 = trim(val)
                    case ('ifile')
                        vcfg%l1b_m = trim(val)
                    case ('geofile')
                        vcfg%geo_m = trim(val)
                    case default
                        cycle
                end select
            else
                cycle
            end if


        end do

        return

    end function load_viirs_config

    !
    ! rm_whitespace()
    !
    !      Removes whitespace from a specified string.
    !
    !      Inputs:              input_string        string to be cleared of whitespace
    !      Outputs:            Returns character array cleared of all whitespace.
    !   Side effects:   None
    !
    !----------------------------------------------------------------------
    function rm_whitespace(input_string)
        implicit none

        character(len=255), intent(in) :: input_string

        character(len=len(input_string)) :: rm_whitespace
        integer                             :: space_index, end_index

        rm_whitespace = input_string

        ! Empty string?
        if (len_trim(rm_whitespace) == 0) then
            return
        end if

        ! Remove frontal spaces.
        do
            space_index = index(rm_whitespace,' ')
            if ( space_index == 1) then
                rm_whitespace = rm_whitespace(2:len(rm_whitespace))
            else
                exit
            end if
        end do

        ! Remove any internal spaces
        do
            space_index = index(rm_whitespace,' ')
            end_index = len_trim(rm_whitespace)
            if (space_index > 0 .AND. space_index < end_index) then
                rm_whitespace = rm_whitespace(1:space_index-1) //           &
                    &    rm_whitespace(space_index+1:end_index)
            else
                exit
            end if
        end do

        return
    end function rm_whitespace

end module deepblue_config
