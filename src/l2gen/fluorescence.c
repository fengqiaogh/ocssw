//---------------------------------------------------------------------
// fluorescence.c -  functions to support chlorophyll fluorescence.    
//                                                                     
// Algorithm Reference:                                                
// Behrenfeld, M.J., T.K. Westberry, E.S. Boss, R.T. O'Malley, D.A. 
// Siegel, J.D. Wiggert, B.A. Franz, C.R. McClain, G.C. Feldman, S.C. 
// Doney, J.K. Moore, G. Dall'Olmo, A. J. Milligan, I. Lima, and N. 
// Mahowald (2009). Satellite-detected fluorescence reveals global 
// physiology of ocean phytoplankton, Biogeosci., 6, 779-794.
//
// Implementation:  B. Franz, Oct 2009.                      
//---------------------------------------------------------------------

#include "l12_proto.h"

#define PARW1 400 /**< Minimum wavelength for PAR */
#define PARW2 700 /**< Maximum wavelength for PAR */
#define PARWN (PARW2-PARW1)+1 /**< Number of wavelengths included in PAR */

static float fqymin = 0.0;
static float fqymax = 0.3;
static float flhmin = 0.0;

/*---------------------------------------------------------------------*/
/* get_fsat - normalized fluorescence line height for each pixel       */
/*             fsat output is in radiance units (mW/cm^2/um/sr)        */
/*---------------------------------------------------------------------*/

/**
 * get_flh - computes normalized fluorescence line height for MODIS, MERIS, OLCI and OCI
 * computation based on nLw measurements (thus the normalized moniker)
 * @param l2rec
 * @param flh
 */
void get_fsat(l2str *l2rec, float flh[]) {
    static int32_t ib665, ib680, ib709;
    static int firstCall = 1;

    int32_t ip, ipb;
    float base;
    float nLw1;
    float nLw2;
    float nLw3;
    float Lf1;
    float Lf2;
    float Lf3;
    float bias = input->flh_offset;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        firstCall = 0;
        ib665 = windex(665., l1file->fwave, l1file->nbands);
        ib680 = windex(680., l1file->fwave, l1file->nbands);
        ib709 = windex(709., l1file->fwave, l1file->nbands);

        if (fabs(l1file->fwave[ib665] - 665) > 2.5){
            printf("No fluorescence algorithm available for this sensor.\n");
            exit(EXIT_FAILURE);
        }
    
        if (fabs(l1file->fwave[ib680] - 680 ) > 2.5){
            printf("No fluorescence algorithm available for this sensor.\n");
            exit(EXIT_FAILURE);

        }
        // special handling for MODIS
        if (fabs(l1file->fwave[ib709] - 709) > 5){
            ib709 = windex(748., l1file->fwave, l1file->nbands);
            if (fabs(l1file->fwave[ib709] - 748) > 5){
                printf("No fluorescence algorithm available for this sensor.\n");
                exit(EXIT_FAILURE);
            }   
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        flh[ip] = BAD_FLT;

        ipb = l1file->nbands * ip;
        /**
         * the nflh algorithm requires nLw values from 665, 680 and 709nm
         */
        nLw1 = l2rec->nLw[ipb + ib665];
        nLw2 = l2rec->nLw[ipb + ib680];
        nLw3 = l2rec->nLw[ipb + ib709];

        Lf1 = l1file->fwave[ib665];
        Lf2 = l1file->fwave[ib680];
        Lf3 = l1file->fwave[ib709];

        /**
         * if the pixel is already masked, or any of the input nLw values are
         * less than -0.01, set the PRODFAIL flag and move along
         */
        if (l1rec->mask[ip] || nLw1 < -0.01 || nLw2 < -0.01 || nLw3 < -0.01) {
            l1rec->flags[ip] |= PRODFAIL;
            continue;

        } else {

            /**
             * fsat (Behrenfeld et al.[2009] equation A2)
             */
            base = nLw3 + (nLw1 - nLw3) * ((Lf3- Lf2) / (Lf3 - Lf1));
            flh[ip] = nLw2 - base;

            /**
             * bias correction as per Behrenfeld et al.
             */
            flh[ip] -= bias;
        
        }
    }
}

/*---------------------------------------------------------------------*/
/* get_flh - fluorescence line height for each pixel in rec            */
/*             flh output is in radiance units (mW/cm^2/um/sr)         */                                               
/*---------------------------------------------------------------------*/

/**
 * get_flh - computes fluorescence line height for MODIS, MERIS, OLCI and OCI
 * computation based on Lw measurements
 * 
 * @param l2rec
 * @param flh
 */
void get_flh(l2str *l2rec, float flh[]) {
    static int32_t ib665, ib680, ib709;
    static int firstCall = 1;

    int32_t ip, ipb;
    float Lw1;
    float Lw2;
    float Lw3;
    float Lf1;
    float Lf2;
    float Lf3;
    float base;
    float bias = input->flh_offset;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        firstCall = 0;
        ib665 = windex(665., l1file->fwave, l1file->nbands);
        ib680 = windex(680., l1file->fwave, l1file->nbands);
        ib709 = windex(709., l1file->fwave, l1file->nbands);
        
        if (fabs(l1file->fwave[ib665] - 665) > 2.5){
            printf("No fluorescence algorithm available for this sensor.\n");
            exit(EXIT_FAILURE);
        }
    
        if (fabs(l1file->fwave[ib680] - 680 ) > 2.5){
            printf("No fluorescence algorithm available for this sensor.\n");
            exit(EXIT_FAILURE);

        }
        // special handling for MODIS
        if (fabs(l1file->fwave[ib709] - 709) > 5){
            ib709 = windex(748., l1file->fwave, l1file->nbands);
            if (fabs(l1file->fwave[ib709] - 748) > 5){
                printf("No fluorescence algorithm available for this sensor.\n");
                exit(EXIT_FAILURE);
            }   
        
        }   
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        flh[ip] = BAD_FLT;

        ipb = l1file->nbands * ip;
        /**
         * the flh algorithm requires nLw values from 667, 678 and 748nm
         */
        Lw1 = l2rec->Lw[ipb + ib665];
        Lw2 = l2rec->Lw[ipb + ib680];
        Lw3 = l2rec->Lw[ipb + ib709];

        Lf1 = l1file->fwave[ib665];
        Lf2 = l1file->fwave[ib680];
        Lf3 = l1file->fwave[ib709];

        // scattering-angle-dependent spectral difference in aerosol model phase 
        // functions are introducing geometric artifacts
        // replace aerosol difference between red bands with simple model
        // rho_a[678] = rho_a[667]*(678/667)^-0.65

        //La  = l2rec->La[ipb+ib678]/l2rec->t_sen[ipb+ib678];   // original aerosol
        //Lap = 0.99*l2rec->La[ipb+ib667]/l2rec->Fo[ib667]*l2rec->Fo[ib678]/l2rec->t_sen[ipb+ib678];
        //Lw2 = Lw2 + La - Lap;
        /**
         * if the pixel is already masked, or any of the input Lw values are
         * less than -0.01, set the PRODFAIL flag and move along
         */
        if (l1rec->mask[ip] || Lw1 < -0.01 || Lw2 < -0.01 || Lw3 < -0.01) {
            l1rec->flags[ip] |= PRODFAIL;
            continue;

        } else {

            /**
             * Lw,f (Behrenfeld et al.[2009] equation A2 and A3)
             */

            base = Lw3 + (Lw1 - Lw3) * ( (Lf3- Lf2) / (Lf3 - (Lf1)));
            flh[ip] = Lw2 - base;

            /**
             *  bias correction as per Behrenfeld et al. [2009]
             */
            flh[ip] -= bias;

            if (flh[ip] < flhmin) {
                flh[ip] = 0.0;
            }
        }
    }
}

/*---------------------------------------------------------------------*/
/* fqy - fluorescence quantum yield                                    */
/*---------------------------------------------------------------------*/

/**
 * get_fqy - function to return fluorescence quantum yield for the MODIS, MERIS, OLCI and OCI instruments
 * internally calls the calc_fsat and get_ipar functions
 *
 * @param l2rec
 * @param fqy
 */
void get_fqy(l2str *l2rec, float fqy[]) {
    static float badval = BAD_FLT;
    static int firstCall = 1;
    static float *ipar;
    static float *fsat;

    int32_t ip;

    l1str *l1rec = l2rec->l1rec;

    if (firstCall) {

        firstCall = 0;

        if ((ipar = (float *) calloc(l1rec->npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d: Unable to allocate space for ipar.\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if ((fsat = (float *) calloc(l1rec->npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d: Unable to allocate space for fsat.\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    // Compute ipar and fsat at all pixels

    get_ipar(l2rec, ipar);
    get_fsat(l2rec, fsat);

    // Compute fluorescence quantum yield for each pixel

    for (ip = 0; ip < l1rec->npix; ip++) {

        fqy[ip] = badval;

        if (!l1rec->mask[ip] && l2rec->chl[ip] > 0.0 && fsat[ip] >= 0.0) {

            // simple form from Behrenfeld, A12
            fqy[ip] = 0.01 * fsat[ip] / (0.0302 * l2rec->chl[ip]) * ipar[ip]
                    * 1e6 / 1590;

            if (!isfinite(fqy[ip])) {
                fqy[ip] = badval;
                l1rec->flags[ip] |= PRODFAIL;
            } else if (fqy[ip] < fqymin || fqy[ip] > fqymax) {
                l1rec->flags[ip] |= PRODWARN;
            }

        } else {
            fqy[ip] = badval;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }
}

// ***************   old stuff *******************

/* we want flhq in same units as arp (Ein/m^2/s) */
/* flh is in mW/cm^2/um/sr / 100 = W/m^2/nm/sr    */
/* hc = 1986.5e-19 J*nm/photon, Lambda 676.7 nm   */
/* radiance->iradiance = 4*pi                     */
/* line height -> area = 43.38 nm                 */
/* 6.023e23 per mole Quanta                       */

/* Notes from modcol (anly8dbl.f90)
 ! 4*$PI - convert radiance to scalar irradiance
 ! energy of 1 photon = hc/lambda
 ! h = 6.6261e-34 [Js]; c = 299,792,458 [m/s] (*10^9 for nm/s)
 ! energy of one photon = (1986.5*10^-19)/Lambda [W s]
 ! where lambda is [nm]
 !
 ! FLHQ units are quanta/m^2/s
 ! to convert FLH into total fluorescence under the fluorescence
 ! curve :
 ! The fluorescence curve is described as a gaussian curve with
 ! center at 683 nm and half maximum emission of 25 nm (from
 ! Collins et al. 1985
 ! If L683 (1 nm width) = 1 then the area under the curve= 26.61
 ! Band14 center= 677 bandwidth= 10 Fluorescence signal read= 8.61
 ! Band13 center= 665 bandwidth= 10 Fluorescence signal read= 2.86
 ! Band15 center= 747 bandwidth= 10 Fluorescence signal read= 1.55e-6
 ! All band readings correspond to a 10 nm bandwidth signal.  For this
 ! reason the signals must be divided by 10 to get the reading per nm
 !
 ! Because band 13 is affected by the fluorescence signal then there
 ! is a contribution of this signal equal to 0.2478 to the baseline
 ! correction.
 !
 ! The conversion factor of FLH into area under the curve can be
 ! computed: area under the curve / (Band14 read - baseline contrib)
 ! factor = 26.61 /(0.861 - 0.2478) = 43.38
 ! rescale: 43.38 (nm) -> 0.04338 (um)
 !
 ! FLHQ = FLH * Lambda/(hc) * (radiance->irradiance)*(FLH->area)
 ! FLHQ is in quanta m-2 s-1 and ARP is in Einstein m-2 s-1
 ! 1 Einstein = 1 mol quanta = 6.023e23 quanta
 */
