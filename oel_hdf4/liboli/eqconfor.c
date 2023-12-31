/*******************************************************************************
NAME                            EQUIDISTANT CONIC 

PURPOSE:	Transforms input longitude and latitude to Easting and Northing
		for the Equidistant Conic projection.  The longitude and
		latitude must be in radians.  The Easting and Northing values
		will be returned in meters.

ALGORITHM REFERENCES

1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
    Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
    State Government Printing Office, Washington D.C., 1987.

2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
    U.S. Geological Survey Professional Paper 1453 , United State Government
    Printing Office, Washington D.C., 1989.
*******************************************************************************/
#include "oli_cproj.h"
#include "oli_local.h"

/* Variables common to all subroutines in this code file
  -----------------------------------------------------*/
static double r_major;		/* major axis 				*/
static double r_minor;		/* minor axis 				*/
static double lon_center;	/* Center longitude (projection center) */
static double e0,e1,e2,e3;	/* eccentricity constants		*/
static double e,es;		/* eccentricity constants		*/
static double ml0;		/* small value m			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/
static double ns;
static double g;
static double rh;


/* Initialize the Equidistant Conic projection
  ------------------------------------------*/
long eqconforint
(
    double r_maj,			/* major axis			*/
    double r_min,			/* minor axis			*/
    double lat1,			/* latitude of standard parallel*/
    double lat2,			/* latitude of standard parallel*/
    double center_lon,		/* center longitude		*/
    double center_lat,		/* center latitude		*/
    double false_east,		/* x offset in meters		*/
    double false_north,		/* y offset in meters		*/
    long mode			/* which format is present A B	*/
)
{
double temp;			/* temporary variable		*/
double sinphi,cosphi;		/* sin and cos values		*/
double ms1,ms2;
double ml1,ml2;

/* Place parameters in static storage for common use
  -------------------------------------------------*/
r_major = r_maj;
r_minor = r_min;
lon_center = center_lon;
false_northing = false_north;
false_easting = false_east;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = gctp_calc_e0(es);
e1 = gctp_calc_e1(es);
e2 = gctp_calc_e2(es);
e3 = gctp_calc_e3(es);

sincos(lat1,&sinphi,&cosphi);
ms1 = gctp_calc_small_radius(e,sinphi,cosphi);
ml1 = gctp_calc_dist_from_equator(e0, e1, e2, e3, lat1);

/* format B
---------*/
if (mode != 0)
   {
   if (fabs(lat1 + lat2) < EPSLN)
      {
      GCTP_PRINT_ERROR("Standard Parallels on opposite sides of equator");
      return(81);
      }
   sincos(lat2,&sinphi,&cosphi);
   ms2 = gctp_calc_small_radius(e,sinphi,cosphi);
   ml2 = gctp_calc_dist_from_equator(e0, e1, e2, e3, lat2);
   if (fabs(lat1 - lat2) >= EPSLN)
      ns = (ms1 - ms2) / (ml2 - ml1);
   else
      ns = sinphi;
   }
else
   ns = sinphi;
g = ml1 + ms1/ns;
ml0 = gctp_calc_dist_from_equator(e0, e1, e2, e3, center_lat);
rh = r_major * (g - ml0);
   

/* Report parameters to the user
  -----------------------------*/
if (mode != 0)
   {
   gctp_print_title("EQUIDISTANT CONIC"); 
   gctp_print_radius2(r_major, r_minor);
   gctp_print_stanparl(lat1,lat2);
   gctp_print_cenlonmer(lon_center);
   gctp_print_origin(center_lat);
   gctp_print_offsetp(false_easting,false_northing);
   }
else 
   {
   gctp_print_title("EQUIDISTANT CONIC"); 
   gctp_print_radius2(r_major, r_minor);
   gctp_print_stparl1(lat1);
   gctp_print_cenlonmer(lon_center);
   gctp_print_origin(center_lat);
   gctp_print_offsetp(false_easting,false_northing);
   }

return(OK);
}


/* Equidistant Conic forward equations--mapping lat,long to x,y
  -----------------------------------------------------------*/
long eqconfor
(
    double lon,			/* (I) Longitude 		*/
    double lat,			/* (I) Latitude 		*/
    double *x,			/* (O) X projection coordinate 	*/
    double *y			/* (O) Y projection coordinate 	*/
)
{
double ml;
double theta;
double rh1;

/* Forward equations
  -----------------*/
ml = gctp_calc_dist_from_equator(e0, e1, e2, e3, lat);
rh1 = r_major * (g - ml);
theta = ns * adjust_lon(lon - lon_center);
*x = false_easting  + rh1 * sin(theta);
*y = false_northing + rh - rh1 * cos(theta);

return(OK);
}
