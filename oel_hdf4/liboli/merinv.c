/*******************************************************************************
NAME                            MERCATOR

PURPOSE:	Transforms input Easting and Northing to longitude and
		latitude for the Mercator projection.  The
		Easting and Northing must be in meters.  The longitude
		and latitude values will be returned in radians.

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
static double lat_origin;	/* center latitude			*/
static double e,es;		/* eccentricity constants		*/
static double m1;		/* small value m			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/


/* Initialize the Mercator projection
  -----------------------------------*/
long merinvint
(
    double r_maj,			/* major axis			*/
    double r_min,			/* minor axis			*/
    double center_lon,		/* center longitude		*/
    double center_lat,		/* center latitude		*/
    double false_east,		/* x offset in meters		*/
    double false_north		/* y offset in meters		*/
)
{
double temp;			/* temporary variable		*/

/* Place parameters in static storage for common use
  -------------------------------------------------*/
r_major = r_maj;
r_minor = r_min;
lon_center = center_lon;
lat_origin = center_lat;
false_northing = false_north;
false_easting = false_east;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
m1 = cos(center_lat)/(sqrt(1.0 - es * sin(center_lat) * sin(center_lat)));

/* Report parameters to the user
  -----------------------------*/
gctp_print_title("MERCATOR"); 
gctp_print_radius2(r_major, r_minor);
gctp_print_cenlonmer(lon_center);
gctp_print_origin(lat_origin);
gctp_print_offsetp(false_easting,false_northing);
return(OK);
}


/* Mercator inverse equations--mapping x,y to lat/long
  --------------------------------------------------*/
long merinv
(
    double x,			/* (O) X projection coordinate 	*/
    double y,			/* (O) Y projection coordinate 	*/
    double *lon,			/* (I) Longitude 		*/
    double *lat			/* (I) Latitude 		*/
)
{
double ts;		/* small t value				*/
long flag;		/* error flag 					*/


/* Inverse equations
  -----------------*/
flag = 0;
x -= false_easting;
y -= false_northing;
ts = exp(-y/(r_major * m1));
flag = gctp_calc_phi2(e,ts,lat);
if (flag != GCTP_SUCCESS)
   return GCTP_ERROR;
*lon = adjust_lon(lon_center + x/(r_major * m1));

return(OK);
}
