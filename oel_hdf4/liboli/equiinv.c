/*******************************************************************************
NAME                             EQUIRECTANGULAR 

PURPOSE:	Transforms input Easting and Northing to longitude and
		latitude for the Equirectangular projection.  The
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
static double lon_center;	/* Center longitude (projection center) */
static double lat_origin;	/* center latitude			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/

/* Initialize the Equirectangular projection
  ----------------------------------------*/
long equiinvint
(
    double r_maj,			/* major axis			*/
    double center_lon,		/* center longitude		*/
    double lat1,			/* latitude of true scale	*/
    double false_east,		/* x offset in meters		*/
    double false_north		/* y offset in meters		*/
)
{

/* Place parameters in static storage for common use
  -------------------------------------------------*/
r_major = r_maj;
lon_center = center_lon;
lat_origin = lat1;
false_northing = false_north;
false_easting = false_east;

/* Report parameters to the user
  -----------------------------*/
gctp_print_title("EQUIRECTANGULAR"); 
gctp_print_radius(r_major);
gctp_print_cenlonmer(lon_center);
gctp_print_origin(lat_origin);
gctp_print_offsetp(false_easting,false_northing);
return(OK);
}


/* Equirectangular inverse equations--mapping x,y to lat/long
  ---------------------------------------------------------*/
long equiinv
(
    double x,			/* (O) X projection coordinate 	*/
    double y,			/* (O) Y projection coordinate 	*/
    double *lon,			/* (I) Longitude 		*/
    double *lat			/* (I) Latitude 		*/
)
{
/* Inverse equations
  -----------------*/
x -= false_easting;
y -= false_northing;
*lat = y / r_major;
if (fabs(*lat) > HALF_PI)
   {
   GCTP_PRINT_ERROR("Input data error");
   return(174);
   }
*lon = adjust_lon(lon_center + x / (r_major * cos(lat_origin)));
return(OK);
}
