/*******************************************************************************
NAME                  LAMBERT AZIMUTHAL EQUAL-AREA
 
PURPOSE:	Transforms input longitude and latitude to Easting and
		Northing for the Lambert Azimuthal Equal-Area projection.  The
		longitude and latitude must be in radians.  The Easting
		and Northing values will be returned in meters.

This function was adapted from the Lambert Azimuthal Equal Area projection
code (FORTRAN) in the General Cartographic Transformation Package software
which is available from the U.S. Geological Survey National Mapping Division.
 
ALGORITHM REFERENCES

1.  "New Equal-Area Map Projections for Noncircular Regions", John P. Snyder,
    The American Cartographer, Vol 15, No. 4, October 1988, pp. 341-355.

2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
    Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
    State Government Printing Office, Washington D.C., 1987.

3.  "Software Documentation for GCTP General Cartographic Transformation
    Package", U.S. Geological Survey National Mapping Division, May 1982.
*******************************************************************************/
#include <stdio.h>
#include "oli_cproj.h"
#include "oli_local.h"

/* Variables common to all subroutines in this code file
  -----------------------------------------------------*/
static double lon_center;	/* Center longitude (projection center) */
static double R;		/* Radius of the earth (sphere)	 	*/
static double sin_lat_o;	/* Sine of the center latitude 		*/
static double cos_lat_o;	/* Cosine of the center latitude 	*/
static double false_easting;	/* x offset in meters			*/
static double false_northing;	/* y offset in meters			*/

/* Initialize the Lambert Azimuthal Equal Area projection
  ------------------------------------------------------*/
long lamazforint
(
    double r, 			/* (I) Radius of the earth (sphere) 	*/
    double center_long,		/* (I) Center longitude 		*/
    double center_lat,		/* (I) Center latitude 			*/
    double false_east,		/* x offset in meters			*/
    double false_north		/* y offset in meters			*/
)
{
/* Place parameters in static storage for common use
  -------------------------------------------------*/
R = r;
lon_center = center_long;
false_easting = false_east;
false_northing = false_north;
sincos(center_lat, &sin_lat_o, &cos_lat_o);

/* Report parameters to the user
  -----------------------------*/
gctp_print_title("LAMBERT AZIMUTHAL EQUAL-AREA"); 
gctp_print_radius(r);
gctp_print_cenlon(center_long);
gctp_print_cenlat(center_lat);
gctp_print_offsetp(false_easting,false_northing);
return(OK);
}

/* Lambert Azimuthal Equal Area forward equations--mapping lat,long to x,y
  -----------------------------------------------------------------------*/
long lamazfor
(
    double lon,			/* (I) Longitude */
    double lat,			/* (I) Latitude */
    double *x,			/* (O) X projection coordinate */
    double *y			/* (O) Y projection coordinate */
)
{
double delta_lon;	/* Delta longitude (Given longitude - center 	*/
double sin_delta_lon;	/* Sine of the delta longitude 			*/
double cos_delta_lon;	/* Cosine of the delta longitude 		*/
double sin_lat;		/* Sine of the given latitude 			*/
double cos_lat;		/* Cosine of the given latitude 		*/
double g;		/* temporary varialbe				*/
double ksp;		/* heigth above elipsiod			*/

/* Forward equations
  -----------------*/
delta_lon = adjust_lon(lon - lon_center);
sincos(lat, &sin_lat, &cos_lat);
sincos(delta_lon, &sin_delta_lon, &cos_delta_lon);
g = sin_lat_o * sin_lat + cos_lat_o * cos_lat * cos_delta_lon;
if (g == -1.0) 
   {
   GCTP_PRINT_ERROR("Point projects to a circle of radius = %f\n", 2.0 * R);
   return(113);
   }
ksp = R * sqrt(2.0 / (1.0 + g));
*x = ksp * cos_lat * sin_delta_lon + false_easting;
*y = ksp * (cos_lat_o * sin_lat - sin_lat_o * cos_lat * cos_delta_lon) + 
	false_northing;
return(OK);
}
