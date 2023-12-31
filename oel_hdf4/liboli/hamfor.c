/*******************************************************************************
NAME                            HAMMER 

PURPOSE:	Transforms input longitude and latitude to Easting and
		Northing for the Hammer projection.  The longitude
		and latitude must be in radians.  The Easting and
		Northing values will be returned in meters.
  
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
#include "oli_cproj.h"
#include "oli_local.h"

/* Variables common to all subroutines in this code file
  -----------------------------------------------------*/
static double lon_center;	/* Center longitude (projection center) */
static double R;		/* Radius of the earth (sphere)	 	*/
static double false_easting;	/* x offset in meters			*/
static double false_northing;	/* y offset in meters			*/

/* Initialize the HAMMER projection
  -------------------------------*/
long hamforint
(
    double r, 			/* (I) Radius of the earth (sphere) 	*/
    double center_long,		/* (I) Center longitude 		*/
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

/* Report parameters to the user
  -----------------------------*/
gctp_print_title("HAMMER"); 
gctp_print_radius(r);
gctp_print_cenlon(center_long);
gctp_print_offsetp(false_easting,false_northing);
return(OK);
}

/* HAMMER forward equations--mapping lat,long to x,y
  ------------------------------------------------------------*/
long hamfor
(
    double lon,			/* (I) Longitude */
    double lat,			/* (I) Latitude */
    double *x,			/* (O) X projection coordinate */
    double *y			/* (O) Y projection coordinate */
)

{
double dlon;
double fac;

/* Forward equations
  -----------------*/
dlon = adjust_lon(lon - lon_center);

fac  = R * 1.414213562 / sqrt(1.0 + cos(lat) * cos(dlon / 2.0));
*x = false_easting + fac * 2.0 * cos(lat) * sin(dlon / 2.0);
*y = false_northing + fac * sin(lat);

return(OK);
}
