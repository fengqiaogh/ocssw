/*******************************************************************************
NAME                 GENERAL VERTICAL NEAR-SIDE PERSPECTIVE 

PURPOSE:	Transforms input longitude and latitude to Easting and
		Northing for the General Vertical Near-Side Perspective
		projection.  The longitude and latitude must be in
		radians.  The Easting and Northing values will be
		returned in meters.

This function was adapted from the General Vertical Near-Side Perspective
projection code (FORTRAN) in the General Cartographic Transformation Package
software which is available from the U.S. Geological Survey National Mapping
Division.
 
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
static double p;		/* Height above sphere			*/
static double sin_p15;		/* Sine of the center latitude 		*/
static double cos_p15;		/* Cosine of the center latitude 	*/
static double false_easting;	/* x offset in meters			*/
static double false_northing;	/* y offset in meters			*/

/* Initialize the General Vertical Near-Side Perspective projection
  ---------------------------------------------------------------*/
long gvnspforint
(
    double r, 			/* (I) Radius of the earth (sphere) 	*/
    double h,			/* height above sphere			*/
    double center_long,		/* (I) Center longitude 		*/
    double center_lat,		/* (I) Center latitude 			*/
    double false_east,		/* x offset in meters			*/
    double false_north		/* y offset in meters			*/
)
{
/* Place parameters in static storage for common use
  -------------------------------------------------*/
R = r;
p = 1.0 + h / R;
lon_center = center_long;
false_easting = false_east;
false_northing = false_north;
sincos(center_lat, &sin_p15, &cos_p15);

/* Report parameters to the user
  -----------------------------*/
gctp_print_title("GENERAL VERTICAL NEAR-SIDE PERSPECTIVE"); 
gctp_print_radius(r);
gctp_print_genrpt(h,"Height of Point Above Surface of Sphere:   ");
gctp_print_cenlon(center_long);
gctp_print_cenlat(center_lat);
gctp_print_offsetp(false_easting,false_northing);
return(OK);
}

/* General Vertical Near-Side Perspective forward equations--mapping 
   lat,long to x,y
  ----------------------------------------------------------------*/
long gvnspfor
(
    double lon,			/* (I) Longitude */
    double lat,			/* (I) Latitude */
    double *x,			/* (O) X projection coordinate */
    double *y			/* (O) Y projection coordinate */
)

{
double dlon;
double sinphi,cosphi;
double coslon;
double g;
double ksp;

/* Forward equations
  -----------------*/
dlon = adjust_lon(lon - lon_center);
sincos(lat,&sinphi,&cosphi);
coslon = cos(dlon);
g = sin_p15 * sinphi + cos_p15 * cosphi * coslon;
if (g < (1.0/ p))
   {
   GCTP_PRINT_ERROR("Point cannot be projected");
   return(153);
   }
ksp = (p - 1.0)/(p - g);
*x = false_easting + R * ksp * cosphi * sin(dlon);
*y = false_northing + R * ksp * (cos_p15 * sinphi - sin_p15 * cosphi * coslon);

return(OK);
}
