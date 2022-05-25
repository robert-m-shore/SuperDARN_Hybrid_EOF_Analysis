#ifndef _CENTROID_CALCULATION_H
#define _CENTROID_CALCULATION_H

char* concat(const char *s1, const char *s2);

void itoa_new(int n, char s[]);

void reverse_new(char s[]);

double slant_range(int frang,int rsep,
		   double rxris,double range_edge,
		   int range_gate);

void geodtgc(int iopt,double *gdlat,double *gdlon,
             double *grho,double *glat,
	     double *glon,double *del);

void fldpnt(double rrho,double rlat,double rlon,double ral,
	    double rel,double r,double *frho,double *flat,
	    double *flon);

void geocnvrt(double gdlat,double gdlon,
	      double xal,double xel,double *ral,double *rel);

void fldpnth_vhm(double gdlat,double gdlon,double psi,double bore,
                 double xh,double r,double *frho,double *flat,
                 double *flon, int proptype);

void RPosGeo_vhm_rms(int center,int bcrd,int rcrd,
                double bmsep,double recrise,
                int maxbeam,
                double geolat,double geolon,double boresite,
                int frang,int rsep,
                int rxrise,double height,
		 double *rho,double *lat,double *lng, int vhmtype);

#endif
