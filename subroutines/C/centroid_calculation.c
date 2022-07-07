#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rmath.h"
#include "centroid_calculation.h"

int main(void)
{

/* Definitions */
int i_radar;

/* Define parameters used in call of RPosGeo */
/* Definitions from /users/gchi/superdarn/coverage/code/rblocs_vhm.c */
int ipbeam, iprange;
int center=1;
int vhmtype=0;
int frang=180;/* to match SuperDARN data files */
/* int frang=165; to match Gareth's existing centroid metadata set. */
int rsep=45;
double height=0.0;
double rho,lat,lng;

/* Define radar acronymns: */
const char *radar_identifiers[36];
radar_identifiers[0] = "ade";
radar_identifiers[1] = "adw";
radar_identifiers[2] = "bks";
radar_identifiers[3] = "cve";
radar_identifiers[4] = "cvw";
radar_identifiers[5] = "cly";
radar_identifiers[6] = "fhe";
radar_identifiers[7] = "fhw";
radar_identifiers[8] = "gbr";
radar_identifiers[9] = "han";
radar_identifiers[10] = "hok";
radar_identifiers[11] = "hkw";
radar_identifiers[12] = "inv";
radar_identifiers[13] = "kap";
radar_identifiers[14] = "ksr";
radar_identifiers[15] = "kod";
radar_identifiers[16] = "pyk";
radar_identifiers[17] = "pgr";
radar_identifiers[18] = "rkn";
radar_identifiers[19] = "sas";
radar_identifiers[20] = "sch";
radar_identifiers[21] = "sto";
radar_identifiers[22] = "wal";
radar_identifiers[23] = "bpk";
radar_identifiers[24] = "dce";
radar_identifiers[25] = "fir";
radar_identifiers[26] = "hal";
radar_identifiers[27] = "ker";
radar_identifiers[28] = "mcm";
radar_identifiers[29] = "san";
radar_identifiers[30] = "sps";
radar_identifiers[31] = "sye";
radar_identifiers[32] = "sys";
radar_identifiers[33] = "tig";
radar_identifiers[34] = "unw";
radar_identifiers[35] = "zho";


for(i_radar = 0; i_radar < 36; i_radar = i_radar + 1){
	/* Load the hardware file for this radar and determine the number of rows of data in it: */
	const char* radar_acronym = radar_identifiers[i_radar];
	char* file_string = concat("/local/users/robore/Code/JH/SHEAR/metadata/hardware_files/hdw.dat.", radar_acronym);
	/* Data check: Display hardware filename: */
	printf("Hardware filename: %s\n",file_string);

	/* Declare and open hardware file for this radar: */
	FILE *hardware_file;
	hardware_file = fopen(file_string,"r");

	/* loop over each line of the hardware file, and compute centroid positions for non-header lines: */
	char line[256];
	int data_line_count = 0;
	while (fgets(line, sizeof(line), hardware_file)) {
		/* note that fgets doesn't strip the terminating \n, checking its presence would allow to handle lines longer that sizeof(line). */

	    /* Data check: Display line(s) of hardware file: */
	    /* printf("Single line of hardware file: %s", line); */
	    char first_letter_of_line[1];
		strncpy(first_letter_of_line, line,1);
		first_letter_of_line[0] = 0; /* null terminate destination */

		/* Data check: display first letter of line: */
		/* printf("First letter of this line: %s\n",first_letter_of_line); */

	    /* If the first letter of the line is not indicative of a header, store the line in an array: */
	    if(strcmp(first_letter_of_line,"#")<0){
			/* Data check: Display line(s) of hardware file which are not headers: */
			/* printf("Non-header line in full: %s", line); */

	    	/* Advance valid-data increment indicator, which controls the output filename: */
	    	data_line_count = data_line_count + 1;

	    	/* Data check: Display the monotonic fiducial counter for the data line in question: */
	    	/* printf("Data line counter: %d\n", data_line_count); */

	    	/* Parse line: */
			int i_part = 0;
			char *line_part = strtok(line, " ");
			char *text_array[19];

			while (line_part != NULL){
			    text_array[i_part++] = line_part;
			    line_part = strtok(NULL, " ");
			}/* end while: loop over each part of the line. */

			/* Data check: Display each part of the parsed line: */
			/* for (i_part = 0; i_part < 19; i_part = i_part + 1){ */
			/*     printf("%s\n", text_array[i_part]); */
			/* }end for: loop over each part of the line and display it. */

			/* Convert (required) text values to numbers: */
			/* int station_ID = atoi(text_array[0]); */ /* 0 means the first element of this array. */
			/* printf("Station ID: %d\n", station_ID);Data check: display this element of the array. */
			int last_valid_year = atoi(text_array[1]);
			int last_valid_second_of_year = atoi(text_array[2]);
			double radar_geodetic_lat = atof(text_array[3]);
			double radar_geodetic_lon = atof(text_array[4]);
			/* double radar_alt = atof(text_array[5]); */
			double scanning_boresite = atof(text_array[6]);
			double beam_separation = atof(text_array[7]);
			/* int doppler_velocity_sign = atoi(text_array[8]); */
			/* double attenuator_step = atof(text_array[9]); */
			/* double t_diff = atof(text_array[10]); */
			/* int phase_sign = atoi(text_array[11]); */
			/* double interferometer_offset_coord_1 = atof(text_array[12]); */
			/* double interferometer_offset_coord_2 = atof(text_array[13]); */
			/* double interferometer_offset_coord_3 = atof(text_array[14]); */
			double rx_rise_time = atof(text_array[15]);
			/* int num_attenuation_stages = atoi(text_array[16]); */
			int number_of_gates = atoi(text_array[17]);
			int number_of_beams = atoi(text_array[18]);

			/* Data check: */
			/* printf("Rx rise time text: %s\n", text_array[15]); */
			/* printf("Rx rise time data: %f\n", rx_rise_time); */

			/* We would define here a stripped-down version of RadarSite used by rblocs_vhm.c and defined in radar.h, */
			/* but we only need the variables bmsep, recrise, maxbeam, geolat, geolon and boresite, so we load them into RPoseGeo directly. */
			/* In order, they are: beam separation (degrees), Analog Rx rise time (microseconds), */
			/* maximum number of beams for the radar, radar latitude (geodetic), radar longitude (geodetic), */
			/* boresite angle (from some reference line). */

            /* Define a filename based on the station acronym, and monotonic fiducial counter for the data line in question: */
			char* outfile_string_1 = concat("/local/users/robore/Code/JH/SHEAR/metadata/centroid_locations/centroid_locations_", radar_acronym);
			char* outfile_string_2 = concat(outfile_string_1, "_");
			char outfile_number[1];
			itoa_new(data_line_count,outfile_number);
			char* outfile_string_3 = concat(outfile_string_2, outfile_number);
			char* outfile_string = concat(outfile_string_3, ".dat");

			/* Data check: Display output filename: */
			printf("Output filename: %s\n", outfile_string);

			/* Open output file: */
			FILE *out_file;
			out_file = fopen(outfile_string,"w");

			/* On first line of output file, write the number of beams and gates, the expiry year and the expiry second-of-year: */
			fprintf(out_file,"%d %d %d %d\n",number_of_beams,number_of_gates,last_valid_year,last_valid_second_of_year);

			/* On second line of output file, write the vhm type (integer), the 'height' variable value put into RPosGeo (units unknown), */
			/* the frang and rsep values (both in km). */
			fprintf(out_file,"%d %f %d %d\n",vhmtype,height,frang,rsep);

			/* Data check: print inputs to RPosGeo_vhm_rms: */
			/* printf("center: %d\n", center); */
			/* printf("beam separation: %f\n", beam_separation); */
			/* printf("rx rise time: %f\n", rx_rise_time); */
			/* printf("number of beams: %d\n", number_of_beams); */
			/* printf("radar geodetic lat: %f\n", radar_geodetic_lat); */
			/* printf("radar geodetic lon: %f\n", radar_geodetic_lon); */
			/* printf("scanning boresite: %f\n", scanning_boresite); */
			/* printf("frang: %d\n", frang); */
			/* printf("rsep: %d\n", rsep); */
			/* printf("height: %f\n", height); */
			/* printf("vhmtype: %d\n", vhmtype); */


			/* Compute all centroid positions for this hardware file line. */
			/* Code from /users/gchi/superdarn/coverage/code/rblocs_vhm.c. */
			for (ipbeam=0;ipbeam<number_of_beams;ipbeam=ipbeam+1){
				for (iprange=0;iprange<number_of_gates;iprange=iprange+1){
					/* Call for cenrtoid position calculator: */
					RPosGeo_vhm_rms(center,ipbeam,iprange,beam_separation,rx_rise_time,number_of_beams,radar_geodetic_lat,radar_geodetic_lon,scanning_boresite,frang,rsep,rx_rise_time,height,&rho,&lat,&lng,vhmtype);

					/* RPosGeo_vhm_rms inputs copied from the above line of code: */
					/* (center,ipbeam,iprange,beam_separation,rx_rise_time,number_of_beams,radar_geodetic_lat,radar_geodetic_lon,scanning_boresite,frang,rsep,rx_rise_time,height,&rho,&lat,&lng,vhmtype) */
					/* and expected values: these are the translation of the input names as seen by the RPosGeo_vhm_rms routine: */
					/* (center,bcrd,  rcrd,   bmsep,          recrise,     maxbeam,        geolat,            geolon,            boresite,         frang,rsep,rxrise,      height,rho, lat, lng, vhmtype) */

					/* save output of this beam/gate combination: */
					fprintf(out_file,"%d %d %g %g\n",ipbeam,iprange+1,lat,lng);
				}/* end for: loop over range gates. */
			}/* end for: loop over radar beams. */

			fclose(out_file);

	    }/* end if: conditional check for whether this line of the hardware file is a header or not, and acting accordingly. */

	}/* end while: loop over all lines in hardware file. */


	/* may check feof here to make a difference between eof and io failure -- network timeout for instance */

	fclose(hardware_file);

	/* free(file_string); deallocate the string */


}/* end for: loop over all radars. */

return 0; /* means program was successful. */

}/* end main */


/***********************************************************************/
/* Concatenation function from 'http://stackoverflow.com/questions/8465006/how-to-concatenate-2-strings-in-c'. */

char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);/* +1 for the zero-terminator */
    /* in real code you would check for errors in malloc here */
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

/***********************************************************************/
/* itoa:  convert n to characters in s */
/* From https://en.wikibooks.org/wiki/C_Programming/C_Reference/stdlib.h/itoa. */
 void itoa_new(int n, char s[])
 {
     int i, sign;

     if ((sign = n) < 0)  /* record sign */
         n = -n;          /* make n positive */
     i = 0;
     do {       /* generate digits in reverse order */
         s[i++] = n % 10 + '0';   /* get next digit */
     } while ((n /= 10) > 0);     /* delete it */
     if (sign < 0)
         s[i++] = '-';
     s[i] = '\0';
     reverse_new(s);
 }

/***********************************************************************/
/* reverse:  reverse string s in place */
/* From https://en.wikibooks.org/wiki/C_Programming/C_Reference/stdlib.h/itoa. */
  void reverse_new(char s[])
  {
      int i, j;
      char c;

      for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
          c = s[i];
          s[i] = s[j];
          s[j] = c;
      }
 }

/************************************************************************/

double slant_range(int frang,int rsep,
		   double rxris,double range_edge,
		   int range_gate) {
   int lagfr,smsep;
   lagfr=frang*20/3;
   smsep=rsep*20/3;
   return (lagfr-rxris+(range_gate-1)*smsep+range_edge)*0.15;
}

/************************************************************************/

void geodtgc(int iopt,double *gdlat,double *gdlon,
             double *grho,double *glat,
			 double *glon,double *del) {

   double a=6378.16;
   double f=1.0/298.25;
   double b,e2;

   b=a*(1.0-f);
   e2=(a*a)/(b*b)-1;
   if (iopt>0) {
     *glat=atand( (b*b)/(a*a)*tand(*gdlat));
     *glon=*gdlon;
     if (*glon > 180) *glon=*glon-360;
   } else {
     *gdlat=atand( (a*a)/(b*b)*tand(*glat));
     *gdlon=*glon;
   }
   *grho=a/sqrt(1.0+e2*sind(*glat)*sind(*glat));
   *del=*gdlat-*glat;
}

/************************************************************************/

void fldpnt(double rrho,double rlat,double rlon,double ral,
			double rel,double r,double *frho,double *flat,
                        double *flon) {

   double rx,ry,rz,sx,sy,sz,tx,ty,tz;
   double sinteta;

   /* convert from global spherical to global cartesian*/

   sinteta=sind(90.0-rlat);
   rx=rrho*sinteta*cosd(rlon);
   ry=rrho*sinteta*sind(rlon);
   rz=rrho*cosd(90.0-rlat);

   sx=-r*cosd(rel)*cosd(ral);
   sy=r*cosd(rel)*sind(ral);
   sz=r*sind(rel);

   tx  =  cosd(90.0-rlat)*sx + sind(90.0-rlat)*sz;
   ty  =  sy;
   tz  = -sind(90.0-rlat)*sx + cosd(90.0-rlat)*sz;
   sx  =  cosd(rlon)*tx - sind(rlon)*ty;
   sy  =  sind(rlon)*tx + cosd(rlon)*ty;
   sz  =  tz;

   tx=rx+sx;
   ty=ry+sy;
   tz=rz+sz;

   /* convert from cartesian back to global spherical*/
   *frho=sqrt((tx*tx)+(ty*ty)+(tz*tz));
   *flat=90.0-acosd(tz/(*frho));
   if ((tx==0) && (ty==0)) *flon=0;
   else *flon=atan2d(ty,tx);
}

/************************************************************************/

void geocnvrt(double gdlat,double gdlon,
			  double xal,double xel,double *ral,double *rel) {

  double kxg,kyg,kzg,kxr,kyr,kzr;
  double rrad,rlat,rlon,del;

  kxg=cosd(xel)*sind(xal);
  kyg=cosd(xel)*cosd(xal);
  kzg=sind(xel);
  geodtgc(1,&gdlat,&gdlon,&rrad,&rlat,&rlon,&del);
  kxr=kxg;
  kyr=kyg*cosd(del)+kzg*sind(del);
  kzr=-kyg*sind(del)+kzg*cosd(del);

  *ral=atan2d(kxr,kyr);
  *rel=atand(kzr/sqrt((kxr*kxr)+(kyr*kyr)));
}

/************************************************************************/

void fldpnth_vhm(double gdlat,double gdlon,double psi,double bore,
  double xh,double r,double *frho,double *flat, double *flon, int proptype) {

  double rrad,rlat,rlon,del;
  double tan_azi,azi,rel,xel,fhx,xal,rrho,ral;
  double dum,dum1,dum2,dum3;
  double frad;
  double phi,beta;

  geodtgc(1,&gdlat,&gdlon,&rrad,&rlat,&rlon,&del);
  rrho=rrad;
  frad=rrad;

  do {
    *frho=frad+xh;
    rel=asind(((*frho**frho) - (rrad*rrad) - (r*r)) / (2.0*rrad*r));

    if (proptype == 0){
      xel=rel;
    }
    else{
      phi = acosd(((frad*frad) + (*frho**frho) - (r*r))/(2.0*frad**frho));
      beta = asind((frad*sind(phi/3.0))/(r/3.0));
      xel = 90.0 - beta - (phi/3.0);
    }

    if (((cosd(psi)*cosd(psi))-(sind(xel)*sind(xel)))<0) tan_azi=1e32;
    else tan_azi=sqrt( (sind(psi)*sind(psi))/((cosd(psi)*cosd(psi))-(sind(xel)*sind(xel))));
    if (psi>0) azi=atand(tan_azi)*1.0;
    else azi=atand(tan_azi)*-1.0;
    xal=azi+bore;
    geocnvrt(gdlat,gdlon,xal,xel,&ral,&dum);

    fldpnt(rrho,rlat,rlon,ral,rel,r,frho,flat,flon);
    geodtgc(-1,&dum1,&dum2,&frad,flat,flon,&dum3);
    fhx=*frho-frad;
  } while(fabs(fhx-xh) > 0.5);
}

/**********************************************************************/

void RPosGeo_vhm_rms(int center,int bcrd,int rcrd,
                double bmsep,double recrise,
                int maxbeam,
                double geolat,double geolon,double boresite,
                int frang,int rsep,
                int rxrise,double height,
		 double *rho,double *lat,double *lng, int vhmtype) {

  int proptype=0;
  double vh=0.0;
  double rx;
  double psi,d;
  double re=6356.779;
  double offset;
  double bm_edge=0;
  double range_edge=0;

  if (center==0) {
    bm_edge=-bmsep*0.5;
    range_edge=-0.5*rsep*20/3;
  }

  if (rxrise==0) rx=recrise;
  else rx=rxrise;
  offset=maxbeam/2.0-0.5;
  psi=bmsep*(bcrd-offset)+bm_edge;
  d=slant_range(frang,rsep,rx,range_edge,rcrd+1);

  if (vhmtype == 0){
    if (d <= 787.5){
      proptype = 0;
      vh = 108.974 + (0.0191271*d) + (6.68283e-5*d*d);
    }
    else if ((d > 787.5) && (d <= 2137.5)){
      proptype = 0;
      vh = 384.416 - (0.178640*d) + (1.81405e-4*d*d);
    }
    else if (d > 2137.5){
      proptype = 1;
      vh = 1098.28 - (0.354557*d) + (9.39961e-5*d*d);
    }
  }
  else if (vhmtype == 1){
    proptype = 0;
    vh = 108.974 + (0.0191271*d) + (6.68283e-5*d*d);
  }
  else if (vhmtype == 2){
    proptype = 0;
    vh = 384.416 - (0.178640*d) + (1.81405e-4*d*d);
  }
  else if (vhmtype == 3){
    proptype = 1;
    vh = 1098.28 - (0.354557*d) + (9.39961e-5*d*d);
  }
  else{
    proptype = 0;
    if (height < 90) height=-re+sqrt((re*re)+(2.0*d*re*sind(height))+(d*d));
    if (height <= 150.0) vh = height;
    else{
      if (d <= 600.0) vh = 115.0;
      else if ((d > 600.0) && (d < 800.0))
        vh = (((d-600.0)/200.0) * (height - 115.0)) + 115.0;
      else vh = height;
    }
    if (d < 150.0) vh = (d/150.0) * 115.0;
  }

  fldpnth_vhm(geolat,geolon,psi,boresite,
	      vh,d,rho,lat,lng,proptype);
}

/**********************************************************************/
