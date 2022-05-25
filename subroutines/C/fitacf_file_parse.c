/********************************************************************
*
* fitacf_file_parse.c
*
*Code to read .fitacf binary files, check parameters of data against the
* assumed parameters used when creating the centroid positions, and write
* out the velocities, times, beam and gate numbers in Matlab-friendly format.
*
*Compile format:
* make -f makefile_fitacf_file_parse
*Call format:
* fitacf_file_parse [C binary fitted filename] [number of beams] [number of range gates] ...
*  [distance to first range gate at speed of light in microseconds] [separation of range gates at speed of light in microseconds] ...
*  [output ascii filename]
* e.g. fitacf_file_parse /users/robore/Research/Code/JH/gbr/20011201.gbr.fit 16 120 1200 300 /data/psdcomplexity/eimf/SuperDARN_Data/Velocities/gbr/gbr_20011201_velocities_GEO.dat
*Author: Rob Shore, 2016-11-03, made a fitacf version on 2017-10-17.
*
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <zlib.h>

/* My header file */
#include "fitacf_file_parse.h"

/* These header files are stored somewhere in \\samba.nerc-bas.ac.uk\superdarn\. */
#include "rtypes.h"
#include "dmap.h"
/* #include "limit.h" */
#include "rprm.h"
#include "fitdata.h"
#include "fitread.h"

int main(int argc, char *argv[]){
    /*  Definitions: */
    char *in_filename = NULL;/* the filename of fitted velocities to read in. */
    char *out_filename = NULL;/* the ascii filename of fitted velocities to save out */
    int N_BEAM;/* number of beams. */
    int N_RANGE;/* number of range gates. */
    int frang_in;/* distance to first range gate (frang) for which the centroid positions were calculated.  Units of km. */
    int rsep_in;/* separation of range gates (rsep).  Units of km. */
    /*  More definitions: */
    struct RadarParm *prm;/*  pointer for radar parameter block. */
    struct FitData *fit;/*  pointer for data structure. */
    FILE *fp;/*  pointer to fitted data file. */
    int fd;/*  possibly a counter for the number of fiducial lines in file. */
    float pwr_min = 3.0;/*  units: dB.  Minimum power required for a reading to be considered valid. */
    int beam;/* Value of beam number in radar data record. */
    int maxrange;/* Value of maximum number of data-containing range gates in radar data record. */
    int range;/* to loop over each range gate. */

    /*  Make radar parameter block and data structure: */
	prm = RadarParmMake();
    fit = FitMake();

    /*  Get the centroid parameters and fit filename details from command line: */
    /*  Note: argv[0] is the program name. */
    in_filename = argv[1];
    N_BEAM = atoi(argv[2]);
    N_RANGE = atoi(argv[3]);
    frang_in = atoi(argv[4]);
    rsep_in = atoi(argv[5]);
    out_filename = argv[6];

    /*  Data check: are the inputs being parsed OK? */
    /* printf("filename: %s\n", in_filename); */
    /* printf("N_BEAM: %d\n", N_BEAM); */
    /* printf("N_RANGE: %d\n", N_RANGE); */
    /* printf("frang_in: %d\n", frang_in); */
    /* printf("rsep_in: %d\n", rsep_in); */
    /* printf("out_filename: %s\n", out_filename); */

    /*  Open SuperDARN fit file and read and analyse data: */
    /*  Check for existence of file: */
    fp = fopen(in_filename,"r");
	if (fp==NULL){
	    fprintf(stderr, "ERROR: File not found\n");
	    exit(0);
	}/* end if: checking for existence of SuperDARN fitted data file. */
	fd = fileno(fp);

    /* Declare and open output file: */
    FILE *out_file;/* output filename of velocities. */
    out_file = fopen(out_filename,"w");

    /* Make a variable that allows me to just run a bit of code below once: */
    int run_once = 1;

    /*  Loop over all file records, store valid parts: */
    while(FitRead(fd,prm,fit) != -1){
        /* Data check: tells you if frang or rsep are why it didn't write out a given file: */
        /* if ((run_once == 1) && (fit.rng[0].qflg == 1) && (fit.rng[0].gsct == 0) && (fit.rng[0].p_l >= pwr_min)){ */
        if ((run_once == 1)){
            run_once = run_once + 1;
            if ((prm->lagfr != frang_in) || (prm->smsep != rsep_in)){
                printf("Some records missing in %s\n", out_filename);
                printf("actual frang of first beam: %d\n", prm->lagfr);
                printf("actual rsep of first beam: %d\n", prm->smsep);
                /*printf("quality flag should be 1: %d\n", fit->rng[0].qflg);*/
                /*printf("ground scatter should be 0: %d\n", fit->rng[0].gsct);*/
                /*printf("min power should exceed 3: %g\n", fit->rng[0].p_l);*/
            }/* end if: prints out frang and rsep values for first beam if they don't match the assumed values. */
            if (prm->nrang > N_RANGE){
                printf("Not saving out %s\n", out_filename);
                printf("First beam max range count: %d\n", prm->nrang);
            }/* end if: prints out maximum rabnge gate count for the first beamif it exceeds the assumed value. */
        }/* end if: I only want it to print out frang and rsep mismatches once, so this checks if it's the first record.*/

        if ((prm->lagfr == frang_in) && (prm->smsep == rsep_in)){
            /* Note: typical values from frang_in and rsep_in are 1200 and 300. That is, those are what I assumed when I made the centroid locations. */

            /* Get the beam and range gate number for this record from the radar parameter block: */
            beam = prm->bmnum;/* Single value, varying from 0 to N_BEAM-1 */
            maxrange = prm->nrang;/* Single value, invariant and always about the same as the maximum number of range gates. */

            /* Data check: print out the values of beam and maxrange for each radar data entry, so I can see which varies: */
            /* printf("beam: %d\n", beam); */
            /* printf("maxrange: %d\n", maxrange); */

            /* Data check: print out the beam epoch: */
            /* printf("year: %d\n", prm.time.yr); */

            /* Check that this beam number is within the range of valid beam numbers, likewise the maximum range gate count: */
            if ((beam >= 0) && (beam < N_BEAM) && (maxrange > 0) && (maxrange <= N_RANGE)){
                /* Loop over each range gate within this beam: */
                for (range = 0; range < maxrange; range = range + 1){
                    /* Check that this data record's scatter flag, quality flag and power level are OK: */
                    if ((fit->rng[range].qflg == 1) && (fit->rng[range].p_l >= pwr_min) && (fit->rng[range].gsct == 0)){
                        /* Note: gsct == 0 means that the signal should be ionospheric. 1 means from ground. */

                        /* Print velocity, beam number, gate number and date to file: */
                        fprintf(out_file,"%d %d %d %d %d %d %d %d %g\n",beam,range+1,prm->time.yr,prm->time.mo,prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,fit->rng[range].v);

                    }/* end if: checks if data quality flag is OK and power is above the generic minimum of 3dB.*/
                }/* end for: loop over each range gate within this beam. */
            }/* end if: only allows the data record to be read if the beam and range gate numbers are within expected ranges.*/
        }/* end if: only allows the data record to be read if the frang and rsep variables match the values assumed for them when the centroids were made.*/
    }/* end while: presumably loops over all records in fit file.*/

    /* Close fit and output files: */
    fclose(fp);
    fclose(out_file);

    return 0;/* 0 means program was successful. */
}/*  end of main program. */


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
