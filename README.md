# SuperDARN_Hybrid_EOF_Analysis

Program name: 
    SuperDARN_Hybrid_EOF_Analysis.m


Program objectives:
    This program is designed to infill SuperDARN line-of-sight plasma velocity measurements on an equal-area grid of spatial locations in magnetic latitude and magnetic local time coordinates, producing a reanalysis dataset of the north and east components of the horizontal plasma velocity vectors in the quasi-dipole coordinate frame, spanning the northern polar region.  In the process, the SuperDARN line-of-sight plasma velocities are decomposed into independent spatial and temporal modes of variability, using a temporal batch size of one calendar month at a time.  The key stages of the program are as follows: 
    1.	Loads in autocorrelation functions of line-of-sight signals backscattered from ionospheric electron density irregularities.
    2.	The autocorrelation functions are converted to line-of-sight Doppler velocities by applying linear fits to the autocorrelation function phase variations using version 4.5 of the radar software toolkit (RST 4.5 and FITACF 2.5).
    3.	The estimated locations of the SuperDARN measurements are converted into the quasi-dipole coordinate frame, and the quasi-dipole basis vectors are used to define a relationship between the line-of-sight measurements and the local quasi-dipole geomagnetic north direction – this will allow the full plasma velocity vector to be inferred from multiple line-of-sight measurements.
    4.	The line-of-sight velocity data are binned in 5-minute median temporal bins, equal-area spatial bins defined in quasi-dipole latitude and quasi-dipole magnetic local time, and line-of-sight angle partition bins defined with respect to quasi-dipole geomagnetic north.
    5.	Data selection is applied to the binned line-of-sight plasma velocities to remove ground scatter.
    6.	The binned data for each calendar month are used in two iterative EOF solutions: one which prioritises good temporal coverage, and another which aims to maximise signal recovery in poorly-sampled spatial locations.
    7.	The two EOF models are compared in terms of how well each recovers the original data variability, and the best-performing model for each spatial location is selected.  The resulting hybrid dataset comprises a continuous record of SuperDARN horizontal plasma velocity vectors at each spatial location and at 5-minute temporal cadence, with no missing information, and without relying on external drivers (like the solar wind) to complete the data coverage.


Note: in the description below, the wildcards 'code_directory', 'data_directory', and 'output_directory' are directory paths on the British Antarctic Survey servers, and these variables are already assigned their appropriate directory paths within program SuperDARN_Hybrid_EOF_Analysis.m.  The user may redefine them within that program if required.  If so, then the following (empty) subdirectories must be created within 'data_directory': [data_directory]/Binned_Velocities, and [data_directory]/Velocities.  These are used to store interstitial files produced during the analysis procedure.


Description of input data requirements: 
    
    SuperDARN measurement files:
        Files of autocorrelation functions of line-of-sight signals backscattered from ionospheric electron density irregularities.  These are available at the BAS SuperDARN data mirror (https://www.bas.ac.uk/project/superdarn). To gain access to the SuperDARN data, email bassuperdarndata@bas.ac.uk to set up a user account for the BAS SuperDARN data mirror.  These files may have a variety of valid filename formats.  Prior to July 2006, the files have the suffix .dat, and following this date, they have the suffix .rawacf.  The files by default contain 2-hours of data each: these files are formatted as follows: 
        $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).${FILETYPE}, or 
        $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).$(CHANNEL).${FILETYPE}, e.g.: 
        '20020918.1243.43.gbr.dat' or
        '20070918.1243.43.gbr.rawacf' or
        '20020918.1243.43.fir.a.dat' or 
        '20070918.1243.43.fir.a.rawacf'.
        where the channel identifier is a single letter, and absent for stations which only have a single channel.
        Alternatively, the files may be already concatenated into daily files, and in this case the format is as follows: 
        $(YEAR)$(MONTH)$(DAY).$(STATION).${FILETYPE}, or 
        $(YEAR)$(MONTH)$(DAY).$(STATION).$(CHANNEL).${FILETYPE}, e.g.: 
        '20021218.gbr.dat' or
        '20101218.gbr.rawacf' or
        '20021218.fir.a.dat' or 
        '20101218.fir.a.rawacf'.
        These files are large, and the user may not wish to store them locally.  They are accessed by program SuperDARN_Hybrid_EOF_Analysis.m via connection to //samba.nerc-bas.ac.uk/data/superdarn (which is at /data/superdarn when on the BAS network).
		
	Files stored locally in [code_directory]/metadata: 
        metadata/bin_coordinates/north_polar_region_equal_area_bins.mat: 
            These are the coordinates of the equal area spatial bins used in this analysis. Their format is described within program SuperDARN_Hybrid_EOF_Analysis.m, along with the (commented out) code used to produce these coordinates, which calls the MATLAB programs in [code_directory]/subroutines/equal_area_sphere_partitions.
        metadata/centroid_locations/centroid_locations_ade_1.dat, and other files in that directory: 
            These are the estimated geographic locations of the SuperDARN measurements, along with some records of the parameters of the radar beam configuration.  The exact format of the values in these files is described in program SuperDARN_Hybrid_EOF_Analysis.m, where the variable centroid_data_bulk is first defined.  These files are created by program [code_directory]/subroutines/C/centroid_calculation.c (described below).  There is one file per configuration setup for each radar.  The configuration setups are date-ordered, and are defined by the entries in the files in metadata/hardware_files.
        metadata/hardware_files/hdw.dat.ade, and similar files: 
            These files were downloaded from http://davit1.ece.vt.edu/hdw/[filename for each radar] on 2021/03/19, with subsequent edits to the whitespace for the files for radars cve, cvw, cly, hkw.  These files describe the hardware parameters that are used by the radar control software and the analysis software.  Each file has a header which more completely describes the file contents.


Description of subroutines, stored in [code_directory]/subroutines: 
    
    subroutines/BASH/fit_script_v6: 
        This is a BASH script which calls the make_fit routine in the RST software package.  It is designed to convert .dat files to .fit format for a given time range.  Source: supplied by Gareth Chisham (gchi@bas.ac.uk).
    
	subroutines/BASH/fitacf_script_v6: 
        This is a BASH script which calls the make_fit routine in the RST software package.  It is designed to convert .rawacf files to .fitacf format for a given time range.  Source: supplied by Gareth Chisham (gchi@bas.ac.uk).
    
	subroutines/C/centroid_calculation.c
        This is a C program which applies the virtual height model described in Chisham et al. 2008 (https://doi.org/10.5194/angeo-26-823-2008), which estimates geolocations for the SuperDARN line-of-sight velocity measurements, based on the radar beam configuration parameters.
    
	subroutines/C/fit_file_parse.c
        This is a C program which reads .fit files (in C binary format).  The program checks the parameters of the binary file data against the assumed parameters used when creating the centroid positions (using program subroutines/C/centroid_calculation.c).  It also applies data selection: the values retained have a minimum power level of 3dB, ground scatter flag inactive, quality flag active, and have a range gate value equal to or within the values specified in program SuperDARN_Hybrid_EOF_Analysis.m (minimum range gate 11, maximum range gate 150, be default).  This program finally writes out the velocities, times, beam and gate numbers in friendlier (ascii) format.
    
	subroutines/C/fitacf_file_parse.c
        This is a C program which reads .fitacf files (in C binary format).  The program checks the parameters of the binary file data against the assumed parameters used when creating the centroid positions (using program subroutines/C/centroid_calculation.c).  It also applies data selection: the values retained have a minimum power level of 3dB, ground scatter flag inactive, quality flag active, and have a range gate value equal to or within the values specified in program SuperDARN_Hybrid_EOF_Analysis.m (minimum range gate 11, maximum range gate 150, be default).  This program finally writes out the velocities, times, beam and gate numbers in friendlier (ascii) format.
    
	subroutines/equal_area_sphere_partitions/: 
        These MATLAB scripts were used to produce the coordinates of the equal area spatial bins used in this analysis, which are saved in the file metadata/bin_coordinates/north_polar_region_equal_area_bins.mat.  All these MATLAB scripts were obtained via personal communication from Ciaran Beggan at the British Geological Survey, and originally from Paul Leopardi.  The subroutines are described in Leopardi's thesis (Leopardi, P. C. (2007). Distributing points on the sphere: partitions, separation, quadrature and energy. Doctoral dissertation, University of New South Wales, Sydney, Australia).  Nearby the load-in of the file north_polar_region_equal_area_bins.mat in program SuperDARN_Hybrid_EOF_Analysis.m, example code which calls these MATLAB scripts (and produces the contents of the north_polar_region_equal_area_bins.mat file) is included in commented-out format.
    
	subroutines/quasi_dipole/: 
        The MATLAB scripts and associated supporting files in this directory were provided by Nils Olsen at DTU Space, and are a MATLAB wrapper for the FORTRAN code developed by John Emmert, described in Emmert, J. T., A. D. Richmond, and D. P. Drob, A computationally compact representation of Magnetic-Apex and Quasi-Dipole coordinates with smooth base vectors, J. Geophys. Res., 115, doi:10.1029/2010JA015326, 2010.  This FORTRAN code is an implementation of the techniques first described by Art Richmond in Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.  Please note that the MATLAB executable for the FORTRAN code will only run on a 64 bit Linux system.


Description of outputs: 
    
    Fitted line-of-sight velocity files in C binary format, saved to [data_directory]/Velocities/[radar name]/: 
        These are files of line-of-sight velocities for a given radar.  They are obtained by running the BASH scripts fit_script_v6 and fitacf_script_v6 (in [code_directory]/subroutines/BASH) on files of autocorrelation functions (which have either .dat or .rawacf as suffix).  The fitted velocity files have either .fit or .fitacf as suffix.  The filenames have the following formats: 
        $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).${FILETYPE}, or 
        $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).$(CHANNEL).${FILETYPE}, e.g.: 
        '20020918.1243.43.gbr.fit' or
        '20070918.1243.43.gbr.fitacf' or
        '20020918.1243.43.fir.a.fit' or 
        '20070918.1243.43.fir.a.fitacf'.
        where the channel identifier is a single letter, and absent for stations which only have a single channel.
        Alternatively, the source autocorrelation files may be already concatenated into daily files, and in this case the format of the fitted files is as follows: 
        $(YEAR)$(MONTH)$(DAY).$(STATION).${FILETYPE}, or 
        $(YEAR)$(MONTH)$(DAY).$(STATION).$(CHANNEL).${FILETYPE}, e.g.: 
        '20021218.gbr.fit' or
        '20101218.gbr.fitacf' or
        '20021218.fir.a.fit' or 
        '20101218.fir.a.fitacf'.
        These files are deleted by program SuperDARN_Hybrid_EOF_Analysis.m after their conversion from C binary format to ascii format.
		
    Log files from the BASH scripts, saved to [code_directory]/subroutines/BASH/: 
        These contain the runtime reports from the BASH scripts fit_script_v6 and fitacf_script_v6.  The filename format is cover_*.log.  These log files are deleted by program SuperDARN_Hybrid_EOF_Analysis.m after processing data from all radars in a given calendar month.
    
    Fitted line-of-sight velocity files in ascii format, saved to [data_directory]/Velocities/[radar name]/: 
        The C binary format files of line-of-sight velocity (described above) which are output by the BASH scripts fit_script_v6 and fitacf_script_v6 are converted to these ascii-formatted files by programs fit_file_parse.c and fitacf_file_parse.c (each of which is in [code_directory]/subroutines/C).  The ascii files can either be of daily or subdaily cadence.  Subdaily files have the following format: 
        [radar name]_[year][month][day]_[hour][minute]_[second]_velocities_GEO_[version identifier].dat, e.g.: gbr_20170928_1848_04_velocities_GEO_v7.dat, where GEO indicates that these line-of-sight velocities are still in geographic coordinates, and [version identifier] is a string indicating the program settings used to process the SuperDARN velocities: it is v7 by default, and can be manually set in program SuperDARN_Hybrid_EOF_Analysis.m.
        Daily files have the following format: 
        [radar name]_[year][month][day]_velocities_GEO_[version identifier].dat, e.g.: gbr_20010201_velocities_GEO_v7.dat.
        These ascii files are deleted by program SuperDARN_Hybrid_EOF_Analysis.m after they are loaded into the MATLAB workspace.
    
    Fitted line-of-sight velocity files in MATLAB binary format, in the quasi-dipole coordinate frame, saved to [data_directory]/Velocities/[radar name]/: 
        These are daily files of line-of-sight velocity values, after the SuperDARN measurement locations and velocities are transformed into the quasi-dipole coordinate frame.  The daily files have the following format: 
        [radar name]_[year][month][day]_velocities_QD_[version identifier].dat, e.g.: gbr_20010201_velocities_QD_v7.mat, where QD indicates that these line-of-sight velocities are transformed to be valid in the quasi-dipole frame, and [version identifier] is a string indicating the program settings used to process the SuperDARN velocities: it is v7 by default, and can be manually set in program SuperDARN_Hybrid_EOF_Analysis.m.
    
    Binned data files of line-of-sight velocity, saved to [data_directory]/Binned_Velocities/: 
        These monthly files contain the line-of-sight velocity data after it is binned in regular temporal, spatial (quasi-dipole magnetic latitude and magnetic local time), and line-of-sight bins.  The variables contained in these files (bin_data_vel, bin_time_centroids, bin_fiducial_indices, and bin_angular_partitions) are described fully in the comments of program SuperDARN_Hybrid_EOF_Analysis.m.  The filenames have the following format: 
        [year]_[month]_[run-identifier string for the binned data].mat, e.g. 2001_02_SD_v10_NPC_AP_UTw5_APw6_rQDv7_BinnedData.mat, where the run-identifier string has the following meanings: SD: SuperDARN.  v10: a string indicating the binning and EOF settings used to process the SuperDARN velocities by program SuperDARN_Hybrid_EOF_Analysis.m.  NPC: indicates that the analysis spans the north polar cap region.  AP: 'all partitions', signifying that all angular partition bins of the radar's look directions have been used when forming the EOF analysis basis in time, space, and look-direction.  UTw5: signifies that the temporal bins used to form the irregularly-sampled SuperDARN data into a regular temporal basis for the EOF analysis are each 5 mins in length.  APw6: signifies that the width of the angular partitions (used to convert the radar look-directions into a discrete basis prior to the EOF analysis) is 6 degrees per partition.  rQDv7: a version identifier for the coordinate system rotation  approach used in program SuperDARN_Hybrid_EOF_Analysis.m.
    
    Files of EOF analysis outputs, saved to [output_directory]: 
        Program SuperDARN_Hybrid_EOF_Analysis.m runs two separate EOF analyses per calendar month, and the outputs from each model are saved.  Each EOF analysis output filename contains a run-identifier string with the following format: SD_[binning and EOF version identifier]_NPC_AP_[EOF analysis mode type]_Lanczos10s1r_10m35x_MRt_UTw5_APw6_rQD[rotation version identifier]_ds08_[model type].  The run identifier strings for the two EOF analyses are SD_v10_NPC_AP_Tmode_Lanczos10s1r_10m35x_MRt_UTw5_APw6_rQDv7_ds08_WtABn, and SD_v10_NPC_AP_Tmode_Lanczos10s1r_10m35x_MRt_UTw5_APw6_rQDv7_ds08_WSiABy, which have the following meanings: 
		    SD: SuperDARN.
			v10: a version identifier for the EOF and binning procedure applied in program SuperDARN_Hybrid_EOF_Analysis.m.
			NPC: indicates that the analysis spans the north polar cap region.
			AP: 'all partitions', signifying that all angular partition bins of the radar's look directions have been used when forming the EOF analysis basis in time, space, and look-direction.
			Tmode: the EOF analysis mode type refers to how the covariance matrix is constructed from a binned data matrix (call it X) which has dimensions [times by spatial parameters].  In an S-mode EOF analysis, we form the covariance matrix with X^T*X, hence the covariance matrix has dimensions [space by space]. In a T-mode EOF analysis, we form the covariance matrix with X*X^T, hence covariance matrix has dimensions [time by time].  The S-mode/T-mode terminology is described in Bjornsson and Venegas, 1997 (Björnsson, H., and S. A. Venegas. "A manual for EOF and SVD analyses of climatic data." CCGCR Report 97.1 (1997): 112-134) and Richman, 1986 (Richman, Michael B. "Rotation of principal components." Journal of climatology 6.3 (1986): 293-335).
			Lanczos10s1r: a Lanczos solver is used to significantly speed up the EOF analysis time: rather than solving for several thousand eigenvectors at each iteration, we can solve for 10 (hence '10s'), and return just 1 eigenvector with the highest eigenvalue (hence '1r'), which is sufficient for the infill to converge in amplitude (within acceptable loss bounds).
			10m35x: we solve for 10 modes overall, and we iterate each one 35 times.
			MRt: Mean removed along temporal direction of binned data prior to EOF analysis.
			UTw5: Universal Time width of 5 minutes: signifies that the temporal bins used to form the irregularly-sampled SuperDARN data into a regular temporal basis for the EOF analysis are each 5 mins in length.
			APw6: signifies that the width of the angular partitions (used to convert the radar look-directions into a discrete basis prior to the EOF analysis) is 6 degrees per partition.
			rQDv7: a version identifier for the coordinate system rotation approach used in this program.
			ds08: shorthand for the post-binning data selection applied to the data prior to EOF analysis.
			Wt/Wsi: Wt are temporal weights, and Wsi are inverse spatial weights, which are applied to the data matrix before EOF analysis and removed thereafter.  The temporal weights have a value of 1 where a given row (i.e. set of spatial parameters) of the binned data matrix is fully-filled (i.e. has been fully measured by SuperDARN), and 0 where that row is empty.  This will act to decrease the relative contribution of poorly-sampled epochs to the EOF solution for the variance of the data.  The inverse spatial weights have a value of 1 where a given column (i.e. time series) of the binned data matrix is empty (i.e. contains no SuperDARN measurements), and 0 where that column is fully filled.  This weights the data to increase the relative contribution from poorly sampled locations.
			ABy/ABn: defines whether amplitude boosting is applied or not.  Amplitude boosting is the name given to the procedure of regressing the (weaker) EOF infill solution onto the (stronger) original data after the infill iterations stages have been performed, in order to determine a multiplicative factor for the final infill values.  We apply amplitude boosting the EOF analysis for which the poorly sampled locations were up-weighted -- i.e. down-weighting the well-sampled locations, because it is likely that in this case the iterative infill (which starts with an initial infill of zero, hence is always of lower amplitude than the data itself) will not otherwise converge in amplitude with the original data in the specified number of iterations.
		There are multiple EOF outputs files saved by program SuperDARN_Hybrid_EOF_Analysis.m, as follows: 
            - Unitless temporal eigenvectors: these are the de-weighted temporal eigenvectors output by the EOF analysis.  Filename format: [output_directory]/[year]-[month]_[run identifier string]M1TemporalEigVecs_Deweighted.mat.
            - Unitless spatial eigenvectors: these are the de-weighted spatial eigenvectors output by the EOF analysis.  Filename format: [output_directory]/[year]-[month]_[run identifier string]M1SpatialEigVecs_Deweighted.mat.
            - Temporal eigenvectors with units of nanotesla: these are the de-weighted temporal eigenvectors which have been altered such that they have units of nanotesla and can be plotted as time series of amplitudes of each mode of variability.  Filename format: [output_directory]/[year]-[month]_[run identifier string]M1TemporalEigVecs_Deweighted_Xunits.mat.
            - Spatial eigenvectors with normalised amplitudes: these are the de-weighted spatial eigenvectors which have been altered such that their maximum amplitude is 1 (either positive or negative).  These can be multiplied with the temporal eigenvectors with units of nanotesla to return the nanotesla contribution of a particular mode of variability to the full dataset at any epoch.  Filename format: [output_directory]/[year]-[month]_[run identifier string]M1SpatialEigVecs_Deweighted_Normed.mat.
            - The means of the binned data: these are the (default temporal) means of the binned data matrix, removed prior to the EOF analysis.  They should be added back on to any reconstruction of the binned data based on the EOF modes.  Filename format: [output_directory]/[year]-[month]_[run identifier string]RemovedMeans.mat
    
    Files of hybrid reanalysis datasets, saved to [output_directory]:
        Program SuperDARN_Hybrid_EOF_Analysis.m applies a sinusoid fitting procedure to the spatial eigenvectors of each of the EOF models in a given month, then determines (for each spatial location) which of the models best recovers the original data variability.
		    The model choice values are saved out for each calendar month, and comprise a vector of size [spatial bins by 1], the values of which are either '1' (indicating that the 'WtABn' model has the best recovery of the original data for this location for this month), or '2' (indicating that the 'WsiABy' model has the best recovery of the original data for this location for this month). The model choice vector filename has the format: [output_directory]/[year]-[month]_HybridModelChoiceIndices.mat.  
		    The best performing models are combined in a hybrid reanalysis, of size [all times in one month by all spatial locations], for the north and east plasma velocity vector components.  The hybrid data matrix filename has the format: [output_directory]/[year]-[month]_HybridModelReanalysisDataset.mat.  The format of the variables saved within this file are described in the code comments of program SuperDARN_Hybrid_EOF_Analysis.m.


Instructions for running the hybrid EOF analysis: 
    
    Program SuperDARN_Hybrid_EOF_Analysis.m is designed to run on the British Antarctic Survey high-performance cluster of workstations, and was developed in MATLAB R2012a (although it appears to run fine in later versions too).  It is set up to run using the same environment, i.e. it is not trivially portable to another machine.  If running on another machine, the SuperDARN Radar Software Toolkit (RST) will need to be installed: RST 4.5 (the version used in this analysis) is available at https://zenodo.org/record/4435297#.YsWlIXbMKUk, along with installation instructions.
    Prior to running program SuperDARN_Hybrid_EOF_Analysis.m, the user should take the following steps:
        Link to the RST libraries by running the following commands. To avoid doing this every time, add these commands to the user's .cshrc file:
            setenv RSTPATH ~superdarn/rst-4.5
            source "$RSTPATH/.profile.tcsh"
        The MATLAB script was developed in bash shell. Before starting MATLAB, change the active shell to be bash (if not already set) by running the following command:
            bash
        If MATLAB has not previously been started on this login and in this shell, you will need to run these two commands to load, and then start, MATLAB:
            module load hpc/matlab/R2021a
            matlab &
        If required, alter the file permissions of the scripts in [code_directory]/subroutines/BASH/ and the C program executables in [code_directory]/subroutines/C/ using the 'chmod +x [filename]' command to allow executable permission.
        The subroutine and metadata subdirectories (described above) for program SuperDARN_Hybrid_EOF_Analysis.m need to be added to the MATLAB path.
        As stated above, the wildcards 'code_directory', 'data_directory', and 'output_directory' are directory paths on the British Antarctic Survey servers, and these variables are already assigned their appropriate directory paths within program SuperDARN_Hybrid_EOF_Analysis.m.  The user may redefine them within that program if required.  If so, then the following (empty) subdirectories must be created within 'data_directory': [data_directory]/Binned_Velocities, and [data_directory]/Velocities.  These are used to store interstitial files produced during the analysis procedure.
    Following these steps, the user should alter the variables in cell “%% Options: state which calendar months to run for, and where the data are stored:” of program SuperDARN_Hybrid_EOF_Analysis.m to specify the calendar months which will be processed, and the directory locations in which the data are stored (on the BAS servers and the SuperDARN Area Network), and where the outputs will be saved to.


@author: robore@bas.ac.uk. Robert Shore, ORCID: orcid.org/0000-0002-8386-1425.






