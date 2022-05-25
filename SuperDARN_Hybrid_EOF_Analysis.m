%% Introduction: 
%This program was created on 2022-05-16 in order to serve as the basis for
% the SuperDARN EOF analysis code handover.
%!!!! input more detailed metadata about inputs and outputs later.

%Datasets to describe: 
% hardware file for each radar
% centroid positions and metadata
% autocorrelation files
% fitted autocorrelation files
% output ascii fitted velocities

%!!!! description needs completing!

%% Options:

%Specify start and end year of span of data to process:
start_year = 2001;
end_year = 2001;

%Specify start-month and end-month of data to process: this span will apply
% uniformly to each year specified above:
start_month = 09;
end_month = 09;

%Define directory names:
%!!! reinstate!
SuperDARN_area_network_directory = '/local/users/robore/shortcut_SD_data';%this is a symlink to //samba.nerc-bas.ac.uk/data/superdarn. Needs no terminating slash! This is where the SuperDARN data are stored at BAS.
data_directory = '/local/users/robore/Data/SuperDARN_Data/';
output_directory = '/local/users/robore/Data/SuperDARN_Results/';
code_directory = '/local/users/robore/Code/JH/SHEAR/';%SuperDARN Hybrid EOF Analysis Repo.

%%!!!! used for testing in Windows:
% data_directory = 'C:/Users/robore/BAS_Data/SuperDARN_Data/';
% code_directory = 'C:/Users/robore/BAS_Files/Research/Code/JH/SuperDARN_Hybrid_EOF_Analysis_repo/';

%Version indicators: 
% Define version of code used to rotate to Quasi-Dipole coordinates:
rotation_version = 'v7';
binning_and_EOF_version = 'v10';

%Specify geometric ratios:
%Degrees to radians:
rad = pi/180;
%Radians to degrees:
deg = 180/pi;

%--------------------------------------------------------------------------
%Define data binning parameters:

%What binning UT width to use when binning the data?
bin_UT_width = 5;%units of minutes.

%Define the angular bin width. All beams from all radars which fall within
% this range of angular deviations from QD north will be considered to
% point in the same direction. It also includes the directions
% antiparallel to this bin, since they entail the same information, with a
% sign change.
partition_angular_width = 6;%units of degrees: this is the value within which the beams will be clumped together

%Which angular partitions will you select when performing the EOF analysis?
%!!!!partition_choice = 'AP';%AP means All Partitions.

%--------------------------------------------------------------------------
%Data selection settings to be applied before data binning:

%What values of range gate will you allow? For reference, the values span 1
% to 'n'.
range_gate_allowed_min = 11;%this number and up will be used.
range_gate_allowed_max = 150;%this number and lower will be used.

%--------------------------------------------------------------------------
%Data selection settings to be applied after data binning:

%What data selection are you applying, in addition to the minimum power
% level of 3dB, ground scatter flag inactive, quality flag active, and
% certain range gates (specified above) removed?
ds = '08';%Shorthand for: bins with centroid latitude of less than 59.375 were 'NaNned', values below 50ms-1 were 'NaNned', polar bin 'NaNed'.

%--------------------------------------------------------------------------
%Define EOF analysis options:

%How many iterated modes to compute:
EOF_options.number_of_iterated_modes = 10;

%How many iterations to apply per mode:
EOF_options.number_of_iterations_per_mode = 35;

%How many Lanczos solver basis functions to compute in each iteration:
EOF_options.number_of_lanczos_basis_functions_solved = 10;

%How many basis functions of those solved should the Lanczos solver return:
EOF_options.number_of_lanczos_eigenvectors_returned = 1;%will only work for 1, since we're not sorting the output by eigenvalue.

%What dimension should the bin data means be removed from in order to
% centre the data matrix? Note: should generally match the analysis mode:
EOF_options.postbin_mean_removal_dimension = 't';

%What analysis mode to use when forming the covariance matrix:
EOF_options.analysis_mode = 'Tmode';

%What exponent to apply to the weights after they have been calculated? Note
% that any factor other than 1 needs a good basis behind it:
EOF_options.additional_weights_exponent = 1;






%!!!! this needs to be iterated over: 
%What dimension of the data matrix should the weights (if used) be applied to?
%EOF_options.weights_dimension = 't';%'s', 'si', 't', or 'none'. 't' means that the weights were calculated based on the count of data over all spatial bins, for each epoch.

%!!!! this needs to be iterated over: 
%Do you apply amplitude boosting?
%EOF_options.amplitude_boosted = 'n';%'y' or 'n'.


%!!!! get these right, and take out the useless parts!
%Concatenate run-identifier string:
ri_string = ['SD_v09_NPC_' partition_choice '_' EOF_options.analysis_mode '_Lanczos' num2str(EOF_options.number_of_lanczos_basis_functions_solved) 's' ...
    num2str(EOF_options.number_of_lanczos_eigenvectors_returned) 'r_' num2str(EOF_options.number_of_iterated_modes) 'm' num2str(EOF_options.number_of_iterations_per_mode) 'x' ...
    '_MR' EOF_options.postbin_mean_removal_dimension '_UTw' num2str(bin_UT_width)...
    '_APw' num2str(partition_angular_width) '_r' IF_rotation_program_version '_ds' ds '_W' EOF_options.weights_dimension 'Exp' num2str(EOF_options.additional_weights_exponent) ...
    '_AB' EOF_options.amplitude_boosted];
%Meaning:
%SuperDARN
%Binning program v09 (in the QD program set).
%North Polar Cap
%Using one orthogonal pair of directions: the first one starting from north.
%T-mode EOF analysis.
%10 modes solved for with the Lanczos solver method, and 1 returned to form the infill basis.
%Mean removed along temporal direction.
%Number of modes iteratively infilled specified via convergence to data RMS amplitudes but with minimum and maximum of 10 modes processed.
%Number of infill iterations specified via convergence to data RMS amplitudes but with minimum and maximum iteration count set to 35.
%UT width of each bin is 5 minutes.
%Angular partition width is 6 degrees.
%Version of program used to rotate to SM is 'v3'.
%Post-binning data selection approach is '##'.
%Weights applied to the '#' dimension, with exponent '#'.

%!!!! get these right, and take out the useless parts!











%% State list of radar acronyms:
%Obtained from http://vt.superdarn.org/tiki-index.php?page=Radar+Overview,
% additionally saved in 'Z:\SuperDARN_Data\Metadata\Radar_locations_and_Acronyms_from_VTechwebsite_Table.txt'.
%The order is by their colatitude (north first), then alphabetically.
%Obtained in 2016.

radar_identifiers = {...
    'ade';...
    'adw';...
    'bks';...
    'cve';...
    'cvw';...
    'cly';...
    'fhe';...
    'fhw';...
    'gbr';...
    'han';...
    'hok';...
    'hkw';...
    'inv';...
    'kap';...
    'ksr';...
    'kod';...
    'pyk';...
    'pgr';...
    'rkn';...
    'sas';...
    'sch';...
    'sto';...
    'wal';...
    'bpk';...
    'dce';...
    'fir';...
    'hal';...
    'ker';...
    'mcm';...
    'san';...
    'sps';...
    'sye';...
    'sys';...
    'tig';...
    'unw';...
    'zho';...
    };

%% Load equal area bin layout to apply to each month looped over:
%We will eventually bin the SuperDARN data in these spatial bins.

load([code_directory 'metadata/bin_coordinates/north_polar_region_equal_area_bins.mat'])
%Brings in:
% bin_coords_colat: upper and lower co-latitude limits of each equal area bin. Size of matrix is [bins by 2].
% bin_coords_colong: upper and lower longitude limits of each equal area bin. Size of matrix is [bins by 2].
% bin_centroids_colat: co-latitude coordinates of the centroid of each equal area bin. Size of matrix is [bins by 1].
% bin_centroids_colong: longitude coordinates of the centroid of each equal area bin. Size of matrix is [bins by 1].

%Create index of bin numbers, for later use in removing empty bins:
index_original_bin_fiducials = (1:size(bin_coords_colat,1))';%size [NPC bins by 1].

%% Check for existence of radar hardware files, download from the SuperDARN website if required:

%Loop over radar acronyms, and if required, download the hardware file for
% a given radar from the SuperDARN website.
%This loop is placed before the 'main' loop over the radars since the C
% program centroid_calculation.c depends on the files being already
% downloaded.
for i_r = 1:size(radar_identifiers,1)
    if(~exist([code_directory 'metadata/hardware_files/hdw.dat.' radar_identifiers{i_r,1}],'file'))
        urlwrite(['http://davit1.ece.vt.edu/hdw/hdw.dat.' radar_identifiers{i_r,1} '.txt'], [code_directory 'metadata/hardware_files/hdw.dat.' radar_identifiers{i_r,1}]);
        disp(['Hardware file for radar ' radar_identifiers{i_r,1} ' downloaded from web: check it for formatting errors.'])
    end%Conditional: check for file existence, and download it if not found.
end%Loop over all radars, save hardware file for each.
%Important: these files were originally downloaded on 2021/03/19,
% and required subsequent manual edits to whitespace and carriage
% returns for the following radars' hardware files: cve, cvw, cly,
% hkw.

%Display progress indicator:
disp('Hardware files have been found or obtained.')

%% Compute centroid metadata for all radars, using externally-run C code and website-sourced hardware file data:
%This cell is outside the loop over radars since it works for all of them
% at once.

%The program 'centroid_calculation.c' does the following:
%Compute and save one set of centroid metadata (linked to the epochs at
% which they should be superseded by more recent information) for each
% line in each radar's hardware file, in geographic geocentric
% coordinates.

%Using radar 'ade' and its first hardware file line as a test-flag, check
% whether the centroid data exists for all radars, and make it for all
% radars if it does not:
if(~exist([code_directory 'metadata/centroid_locations/centroid_locations_ade_1.dat'],'file'))
    %The program that makes the centroid files is called 'U:\Research\Code\JH\C_Programs\centroid_calculation.c'.
    
    %Compile the C code in directory with the code in it:
    system(['cd ' code_directory 'subroutines/C/ ; make -f makefile_centroid_calculation']);
    
    %Run the C program in directory with the code in it:
    system(['cd ' code_directory 'subroutines/C/ ; ./centroid_calculation']);
    %Outputs centroid positions to [code_directory 'metadata/centroid_locations/'].
    
    %State progress:
    disp('Centroid positions synthesised for all radars.')
else
    disp('Found ready-processed centroid locations.')
end%Conditional: if the centroid location files do not exist, make them using the C code stated.

%% Compile the C program used to parse the fitted velocity binary files into ascii format:

%Compile the C code in directory with the code in it:
system(['cd ' code_directory 'subroutines/C/ ; make -f makefile_fit_file_parse']);
system(['cd ' code_directory 'subroutines/C/ ; make -f makefile_fitacf_file_parse']);

%% Loop over all years and months, and for each calendar month, loop over all radars, then process each month's set of files for that radar in turn:

%Loop over all specified years:
for i_year = start_year:1:end_year
    %% Loop over all specified months for this year:
    
    %Loop over all specified months in current year:
    for i_month = start_month:1:end_month
        %State the year and month that you're processing now:
        disp(['Processing ' num2str(i_year) '-' num2str(i_month,'%2.2d') '.'])
        
        %% Check for existence of binned data file, and bin the data if it is not found:
        
        %Concatenate run-identifier string for binned data, based on
        % options specified at the start of the program:
        binned_data_ri_string = ['SD_' binning_and_EOF_version '_NPC_AP_UTw' num2str(bin_UT_width) '_APw' num2str(partition_angular_width) '_rQD' rotation_version '_BinnedData'];
        %E.g. SD_v10_NPC_AP_UTw5_APw6_rQDv7_BinnedData, which means:
        % SD: SuperDARn.
        % v10: a version identifier for the outputs of this program.
        % NPC: north polar cap region.
        % AP: 'all partitions'. This signifies that all angular partitions
        %     of the radar's look directions have been used when forming
        %     the EOF analysis basis in time, space, and look-direction.
        % UTw5: signifies that the temporal bins used to form the
        %       irregularly-sampled SuperDARN data into a regular temporal
        %       basis for the EOF analysis are each 5 mins in length.
        % APw6: signifies that the width of the angular partitions (used to
        %       convert the radar look-directions into a discrete basis
        %       prior to the EOF analysis) is 6 degrees per partition.
        % rQDv7: a version identifier for the coordinate system rotation
        %        approach used in this program.
        
        %Check for existing binned data for this month and read it in if found:
        binned_data_filename = [data_root 'SuperDARN_Data/Binned_Velocities/' num2str(i_year) '_' num2str(i_month,'%2.2d') '_' binned_data_ri_string '.mat'];
        
        if(exist(binned_data_filename,'file'))
            load(binned_data_filename);
            %Source: this program.
            %Brings in: %!!!!
        else
            %% Define epochs for this month:
            %Define the number of the last day in this month to account for
            % leap years and variable month lengths:
            month_end_day = eomday(i_year,i_month);%size [1 by 1].
            
            %Define month-centre epoch in MJD2000:
            month_mean_epoch_MJD2000 = mean([(datenum([i_year i_month 1 0 0 0]) - 730486); (datenum([i_year i_month month_end_day 0 0 0]) - 730486)]);%size [1 by 1].
            
            %% Loop over all radars and for each: fit autocorrelation files, convert C binary outputs to ascii, rotate data to QD frame and save as MATLAB binary files, then store each data file in this month:
            
            %Set the 'previous day's' index for the applicable hardware
            % file line to some null value to force the computation of the
            % offset angles in a later cell. This must be done prior to the
            % loop(s) over the days of the month, just in case any are
            % skipped due to missing data:
            index_last_applicable_hardware_line_PrevDay = 0;
            
            %Loop over all radars, compute angular offsets from QD north, and
            % store the result:
            each_radar_velocity_data_cell = cell(size(radar_identifiers,1),1);%size [radars by 1].
            each_radar_monthly_data_count = NaN(size(radar_identifiers,1),1);%size [radars by 1].
            for i_r = 1:size(radar_identifiers,1)
                disp(['Processing radar ' radar_identifiers{i_r,1} ' for ' num2str(i_year) '-' num2str(i_month)]);
                
                %% Parse the hardware file for this radar:
                %This is somewhat inefficient as the load-in doesn't really
                % need to be repeated every calendar month, but I could not be
                % bothered to set up a matrix structure and index it later: it
                % is safer to just repeat the load-in.
                
                %Load the radar hardware file contents into a single cell:
                fid_hardware = fopen([code_directory 'metadata/hardware_files/hdw.dat.' radar_identifiers{i_r,1}]);
                hardware_data_cell_combined = textscan(fid_hardware, '%s', 'Delimiter', '\n');%size [1 by 1]. This is a cell array in which each element is a cell array that contains a string.
                fclose(fid_hardware);
                
                %Using approach from https://uk.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string,
                % split the single cell into an array of cells, each of which contains
                % a single row-line of the hardware file text:
                hardware_data_cell_divided = [hardware_data_cell_combined{:}];%cell, size [number of hardware file lines (headers included) by 1].
                %Note that the parse function appears clever enough not to return the
                % last blank line in each hardware file as an entry in the cell.
                
                %Using the character '#' as indicative of a header line (i.e. row of
                % the cell array), define an index of the lines which contain this
                % character anywhere in that line os text.  Given the vagaries of cells
                % in Matlab, this index will be a cell contianing 1 where there is a
                % header line and it will be an empty cell where there is a line
                % containing data, so we convert it to an index of fiducials:
                index_header_lines_cell = strfind(hardware_data_cell_divided, '#');%cell, size [number of hardware file lines (headers included) by 1].
                index_header_lines_fid = find(not(cellfun('isempty', index_header_lines_cell)));%converted to fiducials, size [(number of hardware file rows, header lines only) by 1].
                
                %Remove the cells which contain headers:
                hardware_data_cell_divided(index_header_lines_fid) = [];%this approach removes the indexed cells entirely.  New size of cell: [(number of hardware file data rows) by 1].
                
                %The data values have so far been parsed in string format, so convert
                % them to double and concatenate the records:
                hardware_data_lines = NaN(size(hardware_data_cell_divided,1),19);%size [(number of time-based hardware file entries) by (number of hardware data columns, known a priori)].
                for i_hardware_line = 1:size(hardware_data_cell_divided,1)
                    hardware_data_lines(i_hardware_line,:) = str2num(hardware_data_cell_divided{i_hardware_line,1})'; %#ok<ST2NM>
                end%Loop over each line of the hardware file metadata, and concatenate that information in double precision format.
                
                %Format of columns of radar_hardware_data_web, taken from the hardware
                % file header, in which a more detailed description is given:
                % - 1: Station ID (unique numerical value).
                % - 2: Last year that parameter string is valid. (4 digit year).
                % - 3: Last second of year that parameter string is valid (range 0 to
                %       34163999 for non-leap years). The parameter string giving the
                %       current configuration is assumed to be valid until the last
                %       second of 2999.
                % - 4: Geographic latitude of radar site.
                % - 5: Geographic longitude of radar site.
                % - 6: Altitude of the radar site (metres).
                % - 7: Scanning boresight (Direction of the center beam, measured in
                %       degrees relative to geographic north. CCW rotations are
                %       negative).
                % - 8: Beam separation (Angular separation in degrees between adjacent
                %       beams. Normally 3.24 degrees).
                % - 9: Velocity sign according to Doppler velocity. Set to +1 or -1.
                % - 10: Analog Rx attenuator step (dB).
                % - 11: Tdiff (Propagation time from interferometer array antenna to
                %        phasing matrix input minus propagation time from main array
                %        antenna through transmitter to phasing matrix input. If the
                %        signal from the interferometer comes first, then tdiff < 0.
                %        Units are in decimal microseconds).
                % - 12: Phase sign (Cabling errors can lead to a 180 degree shift of
                %        the interferometry phase measurement. +1 indicates that the
                %        sign is correct, -1 indicates that it must be flipped).
                % - 13: Interferometer offset (Displacement of midpoint of
                %        interferometer array from midpoint of main array. This is
                %        given in meters in Cartesian coordinates. X is along the line
                %        of antennas with +X toward higher antenna numbers, Y is along
                %        the array normal direction with +Y in the direction of the
                %        array normal. Z is the altitude difference, +Z up.).
                %        Presumably, this first entry is X.
                % - 14: Interferometer offset, second coordinate.
                % - 15: Interferometer offset, third coordinate.
                % - 16: Analog Rx rise time, given in microseconds. Time delays of less
                %        than ~10 microseconds can be ignored.
                % - 17: Number of analog attenuation stages.
                % - 18: Maximum number of range gates.  Used for allocation of array
                %        storage.
                % - 19: Maximum number of beams.
                
                %Convert temporal records of the hardware data lines to year, month, day
                % format:
                hardware_data_times_year = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                hardware_data_times_month = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                hardware_data_times_day = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                hardware_data_times_hour = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                hardware_data_times_minute = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                hardware_data_times_second = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                for i_hardware_line = 1:size(hardware_data_lines,1)
                    %Use the 'second-of-year' and 'year' records from the hardware file to
                    % figure out the date that it is referring to:
                    [YYYY_temp MM_temp DD_temp hh_temp mm_temp ss_temp] = datevec(hardware_data_lines(i_hardware_line,3)/86400+datenum(hardware_data_lines(i_hardware_line,2),1,1,0,0,0));
                    %Note: the full date string is given in the hardware file header lines,
                    % and the output from this matches, to the second.
                    
                    %Store just the year, month and day of the hardware file line's expiry date:
                    hardware_data_times_year(i_hardware_line,1) = YYYY_temp;
                    hardware_data_times_month(i_hardware_line,1) = MM_temp;
                    hardware_data_times_day(i_hardware_line,1) = DD_temp;
                    hardware_data_times_hour(i_hardware_line,1) = hh_temp;
                    hardware_data_times_minute(i_hardware_line,1) = mm_temp;
                    hardware_data_times_second(i_hardware_line,1) = ss_temp;
                end%Loop over each line of hardware file data, get the temporal records in a more approachable format.
                
                %Format the hardware line expiry dates in Matlab datenum format:
                hardware_data_times_datenum = datenum([hardware_data_times_year  hardware_data_times_month  hardware_data_times_day  hardware_data_times_hour  hardware_data_times_minute  hardware_data_times_second]);%size [(number of time-based hardware file entries) by 1].
                
                %% Load centroid locations for all epochs of the radar's hardware data file:
                
                %Loop over each epoch of the hardware data file for this radar and store
                % the associated geographic-coordinate centroids:
                %Preallocate geographic centroid location data storage for this radar:
                centroid_locations_beam_number = cell(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_locations_range_gate = cell(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_locations_GEO_colatitude = cell(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_locations_GEO_longitude = cell(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_parameters_max_beams = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_parameters_max_gates = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_parameters_frang = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                centroid_parameters_rsep = NaN(size(hardware_data_lines,1),1);%size [(number of time-based hardware file entries) by 1].
                for i_hardware_line = 1:size(hardware_data_lines,1)
                    %Check if the applicable set of centroid locations exists for this
                    % radar, which it should:
                    if(~exist([code_directory 'metadata/centroid_locations/centroid_locations_' radar_identifiers{i_r,1} '_' num2str(i_hardware_line) '.dat'],'file'))
                        disp(['Missing centroid location file for ' radar_identifiers{i_r,1}  '.  Run instructions for centroid_calculation.c given in code of this program.'])
                        return
                        %Run instructions:
                        % The program that makes the centroid files is called 'centroid_calculation.c'.
                        % Compile command: make -f makefile_centroid_calculation
                        % Run command: ./centroid_calculation
                        %The cell 'Compute centroid metadata (...)' should do all this.
                    end%Conditional: if the centroid location file does not exist, make it using the C code stated.
                    
                    %Load the centroid data for the corresponding hardware file line:
                    centroid_data_bulk = importdata([code_directory 'metadata/centroid_locations/centroid_locations_' radar_identifiers{i_r,1} '_' num2str(i_hardware_line) '.dat'],' ');
                    %Source: U:\Research\Code\JH\C_Programs\centroid_calculation.c,
                    % described in an above cell of this program.
                    %Format of centroid_data_bulk:
                    %Top row of 4 columns: The maximum count of beam numbers and range gates, the expiry year and second of year.
                    %Second row of four columns: The vhm type (integer), the 'height' variable value put into RPosGeo (units unknown), the frang and rsep values (both in km).
                    %The other rows of columns are: beam number, range gate number, latitude of centroid and longitude of centroid.
                    %The rows are ordered by beam number, and within that, range gate number.
                    %Store variables for the centroid locations:
                    centroid_locations_beam_number{i_hardware_line,1} = centroid_data_bulk(3:end,1);%stored data is size [centroids by 1]. Range is 0 to 15.
                    centroid_locations_range_gate{i_hardware_line,1} = centroid_data_bulk(3:end,2);%stored data is size [centroids by 1]. Range is 1 to max range gate count.
                    centroid_locations_GEO_colatitude{i_hardware_line,1} = 90 - centroid_data_bulk(3:end,3);%stored data is size [centroids by 1].
                    centroid_locations_GEO_longitude{i_hardware_line,1} = centroid_data_bulk(3:end,4);%stored data is size [centroids by 1].
                    
                    %Store variables for the parameters used to make the centroid
                    % locations:
                    centroid_parameters_max_beams(i_hardware_line,1) = centroid_data_bulk(1,1);%maximum number of beams for this radar. Stored data is size [1 by 1].
                    centroid_parameters_max_gates(i_hardware_line,1) = centroid_data_bulk(1,2);%maximum number of range gates for this radar. Stored data is size [1 by 1].
                    centroid_parameters_frang(i_hardware_line,1) = centroid_data_bulk(2,3);%distance to first range gate, given in km. Stored data is size [1 by 1].
                    centroid_parameters_rsep(i_hardware_line,1) = centroid_data_bulk(2,4);%distance between range gates, given in km. Stored data is size [1 by 1].
                    
                end%Loop over each epoch of the hardare data file for this radar.
                
                %% For the present radar, convert the month of autocorrelation files to fitted velocities format using external BASH script:
                
                %Define which BASH script to run, based on the date at which
                % the .dat files are formatted as .rawacf files:
                if(i_year < 2006)
                    file_format_old_or_new = 'old';
                elseif(i_year > 2006)
                    file_format_old_or_new = 'new';
                elseif(i_year == 2006)
                    if(i_month <= 6)
                        file_format_old_or_new = 'old';
                    elseif(i_month >= 7)
                        file_format_old_or_new = 'new';
                    end%Conditional: file format switches from old to new between months 6 and 7 of 2006.
                end%Conditional: file format switches from old to new in 2006.
                
                %Run a BASH script for the timespan under consideration to
                % convert the files (of which there may be several in a given
                % day) to fitted format:
                %For now, we know that these are .dat files, but later we will have the
                % routine check which BASH script to call:
                if(strcmp(file_format_old_or_new,'old'))
                    %Change to the directory where the SuperDARN data is
                    % expected to be stored and run fit_script_v6:
                    system(['cd ' code_directory 'subroutines/BASH/ ; ' ...
                        './fit_script_v6 ' ...
                        '-sd 01/' num2str(i_month,'%2.2d') '/' num2str(i_year) ' ' ...
                        '-ed ' num2str(month_end_day,'%2.2d') '/'  num2str(i_month,'%2.2d') '/' num2str(i_year) ' ' ...
                        '-s ' radar_identifiers{i_r,1} ' ' ...
                        '-in ' SuperDARN_area_network_directory ' ' ...
                        '-out ' data_directory 'Velocities']);
                    %Only works for Matlab run in UNIX/Linux environment.
                elseif(strcmp(file_format_old_or_new,'new'))
                    %Change to the directory where the SuperDARN data is
                    % expected to be stored and run fitacf_script_v6:
                    system(['cd ' code_directory 'subroutines/BASH/ ; ' ...
                        './fitacf_script_v6 ' ...
                        '-sd 01/' num2str(i_month,'%2.2d') '/' num2str(i_year) ' ' ...
                        '-ed ' num2str(month_end_day,'%2.2d') '/'  num2str(i_month,'%2.2d') '/' num2str(i_year) ' ' ...
                        '-s ' radar_identifiers{i_r,1} ' ' ...
                        '-in ' SuperDARN_area_network_directory ' ' ...
                        '-out ' data_directory 'Velocities']);
                    %Only works for Matlab run in UNIX/Linux environment.
                end%Conditional: apply different BASH routine dependent on data file type.
                
                %% Start loop over the month's daily fitted velocity C binary files for this radar to put them in ascii format:
                %In each month, we loop over the days in the month twice:
                % firstly to put the C binary files into ascii format, and
                % secondly to read the data from those ascii files into the
                % MATLAB workspace. The ascii files are output by a C code
                % routine, so the two loops allow time for a file to be written
                % before the program needs to read it in again.
                
                %Make a checker to see if there are any data files for this
                % month: this will allow us to skip the month later if there are
                % none.
                any_data_in_month = 0;
                for i_day = 1:1:end_day
                    %% Loop over a wildcard-formed list of all the files for this radar for this day:
                    
                    %Get Matlab to wildcard-list all the files in a directory:
                    % this will help you loop through any sub-daily files, the
                    % filenames of which cannot be predicted easily, because
                    % they depend on the format of the pre-fitting '.dat' or
                    % '.rawacf' fles, which have one of these formats:
                    % $(YEAR)$(MONTH)$(DAY).$(STATION).${FILETYPE}
                    % or
                    % $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).${FILETYPE}
                    % or
                    % $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).$(CHANNEL).${FILETYPE}
                    %So here, the first '*' accounts for the fact that these
                    % may be daily or subdaily (generally 2-hourly) files (i.e.
                    % it can either be '.' or '.$(HOUR)$(MINUTE).$(SECOND).').
                    % The second '*' accounts for the tendency of some
                    % filenames to have a single-character channel identifier,
                    % which is absent for stations which only have a single
                    % channel (i.e. it can either be '.' or '.$(CHANNEL).').
                    %We constrain the search to just the files we want by
                    % appending the correct filetype, i.e. '.fit' or '.fitacf':
                    if(strcmp(file_format_old_or_new,'old'))
                        single_day_acf_filenames = dir([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '*' radar_identifiers{i_r,1} '*fit']);
                    elseif(strcmp(file_format_old_or_new,'new'))
                        single_day_acf_filenames = dir([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '*' radar_identifiers{i_r,1} '*fitacf']);
                    end%Conditional: filename differs dependent on old/new content format.
                    
                    %Get the actual filenames from the structure/object that is
                    % the variable 'single_day_acf_filenames':
                    single_day_acf_filenames = {single_day_acf_filenames.name};%cell, size {1 by (number of sub-daily filenames in this day)}.
                    
                    %Check for there being no data to iterate over in this day:
                    if(isempty(single_day_acf_filenames))
                        continue%pertains to loop over i_day.
                    end%Conditional: if there's no file to read in, skip the processing of it.
                    
                    %If there are data in at least one day in this month, you
                    % can (at this point) set the 'is there any data in the
                    % month' checker value to 1. this can be redefined for
                    % other days in the month if required:
                    any_data_in_month = 1;%!!!! are we stil using this in the new program?
                    
                    %Loop over the set of sub-daily filenames:
                    for i_daily_file = 1:max(size(single_day_acf_filenames))
                        %% Retrieve daily/sub-daily C-binary filename:
                        
                        %Create string for filename of single daily or
                        % sub-daily fitted autocorrelation data file:
                        single_acf_filename = single_day_acf_filenames{i_daily_file};%Here, the filename is like '20170918.1243.43.gbr.fitacf', or '20101218.gbr.fitacf'.
                        SD_fitted_binary_filename = [data_directory 'Velocities/' radar_identifiers{i_r,1} '/' single_acf_filename];%here, I have added on the full filepath.
                        %Source of file: 'fit_script_v6' or 'fitacf_script_v6', run in this program.
                        
                        %% Use hardware epoch expiry date to determine which centroid location parameters we should be using for this daily/sub-daily file:
                        
                        %Define the epoch of this file: not sure if this is
                        % the start or the end of the file if its a subdaily
                        % file.
                        %Here, we check the filename length. The SuperDARN
                        % files can be either in 2-hourly format, or
                        % concatenated into daily files (as stored on the BAS
                        % servers). Irritatingly, this changes the filename
                        % format, and since the timestamp of the subdaily files
                        % is unknown prior to getting the file, we need to
                        % extract it from the filename, which means checking
                        % that we're extracting the right characters.
                        if(length(single_acf_filename) == 27 | length(single_acf_filename) == 24 | length(single_acf_filename) == 29 | length(single_acf_filename) == 26)
                            %Presumably, these are 2-hour length files, like
                            % '20170918.1243.43.gbr.fitacf' or
                            % '20170918.1243.43.gbr.fit' or
                            % '20170918.1243.43.fir.a.fitacf' or
                            % '20170918.1243.43.fir.a.fit'. The formats are:
                            % $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).${FILETYPE}
                            % or
                            % $(YEAR)$(MONTH)$(DAY).$(HOUR)$(MINUTE).$(SECOND).$(STATION).$(CHANNEL).${FILETYPE}
                            % Where the channel identifier is a single letter,
                            % and absent for stations which only have a single
                            % channel. Filetype is either .fit or .fitacf,
                            % stemming respectively from .dat or .rawacf files.
                            
                            %From the filename, we extract the date
                            % and the time for comparison against the hardware
                            % file cutoff epochs.
                            single_file_epoch = datenum([i_year  i_month  i_day  str2double(single_acf_filename(10:11))  str2double(single_acf_filename(12:13))  str2double(single_acf_filename(15:16))]);
                            
                            %Create string for an ascii version of the fitted
                            % C binary velocities file, with a matching
                            % filename.
                            SD_fitted_ascii_filename = [data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
                                radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') ...
                                '_' single_acf_filename(10:11)  single_acf_filename(12:13) '_' single_acf_filename(15:16) ...
                                '_velocities_GEO_' rotation_version '.dat'];
                        elseif(length(single_acf_filename) == 19 | length(single_acf_filename) == 16 | length(single_acf_filename) == 21 | length(single_acf_filename) == 18)
                            %Presumably, these are daily files, like
                            % '20101218.gbr.fitacf' or
                            % '20101218.gbr.fit' or
                            % '20101218.fir.a.fitacf' or
                            % '20101218.fir.a.fit'. The format is:
                            % $(YEAR)$(MONTH)$(DAY).$(STATION).${FILETYPE}
                            % or
                            % $(YEAR)$(MONTH)$(DAY).$(STATION).$(CHANNEL).${FILETYPE}
                            % Where the channel identifier is a single letter,
                            % and absent for stations which only have a single
                            % channel. Filetype is either .fit or .fitacf,
                            % stemming respectively from .dat or .rawacf files.
                            
                            %From the filename, we extract just the date
                            % for comparison against the hardware file cutoff
                            % epochs.
                            single_file_epoch = datenum([i_year  i_month  i_day]);
                            
                            %Create string for an ascii version of the fitted
                            % C binary velocities file, with a matching
                            % filename.
                            SD_fitted_ascii_filename = [data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
                                radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') ...
                                '_velocities_GEO_' rotation_version '.dat'];
                        else
                            disp('Warning: the fitted acf filenames output by fit*_script_v* have an unexpected length. This program cannot parse them. Terminating')
                            return
                        end%Conditional: check filename length of the fitted C binary files output by the fit*_script_v* BASH script(s).
                        
                        %For the date entailed in the C binary filename, find the
                        % latest applicable hardware file line, given that the
                        % hardware_data_times_year, etc. variables are expiry dates:
                        index_last_applicable_hardware_line = find(hardware_data_times_datenum > single_file_epoch,1,'first');%size [1 by 1].
                        
                        %% Run external C code routine to put the daily/sub-daily velocity data file in readable (ascii) format:
                        
                        %Run some basic checks against the hardware file information: this was
                        % parsed in C to make the centroids, so it is worth double-checking:
                        if((centroid_parameters_max_beams(index_last_applicable_hardware_line,1) ~= hardware_data_lines(index_last_applicable_hardware_line,19)) | ...
                                (centroid_parameters_max_gates(index_last_applicable_hardware_line,1) ~= hardware_data_lines(index_last_applicable_hardware_line,18))) %#ok<OR2>
                            disp('Error in hardware file line and centroid data locations: terminating.')
                            return
                        end%Conditional: error-check.
                        
                        %The fitted files are C binaries in which the expected values of the
                        % distance to the first range gate and the separation between centroids are
                        % both expected as two-way travel times (at the speed of light), given in
                        % microseconds.  So here we convert the km values to mew seconds:
                        % ... for frang:
                        frang_km = centroid_parameters_frang(index_last_applicable_hardware_line,1);%new attribution for ease of reference. size [1 by 1].
                        frang_km = frang_km * 2;%doubling value to account for the two way travel time. Thus, 180 becomes 360.
                        frang_metres = frang_km * 1e3;%conversion from km to metres.
                        frang_seconds = frang_metres / 299792458;%t = d / s.
                        frang_mew_seconds = frang_seconds * 1e6;%conversion from seconds to microseconds.
                        %Reduce the precision of the value by rounding to the nearest hundred:
                        frang_mew_seconds = frang_mew_seconds / 100;
                        frang_mew_seconds = round(frang_mew_seconds);%rounding to eliminate differences between line of sight and along-ground travel distances.
                        frang_mew_seconds = frang_mew_seconds * 100;%should equal 1200.
                        % ... for rsep:
                        rsep_km = centroid_parameters_rsep(index_last_applicable_hardware_line,1);%new attribution for ease of reference. size [1 by 1].
                        rsep_km = rsep_km * 2;%doubling value to account for the two way travel time. Thus, 45 becomes 90.
                        rsep_metres = rsep_km * 1e3;%conversion from km to metres.
                        rsep_seconds = rsep_metres / 299792458;%t = d / s.
                        rsep_mew_seconds = rsep_seconds * 1e6;%conversion from seconds to microseconds.
                        %Reduce the precision of the value by rounding to the nearest hundred:
                        rsep_mew_seconds = rsep_mew_seconds / 100;
                        rsep_mew_seconds = round(rsep_mew_seconds);%rounding to eliminate differences between line of sight and along-ground travel distances.
                        rsep_mew_seconds = rsep_mew_seconds * 100;%should equal 300.
                        
                        %Parameters which we need to pass to fit_file_parse.c: these were used to make the
                        % centroid positions by 'centroid_calculation.c', so it is important that
                        % they match the radar data file:
                        % - number of beams.
                        % - number of range gates.
                        % - distance to first range gate (frang).
                        % - separation of range gates (rsep).
                        % - the filename to read in
                        % - the filename to save out
                        
                        %Run a C routine which takes in the parameters under which
                        % the centroids were computed, checks them against the
                        % fitted data and outputs the velocities, times, gate and
                        % beam numbers:
                        %Change directory to be in the same location as the C
                        % program and call the C routine for the daily file:
                        if(strcmp(file_format_old_or_new,'old'))
                            system(['cd ' code_directory 'subroutines/C/ ; ' ...
                                './fit_file_parse ' ...
                                SD_fitted_binary_filename ' ' ...
                                num2str(centroid_parameters_max_beams(index_last_applicable_hardware_line,1)) ' ' ...
                                num2str(centroid_parameters_max_gates(index_last_applicable_hardware_line,1)) ' ' ...
                                num2str(frang_mew_seconds) ' ' ...
                                num2str(rsep_mew_seconds) ' ' ...
                                SD_fitted_ascii_filename]);
                        elseif(strcmp(file_format_old_or_new,'new'));
                            system(['cd ' code_directory 'subroutines/C/ ; ' ...
                                './fitacf_file_parse ' ...
                                SD_fitted_binary_filename ' ' ...
                                num2str(centroid_parameters_max_beams(index_last_applicable_hardware_line,1)) ' ' ...
                                num2str(centroid_parameters_max_gates(index_last_applicable_hardware_line,1)) ' ' ...
                                num2str(frang_mew_seconds) ' ' ...
                                num2str(rsep_mew_seconds) ' ' ...
                                SD_fitted_ascii_filename]);
                        end%Conditional: C program to call differs dependent on old/new content format.
                        
                        %% Delete the .fit files as you go:
                        
                        delete(SD_fitted_binary_filename);
                        
                    end%Loop over the daily/sub-daily fitted autocorrelation files contained within this day.
                    
                end%Loop over each day in the present month for this present radar, write the ascii files.
                
                %If there were no data files found in any day for this month,
                % then skip the processing for (the rest of) this month for
                % this radar. This effects the defnintion of the
                % index_last_applicable_hardware_line variable, which retains
                % its value between months (and if you're only running a single
                % month, between radars) if it is not redefined on at least
                % one day for this month for this radar.
                if(any_data_in_month == 0)
                    disp(['No data in this month for radar ' radar_identifiers{i_r,1} ': skipping this month.'])
                    continue%pertains to loop over i_r.
                end%Conditional: skip processing the month of data for this radar if there's no data to process.
                
                %% Make variables for the radar location in GEO coordinates:
                radar_GEO_colatitude = (90 - hardware_data_lines(index_last_applicable_hardware_line,4));%size [1 by 1].
                radar_GEO_colongitude = (hardware_data_lines(index_last_applicable_hardware_line,5) + 360);%converted to co-longitude.
                
                %% Start loop over the month's daily fitted velocity files for this radar to read the ascii format files, %!!!! (save out MAATLAB binary files)?, and aggregate the data for this radar over all days in the month:
                %The reason there are two loops over each day of data within
                % the month is that I do not trust C to finish writing the ascii
                % files before I ask Matlab to read them in,  This way, the
                % first file should be written long before I ask Matlab to read
                % it.
                
                %Preallocate storage for the data aggregated for this radar
                % over all days in the month:
                single_radar_velocity_data_cell = cell(month_end_day,1);%size [number of days in month by 1].
                single_radar_daily_data_count = NaN(month_end_day,1);%size [number of days in month by 1].
                
                %Start (second) loop over the days of the month for this
                % radar:
                for i_day = 1:1:end_day
                    %% Loop over a wildcard-formed list of all the ascii files for this radar for this day:
                    
                    %Get Matlab to wildcard-list all the files in a directory:
                    % this will help you loop through the daily/sub-daily files, the
                    % filenames of which cannot be predicted easily (because
                    % even though I just made them in this program, the hour,
                    % minute, and second elements are not predictable).
                    single_day_ascii_filenames = dir([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
                        radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') ...
                        '*velocities_GEO_' rotation_version '.dat']);
                    
                    %Get the actual filenames from the structure/object that is
                    % the variable 'subdaily_acf_filenames':
                    single_day_ascii_filenames = {single_day_ascii_filenames.name};%cell, size {1 by (number of sub-daily filenames in this day)}.
                    
                    %Check for there being no data to iterate over in this day:
                    if(isempty(single_day_ascii_filenames))
                        continue%pertains to loop over i_day.
                    end%Conditional: if there's no file to read in, skip the processing of it.
                    
                    %Preallocate storage for the daily set of all sub-daily velocity data:
                    velocity_data = [];
                    %Loop over the set of sub-daily filenames:
                    for i_daily_file = 1:max(size(single_day_ascii_filenames))
                        %% Retrieve daily/sub-daily ascii filename:
                        
                        %Create string for filename of single sub-daily fitted
                        % autocorrelation data file:
                        SD_fitted_ascii_filename = single_day_ascii_filenames{i_daily_file};%this is like 'gbr_20170928_1848_04_velocities_GEO_v7.dat', or 'gbr_20170928_velocities_GEO_v7.dat'.
                        %Source of file: 'fit_script_v6'/'fitacf_script_v6', run in this program.
                        
                        %% Use hardware epoch expiry date to determine which centroid location parameters we should be using for this sub-daily file:
                        
                        %Define the epoch of this file: not sure if this is
                        % the start or the end of the file if its a subdaily
                        % file.
                        %Here, we check the filename length. The SuperDARN
                        % files can be either in 2-hourly format, or
                        % concatenated into daily files (as stored on the BAS
                        % servers). Irritatingly, this changes the filename
                        % format, and since the timestamp of the subdaily files
                        % is unknown prior to getting the file, we need to
                        % extract it from the filename, which means checking
                        % that we're extracting the right characters.
                        if(length(SD_fitted_ascii_filename) == 42)
                            %Presumably, these are 2-hour length files, like
                            % 'gbr_20170928_1848_04_velocities_GEO_v7.dat'. We
                            % extract the date and the time for comparison
                            % against the hardware file cutoff epochs.
                            single_file_epoch = datenum([i_year  i_month  i_day  str2double(SD_fitted_ascii_filename(14:15))  str2double(SD_fitted_ascii_filename(16:17))  str2double(SD_fitted_ascii_filename(19:20))]);
                        elseif(length(SD_fitted_ascii_filename) == 34)
                            %Presumably, these are daily files, like
                            % 'gbr_19970618_velocities_GEO_v7.dat'. We
                            % use just the date for comparison against the
                            % hardware file cutoff epochs.
                            single_file_epoch = datenum([i_year  i_month  i_day]);
                        else
                            disp('Warning: the ascii versions (made by this prorgam) of fitted acf filenames output by fit*_script_v* have an unexpected length. This program cannot parse them. Terminating')
                            return
                        end%Conditional: check filename length of the fitted C binary files output by the fit*_script_v* BASH script(s).
                        
                        %For the date entailed in the sub-daily filename, find the
                        % latest applicable hardware file line, given that the
                        % hardware_data_times_year, etc. variables are expiry dates:
                        index_last_applicable_hardware_line = find(hardware_data_times_datenum > single_file_epoch,1,'first');%size [1 by 1].
                        
                        %Extract the pertinent set of centroid positions and metadata:
                        centroid_locations_beam_number_single_day = centroid_locations_beam_number{index_last_applicable_hardware_line,1};%size [centroids by 1]. Range is 0 to 15.
                        centroid_locations_range_gate_single_day = centroid_locations_range_gate{index_last_applicable_hardware_line,1};%size [centroids by 1]. Range is 1 to max range gate count.
                        centroid_locations_GEO_colatitude_single_day = centroid_locations_GEO_colatitude{index_last_applicable_hardware_line,1};%size [centroids by 1].
                        centroid_locations_GEO_longitude_single_day = centroid_locations_GEO_longitude{index_last_applicable_hardware_line,1};%size [centroids by 1].
                        
                        %% Load the velocity data which was just saved in the ascii file:
                        
                        %Convert the ascii filename to a filepath:
                        SD_fitted_ascii_filename = [data_directory 'Velocities/' radar_identifiers{i_r,1} '/' SD_fitted_ascii_filename];%#ok<AGROW> %this is like 'gbr_20170928_1848_04_velocities_GEO_v6.dat', or similar, and here I have added on the full filepath.
                        
                        %For at least radar kod, I have noticed that the ascii file
                        % may exist, but that it may also be completely empty:
                        % check for this here:
                        ascii_file_details = dir(SD_fitted_ascii_filename);
                        if(ascii_file_details.bytes < 10)
                            %In this case, you should take the existence of an
                            % empty file to indicate some error in the processing.
                            continue%pertains to loop over i_subdaily_file.
                        end%Conditional: if the ascii file is empty, skip the processing of it.
                        
                        %Load velocity data:
                        velocity_data_bulk = dlmread(SD_fitted_ascii_filename);%size [(radar data readings in this day) by 9].
                        %Source: an above cell of this program, running fit_file_parse.c.
                        %Format: a number of rows of existing data, each with the
                        % column format:
                        % - beam number (0 to i)
                        % - range gate number (1 to j)
                        % - year of beam
                        % - month of beam
                        % - day of beam
                        % - hour of beam
                        % - minute of beam
                        % - second of beam
                        % - velocity at centroid (units metres per second).
                        %The order of the rows is as given by the C binary files
                        % read into by fit_file_parse.c. The C program reads each
                        % beam at a time, then outputs all range gates within that
                        % beam one by one.  Thus I would assume the ordering is
                        % temporal, then by range gate, but this is not tested.
                        
                        %% Convert the centroid positions to QD using the epochs of the beams:
                        %We can apply the rotation to QD for multiple times at
                        % once, but can only input a single location at once. So we
                        % loop over all the locations present in this day's
                        % velocity data:
                        
                        %Determine the number of unique locations within this day
                        % of data, using all velocity data rows for the beam and
                        % range gate number columns:
                        number_of_unique_locations = size(unique(velocity_data_bulk(:,1:2),'rows'),1);%size [1 by 1];
                        
                        %Specify the beam numbers and range gates of the unique
                        % locations, for reference in the loop over them:
                        unique_locations_beams_gates = unique(velocity_data_bulk(:,1:2),'rows');%size [unique locations by 2].  Column format: beam numbers, range gate numbers.
                        
                        %Data check: there should not be more unique locations than
                        % permitted by the centroid position parameters:
                        if(number_of_unique_locations > ((max(centroid_locations_beam_number_single_day)+1) * max(centroid_locations_range_gate_single_day)))
                            disp('Unexpected number of unique locations: possible error in code: terminating.')
                            return
                        end%Conditional: data check.
                        
                        %Loop over each of these unique locations and perform the
                        % rotation to QD coordinates for all entailed epochs:
                        %Preallocate centroid QD coordinate storage, and offset
                        % angles:
                        velocity_data_centroid_theta_QD = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        velocity_data_centroid_phi_QD = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        velocity_data_centroid_QD_quadrant_span = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        velocity_data_centroid_beam_QD_north_offset = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        %Order of columns is (i.e. MUST be) the same as velocity_data_bulk.
                        %Variables used to test that the beam and range gate order
                        % is the same in velocity_data_centroid_theta_QD (and phi
                        % equivalent) as in velocity_data_bulk:
                        temp_stored_beam_numbers = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        temp_stored_gate_numbers = NaN(size(velocity_data_bulk,1),1);%size [(radar data readings in this day) by 1];
                        %Make an index of velocity data rows to later remove --
                        % this is initially based on errors in the radar measurement
                        % temporal records:
                        index_temporal_errors = [];
                        for i_centroid = 1:number_of_unique_locations
                            %% Find the data which were measured at this single centroid location:
                            %Specify the temporally invariant geographic colatitude
                            % and longitude of this (unique) radar centroid
                            % location by linking up the beam and range gate
                            % numbers between those specified in the ascii (made
                            % from a C binary file) and those specified in the
                            % centroid metadata locations synthesised by the C
                            % program centroid_calculation.c.  This will return a
                            % single index value pertaining to the matrices of
                            % centroid locations:
                            index_locations_centroid = find((centroid_locations_beam_number_single_day == unique_locations_beams_gates(i_centroid,1)) & (centroid_locations_range_gate_single_day == unique_locations_beams_gates(i_centroid,2)));%size [1 by 1].  Pertains to all centroid_locations_*_single_day variables.
                            single_centroid_GEO_colatitude = centroid_locations_GEO_colatitude_single_day(index_locations_centroid,1);%size [1 by 1].
                            single_centroid_GEO_colongitude = centroid_locations_GEO_longitude_single_day(index_locations_centroid,1) + 360;%size [1 by 1]. Defined in -180 to 180 range, converted here to co-longitude.
                            
                            %Specify which rows of velocity_data_bulk this unique
                            % location pertains to.  These rows will be taken from
                            % the velocity_data_bulk matrix, and eventually stored
                            % in velocity_data_centroid_theta_SM (and phi
                            % equivalent) in the same order as in velocity_data_bulk:
                            index_locations_velocity = find((velocity_data_bulk(:,1) == unique_locations_beams_gates(i_centroid,1)) & (velocity_data_bulk(:,2) == unique_locations_beams_gates(i_centroid,2)));%size [(number of measurements taken at this location in one day) by 1]. Pertains to velocity_data_bulk.
                            
                            %% Compute QD basis vectors at the month-centre epoch, then the horizontal angle between them:
                            %We apply two separate conversions to QD: one in which
                            % the offset to local QD north is computed for the
                            % monthly mean epoch and used to assign an angular
                            % partition fiducial to the centroid for this month,
                            % and another in which the QD locations are computed.
                            
                            %Apply QD rotation to this centroid at the month-centre
                            % epoch:
                            [QD_lat, QD_lon, Apex_lat, MLT, F1_YX, F2_YX] = qdipole2(month_mean_epoch_MJD2000, 1, single_centroid_GEO_colatitude.*rad, single_centroid_GEO_colongitude.*rad, [code_directory 'subroutines/quasi_dipole/apexsh_1980-2020.txt']);%t, r, theta, phi, filename of coefficients.
                            %F1 is the mostly eastwards QD basis vector, F2 is the
                            % mostly northwards QD basis vector.  The columns of F#
                            % are each Y, X. These vectors define the QD basis
                            % vectors in terms of GEO component-vectors.
                            
                            %Change F1 and F2 to be expressed in the phi and theta
                            % directions, rather than the Y and X direction.  This
                            % involves reversing the sign of the second column.
                            % Hence, we will obtain an identical F# vector, in the
                            % [R]TP system rather than in the [Z]XY system.
                            F1_PT = [F1_YX(:,1)  (-F1_YX(:,2))];%PT = Phi, Theta columns.
                            F2_PT = [F2_YX(:,1)  (-F2_YX(:,2))];%PT = Phi, Theta columns.
                            
                            %Define vertical-direction unit vector:
                            k = [1;0;0];%R (positive upwards), Theta and Phi both zero.  Size [3 by 1], elements R;T;P.
                            
                            %Compute F according to Laundal and Gjerloev 2014:
                            % it's (F1xF2).k:
                            F = abs(dot(cross([0; F1_PT(:,2); F1_PT(:,1)], [0; F2_PT(:,2); F2_PT(:,1)]), k));%both inputs are of form R,T,P in GEO.
                            
                            %Compute norms of F1 and F2, and scale by F: these are the magnitudes of
                            % the QD basis vectors:
                            F2_norm = sqrt(sum(F2_PT .^ 2));
                            F1_norm = sqrt(sum(F1_PT .^ 2));
                            QDthetaF2BVmagnitude = F2_norm ./ F;%scalar: to be applied later.
                            QDphiF1BVmagnitude =   F1_norm ./ F;
                            
                            
                            %Define 'R,Theta,Phi' versions of the F2 (north) and F1
                            % (east) QD basis vectors, for rotation to 3D
                            % Cartesian so we can get the angle between them:
                            QDthetaF2BV = [0; F2_PT(1,2); F2_PT(1,1)];%size [3 by 1], rows R, Theta, Phi components of horizontal basis vector.
                            QDphiF1BV =   [0; F1_PT(1,2); F1_PT(1,1)];%size [3 by 1], rows R, Theta, Phi components of horizontal basis vector.
                            
                            %Make them unit vectors to aid in the angle
                            % computation:
                            QDthetaF2BV = QDthetaF2BV ./ norm(QDthetaF2BV);
                            QDphiF1BV = QDphiF1BV ./ norm(QDphiF1BV);
                            
                            %Define a spherical polar unit vector pointing
                            % towards the local north pole, from the QD F2 vector:
                            pointing_vector_QD_north_sph = QDthetaF2BV;
                            
                            %Convert to 3D Cartesian components:
                            %Make a rotation matrix: code from Overlap_SegmentsMakerHost_and_Solver_fn_DQI.m
                            % and gg2gm_tmodmulti.m, reference is Langel and Hinze,
                            % p.115. Coordinates used for rotation basis are those of
                            % the present centroid location, which is the origin of the
                            % QD basis vectors.
                            c_t = cos(single_centroid_GEO_colatitude .* rad);
                            s_t = sin(single_centroid_GEO_colatitude .* rad);
                            c_p = cos(single_centroid_GEO_colongitude .* rad);
                            s_p = sin(single_centroid_GEO_colongitude .* rad);
                            R_Sph2Cart = zeros(3,3);
                            R_Sph2Cart(1,:) = [(s_t.*c_p)  (c_t.*c_p)  (-s_p)];
                            R_Sph2Cart(2,:) = [(s_t.*s_p)  (c_t.*s_p)  (+c_p)];
                            R_Sph2Cart(3,:) = [ c_t        (-s_t)       0    ];
                            
                            %Convert the QD local north pointing vector to Cartesian:
                            pointing_vector_QD_north_Cart = R_Sph2Cart * pointing_vector_QD_north_sph;%size is [3 by 1] with columns of x,y,z.
                            
                            %Convert the spherical polar QD basis vectors to Cartesian:
                            QDthetaF2BV_Cart = R_Sph2Cart * QDthetaF2BV;%size is [3 by 1] with columns of x,y,z.
                            QDphiF1BV_Cart = R_Sph2Cart * QDphiF1BV;%size is [3 by 1] with columns of x,y,z.
                            %The above code applies:
                            %Br = QDtheta_F1_x_k(1,1);
                            %Btheta = QDtheta_F1_x_k(2,1);
                            %Bphi = QDtheta_F1_x_k(3,1);
                            %Bx = s_t.*c_p.*Br + c_t.*c_p.*Btheta - s_p.*Bphi
                            %By = s_t.*s_p.*Br + c_t.*s_p.*Btheta + c_p.*Bphi;
                            %Bz = c_t.*Br - s_t.*Btheta;
                            
                            %Normalise the vectors:
                            pointing_vector_north = pointing_vector_QD_north_Cart ./ norm(pointing_vector_QD_north_Cart);%size [3 by 1]. Format [x;y;z] Cartesian components.
                            QDthetaF2BV_Cart = QDthetaF2BV_Cart ./ norm(QDthetaF2BV_Cart);%size [3 by 1]. Format [x;y;z] Cartesian components.
                            QDphiF1BV_Cart = QDphiF1BV_Cart ./ norm(QDphiF1BV_Cart);%size [3 by 1]. Format [x;y;z] Cartesian components.
                            
                            %% Compute beam pointing direction in GEO, in Cartesian frame:
                            %Radar beam angle computed using GEO locations, since
                            % QD basis vectors are defined with respect to GEO
                            % coordinates.
                            
                            %Transformation from spherical to Cartesian coordinates...
                            % ...for centroids:
                            s_t = sin(single_centroid_GEO_colatitude.*rad);%theta coordinate. Size [1 by 1].
                            c_t = cos(single_centroid_GEO_colatitude.*rad);%theta coordinate. Size [1 by 1].
                            s_p = sin(single_centroid_GEO_colongitude.*rad);%phi coordinate. Size [1 by 1].
                            c_p = cos(single_centroid_GEO_colongitude.*rad);%phi coordinate. Size [1 by 1].
                            centroid_Cart_x = 6371.2 .* s_t .* c_p;%size [1 by 1].
                            centroid_Cart_y = 6371.2 .* s_t .* s_p;%size [1 by 1].
                            centroid_Cart_z = 6371.2 .* c_t;%size [1 by 1].
                            
                            % ...for radar location:
                            s_t = sin(radar_GEO_colatitude.*rad);%theta coordinate. Size [1 by 1].
                            c_t = cos(radar_GEO_colatitude.*rad);%theta coordinate. Size [1 by 1].
                            s_p = sin(radar_GEO_colongitude.*rad);%phi coordinate. Size [1 by 1].
                            c_p = cos(radar_GEO_colongitude.*rad);%phi coordinate. Size [1 by 1].
                            radar_Cart_x = 6371.2 .* s_t .* c_p;%size [1 by 1].
                            radar_Cart_y = 6371.2 .* s_t .* s_p;%size [1 by 1].
                            radar_Cart_z = 6371.2 .* c_t;%size [1 by 1].
                            
                            %For ease of reference define two vectors, each from the centre
                            % of the Earth to one of the points required to determine the angle
                            % of the radar beam:
                            %For the radar position in for the present epoch:
                            radial_vector_radar = [radar_Cart_x; radar_Cart_y; radar_Cart_z];%size [3 by 1]. Vector from origin (Earth's centre) to the radar position, in Cartesian coordinates.
                            %For the centroid in question at this point in the loop:
                            radial_vector_centroid = [centroid_Cart_x; centroid_Cart_y; centroid_Cart_z];%size [3 by 1]. Vector from origin (Earth's centre) to the beam/gate centroid position, in Cartesian coordinates.
                            
                            %Defining the pointing direction along the surface of the sphere:
                            %The cross product of the radial vector from Earth's
                            % centre to the centroid location, and the radial
                            % vector from Earth's centre to the radar position,
                            % will result in a vector with origin at the centroid,
                            % directed perpendicular to those two. The subsequent
                            % cross product of the output of the first cross
                            % product, with the radial vector from Earth's centre
                            % to the centroid position, will produce a vector
                            % located at the centroid, directed along the Earth's
                            % surface in the pointing direction of the radar beam
                            % (which is directed in a positive sense towards the
                            % radar, since SuperDARN measures Doppler line of sight
                            % velocities).
                            pointing_vector_beam = cross(cross(radial_vector_centroid,radial_vector_radar),radial_vector_centroid);%size [3 by 1]. Format [x;y;z] Cartesian components.
                            
                            %Normalise the vector:
                            pointing_vector_beam = pointing_vector_beam ./ norm(pointing_vector_beam);%size [3 by 1]. Format [x;y;z] Cartesian components.
                            
                            %% Rotate the 'north' and 'beam' pointing vectors (also the QD basis vectors) to the 'top' of the sphere:
                            
                            %Define rotation matrix about z axis, based on phi coordinate
                            % angle, to be applied in Cartesian coordinates.  This rotates the
                            % x unit vector to point through the Greenwich meridian:
                            c_a = cos(single_centroid_GEO_colongitude.*rad);%size [1 by 1].
                            s_a = sin(single_centroid_GEO_colongitude.*rad);%size [1 by 1].
                            R_z = [[+c_a s_a 0]
                                [-s_a c_a 0]
                                [  0   0  1]];%rotate about z axis.
                            
                            %Define rotation matrix about y axis, based on theta coordinate
                            % angle, to be applied in Cartesian coordinates.  This rotates the
                            % x' unit vector to point along the Cartesian horizontal plane:
                            c_a = cos(single_centroid_GEO_colatitude.*rad);%size [1 by 1].
                            s_a = sin(single_centroid_GEO_colatitude.*rad);%size [1 by 1].
                            R_y = [[c_a   0  -s_a]
                                [ 0    1    0 ]
                                [s_a   0   c_a]];%rotate about y' axis.
                            
                            %Combine the rotation matrices:
                            R = R_y * R_z;%z rotation, followed by y' rotation. Size [3 by 3].
                            
                            %Apply the rotation matrix to the two pointing vectors:
                            % Beam vector:
                            pv_beam_x = R(1,1).*pointing_vector_beam(1,1) + R(1,2).*pointing_vector_beam(2,1)  +R(1,3).*pointing_vector_beam(3,1);%size [1 by 1].
                            pv_beam_y = R(2,1).*pointing_vector_beam(1,1) + R(2,2).*pointing_vector_beam(2,1)  +R(2,3).*pointing_vector_beam(3,1);%size [1 by 1].
                            pv_beam_z = R(3,1).*pointing_vector_beam(1,1) + R(3,2).*pointing_vector_beam(2,1)  +R(3,3).*pointing_vector_beam(3,1);%size [1 by 1].
                            % North vector:
                            pv_north_x = R(1,1).*pointing_vector_north(1,1) + R(1,2).*pointing_vector_north(2,1)  +R(1,3).*pointing_vector_north(3,1);%size [1 by 1].
                            pv_north_y = R(2,1).*pointing_vector_north(1,1) + R(2,2).*pointing_vector_north(2,1)  +R(2,3).*pointing_vector_north(3,1);%size [1 by 1].
                            pv_north_z = R(3,1).*pointing_vector_north(1,1) + R(3,2).*pointing_vector_north(2,1)  +R(3,3).*pointing_vector_north(3,1);%size [1 by 1].
                            
                            % QD south vector:
                            QD_south_x = R(1,1).*QDthetaF2BV_Cart(1,1) + R(1,2).*QDthetaF2BV_Cart(2,1)  +R(1,3).*QDthetaF2BV_Cart(3,1);%size [1 by 1].
                            QD_south_y = R(2,1).*QDthetaF2BV_Cart(1,1) + R(2,2).*QDthetaF2BV_Cart(2,1)  +R(2,3).*QDthetaF2BV_Cart(3,1);%size [1 by 1].
                            QD_south_z = R(3,1).*QDthetaF2BV_Cart(1,1) + R(3,2).*QDthetaF2BV_Cart(2,1)  +R(3,3).*QDthetaF2BV_Cart(3,1);%size [1 by 1].
                            % QD east vector:
                            QD_east_x = R(1,1).*QDphiF1BV_Cart(1,1) + R(1,2).*QDphiF1BV_Cart(2,1)  +R(1,3).*QDphiF1BV_Cart(3,1);%size [1 by 1].
                            QD_east_y = R(2,1).*QDphiF1BV_Cart(1,1) + R(2,2).*QDphiF1BV_Cart(2,1)  +R(2,3).*QDphiF1BV_Cart(3,1);%size [1 by 1].
                            QD_east_z = R(3,1).*QDphiF1BV_Cart(1,1) + R(3,2).*QDphiF1BV_Cart(2,1)  +R(3,3).*QDphiF1BV_Cart(3,1);%size [1 by 1].
                            
                            %Check that the radial component of the rotated vectors is tiny,
                            % thus that they are in the Cartesian horizontal plane:
                            if(abs(pv_beam_z) > 0.0001 | abs(pv_north_z) > 0.0001) %#ok<OR2>
                                disp('Possible error in rotation of beam- or north-pointing vector.  Terminating.')
                                return
                            end%Conditional: error-check.
                            
                            %% Compute QD north offset at single centroid location from the Cartesian unit vectors:
                            
                            %Use the horizontal components of the rotated vectors to get the
                            % angle between them, retaining the quadrant information:
                            %This rotates clockwise from the horizontal components of the north
                            % pointing vector and stops (returning the angular value) when it
                            % first encounters the radar beam horizontal direction:
                            %Approach from 'https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151925', message 77:
                            beam_offset_angle_from_QD_north = 360 - mod(atan2(det([[pv_north_x pv_north_y];[pv_beam_x pv_beam_y]]),dot([pv_north_x pv_north_y],[pv_beam_x pv_beam_y])),2*pi) .* deg;%size [1 by 1].
                            
                            %% Compute angle between QD basis vectors in Cartesian frame:
                            
                            %Compute angle between the two QD vectors using the horizontal components
                            % of the component vectors, retaining the quadrant information:
                            %This rotates clockwise from the horizontal components
                            % of the QD east pointing vector and stops (returning the angular value) when it
                            % first encounters the QD south pointing vector:
                            %Approach from 'https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151925', message 77:
                            QD_offset_angle = 360 - mod(atan2(det([[QD_east_x QD_east_y];[QD_south_x QD_south_y]]),dot([QD_east_x QD_east_y],[QD_south_x QD_south_y])),2*pi) .* deg;%size [1 by 1].
                            
                            %% Angle normalisation and magnitude scaling:
                            
                            %Account for the non-orthogonality of the QD basis vectors:
                            %Here we normalise the QD north
                            % offset angles of the beam by the difference between
                            % the [QDnorth to QDeast angle] and [90 degrees]. Thus
                            % when we bin in an even set of angular partitions and
                            % fit a cosine to get two orthogonal components from
                            % the ~30 angular partitions, we will already have
                            % accounted for the QD non-orthogonality, at each
                            % centroid. Our QD model can then be validly
                            % re-projected back into the GEO frame, by use of the QD
                            % basis vectors. It will also be comparable to the EIMF
                            % model frame.
                            %Directly alter the QD north offset angle by the
                            % difference between the QD north and QD east basis
                            % vectors:
                            beam_offset_angle_from_QD_north = beam_offset_angle_from_QD_north ./ ((360 - QD_offset_angle) ./ 90);
                            
                            %Compute the projection of the line of sight vector
                            % onto the F2 and 'F2 + 90 degrees' vectors
                            beam_projection_onto_F2 = cos(beam_offset_angle_from_QD_north .* rad);
                            beam_projection_onto_F2orth = sin(beam_offset_angle_from_QD_north .* rad);
                            %Note: at this point, norm([beam_projection_onto_F2  beam_projection_onto_F2orth]) = 1.
                            
                            %Multiply the 'amounts of the line of sight vector'
                            % which are in the directions of the F2 and 'F2 + 90
                            % degrees' vectors in order to get a scalar magnitude
                            % with which to scale all the velocity values from this
                            % centroid position, for the entire month:
                            single_centroid_scaling_value = sqrt(((beam_projection_onto_F2 .* QDthetaF2BVmagnitude) .^ 2) + ((beam_projection_onto_F2orth .* QDphiF1BVmagnitude) .^ 2));
                            
                            %Directly alter the velocity values from this centroid
                            % by this amount:
                            velocity_data_bulk(index_locations_velocity,9) = velocity_data_bulk(index_locations_velocity,9) .* single_centroid_scaling_value;
                            
                            %% Assign QD quadrant span and beam offset angles to all instances of this beam/gate centroid for this day:
                            %Note: computed for the month-centre epoch only.
                            
                            velocity_data_centroid_QD_quadrant_span(index_locations_velocity,1) = QD_offset_angle;
                            velocity_data_centroid_beam_QD_north_offset(index_locations_velocity,1) = beam_offset_angle_from_QD_north;
                            
                            %% Create an index of temporally-incorrect records, to later remove from the bulk velocity values:
                            
                            %The requirement for this piece of code came because
                            % the QD rotation was failing for some centroids which
                            % had epoch records which were 20 years incorrect (for
                            % radar han, 2001-07).  Here, I add to an index of
                            % elements of the entire set of velocities to remove
                            % once processing is complete for this radar and for
                            % this, and before saving out the daily data file for
                            % the radar. Also, these values cannot be included in
                            % the QD coordinate tranformation, since it has a
                            % limited working time range. The +-0.5 parts are to
                            % prevent rounding errors removing otherwise-fine
                            % records.
                            %Index the incorrect elements of index_locations_velocity:
                            index_temporal_errors_relative_fiducials_part1 = find(velocity_data_bulk(index_locations_velocity,3) > (i_year + 0.5));%pertains to rows of index_locations_velocity, or matrix-subsets made with it.
                            index_temporal_errors_relative_fiducials_part2 = find(velocity_data_bulk(index_locations_velocity,3) < (i_year - 0.5));%pertains to rows of index_locations_velocity, or matrix-subsets made with it.
                            index_temporal_errors_relative_fiducials = [index_temporal_errors_relative_fiducials_part1;index_temporal_errors_relative_fiducials_part2];%pertains to rows of index_locations_velocity, or matrix-subsets made with it.
                            %Index the incorrect rows of velocity_data_bulk:
                            index_temporal_errors = [index_temporal_errors; index_locations_velocity(index_temporal_errors_relative_fiducials,1)]; %#ok<AGROW> %pertains to rows of velocity_data_bulk.
                            %Note: index_locations_velocity is of size [n by 1], so
                            % we are using its rows here. The values of
                            % index_locations_velocity pertain to the rows of
                            % velocity_data_bulk. The lines of code above find the
                            % elements of [the subset of velocity_data_bulk which is
                            % indexed by index_locations_velocity] which have
                            % incorrect year records. That will return a fiducial
                            % stating which elements of index_locations_velocity
                            % relate to incorrect records in velocity_data_bulk.
                            % We add those offending elements to a storage matrix
                            % called index_temporal_errors, which pertains to the
                            % rows of velocity_data_bulk (these are the elements to
                            % remove from velocity_data later).
                            
                            %% Compute QD positions of centroid data at all epochs in day:
                            
                            %Convert the range of epochs sampled at this single
                            % location to MJD2000 format:
                            single_location_epochs_MJD2000 = datenum(velocity_data_bulk(index_locations_velocity,3:8)) - 730486;%inputting year-to-second records. Output size is [(number of measurements taken at this location in one day) by 1].
                            
                            %If any records for this radar's daily file (and the
                            % measurements from it at this particular beam/gate
                            % centroid) have failed the temporal-error test in the
                            % above cell, alter the MDJ2000 epoch to that of the
                            % month's centroid. This way, the QD rotation will not
                            % break, and the values it returns (for the incorrect
                            % rows only) can be removed later:
                            if(~isempty(index_temporal_errors_relative_fiducials))
                                single_location_epochs_MJD2000(index_temporal_errors_relative_fiducials,1) = month_mean_epoch_MJD2000;
                            end%Conditional: if there are no offending data to alter, don't alter them.
                            
                            %Calculate QD positions for this centroid for all epochs:
                            [QD_lat, QD_lon, Apex_lat, MLT, F1_YX, F2_YX] = qdipole2(single_location_epochs_MJD2000, ones(size(single_location_epochs_MJD2000)), ones(size(single_location_epochs_MJD2000)).*single_centroid_GEO_colatitude.*rad, ones(size(single_location_epochs_MJD2000)).*single_centroid_GEO_colongitude.*rad, [code_directory 'subroutines/quasi_dipole/apexsh_1980-2020.txt']);%t, r, theta, phi, filename of coefficients.
                            
                            %Compute the theta and phi positions for this radar,
                            % for all epochs, in QD MLT and QD colatitude:
                            velocity_data_centroid_theta_QD(index_locations_velocity,1) = 90 - QD_lat;
                            velocity_data_centroid_phi_QD(index_locations_velocity,1) = MLT.*15;
                            
                            %Store the beam and gate numbers for data checks:
                            temp_stored_beam_numbers(index_locations_velocity,1) = velocity_data_bulk(index_locations_velocity,1);
                            temp_stored_gate_numbers(index_locations_velocity,1) = velocity_data_bulk(index_locations_velocity,2);
                            
                        end%Loop over each unique centroid location.
                        
                        %Data check:
                        if((sum(temp_stored_beam_numbers - velocity_data_bulk(:,1)) > 0) || (sum(temp_stored_gate_numbers - velocity_data_bulk(:,2)) > 0))
                            disp('Error in assignation of velocity locations: terminating.')
                            return
                        end%Conditional: data check for order of stored locaitons being the same as that already used for the velocities.
                        
                        %% Concatenate data to usable format:
                        
                        %Concatenate the SM coordinates with the velocities:
                        velocity_data_subdaily_part = [velocity_data_bulk(:,1:2)  velocity_data_centroid_theta_QD  velocity_data_centroid_phi_QD  velocity_data_bulk(:,3:end)  velocity_data_centroid_QD_quadrant_span  velocity_data_centroid_beam_QD_north_offset];
                        %Format of columns:
                        % - Beam number (0 to n).
                        % - Range gate (1 to n).
                        % - SM theta centroid coordinate.
                        % - SM phi centroid coordinate.
                        % - Year.
                        % - Month.
                        % - Day.
                        % - Hour.
                        % - Minute.
                        % - Second.
                        % - Velocity.
                        % - Angle between the month-centre-epoch QD basis vectors at each centroid.
                        % - Angle between beam direction and QD north.
                        
                        %Remove the rows of data which failed the temporal checks:
                        if(~isempty(index_temporal_errors))
                            velocity_data_subdaily_part(index_temporal_errors,:) = [];
                        end%Conditional: if there are no offending data to remove, don't remove them.
                        
                        %% Store this sub-daily part of the day's data:
                        
                        %Add the sub-daily file's contents to a larger matrix:
                        velocity_data = [velocity_data; velocity_data_subdaily_part]; %#ok<AGROW>
                        %Format of 13 columns:
                        % - Beam number (0 to i).
                        % - Range gate (1 to j).
                        % - QD theta centroid coordinate.
                        % - QD phi centroid coordinate.
                        % - Year.
                        % - Month.
                        % - Day.
                        % - Hour.
                        % - Minute.
                        % - Second.
                        % - Velocity.
                        % - Angle between the month-centre-epoch QD basis vectors at each centroid.
                        % - Angle between beam direction and QD north.
                        
                        %% %!!!! not sure if this is required any longer! Clean up the daily/sub-daily ascii file used to contribute to the MATLAB binary file:
                        
                        %Health warning: if you terminate the proecss partway
                        % through a given day, then it's probably best to
                        % delete that day's files and start oover.
                        delete(SD_fitted_ascii_filename);
                        
                    end%Loop over the daily/sub-daily fitted autocorrelation files (parsed to ascii) contained within this day.
                    
                    %% %!!!! not sure if this is required any longer! Save out the daily file of velocities:
                    
                    %Write-out the file:
                    save([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
                        radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '_velocities_QD_' rotation_version '.mat'],'velocity_data');
                    
                    %% Check for missing data:
                    
                    %If the velocity_data matrix is empty, skip this day:
                    if(isempty(velocity_data))
                        disp(['No data for day ' num2str(i_day) ', radar ' radar_identifiers{i_r,1}])
                        continue%pertains to loop over i_day.
                    end%Conditional: skip the processing of this daily file is it's empty.
                    
                    %% Data selection: remove selected range gates:
                    %Removing all records from range gates above or below
                    % the allowed maximum or minimum:
                    velocity_data(velocity_data(:,2) < range_gate_allowed_min,:) = []; %#ok<SAGROW>
                    velocity_data(velocity_data(:,2) > range_gate_allowed_max,:) = []; %#ok<SAGROW>
                    
                    %% Aggregate the daily data in a larger matrix for the whole month:
                    
                    single_radar_velocity_data_cell{i_day,1} = velocity_data;%size of storage: [beam/gate centroids by 13].
                    %Format of 13 columns:
                    % - Beam number (0 to i).
                    % - Range gate (1 to j).
                    % - QD theta centroid coordinate.
                    % - QD phi centroid coordinate.
                    % - Year.
                    % - Month.
                    % - Day.
                    % - Hour.
                    % - Minute.
                    % - Second.
                    % - Velocity.
                    % - Angle between the month-centre-epoch QD basis vectors at each centroid.
                    % - Angle between beam direction and QD north.
                    
                    %Store the number of data in this day for this radar:
                    single_radar_daily_data_count(i_day,1) = size(velocity_data,1);%size of stored value: [1 by 1].
                    
                end%Loop over each day in the present month for this present radar, read the ascii files.
                
                %% Restructure the cell of radar data into a matrix array, and assign centroid-specific angular partition information:
                
                %Make an index of empty days:
                index_missing_days = find(isnan(single_radar_daily_data_count));%size [missing days by 1].
                
                %Check if all days of data are missing, in which case, skip
                % this radar:
                if(size(index_missing_days,1) == month_end_day)
                    %State that there were no data in this month:
                    each_radar_monthly_data_count(i_r,1) = 0;
                    
                    %Skip the part of the loop that stores the data from this radar:
                    continue%pertains to loop over i_r.
                end%Conditional: data check for there being no data for this radar.
                
                %Remove missing days of data:
                single_radar_velocity_data_cell(index_missing_days,:) = [];
                single_radar_daily_data_count(index_missing_days,:) = [];
                
                %Compute and store the number of data in the month:
                each_radar_monthly_data_count(i_r,1) = sum(single_radar_daily_data_count);%size of stored value: [1 by 1].
                
                %Loop over the valid days and concatenate data into a single matrix:
                %Preallocate storage:
                single_radar_velocity_data_matrix = NaN(each_radar_monthly_data_count(i_r,1),10);%size [(all measurements for this radar in one month) by 14].
                %Columns:
                % - QD theta centroid coordinate.
                % - QD phi centroid coordinate.
                % - Year.
                % - Month.
                % - Day.
                % - Hour.
                % - Minute.
                % - Second.
                % - Velocity.
                % - Angle between beam direction and QD north (i.e. angular offset from QD north).
                %
                %Not stored here: numerical radar identifier (not required yet).
                counter = 1;
                for i_existing_day = 1:size(single_radar_velocity_data_cell,1)
                    %Extract and store:
                    single_radar_velocity_data_matrix(counter:((counter + single_radar_daily_data_count(i_existing_day,:)) - 1),:) = [single_radar_velocity_data_cell{i_existing_day,1}(:,3:11)  ...
                        single_radar_velocity_data_cell{i_existing_day,1}(:,13)];%size of extraction is [(all measurements for this radar in one day) by 13].
                    
                    %Advance counter increment:
                    counter = counter + single_radar_daily_data_count(i_existing_day,:);
                end%Loop over all populated days in the data from this radar.
                
                %Store the data from this radar in the set for all radars:
                each_radar_velocity_data_cell{i_r,1} = single_radar_velocity_data_matrix;
                %Format of 10 columns:
                % - QD theta centroid coordinate.
                % - QD phi centroid coordinate.
                % - Year.
                % - Month.
                % - Day.
                % - Hour.
                % - Minute.
                % - Second.
                % - Velocity.
                % - Angle between beam direction and QD north (i.e. angular offset from QD north).
                
            end%Loop over all radars, load data from each.
            
            %After looping over all radars, delete all the log files produced by the BASH scripts:
            delete([code_directory 'subroutines/BASH/cover_*.log'])
            
            %% %!!!! more data binning code to go here........
            
            %% %!!!! code below yet to be integrated: 
            
            %% Restructure the cell of data from all radars into a continuous matrix:
            
            %Make an index of empty radars:
            index_missing_radars = find(each_radar_monthly_data_count == 0);%size [missing radars by 1].
            
            %Check if all radars are missing, in which case, skip this month:
            if(size(index_missing_radars,1) == size(radar_identifiers,1))
                continue%pertains to loop over i_month.
            end%Conditional: data check for there being no data for this month.
            
            %Make a duplicate copy of the radar numerical identifiers,
            % since it's about to be modifed and modifying the original
            % will cause issues if multiple months are processed in the
            % same program run:
            radar_numerical_codes_temp = radar_numerical_codes;%Note that this occurs in the loop over individual months, and is placed after aggregating all rdara data in this month.
            
            %Remove missing radars:
            each_radar_velocity_data_cell(index_missing_radars,:) = [];
            each_radar_monthly_data_count(index_missing_radars,:) = [];
            radar_numerical_codes_temp(index_missing_radars,:) = [];%now pertains to same rows as the cell each_radar_velocity_data_cell.
            
            %Compute total number of data in the month from all radars:
            monthly_data_count_total = sum(each_radar_monthly_data_count);%size of stored value: [1 by 1].
            
            %Loop over the existing radars and concatenate data into a single matrix:
            %Preallocate storage:
            all_radars_velocity_data_matrix = NaN(monthly_data_count_total,11);%size [(all measurements for all radars in one month) by 11].
            %Columns:
            % - QD theta centroid coordinate.
            % - QD phi centroid coordinate.
            % - Year.
            % - Month.
            % - Day.
            % - Hour.
            % - Minute.
            % - Second.
            % - Velocity.
            % - Angular offset of beam from QD north.
            % - Numerical radar identifier.
            %Sorted by radar, then by time (presumably).
            counter = 1;
            for i_existing_radar = 1:size(each_radar_velocity_data_cell,1)
                %Extract and store:
                all_radars_velocity_data_matrix(counter:((counter + each_radar_monthly_data_count(i_existing_radar,:)) - 1),:) = ...
                    [each_radar_velocity_data_cell{i_existing_radar,1}  ...
                    repmat(radar_numerical_codes_temp(i_existing_radar,1),[size(each_radar_velocity_data_cell{i_existing_radar,1},1) 1])];
                
                %Advance counter increment:
                counter = counter + each_radar_monthly_data_count(i_existing_radar,:);
            end%Loop over all populated days in the data from this radar.
            
            %% Sense checks:
            
            if(sum(sum(isnan(all_radars_velocity_data_matrix))) > 0)
                disp('Warning: NaNs found in all_radars_velocity_data_matrix.')
            end%Conditional: data check.
            
            %% Partition the data from all radars into angular-offset bins:
            
            %Define the angular partition ranges:
            % Use the partition_angular_width value to define partitions of
            % similar pointing, wrt local north. We define the bins over two
            % quadrants (i.e. 0 to 180 degrees), and in each angular range we
            % also include the values in the antiparallel direction range:
            angle_partition_delimiting_edges = [((0 - (partition_angular_width / 2)):partition_angular_width:((180 - (partition_angular_width / 2)) - partition_angular_width))'  ...
                ((0 + (partition_angular_width / 2)):partition_angular_width:(180 - (partition_angular_width / 2)))'];%size [number of angular partitions by 2].
            
            %Define the angular partition centroids:
            angle_partition_centroids = mean(angle_partition_delimiting_edges,2);%size [number of angular partitions by 1].
            
            %Loop over each partition, find the data inside it:
            all_radars_velocity_partitioned = cell(size(angle_partition_delimiting_edges,1),1);%size [number of angular partitions by 1].
            for i_p = 1:size(angle_partition_delimiting_edges,1)
                %Find the radar centroids which have this angular separation from
                % north, and the antiparallel direction too:
                if(i_p == 1)
                    index_single_partition_parallel_part1 = find(all_radars_velocity_data_matrix(:,10) >= (360 + angle_partition_delimiting_edges(i_p,1)));
                    index_single_partition_parallel_part2 = find(all_radars_velocity_data_matrix(:,10) < angle_partition_delimiting_edges(i_p,2));%pertains to rows of all_radars_velocity_data_matrix.
                    index_single_partition_parallel = [index_single_partition_parallel_part1;index_single_partition_parallel_part2];%pertains to rows of all_radars_velocity_data_matrix.  Unordered!
                    %Sort the index:
                    index_single_partition_parallel = sort(index_single_partition_parallel);
                else
                    index_single_partition_parallel = find(all_radars_velocity_data_matrix(:,10) >= angle_partition_delimiting_edges(i_p,1) & all_radars_velocity_data_matrix(:,10) < angle_partition_delimiting_edges(i_p,2));%pertains to rows of all_radars_velocity_data_matrix.
                end%Conditional: first partition has a negative (first) edge value, thus needs treating differently.
                index_single_partition_antiparallel = find(all_radars_velocity_data_matrix(:,10) >= (angle_partition_delimiting_edges(i_p,1) + 180) & all_radars_velocity_data_matrix(:,10) < (angle_partition_delimiting_edges(i_p,2) + 180));%pertains to rows of all_radars_velocity_data_matrix.
                
                %Extract the relevant data, combining the parallel and antiparallel
                % direction bins and reversing the signs of the velocities as
                % appropriate:
                %Some logical switches are required here since the output of datenum
                % for a zero-row input is 6 columns (by 0 rows), rather than 1 column,
                % and this may make the concatenation unstable in future releases of
                % Matlab:
                if(isempty(index_single_partition_parallel) & ~isempty(index_single_partition_antiparallel))
                    all_radars_velocity_partitioned{i_p,1} = [...
                        [(datenum(all_radars_velocity_data_matrix(index_single_partition_antiparallel,3:8)) - 730486) ...%antiparallel times (MJD2000).
                        all_radars_velocity_data_matrix(index_single_partition_antiparallel,1:2) ...%antiparallel centroid locations (theta, phi).
                        (-all_radars_velocity_data_matrix(index_single_partition_antiparallel,9)) ...%antiparallel velocity at centroid.
                        all_radars_velocity_data_matrix(index_single_partition_antiparallel,11)]...%antiparallel radar identifier.
                        ];%Size of all_radars_velocity_partitioned is [measurements in partition by 5].
                elseif(isempty(index_single_partition_antiparallel) & ~isempty(index_single_partition_parallel))
                    all_radars_velocity_partitioned{i_p,1} = [...
                        [(datenum(all_radars_velocity_data_matrix(index_single_partition_parallel,3:8)) - 730486) ...%parallel times (MJD2000).
                        all_radars_velocity_data_matrix(index_single_partition_parallel,1:2) ...%parallel centroid locations (theta, phi).
                        all_radars_velocity_data_matrix(index_single_partition_parallel,9) ...%parallel velocity at centroid.
                        all_radars_velocity_data_matrix(index_single_partition_parallel,11)] ...%parallel radar identifier.
                        ];%Size of all_radars_velocity_partitioned is [measurements in partition by 5].
                else
                    %Presumably if both directions are filled or if both are unfilled,
                    % the column counts will match.
                    all_radars_velocity_partitioned{i_p,1} = [...
                        [(datenum(all_radars_velocity_data_matrix(index_single_partition_parallel,3:8)) - 730486) ...%parallel times (MJD2000).
                        all_radars_velocity_data_matrix(index_single_partition_parallel,1:2) ...%parallel centroid locations (theta, phi).
                        all_radars_velocity_data_matrix(index_single_partition_parallel,9) ...%parallel velocity at centroid.
                        all_radars_velocity_data_matrix(index_single_partition_parallel,11)]; ...%parallel radar identifier.
                        [(datenum(all_radars_velocity_data_matrix(index_single_partition_antiparallel,3:8)) - 730486) ...%antiparallel times (MJD2000).
                        all_radars_velocity_data_matrix(index_single_partition_antiparallel,1:2) ...%antiparallel centroid locations (theta, phi).
                        (-all_radars_velocity_data_matrix(index_single_partition_antiparallel,9)) ...%antiparallel velocity at centroid.
                        all_radars_velocity_data_matrix(index_single_partition_antiparallel,11)]...%antiparallel radar identifier.
                        ];%Size of all_radars_velocity_partitioned is [measurements in partition by 5].
                end%Conditional: helping the concatenation not have to combine unmatching column counts.
                %Format of columns of all_radars_velocity_partitioned:
                % - measurement time (MJD2000).
                % - centroid location QD theta.
                % - centroid location QD phi.
                % - velocity at centroid.
                % - numerical radar identifier.
                %Not stored here: direction from north at centroid, since
                % we can get that information (coarsely) from the angular
                % partition.
            end%Loop over the number of angular partitions.
            
            %% Bin the data from each angular partition by temporally discretising measurements from all radars, prior to binning spatially:
            
            %Note: SuperMAG is defined as centered boxcar 1-minute means,
            % of which the EIMF model (which we want this to be compearable
            % to) takes 5-min means.  Hence, the EIMF model's 5-min epochs
            % span -0.5 to 4.5 minutes each, and they have a centroid epoch
            % at 2 mins.  We do the same here.  Note that the way the code
            % is set up, we lose 30 seconds of data at thehstart and end of
            % the month -- this is considered negligible and the mean is
            % taken anyway.
            
            %Bin temporal width parameters:
            %Convert bin UT width to MJD2000:
            bin_UT_width_MJD = bin_UT_width / 60 / 24;%minutes to hours, hours to decimal day fraction.
            %Define monthly start and end times in MJD2000:
            this_month_start_MJD2000 = datenum([i_year  i_month       01 0 0 0]) - 730486;
            next_month_start_MJD2000 = datenum([i_year  (i_month + 1) 01 0 0 0]) - 730486;%'month+1' known to work even for month 12.
            %Subtract 30 seconds from the monthly start and end times:
            this_month_start_MJD2000 = this_month_start_MJD2000 - (30/60/60/24);%converting 30 seconds from minutes to decimal day format.
            next_month_start_MJD2000 = next_month_start_MJD2000 - (30/60/60/24);%converting 30 seconds from minutes to decimal day format.
            
            %Define bracketing epochs of each bin UT width:
            bin_time_bracketing_epochs = [(this_month_start_MJD2000:bin_UT_width_MJD:(next_month_start_MJD2000 - bin_UT_width_MJD))'  ...
                ((this_month_start_MJD2000 + bin_UT_width_MJD):bin_UT_width_MJD:next_month_start_MJD2000)'];%size [bin UT widths by 2]. Columns: start-time of bin UT width, end-time.
            %Note the above code line may be susceptible to rounding errors.
            %Define bin-time centroids:
            bin_time_centroids = mean(bin_time_bracketing_epochs,2);%size [bin UT widths by 1].
            
            %Check for numerical errors:
            if((month_end_day / bin_UT_width_MJD) ~= round(month_end_day / bin_UT_width_MJD))
                disp('Numerical error in bin temporal definition.  Terminating.')
                return
            end%Conditional: data check.
            
            %Preallocate storage for each partition direction:
            bin_data_vel_all_dirs = NaN(size(bin_coords_colat,1),size(bin_time_bracketing_epochs,1),size(angle_partition_centroids,1));%size [bins by times by angular partitions].
            bin_contrib_radars_all_dirs = cell(size(bin_coords_colat,1),size(bin_time_bracketing_epochs,1),size(angle_partition_centroids,1));%size [bins by times by angular partitions].
            bin_contrib_radars_all_dirs(:) = {NaN};%preallocation: the matrix version of this used in v08 (and previous) had NaNs as default.
            filled_bins_count = NaN(size(angle_partition_centroids,1),1);%size [angular partitions by 1].
            %Loop over each angular partition centroid direction:
            for i_angular_partition = 1:size(angle_partition_centroids,1)
                disp(['Processing angular partition ' num2str(i_angular_partition) ' of ' num2str(size(angle_partition_centroids,1))])
                
                %Extract the values for this angular partition direction,
                % from all radars:
                all_radars_velocity_SingleDir = all_radars_velocity_partitioned{i_angular_partition,1};%size [measurements in this angular partition by 5].
                %Format of columns:
                % - measurement time (MJD2000).
                % - centroid location QD theta.
                % - centroid location QD phi.
                % - velocity at centroid.
                % - numerical radar identifier.
                
                %Loop over each bin UT width:
                for i_temporal_bin = 1:size(bin_time_bracketing_epochs,1)
                    if((i_temporal_bin - mod(i_temporal_bin,1000)) == i_temporal_bin)
                        disp(['...processing UT width ' num2str(i_temporal_bin) ' of ' num2str(size(bin_time_bracketing_epochs,1))]);
                    end%Conditional: progress-indicator, every 1000 iterations.
                    
                    %Find all radar records for this direction within this
                    % time span, where 'sUTw' means single UT width: 
                    index_SingleDir_sUTw = find((all_radars_velocity_SingleDir(:,1) >= bin_time_bracketing_epochs(i_temporal_bin,1)) & (all_radars_velocity_SingleDir(:,1) < bin_time_bracketing_epochs(i_temporal_bin,2)));%size [single-epoch measurements by 1]. Pertains to rows of all_radars_velocity_SingleDir.
                    
                    %If there's no data in this epoch, skip it:
                    if(isempty(index_SingleDir_sUTw))
                        disp(['No data in bin width ' num2str(i_temporal_bin)])
                        continue%pertains to loop over i_temporal_bin.
                    end%Conditional: data check.
                    
                    %Extract single bin UT width of radar records for later reference:
                    all_radars_velocity_SingleDir_sUTw = all_radars_velocity_SingleDir(index_SingleDir_sUTw,:);%size [single-epoch measurements by 5]. Columns: time, theta coordinate, phi coordinate, velocity, radar identifier.
                    
                    %Loop over spatial bins:
                    for i_spatial_bin = 1:size(bin_coords_colat,1)
                        
                        %Extract index of all 'direction 1' locations in latitude band of present bin:
                        index_colat_SingleDir = find(all_radars_velocity_SingleDir_sUTw(:,2) >= bin_coords_colat(i_spatial_bin,1) & all_radars_velocity_SingleDir_sUTw(:,2) <= bin_coords_colat(i_spatial_bin,2));%size [n by 1]. Pertains to rows of all_radars_velocity_dir1_sUTw.
                        
                        %Extract index of all 'direction 1' locations in longitude band of present bin:
                        if(bin_coords_colong(i_spatial_bin,2) > 360)
                            index_long_SingleDir_part1 = find(all_radars_velocity_SingleDir_sUTw(:,3) >= bin_coords_colong(i_spatial_bin,1));
                            index_long_SingleDir_part2 = find(all_radars_velocity_SingleDir_sUTw(:,3) <= mod(bin_coords_colong(i_spatial_bin,2),360));
                            %Inputs to ismembc must be sorted, so apply a sort:
                            index_long_SingleDir = sort([index_long_SingleDir_part1;index_long_SingleDir_part2]);
                        else
                            index_long_SingleDir = find(all_radars_velocity_SingleDir_sUTw(:,3) >= bin_coords_colong(i_spatial_bin,1) & all_radars_velocity_SingleDir_sUTw(:,3) <= bin_coords_colong(i_spatial_bin,2));
                        end%Conditional: the bin longitudes occasionally extend over 360 degrees, and this must be accounted for.
                        
                        %Define index of values within this bin:
                        index_SingleUTwidth_SingleBin_SingleDir = index_colat_SingleDir(ismembc(index_colat_SingleDir,index_long_SingleDir));
                        
                        if(isempty(index_SingleUTwidth_SingleBin_SingleDir))%nothing in bin.
                            %Assign dummy value to empty bin:
                            %Code line removed since we start with NaNs.
                        elseif(size(index_SingleUTwidth_SingleBin_SingleDir,1) == 1)%one record in bin.
                            %Assign the single measurement to the bin:
                            bin_data_vel_all_dirs(i_spatial_bin,i_temporal_bin,i_angular_partition) = all_radars_velocity_SingleDir_sUTw(index_SingleUTwidth_SingleBin_SingleDir,4);
                            
                            %Assign the number of the radar used to make the measurement:
                            bin_contrib_radars_all_dirs{i_spatial_bin,i_temporal_bin,i_angular_partition} = all_radars_velocity_SingleDir_sUTw(index_SingleUTwidth_SingleBin_SingleDir,5);
                        elseif(size(index_SingleUTwidth_SingleBin_SingleDir,1) > 1)%more than one record in bin.
                            %Extract this direction's measurements within this bin,
                            % where 'sB' means single bin:
                            all_radars_velocity_SingleDir_sUTw_sB = all_radars_velocity_SingleDir_sUTw(index_SingleUTwidth_SingleBin_SingleDir,:);%size [n by 5]. Columns: time, theta coordinate, phi coordinate, velocity, radar identifier.
                            
                            %Take the median of all measurements in the cell:
                            bin_data_vel_all_dirs(i_spatial_bin,i_temporal_bin,i_angular_partition) = nanmedian(all_radars_velocity_SingleDir_sUTw_sB(:,4));
                            
                            %Store the identifiers of all radars which
                            % contributed to this bin, of which there may
                            % be more than one:
                            bin_contrib_radars_all_dirs{i_spatial_bin,i_temporal_bin,i_angular_partition} = all_radars_velocity_SingleDir_sUTw_sB(:,5);%size [n by 1].
                            
                        end%Conditional: computes and stores bin value from measurements in colat, long, time span of bin.
                    end%Loop over each bin, store the data for this direction.
                end%Loop over each bin UT width.
                
                %Count the number of filled bins in the binned set from
                % this angular partition direction. The 'sum' tells you the
                % count of missing temporal records in each spatial bin.
                % The 'find' function tells you which of these spatial bins
                % has a missing-records count equal to the number of
                % temporal slots (i.e. is completely empty). The 'size'
                % function thus tells you how many empty bins there are for
                % this angular partition, and we take the difference
                % between that and the number of NPC bins to get the number
                % of filled bins:
                filled_bins_count(i_angular_partition,1) = size(bin_data_vel_all_dirs,1) - size(find(sum(isnan(bin_data_vel_all_dirs(:,:,i_angular_partition)),2) == size(bin_time_bracketing_epochs,1)),1);%size of stored value is [1 by 1].
                
            end%Loop over each angular partition centroid direction, assign the entailed data to a full set of [space by time] bins.
            
            %Check that there are no NaNs in the filled bin count, or the
            % following bits of code might not work:
            if(any(isnan(filled_bins_count)))
                disp('Warning: NaNs encountered in filled bin counter: terminating.')
                return
            end%Conditional: data check.
            
            %% Remove bins with no contributions, index remaining bins per partition, concatenate for storage:
            
            %Loop over each angular partition:
            %Preallocate storage:
            bin_data_vel = NaN(size(bin_time_bracketing_epochs,1),sum(filled_bins_count));%size [(bin UT width times) by (all bins which had a contribution for all angular partitions)].
            bin_contrib_radars = cell(size(bin_time_bracketing_epochs,1),sum(filled_bins_count));%size [(bin UT width times) by (all bins which had a contribution for all angular partitions)].
            bin_contrib_radars(:) = {NaN};%preallocation: NaNs as default.
            bin_fiducial_indices = NaN(1,sum(filled_bins_count,1));%size [1 by (all bins which had a contribution for all angular partitions)].
            bin_angular_partitions = NaN(1,sum(filled_bins_count,1));%size [1 by (all bins which had a contribution for all angular partitions)].
            storage_counter = 1;
            for i_angular_partition = 1:size(angle_partition_centroids,1)
                %Extract temporary versions of the data of interest:
                bin_data_vel_all_dirs_SinglePartition = bin_data_vel_all_dirs(:,:,i_angular_partition)';%size of extraction after transpose is [UT widths by NPC bins].
                bin_contrib_radars_all_dirs_SinglePartition = bin_contrib_radars_all_dirs(:,:,i_angular_partition)';%size of extraction after transpose is [UT widths by NPC bins]. Cell format!
                index_original_bin_fiducials_SinglePartition = index_original_bin_fiducials';%size [1 by NPC bins] after transpose.
                
                %Crop the data to just the bins that were contributed to:
                index_missing_bins = find(sum(isnan(bin_data_vel_all_dirs_SinglePartition),1) == size(bin_time_bracketing_epochs,1));%pertains to columns of bin_data_vel_all_dirs_SinglePartition.
                bin_data_vel_all_dirs_SinglePartition(:,index_missing_bins) = [];%now of size [UT widths by bins with a contributing value].
                bin_contrib_radars_all_dirs_SinglePartition(:,index_missing_bins) = [];%now of size [UT widths by bins with a contributing value].
                index_original_bin_fiducials_SinglePartition(:,index_missing_bins) = [];%now of size [1 by bins with a contributing value]. Tells you which bins the cropped values pertain to, so it's very important.
                
                %Store the data:
                bin_data_vel(:,storage_counter:(storage_counter + filled_bins_count(i_angular_partition,1) - 1)) = bin_data_vel_all_dirs_SinglePartition;
                bin_contrib_radars(:,storage_counter:(storage_counter + filled_bins_count(i_angular_partition,1) - 1)) = bin_contrib_radars_all_dirs_SinglePartition;%cell format!
                bin_fiducial_indices(:,storage_counter:(storage_counter + filled_bins_count(i_angular_partition,1) - 1)) = index_original_bin_fiducials_SinglePartition;
                
                %Also store the angular partition that these data were
                % binned from, since there is no other record of this
                % information:
                bin_angular_partitions(:,storage_counter:(storage_counter + filled_bins_count(i_angular_partition,1) - 1)) = angle_partition_centroids(i_angular_partition,1);%size of stored extraction is [1 by 1].
                
                %Advance storage counter:
                storage_counter = storage_counter + filled_bins_count(i_angular_partition,1);
            end%Loop over each angular partition centroid direction, concatenate data files intoa  binned data matrix.
            
            %% Save binned data:
            save(binned_data_filename,'bin_data_vel', 'bin_time_centroids', 'bin_fiducial_indices', 'bin_angular_partitions');
            %save(bin_contrib_radars_filename, 'bin_contrib_radars', version,'-v7.3')%Removes compression, allows saving out of the cell variable.
            
        end%Conditional: if binned data file exists, read it in.  Otherwise, bin the data.
        
        %% Data selection, post-binning:
        %Note: I am aware that this is slightly 'hacky' but this is the
        % only way I could think of keeping the ri_string small whilst also
        % retaining some record of what selection has been applied to the
        % data.
        
        if(strcmp(ds,'08'))
            disp(['Data selection version ' ds ' applied.'])
            
            %Bins with centroid latitude of 59.4 (rounded) are known to
            % provide poor post-EOF infill performance due to their
            % typically inadequate coverage in years prior to the
            % introduction of midlatitude radars. Hence, any values
            % equatorward of that latitude are removed by removing bins
            % with centroid colatitudes greater than 30:
            index_removal_bins = find(bin_centroids_colat > 30);
            fiducials_to_remove = index_original_bin_fiducials(index_removal_bins,1);
            
            for i_removal = 1:size(fiducials_to_remove,1)
                bins_to_NaN = find(bin_fiducial_indices == fiducials_to_remove(i_removal,1));
                bin_data_vel(:,bins_to_NaN) = NaN;
            end%Loop over removal indices.
            
            
            %Remove binned values below 50ms-1, to reduce the impact of
            % ground scatter on the EOF solution:
            index_data_removal = find((bin_data_vel > -50) & (bin_data_vel < 50));
            bin_data_vel(index_data_removal) = NaN;
            
            
            %Remove all contributions in polar bin:
            index_removal_bins = find(bin_centroids_colat == 0);
            fiducials_to_remove = index_original_bin_fiducials(index_removal_bins,1);
            
            for i_removal = 1:size(fiducials_to_remove,1)
                bins_to_NaN = find(bin_fiducial_indices == fiducials_to_remove(i_removal,1));
                bin_data_vel(:,bins_to_NaN) = NaN;
            end%Loop over removal indices.
            
        else
            disp('Incorrect data selection applied: check settings.')
            return
        end%Conditional: what data selection to apply?
        
        %% Check data matrix dimensions, reset them to S-mode format if required:
        
        %Data check: this program will only work if bin_data_vel is of size
        % [times by bins].
        if(size(bin_data_vel,1) ~= size(bin_time_centroids,1))
            disp('The bin data was size [bins by times], but it has been transposed, and now it is [times by bins].')
            bin_data_vel = bin_data_vel';%Now of size [times by bins].
        elseif(size(bin_data_vel,1) == size(bin_time_centroids,1))
            disp('The bin data is already size [times by bins].')
        else
            disp('Warning: bin data not of expected dimensions or size!  Terminating.')
            return
        end%Conditional: is the length of the first dimension of bin_data_B equal to the count of bins?
        
        %!!!! remove this once the check has confirmed the data matrix size!
        
        %% Calculate weighting matrices:
        
        %We compute two weighting matrices. Each has only diagonal values,
        % with the off-diagonal elements being zero.
        % - The 'temporal' weights have a value of 1 where a given row of
        %   bin_data_vel is fully-filled, and 0 where that row is empty.
        %   This will act to decrease the relative contribution of
        %   poorly-sampled epochs to the EOF solution for the variance of
        %   the data.
        % - The 'inverse spatial' weights have a value of 1 where a given
        %   column of bin_data_vel is empty, and 0 where that column is
        %   fully filled. This weights the data to increase the relative
        %   contribution from poorly sampled locations.
        
        %Compute temporal weights:
        %Calculate the missing values in all 'spatial' dimension elements, for
        % each 'epoch':
        W_temporal = eye(size(bin_data_vel,1));%size [times by times].
        
        %Loop over each diagonal element of the weighting matrix:
        for i_t = 1:size(bin_data_vel,1)
            %Calculate weight value for this epoch from the number of filled
            % spatial data points divided by the number of possible spatial
            % data points, thus making a weight of 1 for fully-filled, and zero
            % for empty:
            W_temporal(i_t,i_t) = (numel(bin_data_vel(i_t,:)) - sum(isnan(bin_data_vel(i_t,:)))) ./ size(bin_data_vel,2);
        end%Loop over each diagonal element of the weighting matrix.
        
        %Compute inverse spatial weights:
        %Calculate the missing values in all temporal dimension elements, for
        % each 'spatial' bin:
        W_inverse_spatial = eye(size(bin_data_vel,2));%size [spatial elements by spatial elements].
        
        %Loop over each diagonal element of the weighting matrix:
        for i_s = 1:size(bin_data_vel,2)
            %Calculate weight value for this epoch from the number of filled
            % spatial data points divided by the number of possible spatial
            % data points, thus making a weight of 1 for fully-filled, and zero
            % for empty:
            W_inverse_spatial(i_s,i_s) = (numel(bin_data_vel(:,i_s)) - sum(isnan(bin_data_vel(:,i_s)))) ./ size(bin_data_vel,1);
        end%Loop over each diagonal element of the weighting matrix.
        
        %Extract the computed weights:
        W_inverse_spatial_diag = diag(W_inverse_spatial);
        
        %Loop over each diagonal element of the weighting matrix:
        for i_s = 1:size(bin_data_vel,2)
            W_inverse_spatial(i_s,i_s) = 1 ./ W_inverse_spatial_diag(i_s);
        end%Loop over each diagonal element of the weighting matrix.
        
        %% Calculate mean of binned data along specified dimension of data matrix, for removal prior to EOF analysis:
        
        %Calculate the bin data means.  Note that the data matrix is
        % S-mode-format at this stage, so it is of size [times by 
        % bins], so the first dimension is times.  
        if(strcmp(EOF_options.postbin_mean_removal_dimension, 's'))
            removed_bin_data_mean = nanmean(bin_data_vel,2);%Mean calculated along spatial dimension. Size [times by 1].
        elseif(strcmp(EOF_options.postbin_mean_removal_dimension, 't'))
            removed_bin_data_mean = nanmean(bin_data_vel,1);%Mean calculated along temporal dimension. Size [1 by bins].
        else
            disp('Post-binning mean removal option incorrectly specified.  Terminating.')
            return
        end%Conditional: which dimension to calculate the means across, and wehich to replicate them across?
        
        %% Perform standard EOF analysis of binned data matrix with temporal weights:
        %We perform two EOF analyses: one with amplitude boosting applied,
        % the other without amplitude boosting. The meaning of 'amplitude
        % boosting' is explained later.
        
        %Input format notes:
        %
        %The input data matrix must be in the format [times by bins] (which I refer
        % to as 'S-mode'). This is because the covariance matrix equations take the
        % [times by bins] layout as an assumed starting point, and it is the
        % covariance matrix which controls whether the analysis is S-mode or
        % T-mode.
        %
        %The input data matrix must also contain NaNs for any un-sampled bin
        % elements.
        
        %Display EOF analysis parameters:
        disp(['Starting ' EOF_options.analysis_mode ' EOF analysis: '])
        disp(['... data matrix size [t by s] = ' num2str(size(bin_data_vel,1)) ' by ' num2str(size(bin_data_vel,2))])
        disp(['... mean removal dimension: ' EOF_options.postbin_mean_removal_dimension])
        disp(['... weights of size ' num2str(size(W_temporal,1)) ' with temporal basis (decreasing the relative contribution of poorly-sampled epochs)'])
        disp(['... iterating for ' num2str(EOF_options.number_of_iterated_modes) ' modes ' num2str(EOF_options.number_of_iterations_per_mode) ' times each'])
        disp(['... using Lanczos solver, solving ' num2str(EOF_options.number_of_lanczos_basis_functions_solved) ' eigenvectors per iteration and returning ' num2str(EOF_options.number_of_lanczos_eigenvectors_returned) '.'])
        
        %Define run-identifier string, listing the EOF analysis parameters
        % for a later filename:
        EOF_WtABn_ri_string = ['SD_' binning_and_EOF_version '_NPC_AP_' EOF_options.analysis_mode '_Lanczos' num2str(EOF_options.number_of_lanczos_basis_functions_solved) 's' ...
            num2str(EOF_options.number_of_lanczos_eigenvectors_returned) 'r_' num2str(EOF_options.number_of_iterated_modes) 'm' num2str(EOF_options.number_of_iterations_per_mode) 'x' ...
            '_MR' EOF_options.postbin_mean_removal_dimension '_UTw' num2str(bin_UT_width)...
            '_APw' num2str(partition_angular_width) '_rQD' rotation_version '_ds' ds '_WtABn'];
        %i.e 'SD_v10_NPC_AP_Tmode_Lanczos10s1r_10m35x_MRt_UTw5_APw6_rQDv7_ds08_WtABn'.
        %Meaning:
        % 'SD': SuperDARN.
        % 'v10': a version identifier for the outputs of this program.
        % 'NPC': North Polar Cap region selected for analysis.
        % 'AP': 'all partitions'. This signifies that all angular partitions
        %       of the radar's look directions have been used when forming
        %       the EOF analysis basis in time, space, and look-direction.
        % 'Tmode': at this [point in the code, the data matrix
        %          'bin_data_vel' (call it X) has dimensions of [times by
        %          bins]. The S-mode/T-mode terminology is described in
        %          Bjornsson and Venegas (1997) and Richman (1986). In an
        %          S-mode EOF analysis, we form the covariance matrix with
        %          X^T*X, hence the covariance matrix has dimensions (space
        %          by space). In a T-mode EOF analysis, we form the
        %          covariance matrix with X*X^T, hence covariance matrix
        %          has dimensions (time by time). Here, T-mode is selected
        %          since there are typically more SuperDARN spatial bins
        %          than temporal elements for a given monthly analysis.
        % 'Lanczos10s1r': a Lanczos solver is used to significantly speed
        %                 up the EOF analysis time: rather than solving for
        %                 several thousand eigenvectors at each iteration,
        %                 we can solve for 10, and return just 1 
        %                 eigenvector (with the highest eigenvalue), which
        %                 is sufficient for the infill to converge in
        %                 amplitude (within acceptable loss bounds). Note
        %                 that the default MATLAB specification for a
        %                 Lanczos solution is to solve for fewer than 10
        %                 modes, which does not provide an acceptable
        %                 infill solution.
        % '10m35x': we solve for 10 moeds overall, and we iterate each one
        %           35 times.
        % 'MRt': Mean removed along temporal direction.
        % 'UTw5': signifies that the temporal bins used to form the
        %         irregularly-sampled SuperDARN data into a regular temporal
        %         basis for the EOF analysis are each 5 mins in length.
        % 'APw6': signifies that the width of the angular partitions (used to
        %         convert the radar look-directions into a discrete basis
        %         prior to the EOF analysis) is 6 degrees per partition.
        % 'rQDv7': a version identifier for the coordinate system rotation
        %          approach used in this program.
        % 'ds08': shorthand for the post-binning data seelction applied to
        %         the data prior to EOF analysis.
        % 'Wt': using 'temporal' weights, defined above.
        % 'ABn': amplitude boosting option set to 'no'. Amplitude boosting
        %        is defined later on.
        
        %Define index of gaps in the data: these will be infilled in the iteration
        % process.
        NaNs_index = find(isnan(bin_data_vel));%Important: returns a 1D vector: defined for S-mode dimensions (i.e. [times by space]), MUST be used ONLY for matrices with S-mode layouts!
        
        %Calculate the bin data means and centre the data matrix.  Note that the
        % data matrix is S-mode-format at this stage, so it is of size [times by
        % bins], so the first dimension is times.  
        if(strcmp(EOF_options.postbin_mean_removal_dimension, 's'))
            %Display warning if the combination of mean removal and analysis mode
            % will give strange results:
            if(strcmp(EOF_options.analysis_mode,'Tmode'))
                disp('Warning: covariance matrix will be T-mode, and spatial-dimension mean removal is not recommended.')
            end%Conditional: sense-check.
            
            %Remove previously calculated mean:
            X_c = bin_data_vel - repmat(removed_bin_data_mean, [1 size(bin_data_vel,2)]);%The EOFs should now be zero-mean along the dimension the means were removed.
        elseif(strcmp(EOF_options.postbin_mean_removal_dimension, 't'))
            %Display warning if the combination of mean removal and analysis mode
            % will give strange results:
            if(strcmp(EOF_options.analysis_mode,'Smode'))
                disp('Warning: covariance matrix will be S-mode, and temporal-dimension mean removal is not recommended.')
            end%Conditional: sense-check.
            
            %Remove previously calculated mean:
            X_c = bin_data_vel - repmat(removed_bin_data_mean, [size(bin_data_vel,1) 1]);%The EOFs should now be zero-mean along the dimension the means were removed.
        else
            disp('Post-binning mean removal option incorrectly specified.  Terminating.')
            return
        end%Conditional: which dimension to calculate the means across, and wehich to replicate them across?
        disp('Centred data matrix prior to iteration.')
        %The EOF analysis will now pertain to this datum, and there is no
        % need to remove any further mean in the subsequent iterative EOF
        % analyses.
        
        %The weights can only be applied to non-NaN values of the centered
        % data matrix, so we replace them with zeros, weight, then
        % re-insert the NaNs in case the zeros were altered in the
        % weighting process. Here we index the NaNs.
        X_c(NaNs_index) = 0;
        
        %Apply the weights to the data matrix, assuming that the data matrix is
        % dimensioned [times by space]. For why we apply the extra exponent of 0.5, refer to Jolliffe, 2002. Another reason is that it makes the reconstruction algebra easier.
        X_c_w = (W_temporal .^ 0.5) * X_c;
        
        %Replace the NaN elements:
        X_c_w(NaNs_index) = NaN;
        
        %Compute de-weighting matrix, from the inverse of the weighting matrix,
        % accounting for any exponents applied to the weights before it was applied
        % to the data matrix. The inverse weighting matrix is designed for
        % application to the reconstructed data matrix, or (if there are no
        % off-diagonal elements) to the eigenvector with the same dimension as the
        % square weights matrix:
        W_removal = W .^ (-1 .* 0.5);%the factor of a half is because we apply it to the data matrix, not the covariance matrix.
        W_removal(W_removal == Inf) = 0;
        
        %Set eigs options structure to the value of the number of desired
        % Lanczos eigenvectors for which to solve:
        opts.p = EOF_options.number_of_lanczos_basis_functions_solved;
        
        %Preallocate storage for the 'correction' ensemble matrix, which starts off
        % as a matrix of zeros. At each correction stage, we reconstruct the data
        % matrix based on the fully-iterated mode 1, and add this to the ensemble
        % correction matrix. This matrix is then used to remove (from the original
        % data matrix) the combined signal of all fully-iterated modes which have
        % larger eigenvalues than the one we're presently iterating for.
        X_c_w_correction_ensemble = zeros(size(X_c_w));%initial value of zero. This is the ensemble recon of modes of more importance than the one we're iterating for.
        %Storage preallocation for each fully-iterated mode, in weighted format:
        WtABn_T_w_iterN_mAll = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes]. Stores PC of iteration-N mode 1 only.
        WtABn_V_w_iterN_mAll = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes]. Stores eigenvector of iteration-N mode 1 only.
        WtABn_L_w_iterN_mAll = NaN(EOF_options.number_of_lanczos_basis_functions_solved,EOF_options.number_of_iterated_modes);%size [number of returned eigenvectors per iteration by number of iterated modes].
        %Storage preallocation for each fully-iterated mode, in de-weighted format:
        WtABn_T_dw_iterN_mAll = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes].
        WtABn_V_dw_iterN_mAll = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes].
        %Storage preallocation for each fully-iterated mode, in de-weighted format
        % with the temporal eigenvector in the units of the input data, and the
        % spatial eigenvector normalised (sign-preserved, so these can still be
        % used for reconstruction).
        WtABn_T_dw_iterN_mAll_Xunits = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes].
        WtABn_V_dw_iterN_mAll_UnitNorm = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes].
        %Storage preallocation for RMS of residual of [mode 1 prediction (at each
        % iteration stage)] minus [level-corrected data matrix], ignoring elements
        % that were originally NaN. In short, an estimate of how well that mode
        % represents the data for its correction level, and how this varies with
        % iteration:
        WtABn_mode_1_prediction_residual_rms_mAll = NaN(EOF_options.number_of_iterations_per_mode,EOF_options.number_of_iterated_modes);%size [number of iterations by number of modes].
        for i_mode = 1:(EOF_options.number_of_iterated_modes)
            %Iterative EOF analysis:
            disp(['WtABn: Iterating mode 1 for correction level ' num2str(i_mode)])
            
            %Correct the input data for the signal of the ensemble of previous
            % iterated modes, and then iterate for the present mode:
            X_c_w_SingleCorrectionLevel = X_c_w - X_c_w_correction_ensemble;%has no infill at present, but it will do shortly.
            
            %Iterating and infilling:
            X_c_w_infill = zeros(size(X_c_w_SingleCorrectionLevel));%Initial allocation only. The first EOF analysis must infill the gaps with zeros.
            for i_iter = 1:(EOF_options.number_of_iterations_per_mode)
                disp(['... iteration ' num2str(i_iter)]);
                
                %Infill the missing values in this stage of the iteration of the
                % data:
                X_c_w_SingleCorrectionLevel(NaNs_index) = X_c_w_infill(NaNs_index);
                
                %Solve for eigenvectors, calculate projection of these onto the data matrix:
                %Note that for both S- and T-mode, T are time-supported, and V are
                % space-supported. Which of these is the projection-derived
                % eigenvector set and which was solved for directly will vary
                % between the analysis 'modes'. The solved-for eigenvector (which
                % in T-mode is a temporal pattern and in S-mode is a spatial
                % pattern) is unitless, whilst the projected eigenvector has a set
                % of amplitudes which scale according to eigenvalue, the units of
                % which are the units of the input data scaled by the amplitude of
                % the associated solved-for eigenvector). In each S-/T-mode case,
                % the same reconstruction equations (X = T*V') still work fine.
                if(strcmp(EOF_options.analysis_mode, 'Smode'))
                    %Create covariance matrix:
                    tic
                    cov_w = (1/(size(X_c_w_SingleCorrectionLevel,1)-1)) .* (X_c_w_SingleCorrectionLevel' * X_c_w_SingleCorrectionLevel);%X^T*X, or 'S-mode' (Richman, 1986).
                    run_time = toc;
                    disp(['Construction of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Solve for n temporal eigenvectors with largest magnitude
                    % eigenvalues using Lanczos solver:
                    tic
                    [V_w_iterN_mi, L_w_iterN_mi] = eigs(cov_w,EOF_options.number_of_lanczos_eigenvectors_returned,'lm',opts);%outputs are weighted-version.
                    run_time = toc;
                    disp(['Lanczos solution of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Create (weighted) temporal expansion coefficient(s) by projecting
                    % the weighted spatial eigenvector(s) through the infilled data
                    % matrix:
                    T_w_iterN_mi = X_c_w_SingleCorrectionLevel * V_w_iterN_mi;%T is the matrix of temporal coefficients relating to each EOF: AKA the Principal Components or Expansion Coefficients.
                    
                elseif(strcmp(EOF_options.analysis_mode, 'Tmode'))
                    %Create covariance matrix:
                    tic
                    cov_w = (1/(size(X_c_w_SingleCorrectionLevel,1)-1)) .* (X_c_w_SingleCorrectionLevel * X_c_w_SingleCorrectionLevel');%X*X^T, or 'T-mode' (Richman, 1986).
                    run_time = toc;
                    disp(['Construction of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Solve for n temporal eigenvectors with largest magnitude
                    % eigenvalues using Lanczos solver:
                    tic
                    [T_w_iterN_mi, L_w_iterN_mi] = eigs(cov_w,EOF_options.number_of_lanczos_eigenvectors_returned,'lm',opts);%outputs are weighted-version.
                    run_time = toc;
                    disp(['Lanczos solution of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Create (weighted) spatial expansion coefficient(s) by projecting
                    % the weighted temporal eigenvector(s) through the infilled data
                    % matrix:
                    V_w_iterN_mi = X_c_w_SingleCorrectionLevel' * T_w_iterN_mi;%V is the matrix of spatial coefficients relating to each temporal T-mode EOF eigenvector. They are Expansion Coefficients, but they're also eigenvectors of their own.
                    
                end%Conditional: are the solved-for eigenvectors spatially supported (S-mode) or temporally supported (T-mode)
                
                %Create (weighted) reconstruction based on the first mode only,
                % replace the infill matrix values with it: assumes that
                % the standard MATLAB reverse-column-order of eigenvalue
                % size applies, such that the eigenvector with the largest
                % eigenvalue is in the last column:
                X_c_w_infill = T_w_iterN_mi(:,end) * V_w_iterN_mi(:,end)';
                
                %Compute RMS of residual of de-weighted prediction of mode 1 at
                % de-weighted original data, based on infilled data matrix, but
                % removing infilled elements prior to calculation of RMS:
                mode_1_prediction_residual = (W_removal * X_c_w_SingleCorrectionLevel) - (W_removal * X_c_w_infill);%De-weighting formula from Joliffe 2002.
                
                %Re-form residual into NaN-free vector:
                mode_1_prediction_residual(NaNs_index) = [];
                
                %Compute and store RMS:
                WtABn_mode_1_prediction_residual_rms_mAll(i_iter,i_mode) = rms(mode_1_prediction_residual);
                
            end%Loop over several EOF analyses, use first mode to infill data in each case.
            
            %Now update the ensemble of the mode-based reconstructed data matrix with
            % the fully-iterated reconstruction from the latest mode:
            X_c_w_correction_ensemble = X_c_w_correction_ensemble + X_c_w_infill;%now ready to remove from the next iteration of data.
            
            %After iterating, produce de-weighted and re-scaled versions of
            % each eigenvector for the current correction-level.
            %If the weighting was applied to the spatial dimension, we
            % de-weight the spatial eigenvector, and if it was applied to
            % the temporal dimension, we de-weight the temporal
            % eigenvector. here, the weights were applied to the temporal
            % dimension.
            %Make a de-weighted version of the temporal eigenvector. Note that its
            % amplitude is not scaled to the units of the data matrix:
            T_dw_iterN_mi = W_removal * T_w_iterN_mi;
            
            %Now that we have de-weighted the temporal eigenvector, all the
            % re-scaling required for a reconstruction has already been applied and
            % the de-weighted spatial eigenvector is the same as the one that came
            % out of the weighted analysis.  Note that this is only the case for
            % weights applied to the temporal dimension, with no off-diagonal values:
            V_dw_iterN_mi = V_w_iterN_mi;
            
            %To obtain a version of the de-weighted temporal eigenvector in the
            % units of the data matrix, we apply:
            T_dw_iterN_mi_Xunits = T_dw_iterN_mi .* abs(V_dw_iterN_mi(find(abs(V_dw_iterN_mi) == max(abs(V_dw_iterN_mi)),1,'first'),:));
            
            %The corrolary to the X-unit temporal eigenvector is a unit-norm
            % spatial eigenvector, de-weighted:
            V_dw_iterN_mi_UnitNorm = V_dw_iterN_mi ./ max(abs(V_dw_iterN_mi));
            
            %Store fully iterated mode 1 for this correction level of the data
            % matrix:
            WtABn_T_w_iterN_mAll(:,i_mode) = T_w_iterN_mi(:,end);%weighted-version temporal eigenvector
            WtABn_V_w_iterN_mAll(:,i_mode) = V_w_iterN_mi(:,end);%weighted-version spatial eigenvector
            WtABn_L_w_iterN_mAll(:,i_mode) = diag(L_w_iterN_mi);%size of stored part is [number of returned modes in the EOF solution by 1].
            WtABn_T_dw_iterN_mAll(:,i_mode) = T_dw_iterN_mi;
            WtABn_V_dw_iterN_mAll(:,i_mode) = V_dw_iterN_mi;
            WtABn_T_dw_iterN_mAll_Xunits(:,i_mode) = T_dw_iterN_mi_Xunits;%de-weighted X-unit temporal eigenvector
            WtABn_V_dw_iterN_mAll_UnitNorm(:,i_mode) = V_dw_iterN_mi_UnitNorm;%de-weighted normalised spatial eigenvector (due to normalisation, de-weighting only has an effect on this if the weights were applied to the spatial dimension).
            
        end%Loop over each iteration process for the next mode in the line.
        
        %% Perform amplitude-boosted EOF analysis of binned data matrix with inverse spatial weights:
        
        %Input format notes:
        %
        %The input data matrix must be in the format [times by bins] (which I refer
        % to as 'S-mode'). This is because the covariance matrix equations take the
        % [times by bins] layout as an assumed starting point, and it is the
        % covariance matrix which controls whether the analysis is S-mode or
        % T-mode.
        %
        %The input data matrix must also contain NaNs for any un-sampled bin
        % elements.
        
        %Display EOF analysis parameters:
        disp(['Starting ' EOF_options.analysis_mode ' EOF analysis: '])
        disp(['... data matrix size [t by s] = ' num2str(size(bin_data_vel,1)) ' by ' num2str(size(bin_data_vel,2))])
        disp(['... mean removal dimension: ' EOF_options.postbin_mean_removal_dimension])
        disp(['... weights of size ' num2str(size(W_temporal,1)) ' with temporal basis (decreasing the relative contribution of poorly-sampled epochs)'])
        disp(['... iterating for ' num2str(EOF_options.number_of_iterated_modes) ' modes ' num2str(EOF_options.number_of_iterations_per_mode) ' times each'])
        disp('... with aplitude boosting applied')
        disp(['... using Lanczos solver, solving ' num2str(EOF_options.number_of_lanczos_basis_functions_solved) ' eigenvectors per iteration and returning ' num2str(EOF_options.number_of_lanczos_eigenvectors_returned) '.'])
        
        %Define run-identifier string, listing the EOF analysis parameters
        % for a later filename:
        EOF_WsiABy_ri_string = ['SD_' binning_and_EOF_version '_NPC_AP_' EOF_options.analysis_mode '_Lanczos' num2str(EOF_options.number_of_lanczos_basis_functions_solved) 's' ...
            num2str(EOF_options.number_of_lanczos_eigenvectors_returned) 'r_' num2str(EOF_options.number_of_iterated_modes) 'm' num2str(EOF_options.number_of_iterations_per_mode) 'x' ...
            '_MR' EOF_options.postbin_mean_removal_dimension '_UTw' num2str(bin_UT_width)...
            '_APw' num2str(partition_angular_width) '_rQD' rotation_version '_ds' ds '_WsiABy'];
        %i.e 'SD_v10_NPC_AP_Tmode_Lanczos10s1r_10m35x_MRt_UTw5_APw6_rQDv7_ds08_WsiABy'.
        %Meaning:
        % 'SD': SuperDARN.
        % 'v10': a version identifier for the outputs of this program.
        % 'NPC': North Polar Cap region selected for analysis.
        % 'AP': 'all partitions'. This signifies that all angular partitions
        %       of the radar's look directions have been used when forming
        %       the EOF analysis basis in time, space, and look-direction.
        % 'Tmode': at this [point in the code, the data matrix
        %          'bin_data_vel' (call it X) has dimensions of [times by
        %          bins]. The S-mode/T-mode terminology is described in
        %          Bjornsson and Venegas (1997) and Richman (1986). In an
        %          S-mode EOF analysis, we form the covariance matrix with
        %          X^T*X, hence the covariance matrix has dimensions (space
        %          by space). In a T-mode EOF analysis, we form the
        %          covariance matrix with X*X^T, hence covariance matrix
        %          has dimensions (time by time). Here, T-mode is selected
        %          since there are typically more SuperDARN spatial bins
        %          than temporal elements for a given monthly analysis.
        % 'Lanczos10s1r': a Lanczos solver is used to significantly speed
        %                 up the EOF analysis time: rather than solving for
        %                 several thousand eigenvectors at each iteration,
        %                 we can solve for 10, and return just 1 
        %                 eigenvector (with the highest eigenvalue), which
        %                 is sufficient for the infill to converge in
        %                 amplitude (within acceptable loss bounds). Note
        %                 that the default MATLAB specification for a
        %                 Lanczos solution is to solve for fewer than 10
        %                 modes, which does not provide an acceptable
        %                 infill solution.
        % '10m35x': we solve for 10 moeds overall, and we iterate each one
        %           35 times.
        % 'MRt': Mean removed along temporal direction.
        % 'UTw5': signifies that the temporal bins used to form the
        %         irregularly-sampled SuperDARN data into a regular temporal
        %         basis for the EOF analysis are each 5 mins in length.
        % 'APw6': signifies that the width of the angular partitions (used to
        %         convert the radar look-directions into a discrete basis
        %         prior to the EOF analysis) is 6 degrees per partition.
        % 'rQDv7': a version identifier for the coordinate system rotation
        %          approach used in this program.
        % 'ds08': shorthand for the post-binning data seelction applied to
        %         the data prior to EOF analysis.
        % 'Wsi': using 'inverse spatial' weights, defined above.
        % 'ABy': amplitude boosting option set to 'yes'. Amplitude boosting
        %        is defined later on.
        
        %Define index of gaps in the data: these will be infilled in the iteration
        % process.
        NaNs_index = find(isnan(bin_data_vel));%Important: returns a 1D vector: defined for S-mode dinemsions (i.e. [times by space]), MUST be used ONLY for matrices with S-mode layouts!
        
        %Define index of present values in the data:
        not_NaNs_index = find(~isnan(bin_data_vel));%Important: returns a 1D vector: defined for S-mode dinemsions (i.e. [times by space]), MUST be used ONLY for matrices with S-mode layouts!
        
        %Calculate the bin data means and centre the data matrix.  Note that the
        % data matrix is S-mode-format at this stage, so it is of size [times by
        % bins*2], so the first dimension is times, etc.  The EOF analysis
        % will now pertain to this datum, and there is no need to remove
        % any further mean.
        if(strcmp(EOF_options.postbin_mean_removal_dimension, 's'))
            %Display warning if the combination of mean removal and analysis mode
            % will give strange results:
            if(strcmp(EOF_options.analysis_mode,'Tmode'))
                disp('Warning: covariance matrix will be T-mode, and spatial-dimension mean removal is not recommended.')
            end%Conditional: sense-check.
            
            %Remove previously calculated mean:
            X_c = bin_data_vel - repmat(removed_bin_data_mean, [1 size(bin_data_vel,2)]);%The EOFs should now be zero-mean along the dimension the means were removed.
        elseif(strcmp(EOF_options.postbin_mean_removal_dimension, 't'))
            %Display warning if the combination of mean removal and analysis mode
            % will give strange results:
            if(strcmp(EOF_options.analysis_mode,'Smode'))
                disp('Warning: covariance matrix will be S-mode, and temporal-dimension mean removal is not recommended.')
            end%Conditional: sense-check.
            
            %Calculate and remove mean:
            X_c = bin_data_vel - repmat(removed_bin_data_mean, [size(bin_data_vel,1) 1]);%The EOFs should now be zero-mean along the dimension the means were removed.
        else
            disp('Post-binning mean removal option incorrectly specified.  Terminating.')
            return
        end%Conditional: which dimension to calculate the means across, and wehich to replicate them across?
        disp('Centred data matrix prior to iteration.')
        
        %The weights can only be applied to non-NaN values of the centered
        % data matrix, so we replace them with zeros, weight, then
        % re-insert the NaNs in case the zeros were altered in the
        % weighting process. Here we index the NaNs.
        X_c(NaNs_index) = 0;
        
        %Apply the weights to the data matrix, assuming that the data matrix is
        % dimensioned [times by space]. For why we apply the extra exponent
        % of 0.5, refer to Jolliffe, 2002. Another reason is that it makes
        % the reconstruction algebra easier.
        X_c_w = X_c * (W_inverse_spatial .^ 0.5);
        
        %Replace the NaN elements:
        X_c_w(NaNs_index) = NaN;
        
        %Compute de-weighting matrix, from the inverse of the weighting matrix,
        % accounting for any exponents applied to the weights before it was applied
        % to the data matrix. The inverse weighting matrix is designed for
        % application to the reconstructed data matrix, or (if there are no
        % off-diagonal elements) to the eigenvector with the same dimension as the
        % square weights matrix:
        W_removal = W_inverse_spatial .^ (-1 .* 0.5);%the factor of a half is because we apply it to the data matrix, not the covariance matrix.
        W_removal(W_removal == Inf) = 0;
        
        %Set eigs options structure to the value of the number of desired
        % Lanczos eigenvectors for which to solve:
        opts.p = EOF_options.number_of_lanczos_basis_functions_solved;
        
        %Preallocate storage for the 'correction' ensemble matrix, which starts off
        % as a matrix of zeros. At each correction stage, we reconstruct the data
        % matrix based on the fully-iterated mode 1, and add this to the ensemble
        % correction matrix. This matrix is then used to remove (from the original
        % data matrix) the combined signal of all fully-iterated modes which have
        % larger eigenvalues than the one we're presently iterating for.
        X_c_w_correction_ensemble = zeros(size(X_c_w));%initial value of zero. This is the ensemble recon of modes of more importance than the one we're iterating for.
        %Storage preallocation for each fully-iterated mode, in weighted format:
        WsiABy_T_w_iterN_mAll = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes]. Stores PC of iteration-N mode 1 only.
        WsiABy_V_w_iterN_mAll = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes]. Stores eigenvector of iteration-N mode 1 only.
        WsiABy_L_w_iterN_mAll = NaN(EOF_options.number_of_lanczos_basis_functions_solved,EOF_options.number_of_iterated_modes);%size [number of returned eigenvectors per iteration by number of iterated modes].
        %Storage preallocation for each fully-iterated mode, in de-weighted format:
        WsiABy_T_dw_iterN_mAll = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes].
        WsiABy_V_dw_iterN_mAll = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes].
        %Storage preallocation for each fully-iterated mode, in de-weighted format
        % with the temporal eigenvector in the units of the input data, and the
        % spatial eigenvector normalised (sign-preserved, so tehse can still be
        % used for reconstruction).
        WsiABy_T_dw_iterN_mAll_Xunits = NaN(size(X_c_w,1),EOF_options.number_of_iterated_modes);%Size [times by number of iterated modes].
        WsiABy_V_dw_iterN_mAll_UnitNorm = NaN(size(X_c_w,2),EOF_options.number_of_iterated_modes);%Size [bins by number of iterated modes].
        %Storage preallocation for RMS of residual of [mode 1 prediction (at each
        % iteration stage)] minus [level-corrected data matrix], ignoring elements
        % that were originally NaN. In short, an estimate of how well that mode
        % represents the data for its correction level, and how this varies with
        % iteration:
        WsiABy_mode_1_prediction_residual_rms_mAll = NaN(EOF_options.number_of_iterations_per_mode,EOF_options.number_of_iterated_modes);%size [number of iterations by number of modes].
        for i_mode = 1:(EOF_options.number_of_iterated_modes)
            disp(['WsiABy: Iterating mode 1 for correction level ' num2str(i_mode)])
            
            %Correct the input data for the signal of the ensemble of previous
            % iterated modes, and then iterate for the present mode:
            X_c_w_SingleCorrectionLevel = X_c_w - X_c_w_correction_ensemble;%has no infill at present, but it will do shortly.
            
            %Iterating and infilling:
            X_c_w_infill = zeros(size(X_c_w_SingleCorrectionLevel));%Initial allocation only. The first EOF analysis must infill the gaps with zeros.
            for i_iter = 1:(EOF_options.number_of_iterations_per_mode)
                disp(['... iteration ' num2str(i_iter)]);
                
                %If it is the last iteration stage, boost the amplitude of the
                % infill before it is used for one last EOF analysis. After
                % this procedure, a final EOF analysis will be performed for this
                % corrective stage, which should rearrange the subsequent mode
                % structure to better match the real data.
                if(i_iter == (EOF_options.number_of_iterations_per_mode))
                    %For clarity, make variables for the infill, and the real data,
                    % where they co-exist:
                    infill_temp = X_c_w_infill(not_NaNs_index);%size [(number of existent data points) by 1], because the index is a 1D vector.
                    real_data_temp = X_c_w_SingleCorrectionLevel(not_NaNs_index);%size [(number of existent data points) by 1], because the index is a 1D vector.
                    
                    %Regress the (presumed weaker) infill onto the (presumed
                    % stronger) original data to determine a multiplicative factor
                    % for the infill:
                    multiplicative_factor = ((infill_temp' * infill_temp) \ (infill_temp' * real_data_temp));%size [1 by 1].
                    
                    %Directly increase the amplitude of the (almost fully iterated)
                    % original data using the multiplicative factor you just
                    % determined:
                    X_c_w_infill = X_c_w_infill .* multiplicative_factor;
                    
                end%Conditional: if it is the last iteration, boost the infill amplitude.
                
                %Infill the missing values in this stage of the iteration of the
                % data:
                X_c_w_SingleCorrectionLevel(NaNs_index) = X_c_w_infill(NaNs_index);
                
                %Solve for eigenvectors, calculate projection of these onto the data matrix:
                %Note that for both S- and T-mode, T are time-supported, and V are
                % space-supported. Which of these is the projection-derived
                % eigenvector set and which was solved for directly will vary
                % between the analysis 'modes'. The solved-for eigenvector (which
                % in T-mode is a temporal pattern and in S-mode is a spatial
                % pattern) is unitless, whilst the projected eigenvector has a set
                % of amplitudes which scale according to eigenvalue, the units of
                % which are the units of the input data scaled by the amplitude of
                % the associated solved-for eigenvector). In each S-/T-mode case,
                % the same recon equations (X = T*V') still work fine.
                if(strcmp(EOF_options.analysis_mode, 'Smode'))
                    %Create covariance matrix:
                    tic
                    cov_w = (1/(size(X_c_w_SingleCorrectionLevel,1)-1)) .* (X_c_w_SingleCorrectionLevel' * X_c_w_SingleCorrectionLevel);%X^T*X, or 'S-mode' (Richman, 1986).
                    run_time = toc;
                    disp(['Construction of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Solve for n temporal eigenvectors with largest magnitude
                    % eigenvalues using Lanczos solver:
                    tic
                    [V_w_iterN_mi, L_w_iterN_mi] = eigs(cov_w,EOF_options.number_of_lanczos_eigenvectors_returned,'lm',opts);%outputs are weighted-version.
                    run_time = toc;
                    disp(['Lanczos solution of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Create (weighted) temporal expansion coefficient(s) by projecting
                    % the weighted spatial eigenvector(s) through the infilled data
                    % matrix:
                    T_w_iterN_mi = X_c_w_SingleCorrectionLevel * V_w_iterN_mi;%T is the matrix of temporal coefficients relating to each EOF: AKA the Principal Components or Expansion Coefficients.
                    
                elseif(strcmp(EOF_options.analysis_mode, 'Tmode'))
                    %Create covariance matrix:
                    tic
                    cov_w = (1/(size(X_c_w_SingleCorrectionLevel,1)-1)) .* (X_c_w_SingleCorrectionLevel * X_c_w_SingleCorrectionLevel');%X*X^T, or 'T-mode' (Richman, 1986).
                    run_time = toc;
                    disp(['Construction of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Solve for n temporal eigenvectors with largest magnitude
                    % eigenvalues using Lanczos solver:
                    tic
                    [T_w_iterN_mi, L_w_iterN_mi] = eigs(cov_w,EOF_options.number_of_lanczos_eigenvectors_returned,'lm',opts);%outputs are weighted-version.
                    run_time = toc;
                    disp(['Lanczos solution of [' num2str(size(cov_w,1)) ' by ' num2str(size(cov_w,2)) '] covariance matrix took ' num2str(run_time) ' seconds to run.'])
                    
                    %Create (weighted) spatial expansion coefficient(s) by projecting
                    % the weighted temporal eigenvector(s) through the infilled data
                    % matrix:
                    V_w_iterN_mi = X_c_w_SingleCorrectionLevel' * T_w_iterN_mi;%V is the matrix of spatial coefficients relating to each temporal T-mode EOF eigenvector. They are Expansion Coefficients, but they're also eigenvectors of their own.
                    
                end%Conditional: are the solved-for eigenvectors spatially supported (S-mode) or temporally supported (T-mode)
                
                %Create (weighted) reconstruction based on the first mode only,
                % replace the infill matrix values with it:
                X_c_w_infill = T_w_iterN_mi(:,end) * V_w_iterN_mi(:,end)';
                
                %Compute RMS of residual of de-weighted prediction of mode 1 at
                % de-weighted original data, based on infilled data matrix, but
                % removing infilled elements prior to calculation of RMS:
                mode_1_prediction_residual = (X_c_w_SingleCorrectionLevel * W_removal) - (X_c_w_infill * W_removal);%De-weighting formula from Joliffe 2002.
                
                %Re-form residual into NaN-free vector:
                mode_1_prediction_residual(NaNs_index) = [];
                
                %Compute and store RMS:
                WsiABy_mode_1_prediction_residual_rms_mAll(i_iter,i_mode) = rms(mode_1_prediction_residual);
                
            end%Loop over several EOF analyses, use first mode to infill data in each case.
            
            %Now update the ensemble of the mode-based reconstructed data matrix with
            % the fully-iterated reconstruction from the latest mode:
            X_c_w_correction_ensemble = X_c_w_correction_ensemble + X_c_w_infill;%now ready to remove from the next iteration of data.
            
            %After iterating, produce de-weighted and re-scaled versions of
            % each eigenvector for the current correction-level:
            
            %If the weighting was applied to the spatial dimension, we de-weight the
            % spatial eigenvector, and if it was applied to the temporal dimension,
            % we de-weight the temporal eigenvector:
            %Make a de-weighted version of the spatial eigenvector:
            V_dw_iterN_mi = W_removal * V_w_iterN_mi;
            
            %Now that we have de-weighted the spatial eigenvector, all the
            % re-scaling required for a reconstruction has already been applied and
            % the de-weighted temporal eigenvector is the same as the one that came
            % out of the weighted analysis.  Note that this is only the case for
            % weights applied to the spatial dimension, with no off-diagonal values:
            T_dw_iterN_mi = T_w_iterN_mi;
            
            %To obtain a version of the de-weighted temporal eigenvector in the
            % units of the data matrix, we apply:
            T_dw_iterN_mi_Xunits = T_dw_iterN_mi .* abs(V_dw_iterN_mi(find(abs(V_dw_iterN_mi) == max(abs(V_dw_iterN_mi)),1,'first'),:));
            
            %The corrolary to the X-unit temporal eigenvector is a unit-norm
            % spatial eigenvector, de-weighted:
            V_dw_iterN_mi_UnitNorm = V_dw_iterN_mi ./ max(abs(V_dw_iterN_mi));
            
            %Store fully iterated mode 1 for this correction level of the data
            % matrix:
            WsiABy_T_w_iterN_mAll(:,i_mode) = T_w_iterN_mi(:,end);
            WsiABy_V_w_iterN_mAll(:,i_mode) = V_w_iterN_mi(:,end);
            WsiABy_L_w_iterN_mAll(:,i_mode) = diag(L_w_iterN_mi);%size of stored part is [number of returned modes in the EOF solution by 1].
            WsiABy_T_dw_iterN_mAll(:,i_mode) = T_dw_iterN_mi;
            WsiABy_V_dw_iterN_mAll(:,i_mode) = V_dw_iterN_mi;
            WsiABy_T_dw_iterN_mAll_Xunits(:,i_mode) = T_dw_iterN_mi_Xunits;
            WsiABy_V_dw_iterN_mAll_UnitNorm(:,i_mode) = V_dw_iterN_mi_UnitNorm;
            
        end%Loop over each iteration process for the next mode in the line.
        
        %% Save EOF analysis outputs: 
        
        %Save data-matrix-unit-scaled leading temporal eigenvectors -- used for correlations and time series plots:
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WtABn_ri_string 'M1TemporalEigVecs_Deweighted_Xunits.mat'], 'WtABn_T_dw_iterN_mAll_Xunits');
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WsiABy_ri_string 'M1TemporalEigVecs_Deweighted_Xunits.mat'], 'WsiABy_T_dw_iterN_mAll_Xunits');
        
        %Save normalised leading spatial eigenvectors -- used for reconstructions using the data-matrix-unit temporal eigenvectors:
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WtABn_ri_string 'M1SpatialEigVecs_Deweighted_Normed.mat'], 'WtABn_V_dw_iterN_mAll_UnitNorm');
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WsiABy_ri_string 'M1SpatialEigVecs_Deweighted_Normed.mat'], 'WsiABy_V_dw_iterN_mAll_UnitNorm');
        
        %Save de-weighted leading temporal eigenvectors -- used for reconstructions of the full data matrix:
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WtABn_ri_string 'M1TemporalEigVecs_Deweighted.mat'], 'WtABn_T_dw_iterN_mAll');
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WsiABy_ri_string 'M1TemporalEigVecs_Deweighted.mat'], 'WsiABy_T_dw_iterN_mAll');
        
        %Save de-weighted leading spatial eigenvectors -- used for reconstructions of the full data matrix:
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WtABn_ri_string 'M1SpatialEigVecs_Deweighted.mat'], 'WtABn_V_dw_iterN_mAll');
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WsiABy_ri_string 'M1SpatialEigVecs_Deweighted.mat'], 'WsiABy_V_dw_iterN_mAll');
        
        %Save removed means -- used for reconstructions of the full data matrix:
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WtABn_ri_string 'RemovedMeans.mat'], 'removed_bin_data_mean');
        save([output_directory  num2str(i_year) '-' num2str(i_month,'%2.2d') '_' EOF_WsiABy_ri_string 'RemovedMeans.mat'], 'removed_bin_data_mean');
        
        %% Fit sinusoid functions to the discretised line of sight velocities of the EOF modes (and the mean removed prior to EOF analysis), to get cardinal (inertial) velocity vector components:
        %Source: program 'THeMES_AccuracyPaper_AllPlots_draft02.m'. Mean
        % fitting from program 'THeMES_MethodPaper_AllPlots_Draft06.m'.
        
        %Fit sinusoid functions to WtABn EOF eigenvectors:
        WtABn_fitted_spatial_eigvecs_north = NaN(size(bin_coords_colat,1),EOF_options.number_of_iterated_modes);%size [NPC bins by modes].
        WtABn_fitted_spatial_eigvecs_east = NaN(size(bin_coords_colat,1),EOF_options.number_of_iterated_modes);%size [NPC bins by modes].
        for i_mode = 1:EOF_options.number_of_iterated_modes
            %Start loop over bins, for cosine fitting:
            for i_bin = 1:size(bin_coords_colat,1)
                %Use the indices made when binning to find all the spatial
                % records that exist for this bin:
                index_SingleBin = find(bin_fiducial_indices == index_original_bin_fiducials(i_bin,1));%index pertains to columns of orig_data_mean, bin_data_vel, bin_contrib_radars, bin_angular_partitions and bin_fiducial_indices, and rows of V_uw_iterN_mAll.
                
                %Data check for empty bin:
                if(isempty(index_SingleBin))
                    continue%pertains to loop over i_bin.
                end%Conditional: only both performing the fitting if there were any data in this bin.
                
                
                %For ease of reference, extract the velocity eigenvector values
                % and the angular partitions that were used as inputs for this
                % spatial bin:
                ang_partitions_SingleBin = bin_angular_partitions(1,index_SingleBin)';%size after transpose is [number of measurements from all partitions in bin by 1].
                spatial_eigvec_SingleBin = WtABn_V_dw_iterN_mAll_UnitNorm(index_SingleBin,i_mode);
                
                % Function to calculate the sum of residuals for given inputs 'coeffs', for
                % which coeffs(1) = amplitude, coeffs(2) = phase. Note use of cosd instead
                % of cos, since input will be in degrees:
                my_cosfit_fun = @(coeffs) sum((spatial_eigvec_SingleBin - (coeffs(1) .* cosd(ang_partitions_SingleBin - coeffs(2)))).^2);
                
                % Multidimensional unconstrained nonlinear minimization to find
                % best-fitting amplitude and phase of a cosine function:
                guess_values = [1,1];%amplitude, phase.
                [coeffs,fminres] = fminsearch(my_cosfit_fun,guess_values);
                %Format: coeffs(1) = amplitude, coeffs(2) = phase.
                
                %Use phase to determine direction of velocity vector, and
                % calculate north and east components.  The phase is defined
                % clockwise from north, so with trigonometry:
                WtABn_fitted_spatial_eigvecs_north(i_bin,i_mode) = coeffs(1) .* cosd(coeffs(2));%evaluation of the cosine curve at 0 degrees.
                WtABn_fitted_spatial_eigvecs_east(i_bin,i_mode) =  coeffs(1) .* sind(coeffs(2));%evaluation of the cosine curve at 90 degrees.
                
            end%Loop over all NPC bins.
        end%Loop over each EOF mode.
        
        
        %Fit sinusoid functions to WsiABy EOF eigenvectors:
        WsiABy_fitted_spatial_eigvecs_north = NaN(size(bin_coords_colat,1),EOF_options.number_of_iterated_modes);%size [NPC bins by modes].
        WsiABy_fitted_spatial_eigvecs_east = NaN(size(bin_coords_colat,1),EOF_options.number_of_iterated_modes);%size [NPC bins by modes].
        for i_mode = 1:EOF_options.number_of_iterated_modes
            %Start loop over bins, for cosine fitting:
            for i_bin = 1:size(bin_coords_colat,1)
                %Use the indices made when binning to find all the spatial
                % records that exist for this bin:
                index_SingleBin = find(bin_fiducial_indices == index_original_bin_fiducials(i_bin,1));%index pertains to columns of orig_data_mean, bin_data_vel, bin_contrib_radars, bin_angular_partitions and bin_fiducial_indices, and rows of V_uw_iterN_mAll.
                
                %Data check for empty bin:
                if(isempty(index_SingleBin))
                    continue%pertains to loop over i_bin.
                end%Conditional: only both performing the fitting if there were any data in this bin.
                
                
                %For ease of reference, extract the velocity eigenvector values
                % and the angular partitions that were used as inputs for this
                % spatial bin:
                ang_partitions_SingleBin = bin_angular_partitions(1,index_SingleBin)';%size after transpose is [number of measurements from all partitions in bin by 1].
                spatial_eigvec_SingleBin = WsiABy_V_dw_iterN_mAll_UnitNorm(index_SingleBin,i_mode);
                
                % Function to calculate the sum of residuals for given inputs 'coeffs', for
                % which coeffs(1) = amplitude, coeffs(2) = phase. Note use of cosd instead
                % of cos, since input will be in degrees:
                my_cosfit_fun = @(coeffs) sum((spatial_eigvec_SingleBin - (coeffs(1) .* cosd(ang_partitions_SingleBin - coeffs(2)))).^2);
                
                % Multidimensional unconstrained nonlinear minimization to find
                % best-fitting amplitude and phase of a cosine function:
                guess_values = [1,1];%amplitude, phase.
                [coeffs,fminres] = fminsearch(my_cosfit_fun,guess_values);
                %Format: coeffs(1) = amplitude, coeffs(2) = phase.
                
                %Use phase to determine direction of velocity vector, and
                % calculate north and east components.  The phase is defined
                % clockwise from north, so with trigonometry:
                WsiABy_fitted_spatial_eigvecs_north(i_bin,i_mode) = coeffs(1) .* cosd(coeffs(2));%evaluation of the cosine curve at 0 degrees.
                WsiABy_fitted_spatial_eigvecs_east(i_bin,i_mode) =  coeffs(1) .* sind(coeffs(2));%evaluation of the cosine curve at 90 degrees.
                
            end%Loop over all NPC bins.
        end%Loop over each EOF mode.
        
        
        %Fit sinusoid functions to mean removed prior to EOF analysis:
        fitted_removed_bin_data_mean_north = NaN(size(bin_coords_colat,1),1);%size [NPC bins by 1].
        fitted_removed_bin_data_mean_east = NaN(size(bin_coords_colat,1),1);
        fitted_removed_bin_data_mean_ssr_misfit = NaN(size(bin_coords_colat,1),1);
        fitted_removed_bin_data_mean_rms_misfit = NaN(size(bin_coords_colat,1),1);
        for i_bin = 1:559
            %Use the indices made when binning to find all the spatial
            % records that exist for this bin:
            index_SingleBin = find(bin_fiducial_indices == index_original_bin_fiducials(i_bin,1));%index pertains to columns of orig_data_mean, bin_data_vel, bin_contrib_radars, bin_angular_partitions and bin_fiducial_indices, and rows of V_uw_iterN_mAll.
            
            %Data check for empty bin:
            if(isempty(index_SingleBin))
                continue%pertains to loop over i_bin.
            end%Conditional: only both performing the fitting if there were any data in this bin.
            
            %For ease of reference, extract the velocity eigenvector values
            % and the angular partitions that were used as inputs for this
            % spatial bin:
            ang_partitions_SingleBin = bin_angular_partitions(1,index_SingleBin)';%size after transpose is [number of measurements from all partitions in bin by 1].
            removed_mean_SingleBin = removed_bin_data_mean(1,index_SingleBin)';%size after transpose is [number of measurements from all partitions in bin by 1].
            
            %Remove any NaN entries:
            index_NaNs_to_remove = find(isnan(removed_mean_SingleBin));%pertains to rows of ang_partitions_SingleBin and removed_mean_SingleBin;
            ang_partitions_SingleBin(index_NaNs_to_remove) = [];
            removed_mean_SingleBin(index_NaNs_to_remove) = [];
            
            %Note: the partition quadrant information was required in the binning
            % stage to assign the appropriate sign to the binned velocity.
            % After that, negative north became the same as positive south,
            % so the angular partition quadrant information was not stored.
            % Thus, here we fit the cosine to 180 degrees of information:
            % the results from applying to the full 360 degrees of
            % partitions is no different.
            
            % Function to calculate the sum of residuals for given inputs 'coeffs', for
            % which coeffs(1) = amplitude, coeffs(2) = phase. Note use of cosd instead
            % of cos, since input will be in degrees:
            my_cosfit_fun = @(coeffs) sum((removed_mean_SingleBin - (coeffs(1) .* cosd(ang_partitions_SingleBin - coeffs(2)))).^2);
            
            % Multidimensional unconstrained nonlinear minimization to find
            % best-fitting amplitude and phase of a cosine function:
            guess_values = [1,1];%amplitude, phase.
            [coeffs,fminres] = fminsearch(my_cosfit_fun,guess_values);
            %coeffs(1) = amplitude, coeffs(2) = phase.
            
            %Store the solution misfit: this scalar value is equal to the
            % evaluation of the 'my_cosfit_fun' equation with the
            % solution's values for 'coeffs', i.e. it is the sum of squared
            % residuals:
            fitted_removed_bin_data_mean_ssr_misfit(i_bin,1) = fminres;%size [1 by 1].
            
            %Also compute and store the RMS misfit:
            fitted_removed_bin_data_mean_rms_misfit(i_bin,1) = sqrt(sum((removed_mean_SingleBin - (coeffs(1) .* cosd(ang_partitions_SingleBin - coeffs(2)))).^2) ./ length(ang_partitions_SingleBin));%size [1 by 1].
            
            %Remove data if the misfit (for the mean field) is too high:
            if(fminres > (2e6))
                fitted_removed_bin_data_mean_north(i_bin,1) = NaN;
                fitted_removed_bin_data_mean_east(i_bin,1) = NaN;
                continue%pertains to loop over i_bin
            end%Conditional: remove poor fits.
            
            %Use phase to determine direction of velocity vector, and
            % calculate north and east components.  The phase is defined
            % clockwise from north, so via trigonometry:
            fitted_removed_bin_data_mean_north(i_bin,1) = coeffs(1) .* cosd(coeffs(2));%evaluation of the cosine curve at 0 degrees.
            fitted_removed_bin_data_mean_east(i_bin,1) =  coeffs(1) .* sind(coeffs(2));%evaluation of the cosine curve at 90 degrees.
            
        end%Loop over each NPC bin locaiton.
        
        
        
        
        
        
        
        
        
        
    end%Loop over each month in the year.
    
end%Loop over all specified years for which you want the velocities to be processed.

%% 

disp('Program run complete.')






