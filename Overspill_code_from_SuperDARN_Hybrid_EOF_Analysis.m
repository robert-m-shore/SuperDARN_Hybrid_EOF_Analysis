%% Information: saving old code from program SuperDARN_Hybrid_EOF_Analysis.m.


%% %!!!! v5: Start loop over the month's daily fitted velocity files for this radar to read the ascii format files:
%The reason there are two loops over each day of data within
% the month is that I do not trust C to finish writing the ascii
% files before I ask Matlab to read them in,  This way, the
% first file should be written long before I ask Matlab to read
% it.
for i_day = 1:1:month_end_day
    %% Use hardware epoch expiry date to determine which centroid locations we should be using for this day:
    
    %For the date entailed in the loop over year, month and day, find the
    % latest applicable hardware file line, given that the
    % hardware_data_times_year, etc. variables are expiry dates:
    index_last_applicable_hardware_line = find(hardware_data_times_datenum > datenum([i_year  i_month  i_day]),1,'first');%size [1 by 1].
    
    %Extract the pertinent set of centroid positions and metadata:
    centroid_locations_beam_number_single_day = centroid_locations_beam_number{index_last_applicable_hardware_line,1};%size [centroids by 1]. Range is 0 to 15.
    centroid_locations_range_gate_single_day = centroid_locations_range_gate{index_last_applicable_hardware_line,1};%size [centroids by 1]. Range is 1 to max range gate count.
    centroid_locations_GEO_colatitude_single_day = centroid_locations_GEO_colatitude{index_last_applicable_hardware_line,1};%size [centroids by 1].
    centroid_locations_GEO_longitude_single_day = centroid_locations_GEO_longitude{index_last_applicable_hardware_line,1};%size [centroids by 1].
    
    %% Load the velocity data which was just saved in the ascii file:
    
    %Re-specify the filename of ascii fitted velocities:
    SD_fitted_ascii_filename = [data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
        radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '_velocities_GEO_' rotation_version '.dat'];
    
    %Check for the existence of the ascii file of fitted
    % velocities for this radar and this day, and skip the loop
    % if it does not exist:
    if(~exist(SD_fitted_ascii_filename,'file'))
        continue%pertains to loop over i_day.
    end%Conditional: if there's no file to read in, skip the processing of it.
    
    %For at least radar kod, I have noticed that the ascii file
    % may exist, but that it may also be completely empty:
    % check for this here:
    ascii_file_details = dir(SD_fitted_ascii_filename);
    if(ascii_file_details.bytes < 10)
        %In this case, you should take the existence of an
        % empty file to indicate some error in the processing.
        continue%pertains to loop over i_day.
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
        % (which is directed in a positve sense towards the
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
        
        %This is applied in the data binning and EOF program,
        % but we need to apply it here too, even though the
        % results are not saved:
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
        beam_offset_angle_from_QD_north_normalised = beam_offset_angle_from_QD_north ./ ((360 - QD_offset_angle) ./ 90);
        
        %Compute the projection of the line of sight vector
        % onto the F2 and 'F2 + 90 degrees' vectors
        beam_projection_onto_F2 = cos(beam_offset_angle_from_QD_north_normalised .* rad);
        beam_projection_onto_F2orth = sin(beam_offset_angle_from_QD_north_normalised .* rad);
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
    
    %% Concatenate data to usable format and save it to file:
    
    %Concatenate the SM coordinates with the velocities:
    velocity_data = [velocity_data_bulk(:,1:2)  velocity_data_centroid_theta_QD  velocity_data_centroid_phi_QD  velocity_data_bulk(:,3:end)  velocity_data_centroid_QD_quadrant_span  velocity_data_centroid_beam_QD_north_offset];
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
        velocity_data(index_temporal_errors,:) = [];
    end%Conditional: if there are no offending data to remove, don't remove them.
    
    %Write-out the file:
    save([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' ...
        radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '_velocities_QD_' rotation_version '.mat'],'velocity_data');
    
    %% Clean up the files used to produce the matlab binary file:
    
    delete(SD_fitted_ascii_filename);
    
end%Loop over each day in the present month for this present radar, read the ascii files.


%% Check for the existence of any processed files for this radar in this month, and skip the rotation process for this month if they already exist:
                
                %Loop over each day in the month, check for files of type
                % 'radar_YYYYMMDD_velocities_QD_v7.mat', set the matrix you've
                % made to a logical non-number: use NaN to indicate file
                % existence.
                single_radar_single_month_daily_file_existence = zeros(month_end_day,1);
                for i_day = 1:month_end_day
                    %Check for the existence of processed QD-frame files:
                    if(exist([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '_velocities_QD_' rotation_version '.mat'],'file'))
                        single_radar_single_month_daily_file_existence(i_day,1) = NaN;
                    end%Conditional: if the file exists, set this daily element of the existence-checking vector to NaN.
                    
                    %Check for the existence of empty GEO-frame ascii files:
                    % these signify that the program has processed the data for
                    % radar and that day, but found inconsistencies in the
                    % records which led to the data not being put into the
                    % ascii file. Since the second loop over i_day in this
                    % program skips the .mat file processing for empty GEO.dat
                    % files, it also fails to delete them (refer to the line of
                    % code 'if(ascii_file_details.bytes < 10)').  Thus the
                    % empty GEO.dat files remain as a record of processed
                    % radar-days.  These take time to reprocess and would
                    % presumably still result in no extra data, so we skip them
                    % in future processing by enforcing the same 'skip this
                    % iteration of the loop' effect, as if the processed QD.mat
                    % file existed:
                    if(exist([data_directory 'Velocities/' radar_identifiers{i_r,1} '/' radar_identifiers{i_r,1} '_' num2str(i_year) num2str(i_month,'%2.2d') num2str(i_day,'%2.2d') '_velocities_GEO_' rotation_version '.dat'],'file'))
                        single_radar_single_month_daily_file_existence(i_day,1) = NaN;
                    end%Conditional: if the file exists, set this daily element of the existence-checking vector to NaN.
                end%Loop over day of month.
                
                %If any of this month's data for this radar have already been
                % rotated to QD and thus now exist in Matlab binary file
                % format, skip the entire processing for this month for this
                % radar:
                if(any(isnan(single_radar_single_month_daily_file_existence)))
                    disp(['Processed files found for ' radar_identifiers{i_r,1} ' in ' num2str(i_year) '-' num2str(i_month,'%2.2d') ', skipping this month for this radar.'])
                    continue%pertains to loop over i_r.
                    %This is fine since all the process here happens on a
                    % monthly repeat scale, and this allows us to effectively
                    % prevent the re-processing and over-writing of months that
                    % have already been processed. Hard-coded for this program
                    % version.
                end%Conditional: skip this month of processing if it's already been done.
                

