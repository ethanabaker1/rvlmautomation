%Light Study Data Automation Tool
%Ethan Baker, eab109(at)pitt.edu
%Yates Lab, University of Pittsburgh
%Last Updated: 11/14/14 @ 9:40AM;
%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT: This software requires the following dependencies: Signal
%Processing Toolbox, risingEdge.m, risingEdge2.m. M files must be loaded
%into the default MATLAB dir. All SPIKE2 files must be exported as .MAT
%files with the following naming scheme: R(1-6).MAT. HR Channel must be
%generated in SPIKE2 before exporting to MATLAB. 
%%%%%%%%%%%%%%%%%%%%%%

%Get files
path = uigetdir('Select directory with .mat files'); %GUI select dir w/ matlab files, can have other files in dir
s=what(path); %Gets info about selected dir
matfiles=s.mat; % gets .MAT file names
cd(path);
for a=1:numel(matfiles); %Load .MAT files
    load((char(matfiles(a))));
end

%Run Calculations on R1, if it exists.
if exist('R1_Ch14')==1
    %Generate light-time matrix
    light_time_start = R1_Ch14.start; 
    light_time_interval = R1_Ch14.interval;
    light_time_length = R1_Ch14.length;
    light_values = R1_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R1_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R1_Ch12.start;
    table_time_interval = R1_Ch12.interval;
    table_time_length = R1_Ch12.length;
    table_values = R1_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R1_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R1_Ch133.start;
    carotid_time_interval = R1_Ch13.interval;
    carotid_time_length = R1_Ch13.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R1_Ch13.values);

    %Declare HR times
    HR_times = R1_Ch4.times; %Renames time

    %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks); %Determines number of elements in lightpeaks
    
   %The following if statements create the bounds of regions of interest
   %depednent on how many light peaks exist. 
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts); %Finds size of tilt_starts array
    %Create regions dependent on the number of tilt onsets in the file
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R1_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R1_1_table = array2table(R1_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2'); %Check if a second tilt exists in the file
    b = exist('tilt_onset_2'); %Check if a onset exists for the second tilt
    if a+b==2 %If both exist do the following:
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R1_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R1_2_table = array2table(R1_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3'); %Check if 3rd tilt exists
    b = exist('tilt_onset_3'); %Check if 3rd onset time exists
    if a+b==2 %If both exist, do the calculations:
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R1_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R1_3_table = array2table(R1_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end
    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4'); %Check if 4th tilt exists
    b = exist('tilt_onset_4');
    if a+b==2 %If tilt exists & onset time exists, do the calculations:
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R1_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R1_4_table = array2table(R1_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R1_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R1_5_table = array2table(R1_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end
end

%Run Calculations on R2, if it exists.
if exist('R2_Ch14')==1
    %Generate light-time matrix
    light_time_start = R2_Ch14.start; 
    light_time_interval = R2_Ch14.interval;
    light_time_length = R2_Ch14.length;
    light_values = R2_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R2_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R2_Ch12.start;
    table_time_interval = R2_Ch12.interval;
    table_time_length = R2_Ch12.length;
    table_values = R2_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R2_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R2_Ch12.start;
    carotid_time_interval = R2_Ch12.interval;
    carotid_time_length = R2_Ch12.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R2_Ch13.values);

    %Declare HR times
    HR_times = R2_Ch4.times; %Renames time

    %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks);
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts);
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R2_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R2_1_table = array2table(R2_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})
   
    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2');
    b = exist('tilt_onset_2');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R2_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R2_2_table = array2table(R2_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3');
    b = exist('tilt_onset_3');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R2_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R2_3_table = array2table(R2_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4');
    b = exist('tilt_onset_4');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R2_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R2_4_table = array2table(R2_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R2_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R2_5_table = array2table(R2_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end
    
end

%Run Calculations on R3, if it exists.
if exist('R3_Ch14')==1
    %Generate light-time matrix
    light_time_start = R3_Ch14.start; 
    light_time_interval = R3_Ch14.interval;
    light_time_length = R3_Ch14.length;
    light_values = R3_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R3_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R3_Ch12.start;
    table_time_interval = R3_Ch12.interval;
    table_time_length = R3_Ch12.length;
    table_values = R3_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R3_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R3_Ch13.start;
    carotid_time_interval = R3_Ch13.interval;
    carotid_time_length = R3_Ch13.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R3_Ch13.values);

    %Declare HR times
    HR_times = R3_Ch4.times; %Renames time

   %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks);
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts);
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R3_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R3_1_table = array2table(R3_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2');
    b = exist('tilt_onset_2');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R3_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R3_2_table = array2table(R3_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3');
    b = exist('tilt_onset_3');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R3_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R3_3_table = array2table(R3_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4');
    b = exist('tilt_onset_4');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R3_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R3_4_table = array2table(R3_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R3_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R3_5_table = array2table(R3_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end
end

%Run Calculations on R4, if it exists.
if exist('R4_Ch14') ==1
    %Generate light-time matrix
    light_time_start = R4_Ch14.start; 
    light_time_interval = R4_Ch14.interval;
    light_time_length = R4_Ch14.length;
    light_values = R4_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R4_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R4_Ch12.start;
    table_time_interval = R4_Ch12.interval;
    table_time_length = R4_Ch12.length;
    table_values = R4_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R4_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R4_Ch13.start;
    carotid_time_interval = R4_Ch13.interval;
    carotid_time_length = R4_Ch13.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R4_Ch13.values);

    %Declare HR times
    HR_times = R4_Ch4.times; %Renames time

    %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks);
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts);
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R4_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R4_1_table = array2table(R4_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2');
    b = exist('tilt_onset_2');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R4_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R4_2_table = array2table(R4_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3');
    b = exist('tilt_onset_3');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R4_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R4_3_table = array2table(R4_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4');
    b = exist('tilt_onset_4');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R4_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R4_4_table = array2table(R4_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R4_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R4_5_table = array2table(R4_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end
end

%Run Calculations on R5, if it exists.
if exist('R5_Ch14')==1
    %Generate light-time matrix
    light_time_start = R5_Ch14.start; 
    light_time_interval = R5_Ch14.interval;
    light_time_length = R5_Ch14.length;
    light_values = R5_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R5_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R5_Ch12.start;
    table_time_interval = R5_Ch12.interval;
    table_time_length = R5_Ch12.length;
    table_values = R5_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R5_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R5_Ch13.start;
    carotid_time_interval = R5_Ch13.interval;
    carotid_time_length = R5_Ch13.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R5_Ch13.values);

    %Declare HR times
    HR_times = R5_Ch4.times; %Renames time

    %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks);
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts);
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R5_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R5_1_table = array2table(R5_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2');
    b = exist('tilt_onset_2');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R5_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R5_2_table = array2table(R5_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3');
    b = exist('tilt_onset_3');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R5_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R5_3_table = array2table(R5_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'})

    end

    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4');
    b = exist('tilt_onset_4');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R5_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R5_4_table = array2table(R5_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});
    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R5_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R5_5_table = array2table(R5_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});
    end
end

%Run Calculations on R6, if it exists.
if exist('R6_Ch14')==1
    %Generate light-time matrix
    light_time_start = R6_Ch14.start; 
    light_time_interval = R6_Ch14.interval;
    light_time_length = R6_Ch14.length;
    light_values = R6_Ch14.values;
    light_end_time = light_time_length/(1/light_time_interval);
    light_timematrix = (light_time_start:light_time_interval:light_end_time); %Generates an array with time based on start time, end time, recording interval.
    light_timematrix = light_timematrix' ; %Transposes light array.
    light_time_data = horzcat(light_timematrix,R6_Ch14.values); %Constructs 2 column array (Time|Value)

    %Generate Table-Time matrix
    table_time_start = R6_Ch12.start;
    table_time_interval = R6_Ch12.interval;
    table_time_length = R6_Ch12.length;
    table_values = R6_Ch12.values;
    table_end_time = table_time_length/(1/table_time_interval);
    table_timematrix = (table_time_start:table_time_interval:table_end_time);
    table_timematrix = table_timematrix' ;
    table_time_data = horzcat(table_timematrix,R6_Ch12.values);

    %Generate Carotid-Time Matricies
    carotid_time_start = R6_Ch13.start;
    carotid_time_interval = R6_Ch13.interval;
    carotid_time_length = R6_Ch13.length;
    carotid_end_time = carotid_time_length/(1/carotid_time_interval);
    carotid_timematrix = (carotid_time_start:carotid_time_interval:carotid_end_time);
    carotid_timematrix = carotid_timematrix' ;
    carotid_time_data = horzcat(carotid_timematrix,R6_Ch13.values);

    %Declare HR times
    HR_times = R6_Ch4.times; %Renames time

    %Find Peaks in Light Waveform
    smoothData = sgolayfilt(light_values,7,21);
    invSmoothData = -smoothData;
    %%%%%%%
    threshold = invSmoothData(57)+ 0.0185 ;
    data=invSmoothData;
    offsetData = [data(2:end); NaN];
    lightpeaks = find(data < threshold & offsetData > threshold);
    %lightpeaks = risingEdge2(light_timematrix,invSmoothData);
    [m,n] = size(lightpeaks);
    if m==10
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1), lightpeaks(9,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T5 = lightonsets(5,1);
        T5 = light_timematrix(T5,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
        T5_off = T5+2.03;
    elseif m==8
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1), lightpeaks(7,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T4 = lightonsets(4,1);
        T4 =light_timematrix(T4,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
        T4_off = T4+2.03;
    elseif m==6
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1),lightpeaks(5,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T3 = lightonsets(3,1);
        T3 = light_timematrix(T3,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
        T3_off = T3+2.03;
    elseif m==4
        lightonsets = [lightpeaks(1,1),lightpeaks(3,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T2 = lightonsets(2,1);
        T2 = light_timematrix(T2,1);
        T1_off = T1+2.03;
        T2_off = T2+2.03;
    elseif m==2
        lightonsets = [lightpeaks(1,1)];
        lightonsets = lightonsets';
        T1 = lightonsets(1,1);
        T1 = light_timematrix(T1,1);
        T1_off = T1+2.03;
    end

    %Find Table Onset
    tilt_starts = risingEdge(table_timematrix,table_values); %Uses the risingEdge function to find table tilt onset. 
    tilt_starts = tilt_starts/(1/table_time_interval); %Converts to time in s
    [m,n] = size(tilt_starts);
    if m==1
        tilt_onset_1 =tilt_starts(1,1); %Gets 1st table tilt onset time and stores as tilt_onset_1
    elseif m==2
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
    elseif m==3
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
    elseif m==4
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
    elseif m==5
        tilt_onset_1 =tilt_starts(1,1);
        tilt_onset_2 = tilt_starts(2,1);
        tilt_onset_3 = tilt_starts(3,1);
        tilt_onset_4 = tilt_starts(4,1);
        tilt_onset_5 = tilt_starts(5,1);
    end

    %Calculate Regions of Interest 

    %Segment 1 - 5 seconds before Light Onset
    seg1_start = T1 - 5;
    seg1_end = T1;
    x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg1_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
    beat_num = numel(z); 
    HR_interval= seg1_end - seg1_start;
    seg1_avgHR = beat_num/HR_interval;

    %Segment 2a - 5 Seconds After Light Offset
    seg2a_start = T1_off;
    seg2a_end = T1_off + 5;
    x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2a_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
    beat_num = numel(z);
    HR_interval = seg2a_end - seg2a_start;
    seg2a_avgHR = beat_num/HR_interval;

    %Segment 2 - 5seconds before tilt onset
    seg2_start = tilt_onset_1 - 5;
    seg2_end = tilt_onset_1;
    x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg2_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
    beat_num = numel(z);
    HR_interval = seg2_end - seg2_start;
    seg2_avgHR = beat_num/HR_interval;

    %Segment 3 - 0-5s after tilt onset
    seg3_start = tilt_onset_1;
    seg3_end = tilt_onset_1 + 5;
    x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg3_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
    beat_num = numel(z);
    HR_interval = seg3_end - seg3_start;
    seg3_avgHR = beat_num/HR_interval;

    %Segment 4 - 5-10s after tilt onset
    seg4_start = tilt_onset_1 + 5;
    seg4_end = tilt_onset_1 + 10;
    x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
    y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
    seg4_mean_flow = mean(y);
    z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
    beat_num = numel(z);
    HR_interval = seg4_end - seg4_start;
    seg4_avgHR = beat_num/HR_interval;

    %Generate Output Matrix
    R6_Output1 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
    R6_1_table = array2table(R6_Output1,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});

    %Calculate Regions of Interest for second trial, if it exists.

    a = exist('T2');
    b = exist('tilt_onset_2');
    c=8;
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T2 - 5;
        seg1_end = T2;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T2_off;
        seg2a_end = T2_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_2 - 5;
        seg2_end = tilt_onset_2;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_2;
        seg3_end = tilt_onset_2 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_2 + 5;
        seg4_end = tilt_onset_2 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R6_Output2 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R6_2_table = array2table(R6_Output2,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});

    end

    %Calculate Regions of Interest for third trial, if it exists.

    a = exist('T3');
    b = exist('tilt_onset_3');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T3 - 5;
        seg1_end = T3;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T3_off;
        seg2a_end = T3_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_3 - 5;
        seg2_end = tilt_onset_3;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_3;
        seg3_end = tilt_onset_3 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_3 + 5;
        seg4_end = tilt_onset_3 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R6_Output3 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R6_3_table = array2table(R6_Output3,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});

    end

    %Calculate Regions of Interest for 4th trial, if it exists.

    a = exist('T4');
    b = exist('tilt_onset_4');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T4 - 5;
        seg1_end = T4;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T4_off;
        seg2a_end = T4_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_4 - 5;
        seg2_end = tilt_onset_4;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_4;
        seg3_end = tilt_onset_4 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_4 + 5;
        seg4_end = tilt_onset_4 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R6_Output4 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R6_4_table = array2table(R6_Output4,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});

    end

    %Calculate Regions of Interest for 5th trial, if it exists.

    a = exist('T5');
    b = exist('tilt_onset_5');
    if a+b==2
        %Segment 1 - 5 seconds before Light Onset
        seg1_start = T5 - 5;
        seg1_end = T5;
        x = find((carotid_time_data(:,1)>=seg1_start) & (carotid_time_data(:,1)<=seg1_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg1_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg1_start) & (HR_times(:,1)<=seg1_end));
        beat_num = numel(z); 
        HR_interval= seg1_end - seg1_start;
        seg1_avgHR = beat_num/HR_interval;

        %Segment 2a - 5 Seconds After Light Offset
        seg2a_start = T5_off;
        seg2a_end = T5_off + 5;
        x = find((carotid_time_data(:,1)>=seg2a_start) & (carotid_time_data(:,1)<=seg2a_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2a_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2a_start) & (HR_times(:,1)<=seg2a_end));
        beat_num = numel(z);
        HR_interval = seg2a_end - seg2a_start;
        seg2a_avgHR = beat_num/HR_interval;

        %Segment 2 - 5seconds before tilt onset
        seg2_start = tilt_onset_5 - 5;
        seg2_end = tilt_onset_5;
        x = find((carotid_time_data(:,1)>=seg2_start) & (carotid_time_data(:,1)<=seg2_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg2_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg2_start) & (HR_times(:,1)<=seg2_end));
        beat_num = numel(z);
        HR_interval = seg2_end - seg2_start;
        seg2_avgHR = beat_num/HR_interval;

        %Segment 3 - 0-5s after tilt onset
        seg3_start = tilt_onset_5;
        seg3_end = tilt_onset_5 + 5;
        x = find((carotid_time_data(:,1)>=seg3_start) & (carotid_time_data(:,1)<=seg3_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg3_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg3_start) & (HR_times(:,1)<=seg3_end));
        beat_num = numel(z);
        HR_interval = seg3_end - seg3_start;
        seg3_avgHR = beat_num/HR_interval;

        %Segment 4 - 5-10s after tilt onset
        seg4_start = tilt_onset_5 + 5;
        seg4_end = tilt_onset_5 + 10;
        x = find((carotid_time_data(:,1)>=seg4_start) & (carotid_time_data(:,1)<=seg4_end)); %Returns index values of desired blood flow points
        y = carotid_time_data(x,2); %Returns the 2nd column in datatable corresponding to indicies(x)
        seg4_mean_flow = mean(y);
        z = find((HR_times(:,1)>=seg4_start) & (HR_times(:,1)<=seg4_end));
        beat_num = numel(z);
        HR_interval = seg4_end - seg4_start;
        seg4_avgHR = beat_num/HR_interval;

        %Generate Output Matrix
        R6_Output5 = [seg1_avgHR,seg2a_avgHR,seg2_avgHR,seg3_avgHR,seg4_avgHR;seg1_mean_flow,seg2a_mean_flow,seg2_mean_flow,seg3_mean_flow,seg4_mean_flow];
        R6_5_table = array2table(R6_Output5,'VariableNames',{'Seg1' 'Seg2a' 'Seg2' 'Seg3' 'Seg4'}, 'RowNames',{'HR' 'Flow'});

    end
end

%Exports to Excel ----Doesn't quite work properly (yet)
%Another Note: MATLAB only supports writing to Excel files on Windows. 
%filename = 'output.xlsx';
%sheet = 1;
%xlRange = 'A1';
%A = R1_Output1;
%xlswrite(filename,A);
