function extract_4dmri_data_and_save_to_file(json)

jstruct = jsondecode(json);

output_folder = jstruct.output_folder;


phase1_folder = jstruct.phase1_folder;
phase2_folder = jstruct.phase2_folder;
phase3_folder = jstruct.phase3_folder;
mag_folder = jstruct.mag_folder;


% the sequence names associate phase 1,2,3 with directions "in", "ap" and "fh"
% all other relevant metadata for the 4dflow ("WIP") is assumed to be common to the phase1,
% phase2, phase3 and mag folders. We choose one arbitrarily here.
metadata = extract_relevant_dicom_metadata(mag_folder);

sequence_name_phase1 = extract_relevant_dicom_metadata(phase1_folder).sequence_name;
sequence_name_phase2 = extract_relevant_dicom_metadata(phase2_folder).sequence_name;
sequence_name_phase3 = extract_relevant_dicom_metadata(phase3_folder).sequence_name;
[speed1, suffix1] = extract_speed_and_suffix_from_sequence_name(sequence_name_phase1);
[speed2, suffix2] = extract_speed_and_suffix_from_sequence_name(sequence_name_phase2);
[speed3, suffix3] = extract_speed_and_suffix_from_sequence_name(sequence_name_phase3);

if ~(strcmp(speed1,speed2) && strcmp(speed2,speed3))
    error('encoding speed not consistent among phases')
end

metadata.encoding_speed = speed1;

if (suffix1 == "in" && suffix2 == "ap" && suffix3 == "fh")
    phasex_folder = phase2_folder;
    phasey_folder = phase3_folder;
    phasez_folder = phase1_folder;
elseif (suffix1 == "in" && suffix2 == "fh" && suffix3 == "ap")
    phasex_folder = phase3_folder;
    phasey_folder = phase2_folder;
    phasez_folder = phase1_folder;
elseif (suffix1 == "ap" && suffix2 == "in" && suffix3 == "fh")
    phasex_folder = phase1_folder;
    phasey_folder = phase3_folder;
    phasez_folder = phase2_folder;
elseif (suffix1 == "ap" && suffix2 == "fh" && suffix3 == "in")
    phasex_folder = phase1_folder;
    phasey_folder = phase2_folder;
    phasez_folder = phase3_folder;
elseif (suffix1 == "fh" && suffix2 == "in" && suffix3 == "ap")
    phasex_folder = phase3_folder;
    phasey_folder = phase1_folder;
    phasez_folder = phase2_folder;
elseif (suffix1 == "fh" && suffix2 == "ap" && suffix3 == "in")
    phasex_folder = phase2_folder;
    phasey_folder = phase1_folder;
    phasez_folder = phase3_folder;
else
    error('sequence name suffixes should be some permutation of in, ap and fh')
end

disp('writing relevant DICOM metadata to file');
j_wip = jsonencode(metadata, 'PrettyPrint', true);
output_filename = strcat('../', output_folder, 'metadata.json');
fileID = fopen(output_filename, 'w');
fwrite(fileID, j_wip);
fclose(fileID);
disp('finished writing relevant DICOM metadata to file');

disp('extracting dicom data from file and writing to file');
write_phase_data_to_file(metadata, phasex_folder, strcat('../', output_folder, 'phasex_data/'));
write_phase_data_to_file(metadata, phasey_folder, strcat('../', output_folder, 'phasey_data/'));
write_phase_data_to_file(metadata, phasez_folder, strcat('../', output_folder, 'phasez_data/'));
write_phase_data_to_file(metadata, mag_folder, strcat('../', output_folder, 'mag_data/'));
disp('finished extracting dicom data from file and writing to file');


end
