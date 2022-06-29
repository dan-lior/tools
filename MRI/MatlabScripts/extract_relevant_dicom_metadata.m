function [metadata] = extract_relevant_dicom_metadata(input_folder)

    directory_info = struct2cell(dir(input_folder));
    metadata.num_total_images = size(directory_info,2)-2;

    name_of_first_file = cell2mat(directory_info(1,3)); % 1 indicates the name field, 3 indicates the 3rd entry of the directory (first and second being '.' and '..')
    dicom_info_first_file = dicominfo(strcat(input_folder,name_of_first_file));
    % assume that all the dicom files in this directory share this dicom information

    metadata.height = dicom_info_first_file.Height;
    metadata.width = dicom_info_first_file.Width;

    metadata.sequence_name = dicom_info_first_file.SequenceName;

    if isfield(dicom_info_first_file, 'CardiacNumberOfImages') %number of images in a cardiac cycle
        metadata.num_timesteps = dicom_info_first_file.CardiacNumberOfImages;
        assert (rem(metadata.num_total_images, metadata.num_timesteps) == 0, 'error: expected num images to be divisible by num timesteps');
        metadata.depth = metadata.num_total_images/metadata.num_timesteps;    
    else

        metadata.num_acquisitions = 6;
        metadata.acquisition_of_interest = 4;
        assert (rem(metadata.num_total_images, metadata.num_acquisitions) == 0, 'error: expected num images to be divisible by num timesteps');
        metadata.depth = metadata.num_total_images/metadata.num_acquisitions;    
    end
    
    % geometry
    metadata.slice_thickness = dicom_info_first_file.SliceThickness;

    pixel_spacing.x = dicom_info_first_file.PixelSpacing(1);
    pixel_spacing.y = dicom_info_first_file.PixelSpacing(2);
    metadata.pixel_spacing = pixel_spacing;

    pos.x = dicom_info_first_file.ImagePositionPatient(1);
    pos.y = dicom_info_first_file.ImagePositionPatient(2);
    pos.z = dicom_info_first_file.ImagePositionPatient(3);
    metadata.image_position = pos;

    orientation.Xx = dicom_info_first_file.ImageOrientationPatient(1);
    orientation.Xy = dicom_info_first_file.ImageOrientationPatient(2);
    orientation.Xz = dicom_info_first_file.ImageOrientationPatient(3);
    orientation.Yx = dicom_info_first_file.ImageOrientationPatient(4);
    orientation.Yy = dicom_info_first_file.ImageOrientationPatient(5);
    orientation.Yz = dicom_info_first_file.ImageOrientationPatient(6);
    metadata.image_orientation = orientation;

    % data representation    
    sample = dicomread(dicom_info_first_file);
    metadata.data_type = class(sample);
    metadata.bit_depth = dicom_info_first_file.BitDepth;

end