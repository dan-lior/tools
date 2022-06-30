function write_intensity_data_to_file(metadata, input_folder, out_filename)

%    intensity_store = zeros(metadata.height, metadata.width, metadata.depth, metadata.data_type);
    intensity_store = zeros(metadata.width, metadata.height, metadata.depth, metadata.data_type);
    
    slice_map = [];
    
    directory_info = struct2cell(dir(input_folder));

    count = 0;

    for ii = 1:metadata.num_total_images
    
        filename = cell2mat(directory_info(1,ii+2));
        dicom_info = dicominfo(strcat(input_folder,filename));
        
        if (dicom_info.AcquisitionNumber == 4)
            count = count + 1;

            slice_ind = find(slice_map == dicom_info.SliceLocation);
            if isempty(slice_ind) 
                slice_map = [slice_map dicom_info.SliceLocation]; 
                slice_ind = length(slice_map); 
            end
              
            X_p1 = dicomread(dicom_info);
            intensity_store(:,:,count) = X_p1';    % deliberate transpose!
        end
    end

    % sort the times and slices
    [~, sorted_slice_ind] = sort(slice_map);
    
    intensity_store = intensity_store(:,:,sorted_slice_ind);

    write_3d_to_file(intensity_store, out_filename);

end



