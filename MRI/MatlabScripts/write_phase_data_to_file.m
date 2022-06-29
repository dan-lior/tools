function write_phase_data_to_file(metadata, input_folder, output_folder)

    phase1_store = zeros(metadata.num_timesteps, metadata.width,metadata.height,  metadata.depth, metadata.data_type);
    
    slice_map = [];
    time_map = [];
    
    directory_info = struct2cell(dir(input_folder));

    for ii = 1:metadata.num_total_images
       fprintf('%s', '.');
%        disp(ii);
    
      filename = cell2mat(directory_info(1,ii+2));
      dicom_info = dicominfo(strcat(input_folder,filename));
    
%       fprintf('%s %f\n', 'trigger time = ', dicom_info.TriggerTime);
%       fprintf('%s %f\n\n', 'slice location = ', dicom_info.SliceLocation);
      
      % figure out slice and time indices
      slice_ind = find(slice_map == dicom_info.SliceLocation);
      if isempty(slice_ind) 
          slice_map = [slice_map dicom_info.SliceLocation]; 
          slice_ind = length(slice_map); 
      end
      
      time_ind = find(time_map == dicom_info.TriggerTime); 
      if isempty(time_ind) 
          time_map = [time_map dicom_info.TriggerTime]; 
          time_ind = length(time_map); 
      end
    
      X_p1 = dicomread(dicom_info);
      phase1_store(time_ind,:,:,slice_ind) = X_p1'; % deliberate transpose!

    end
    
    % sort the times and slices
    [~, sorted_time_ind] = sort(time_map);
    [~, sorted_slice_ind] = sort(slice_map);
    
    phase1_store = phase1_store(sorted_time_ind,:,:,sorted_slice_ind);
    
    output_subfolder = strcat(output_folder, 'time_slices/');
    mkdir(output_subfolder);    
    
    write_4d_to_file(phase1_store, output_subfolder);

end



