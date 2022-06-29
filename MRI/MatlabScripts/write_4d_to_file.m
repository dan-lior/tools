% write 4d array to file as sequence of 3d projections along its first index
function write_4d_to_file(four_d_array, output_folder)    
    count = 0;
    for it = 1:size(four_d_array, 1)
        fileID = fopen(strcat(output_folder,int2str(it)), 'w');
        temp = fwrite(fileID, four_d_array(it,:,:,:), class(four_d_array));
        count = count + temp;
        fclose(fileID);
    end
    assert(count == numel(four_d_array))
end
