% write 3d array to file
function write_3d_to_file(three_d_array, filename)    
    fileID = fopen(filename, 'w');
    fwrite(fileID, three_d_array, class(three_d_array)); % first dimension is the fastest
    fclose(fileID);
end
