function  [speed, suffix] = extract_speed_and_suffix_from_sequence_name(sequence_name)
% for sequence name of the form .._vxxxyy.., xxx is the encoding speed and the suffix is yy
  
    n = size(sequence_name, 2);
    suffix = sequence_name(n-1:n);
    speed = sequence_name(n-4:n-2);

end
