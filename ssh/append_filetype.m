
function file_list = append_filetype(name_array, append_string)

L = length(name_array);

for l = 1:L
     file_list{l} = [name_array{l} append_string];
end


end