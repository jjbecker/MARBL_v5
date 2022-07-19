function [my_cell] = river_flux_names(lciso_on)
%river_flux_names MARBL "short name" of all riverine flux

my_cell = tracer_names(lciso_on);
for idx = 1:length(my_cell)
    my_cell{idx} = strcat(my_cell{idx}, '_RIV_FLUX');
end

end

