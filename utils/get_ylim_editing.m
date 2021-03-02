function yL = get_ylim_editing(data, ppm_vec)

f=figure;
plot(ppm_vec, data)
xlim([2.2, 4.2]);
yL = get(gca, 'YLim');
close(f);
end

