function long_electrical_noise(data_dir)



if nargin<1
    data_dir = pwd;
end
d = dir(data_dir);


n=1;
for ii=1:length(d)
    tmp = d(ii);
    if contains(tmp.name,'electrical_noise')
        out(n) = generic_generator_template(tmp);
        out(n).type = 'electrical_noise';
        out(n).plotting_func = @mpqc.plot.electrical_noise;
        n=n+1;
    end
end
end