function [d,out] = long_electrical_noise(data_dir)



if nargin<1
    data_dir = pwd;
end
% d = dir(fullfile(data_dir,'**\*.*')%,contains('electrical_noise'));
d = dir(fullfile(data_dir))%,'**\*.*')%,contains('electrical_noise'));
d(~[d.isdir]) = []; % removes the non directories

subD = cell(length(d));
for ii=1:length(d)
    % tmp = d(ii);
    subDir = d(ii).name;
    subD{ii} = dir(fullfile(subDir));%, contains('electrical_noise'));
    if contains(subD{ii}.name, 'electrical_noise')
        out(n).type = 'electrical_noise';
        out(n).plotting_func = @mpqc.plot.electrical_noise;
    end
    % dir(data_dir,d.name)
% if contains(d.name,'electrical_noise')
%     disp(d.name)
% end
end
% n =1;
% if contains(d.name, 'electrical_noise')
%     out(n).type = 'electrical_noise';
%     out(n).plotting_func = @mpqc.plot.electrical_noise;
% end

n=1;
% for ii=1:length(d)
%     tmp = d(ii);
%     if contains(tmp.name,'electrical_noise')
%         % out(n) = generic_generator_template(tmp);
%         out(n).type = 'electrical_noise';
%         out(n).plotting_func = @mpqc.plot.electrical_noise;
%         n=n+1;
%     else
%         disp('No files')
%     end
% end

% t_im = single(imstack(:,:,ii));
%         [n,x] = hist(t_im(:),100);
%         a=area(x,n);
% 
% 
%         a.EdgeColor=[0,0,0.75];
%         a.FaceColor=[0.5,0.5,1];
%         a.LineWidth=2;

end