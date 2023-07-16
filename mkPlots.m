myDir = fileparts(mfilename("fullpath"));

fig = 1;
for item = dir(fullfile(myDir, "naga*.*.*.nc"))'
    fig = plotit(fig, fullfile(item.folder, item.name));
end % for
% fig = plotit(fig, fullfile(myDir, "naga1.west.90.nc"), 109.43, "NAGA1", "Westward 0.1m/s 90deg");
% fig = plotit(fig, fullfile(myDir, "naga2.west.90.nc"), 109.62, "NAGA2", "Westward 0.1m/s 90deg");
% fig = plotit(fig, fullfile(myDir, "naga2.east.270.nc"), 111, "NAGA2", "Eastward 0.1m/s 270deg");
% 
% fig = plotit(fig, fullfile(myDir, "naga1.east.toLine.nc"), 111, "NAGA1", "Eastward 0.1m/s To Line");
% 
% fig = plotit(fig, fullfile(myDir, "naga1.west.CC.nc"), 109.43, "NAGA1", "Westward 0.1m/s Current Corr");
% fig = plotit(fig, fullfile(myDir, "naga2.west.CC.nc"), 109.43, "NAGA2", "Westward 0.1m/s Current Corr");

% fig = plotit(fig, fullfile(myDir, "tkw.nc"), 113, "NAGA2", "Eastward 0.1m/s ToLine");

function fig = plotit(fig, fn)
[dirname, basename] = fileparts(fn);
a = osgl_get_netCDF(fn);
tbl = struct();
for name = ["x", "y", "lat", "lon", "stime"]
    tbl.(name) = a.(name);
end % for
tbl = struct2table(tbl);
tbl.grp = findgroups(tbl.stime);
b = rowfun(@numel, tbl, "InputVariables", "stime", "GroupingVariables", "grp", "OutputVariableNames", "n");
fprintf("n %d %d %d %d %d %s\n", quantile(b.n, [0,0.25,0.5,0.75,1]), basename);
y = rowfun(@(x) x, tbl, "InputVariables", "y", "GroupingVariables", "grp", "OutputFormat", "cell");
x = rowfun(@(x) x, tbl, "InputVariables", "x", "GroupingVariables", "grp", "OutputFormat", "cell");
tracks = cellfun(@(x,y) horzcat(x, y), x, y, "UniformOutput", false);

lat0 = tbl.lat(1);
lineName = "NAGA1";
dirName = "Eastward";

if tbl.lat(1) < 15, lineName = "NAGA2"; end
if tbl.lon(1) > 110, dirName = "Westward"; end

if a.current
    suffix = "Current Corr";
elseif a.toLine
    suffix = "To Line";
else
    suffix = sprintf("heading %g", a.heading);
end % if

tit = sprintf("%s %s %gm/s %s", lineName, dirName, a.dzdt, suffix);
xtit = sprintf("Distance from %g (km)", tbl.lon(1));
ytit = sprintf("Distance north from %g (km)", tbl.lat(1));

figure(fig);
fig = fig + 1;
clf;
hold on;
cellfun(@(x) plot(x(:,1)/1000, x(:,2)/1000, "-"), tracks, "UniformOutput", false);
hold off;
grid on;
xlabel(xtit);
ylabel(ytit);
title(tit);

print(fullfile(dirname, append(basename, ".png")), "-dpng");
end % plotit