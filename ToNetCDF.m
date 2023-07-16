% Convert Javier's ROMS Matlab files to NetCDF files
%
% July-2023, Pat Welch, pat@mousebrains.com

myDir = fileparts(mfilename("fullpath"));
items = dir(fullfile(myDir, "roms*.mat"));
naga1 = {};
naga2 = {};

for index = 1:numel(items)
    item = items(index);
    fn = fullfile(item.folder, item.name);
    fprintf("Loading %s\n", fn);
    a = load(fn);
    name = string(fieldnames(a));
    a = a.(name);
    if endsWith(name, "1")
        naga1{end+1} = a;
    else
        naga2{end+1} = a;
    end % if
end % for

toNetCDF(fullfile(myDir, "naga1.nc"), naga1)
toNetCDF(fullfile(myDir, "naga2.nc"), naga2)

function toNetCDF(fn, a)
arguments
    fn string
    a cell
end % arguments
if exist(fn, "file"), delete(fn); end
ncid = netcdf.create(fn, bitor(netcdf.getConstant("CLOBBER"), netcdf.getConstant("NETCDF4")));
tDim = netcdf.defDim(ncid, "t", numel(a{1}.time));
dDim = netcdf.defDim(ncid, "depth", numel(a));
lDim = netcdf.defDim(ncid, "lon", numel(a{1}.lon));
tID = createVar(ncid, "t", tDim);
dID = createVar(ncid, "depth", dDim);
lonID = createVar(ncid, "lon", lDim);
latID = createVar(ncid, "lat", lDim);
uID = createVar(ncid, "u", [dDim, tDim, lDim]);
vID = createVar(ncid, "v", [dDim, tDim, lDim]);
netcdf.endDef(ncid);

u = nan(numel(a), size(a{1}.u,1), size(a{1}.u,2));
v = nan(size(u));

depth = nan(size(a));

for i = 1:numel(a)
    b = a{i};
    if i == 1
        netcdf.putVar(ncid, tID, posixtime(datetime(b.time, "ConvertFrom", "datenum")));
        netcdf.putVar(ncid, lonID, b.lon);
        netcdf.putVar(ncid, latID, b.lat);
    end % if 1
    depth(i) = b.depth;
    u(i,:,:) = b.u;
    v(i,:,:) = b.v;
end % for i
[depth,ix] = unique(depth);
netcdf.putVar(ncid, dID, depth);
netcdf.putVar(ncid, uID, u(ix,:,:));
netcdf.putVar(ncid, vID, v(ix,:,:));
netcdf.close(ncid);
end % toNetCDF

%%

function ident = createVar(ncid, name, dims)
xtype = netcdf.getConstant("NC_DOUBLE");
ident = netcdf.defVar(ncid, name, xtype, dims);
netcdf.defVarDeflate(ncid, ident, false, true, 5);
if name == "t"
    netcdf.putAtt(ncid, ident, "units", "seconds since 1970-01-01");
    netcdf.putAtt(ncid, ident, "calendar", "proleptic_gregorian");
end % if name
end % createVar