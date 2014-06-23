function result = read_kerner_timing(rundir)

output_files = dir([rundir, '/OUTPUT_0*']);

if strcmp(output_files(1).name, 'OUTPUT_000')
    output_files = output_files(2:end);
end

nproc = length(output_files);
ntimers = 14;

result.ncalls      = zeros(nproc, ntimers);
result.timepercall = zeros(nproc, ntimers);
result.timetotal   = zeros(nproc, ntimers);
result.timeratio   = zeros(nproc, ntimers);

formatSpec = '%*37s%11f%14f%14f%f%[^\n\r]';

for iproc = 1:nproc
    system(['tail -n 18 ', rundir, '/', output_files(iproc).name, ' > OUTPUT_temp']);

    fid = fopen('OUTPUT_temp', 'r');
    dataArray = textscan(fid, formatSpec, ntimers, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
    fclose(fid);

    result.ncalls(iproc, :)      = dataArray{:, 1};
    result.timepercall(iproc, :) = dataArray{:, 2};
    result.timetotal(iproc, :)   = dataArray{:, 3};
    result.timeratio(iproc, :)   = dataArray{:, 4};

end
    
end