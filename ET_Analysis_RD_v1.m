%% ============================================================
% ETa Processing from Landsat + ETo
% Author: INIFAP-CENID RASPA SIJJ
% Description:
%   - Reads ETa GeoTIFF images
%   - Matches ETo data
%   - Computes Kc
%   - Interpolates Kc to daily scale
%   - Computes daily ETa and cumulative ETa
% ============================================================

clear; close all; clc;

%% ===================== 1. INPUTS ============================
files = dir('./ET_data_/AgrET_*.tif');
fileETo = 'ET0_data.xlsx';

FechaInicioETo   = datetime(2025,01,01);
FechaInicioETa   = datetime(2025,02,26);
FechaSiembra     = datetime(2025,03,08);
FechaFinal       = datetime(2025,06,28);

RD = 1; % 1 = reanalysis | 0 = measured

Kcmin = 0.3; %Maize
Kcmax = 1.2; %Maize

%% ===================== 2. MATCH ETa - ETo ===================
n = length(files);
fechasETa = NaT(n,1);

for i = 1:n
    nombre = files(i).name;
    tokens = regexp(nombre, 'AgrET_(\d+)_(\d+)_(\d+)', 'tokens');
    fecha_num = str2double(tokens{1});
    fechasETa(i) = datetime(fecha_num(1), fecha_num(2), fecha_num(3));
end

% Leer ETo
tablaETo = readtable(fileETo);
tablaETo.Date = datetime(tablaETo.Date, ...
    'InputFormat','dd-MMM-yy','Locale','es_ES');

% Selección de columna
if RD == 1
    ETo_col = tablaETo.ET0R_mmDay_;
else
    ETo_col = tablaETo.ET0_mmDay_;
end

% Match fechas
[lia, loc] = ismember(fechasETa, tablaETo.Date);
ETo_asignado = NaN(size(fechasETa));
ETo_asignado(lia) = ETo_col(loc(lia));

% Ordenar
[fechasETa, idx] = sort(fechasETa);
files = files(idx);
ETo_asignado = ETo_asignado(idx);

DOY = day(fechasETa,'dayofyear');
Distance = [NaN; days(diff(fechasETa))];

resultado = table({files.name}', fechasETa, DOY, Distance, ETo_asignado, ...
    'VariableNames', {'Nombre','Fecha','DOY','DeltaDias','ETo'});

%% ===================== 3. READ ETa IMAGES ===================
E = [];

for i = 1:n
    filePath = fullfile(files(i).folder, files(i).name);
    A = double(geotiffread(filePath));
    A(A == 0) = NaN;
    E = cat(3, E, A);

    fprintf('Leyendo imágenes: %.1f %%\n', i*100/n);
end

%% ===================== 4. COMPUTE Kc ========================
[h,k,m] = size(E);

ET0 = zeros(h,k,m);

for i = 1:m
    ET0(:,:,i) = resultado.ETo(i);
end

Kc = E ./ ET0;

%% ===================== 5. INTERPOLATE Kc ====================
distanceETa = resultado.DeltaDias(2:end);
kc_interp = [];

for j = 1:m-1

    dist = distanceETa(j);
    f1 = Kc(:,:,j);
    f2 = Kc(:,:,j+1);

    stepKc = (f2 - f1) / dist;

    for z = 0:dist-1
        temp = f1 + z * stepKc;

        % límites físicos
        mask = ~isnan(temp);
        temp(mask) = max(min(temp(mask),Kcmax),Kcmin);

        kc_interp = cat(3, kc_interp, temp);
    end

    fprintf('Interpolación Kc: %.1f %%\n', j*100/(m-1));
end

Kc_daily = kc_interp;

%% ===================== 6. DAILY ETo =========================
idx = tablaETo.Date >= FechaInicioETa;
tablaETo_filtrada = tablaETo(idx,:);

if RD == 1
    ETkc = tablaETo_filtrada.ET0R_mmDay_;
else
    ETkc = tablaETo_filtrada.ET0_mmDay_;
end

dtotal = length(ETkc);
ET0_daily = zeros(h,k,dtotal);

for i = 1:dtotal
    ET0_daily(:,:,i) = ETkc(i);
end

%% ===================== 7. DAILY ETa =========================
ETa_daily = Kc_daily(:,:,1:dtotal) .* ET0_daily;

%% ===================== 8. TIME SERIES =======================
MediaE = zeros(dtotal,1);
DesviaE = zeros(dtotal,1);

for i = 1:dtotal
    temp = ETa_daily(:,:,i);
    temp(temp == 0) = NaN;
    MediaE(i) = prctile(temp(~isnan(temp)),75);
    DesviaE(i) = nanstd(temp(:));
end

t = FechaInicioETa + days(0:dtotal-1);

figure;
plot(t, MediaE,'b','LineWidth',2); hold on;
plot(t, MediaE+DesviaE,'--r');
plot(t, MediaE-DesviaE,'--r');
grid on;
title('ETa promedio');
ylabel('mm/día');

%% ===================== 9. CUMULATIVE ETa ====================
idx_ini = days(FechaSiembra - FechaInicioETa) + 1;
idx_fin = idx_ini + days(FechaFinal - FechaSiembra);

cumETa = nansum(ETa_daily(:,:,idx_ini:idx_fin),3);
cumETa(cumETa == 0) = NaN;

figure;
imagesc(cumETa); colorbar;
title('ETa acumulada');

%% ===================== 10. EXPORT ===========================
[~, R] = readgeoraster(fullfile(files(1).folder, files(1).name));

epsgCode = 32613;

if RD == 1
    tipo = 'rean';
else
    tipo = 'med';
end

filename = sprintf('cumETa_%s.tif', tipo);

geotiffwrite(filename, cumETa, R, 'CoordRefSysCode', epsgCode);

fprintf('Exportado: %s\n', filename);