function [includePromoters, includeRBS, includeORIs] = get_free_params_for_library(libName)
%GET_FREE_PARAMS_FOR_LIBRARY Devuelve las biopartes libres a estimar según la librería.
%
%   [includePromoters, includeRBS, includeORIs] = get_free_params_for_library(libName)
%
%   libName (string):
%       'L24'  : Biblioteca base (L24)
%       'L30'  : Biblioteca extendida (L30 = L24 + RBS B0034)
%       'L6'   : Sublibrería centrada en B0034 (3 promotores x 2 ORIs)
%       'L4'   : Sublibrería J23100 con 4 RBS (sin el TU B0034)
%       'L5'   : Sublibrería J23100 con 5 RBS (incluye el TU con B0034)
%
%   Salidas:
%       includePromoters : cell array con nombres de promotores cuyas omega_A se estiman
%       includeRBS       : cell array con nombres de RBS cuyos parámetros de traducción se estiman
%       includeORIs      : cell array con nombres de ORIs cuyos copy numbers N_A se estiman
%
%   NOTA:
%   - Los parámetros que NO aparecen en estos arrays se asumen fijados (conocidos a priori
%     o heredados de otra librería, p.ej. N_pSC101 = 5 en L24/L30).
%   - Esto evita incluir columnas de parámetros fijados en el Jacobiano/FIM y garantiza
%     que el rango refleje el problema real de estimación para cada librería.

switch upper(libName)

    case 'L24'
        % ------------------------------------------------------------
        % L24: librería base
        % Estimaste:
        %   - 3 promotores: J23106, J23102, J23101  (omega_A)
        %   - 4 RBS: B0030, B0032, J61100, J61101   (K0, sigma0)
        %   - N_pGreen
        %   - N_pSC101 se fijó a 5 (Thompson 2018) -> NO se incluye
        % ------------------------------------------------------------
        includePromoters = {'J23106','J23102','J23101'};
        includeRBS       = {'B0030','B0032','J61100','J61101'};
        includeORIs      = {'pGreen'};   % sólo N_pGreen libre

    case 'L30'
        % ------------------------------------------------------------
        % L30: librería extendida (añade RBS B0034)
        % Estimaste:
        %   - mismos 3 promotores: J23106, J23102, J23101
        %   - 5 RBS: B0030, B0032, J61100, J61101, B0034
        %   - N_pGreen
        %   - N_pSC101 fijado a 5 -> NO se incluye
        %   - J23100 NO entra aún como promotor libre (se usa luego en L4/L5)
        % ------------------------------------------------------------
        includePromoters = {'J23106','J23102','J23101'};
        includeRBS       = {'B0030','B0032','J61100','J61101','B0034'};
        includeORIs      = {'pGreen'};

    case 'L6'
        % ------------------------------------------------------------
        % L6: sublibrería centrada en B0034 (3 promotores x 2 ORIs, todos con B0034)
        % En la práctica:
        %   - Usas N y omega heredados de L30 (no los reestimas aquí)
        %   - Sólo refinas el parámetro efectivo del RBS B0034
        % Por tanto:
        %   - includePromoters = {}  (todas las omega fijadas)
        %   - includeORIs      = {}  (todos los N fijados)
        %   - includeRBS       = {'B0034'}  (único parámetro nuevo libre)
        % ------------------------------------------------------------
        includePromoters = {};
        includeRBS       = {'B0034'};
        includeORIs      = {};

    case 'L4'
        % ------------------------------------------------------------
        % L4: sublibrería J23100 (pGreen) con 4 RBS (sin B0034)
        %     -> L5 menos la TU con B0034.
        % Estrategia:
        %   - RBS y N se heredan de L30 (fijos)
        %   - Sólo se estima la omega de J23100
        % Por tanto:
        %   - includePromoters = {'J23100'}
        %   - includeRBS       = {}
        %   - includeORIs      = {}
        % ------------------------------------------------------------
        includePromoters = {'J23100'};
        includeRBS       = {};
        includeORIs      = {};

    case 'L5'
        % ------------------------------------------------------------
        % L5: sublibrería J23100 (pGreen) con 5 RBS (incluyendo B0034)
        % Mismo criterio que L4:
        %   - RBS y N se heredan de L30
        %   - Sólo se estima la omega de J23100
        % ------------------------------------------------------------
        includePromoters = {'J23100'};
        includeRBS       = {};
        includeORIs      = {};

    otherwise
        error('get_free_params_for_library: librería "%s" no reconocida.', libName);
end

end
