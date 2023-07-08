function [air_density] = fun_air(temp, pressure)

% pressure in Pa, temp in K

density = (6.022e23.*pressure)./(8.31.*temp); ... air density, molec/cubic m
density = density.*1e-6; ... air density, molec/cubic cm    
air_density = density;