import h5py
import numpy as np
from scipy import interpolate

from gas_corrections_cpp.bin import gas_transmittance

class Gas_Correction_Manager:
    def __init__(self):
        self.l1_data = None
        self.l1_filename = None
        self.gas_transmittance_table = None
        self.gas_transmittance_table_filename = None

    def read_PACE_data(self, l1_filename, **kwargs):
        if self.l1_data is not None and l1_filename == self.l1_filename:
            # Avoid redundant reading from a file we already read data from
            return self.l1_data
        else:
            l1_data = gas_transmittance.L1_Data()
            start_line = kwargs.get("start_line", 0)
            end_line = kwargs.get("end_line", None)
            start_pixel = kwargs.get("start_pixel", 0)
            end_pixel = kwargs.get("end_pixel", None)

            with h5py.File(l1_filename, 'r') as f:
                solar_zenith = 0.01*np.array(f['/geolocation_data/solar_zenith'][start_line:end_line, start_pixel:end_pixel])
                sensor_zenith = 0.01*np.array(f['/geolocation_data/sensor_zenith'][start_line:end_line, start_pixel:end_pixel])
                l1_data.cos_solar_zenith = np.cos(np.deg2rad(solar_zenith))
                l1_data.cos_sensor_zenith = np.cos(np.deg2rad(sensor_zenith))
                blue_wavelengths = np.array(f['/sensor_band_parameters/blue_wavelength'][1:])
                red_wavelengths = np.array(f['/sensor_band_parameters/red_wavelength'][3:])
                sensor_wavelengths = np.zeros(len(blue_wavelengths) + len(red_wavelengths))
                sensor_wavelengths[0:len(blue_wavelengths)] = blue_wavelengths
                sensor_wavelengths[len(blue_wavelengths):] = red_wavelengths

                blue_rhot = np.array(f['/observation_data/rhot_blue'][1:, start_line:end_line, start_pixel:end_pixel])
                red_rhot = np.array(f['/observation_data/rhot_red'][3:, start_line:end_line, start_pixel:end_pixel])
                assert(blue_rhot.shape[1] == red_rhot.shape[1])
                assert(blue_rhot.shape[2] == red_rhot.shape[2])

                rhot = np.zeros((blue_rhot.shape[0] + red_rhot.shape[0], blue_rhot.shape[1], blue_rhot.shape[2]))
                rhot[0:len(blue_wavelengths), :, :] = blue_rhot
                rhot[len(blue_wavelengths):, :, :] = red_rhot
                rhot = np.rollaxis(rhot, 0, 3)
                l1_data.reflectance = rhot/np.pi*l1_data.cos_sensor_zenith[:, :, None]

                assert(len(sensor_wavelengths) == l1_data.reflectance.shape[2])

                l1_data.wavelengths = sensor_wavelengths
                l1_data.num_pixels = len(l1_data.cos_solar_zenith.flatten())
                l1_data.num_wavelengths = len(l1_data.wavelengths)

            self.l1_filename = l1_filename
            self.l1_data = l1_data
            return self.l1_data


    def read_gas_transmittance_table(self, gas_transmittance_table_filename):
        if self.gas_transmittance_table is not None and gas_transmittance_table_filename == self.gas_transmittance_table_filename:
            # Avoid redundant reading from a file we already read data from
            return self.gas_transmittance_table
        else:
            gas_transmittance_table = gas_transmittance.Gas_Transmittance_Lookup_Table()

            netcdf_variables = {
                    'air_mass_factor_mixed' : "air_mass_factor_mixed_gases",
                    'air_mass_factor_wv' : "air_mass_factor_water_vapor",
                    'carbon_dioxide_transmittance' : "co2_transmittance",
                    'carbon_monoxide_transmittance' : "co_transmittance",
                    'methane_transmittance' : "ch4_transmittance",
                    'nitrous_oxide_transmittance' : "n2o_transmittance",
                    'oxygen_transmittance' : "o2_transmittance",
                    'water_vapor' : "water_vapor_concentration",
                    'water_vapor_transmittance' : "h2o_transmittance",
                    'wavelength' : "wavelengths",
                }
            
            netcdf_dimensions = {
                'n_air_mass_factor' : "num_amf_grid_points",
                'n_water_vapor' : "num_water_vapor_concentrations",
                'nmodels' : 'num_models',
                'nwavelengths' : 'num_wavelengths',
            }

            model_map = {
                "Tropical": 0,
                "MidLatSummer" : 1,
                "MidLatWinter" : 2,
                "SubarcticSummer" : 3,
                "SubarcticWinter" : 4,
                "USstandard62" : 5
            }
                
            with h5py.File(gas_transmittance_table_filename, 'r') as f:
                for netcdf_node, var_name in netcdf_variables.items():
                    if netcdf_node in f:
                        var_value = np.array(f[netcdf_node])
                        setattr(gas_transmittance_table, var_name, var_value)

                for netcdf_node, var_name in netcdf_dimensions.items():
                    if netcdf_node in f:
                        var_value = np.array(f[netcdf_node])
                        setattr(gas_transmittance_table, var_name, len(var_value))

                setattr(gas_transmittance_table, 'model', model_map["USstandard62"])

            self.gas_transmittance_table_filename = gas_transmittance_table_filename
            self.gas_transmittance_table = gas_transmittance_table
            return self.gas_transmittance_table


    def ozone_transmittance(self, ancillary_data, use_gas_transmittance_lookup_table=False):
        return gas_transmittance.ozone_transmittance(self.l1_data, ancillary_data, use_gas_transmittance_lookup_table)


    def co2_transmittance(self, use_gas_transmittance_lookup_table=True):
        return gas_transmittance.co2_transmittance(self.l1_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table)
    

    def co_transmittance(self, use_gas_transmittance_lookup_table=True):
        return gas_transmittance.co_transmittance(self.l1_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table)
    

    def ch4_transmittance(self, use_gas_transmittance_lookup_table=False):
        return gas_transmittance.ch4_transmittance(self.l1_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table)
    

    def n2o_transmittance(self, use_gas_transmittance_lookup_table=True):
        return gas_transmittance.n2o_transmittance(self.l1_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table)
    

    def no2_transmittance(self, ancillary_data, use_gas_transmittance_lookup_table=False):
        return gas_transmittance.no2_transmittance(self.l1_data, ancillary_data, use_gas_transmittance_lookup_table)
    

    def o2_transmittance(self, use_gas_transmittance_lookup_table=True, oxygen_A_band_option=gas_transmittance.Oxygen_A_Band_Option.TRANSMITTANCE_TABLE):
        f = interpolate.interp1d(self.gas_transmittance_table.wavelengths, self.gas_transmittance_table.o2_transmittance, axis = 0)
        o2_transmittance_sensor_wavelengths = f(self.l1_data.wavelengths)

        self.gas_transmittance_table.o2_transmittance = o2_transmittance_sensor_wavelengths

        return gas_transmittance.o2_transmittance(self.l1_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table, oxygen_A_band_option)
    

    def h2o_transmittance(self, use_gas_transmittance_lookup_table=True):
        ancillary_data = gas_transmittance.Ancillary_Data()

        ancillary_data.precipitable_water = np.zeros(self.l1_data.cos_solar_zenith.size)
        ancillary_data.water_vapor_bands = np.array([782, 817, 857], dtype=np.float64)
        ancillary_data.num_water_vapor_bands = ancillary_data.water_vapor_bands.size

        f = interpolate.interp1d(self.gas_transmittance_table.wavelengths, self.gas_transmittance_table.h2o_transmittance, axis = 1)
        h2o_transmittance_at_sensor_wavelengths = f(self.l1_data.wavelengths)

        self.gas_transmittance_table.num_wavelengths = len(self.l1_data.wavelengths)
        self.gas_transmittance_table.h2o_transmittance = h2o_transmittance_at_sensor_wavelengths

        # This option is not needed, will remove it from the API soon
        use_gas_transmittance_table = True

        return gas_transmittance.h2o_transmittance(self.l1_data, ancillary_data, self.gas_transmittance_table, use_gas_transmittance_lookup_table, use_gas_transmittance_table)