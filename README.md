
# Atmospheric Gas Correction Library

This library computes atmospheric gas transmittance values for each pixel and each wavelength in an L1 data file that is provided as input. The gas correction algorithms implemented by this library are designed to reproduce the results from NASA's [OCSSW](https://oceandata.sci.gsfc.nasa.gov/ocssw) as a standalone, C++/OpenMP-accelerated library. A Python interface is provided (via pybind11), so it can be called directly from Python workflows such as [OC-SMART]([https://oceancolor.gsfc.nasa.gov/](https://github.com/Light-and-Life-Lab/OC-SMART)).

Transmittance values are computed for eight molecular species:

-   O3
-   NO2
-   CO
-   CO2
-   CH4
-   N2O
-   O2
-   Water Vapor (H2O)

## Installation
### Linux and MacOS

Pre-built wheels for **Linux** and **macOS** are published on the [Releases page](https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library/releases) for each tagged version. 

#### 1. Find the right wheel

Each release lists several `.whl` files, one per combination of Python version and platform. Each wheel has the format {distribution}-{version}-{python tag}-{abi tag}-{platform tag}.whl. Identify the wheel that matches your machine, for example if you are running with Python 3.10 on MacOS 11, the wheel name will look like:

lll_gas_corr_lib-0.1.4-cp310-cp310-macosx_11_0_arm64.whl

Similarly, Python 3.10 on Linux would look like:

lll_gas_corr_lib-0.1.4-cp310-cp310-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl

If you are unsure of your Python version, run:

```bash
python3 --version

```

The `cp3XX` in the filename should match your major/minor version (e.g. Python 3.11 -> `cp311`).

#### 2. Install directly from the release URL

Right-click the appropriate file on the [Releases page](https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library/releases) to copy its link, then run the following from the terminal:

```bash
pip install <paste-the-wheel-url-here>

```

For example:

```bash
pip install https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library/releases/download/v0.1.7/lll_gas_corr_lib-0.1.7-cp310-cp310-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl

```

This installs `gas_corrections_lib` along with its core dependencies (`numpy`, `scipy`, `h5py`).

### Windows

Native Windows wheels are not built or supported. Windows users should install and run this library from within **[WSL (Windows Subsystem for Linux)](https://learn.microsoft.com/en-us/windows/wsl/install)**, then follow the Linux installation instructions above from inside your WSL environment.

## Requirements

-   Python >= 3.10
-   pip (conda not required, but is also supported)

## Usage

Every transmittance function takes an `L1_Data` object, which must be populated with sensor/solar zenith angle arrays (each element of these arrays corresponds to a pixel in the L1B image) as well as the number of pixels and number of wavelengths. In addition, each function requires either a transmittance lookup table (for CH4, CO, CO2, N2O, O2, H2O) or ancillary data (for Ozone, NO2), and each returns a `Gas_Transmittances` object with `sensor_zenith`, `solar_zenith`, and `total` transmittance matrices (each of shape `num_pixels x num_wavelengths`).

### Inputs

Each transmittance function in the library must be called with several arguments. These arguments are dataclasses which contain the data needed to perform the transmittance calculations. Every transmittance function requires an L1 Data dataclass which contains data from the L1 file (such as the sensor and solar zenith angles of each pixel and the sensor wavelengths). 

In addition to the L1 data, the Ozone and NO2 functions require data from ancillary files, such as the concentrations of these gases, and their absorption coefficients. The data from these files must be assigned to an Ancillary Data dataclass before it is passed into the respective transmittance function. 

The other six transmittance functions do not require ancillary data, but they do require a lookup table of transmittance values. Each sensor has an associated lookup table, and these files are included with the library. The Gas Correction Manager has a helper function called read gas transmittance table that takes the name of the lookup table file. This helper function reads the data from the provided file and returns it, which may then be passed into the appropriate transmittance function. 

When manually populating dataclasses, one must be careful to ensure that the dimensions of certain matrices match, and that the num pixels and num wavelengths members of the L1 Data dataclass match those of the matrices. If these assumptions are violated, then the transmittance calculations will not be correct. For this reason, these assumptions are enforced through assertions in the Python API, which will raise an error if any of the assumptions are violated

### Outputs

Upon completion, each transmittance function returns a `Gas_Transmittances` dataclass, which contains three numpy arrays. Each of these arrays has dimensions `np × nw`, where `np` is the number of pixels in the input L1 data (i.e. the `num_pixels` member of the `L1_Data` dataclass passed to the function). Similarly `nw` is the number of wavelengths (`num_wavelengths` member of the `L1_Data` dataclass passed to the function).

-   **`sensor_zenith`**: The _ij_<sup>th</sup> element of this matrix corresponds to the transmittance along the slant path defined by the sensor zenith angle for the _i_<sup>th</sup> pixel and _j_<sup>th</sup> wavelength.
-   **`solar_zenith`**: The _ij_<sup>th</sup> element of this matrix corresponds to the transmittance along the slant path defined by the solar zenith angle for the _i_<sup>th</sup> pixel and _j_<sup>th</sup> wavelength.
-   **`total`**: The _ij_<sup>th</sup> element of this matrix corresponds to the two-way transmittance along the slant paths defined by the solar and sensor zenith angles for the _i_<sup>th</sup> pixel and _j_<sup>th</sup> wavelength.

For CH4, CO, CO2, N2O, O2, and H2O, the transmittance functions use values from lookup tables that are packaged with the library. These lookup tables are exactly the same as those used by OCSSW to perform gas corrections. For Ozone and NO2, no transmittance lookup tables are available, so the following equations are used:

```
T_sen = exp( -τ / cos(θ_sen) )                                        (1)

T_sol = exp( -τ / cos(θ_sol) )                                        (2)

T_tot = T_sen * T_sol = exp( -τ * (1/cos(θ_sen) + 1/cos(θ_sol)) )     (3)
```

where `τ` is the optical depth, `θ_sen` is the sensor zenith angle, and `θ_sol` is the solar zenith angle.

Note that internally, the gas corrections library operates on a 2D matrix of dimensions `np × nw`, but it is often convenient to treat the L1 image data using two dimensions for the pixels and a third for the wavelengths. For example, suppose we wish to compute the transmittance values for a PACE image with 1710 × 1272 pixels at 197 wavelengths. Internally, the library will operate on a matrix with dimensions 2,175,120 × 197, and each of the resulting output transmittance matrices will have these dimensions. After receiving the output from the library, it may be useful to reshape this matrix to one with dimensions 1710 × 1272 × 197. This reshaping can be handled by standard numpy operations.

**Example: CO2 (Uses a transmittance lookup table)**

In this example it is assumed that the user has already read the solar zenith, sensor zenith, and sensor wavelength values from the L1B file and stored them in corresponding variables named `l1b_solz`, `l1b_senz`, and `sensor_wavelengths` respectively. This example then shows how a the CO2 transmittances may be computed, and in this case gas_transmittances.solar_zenith and gas_transmittances.sensor_zenith are both returned.

```python
import numpy as np
from gas_corrections_lib.src import gas_corrections

gas_transmittance_manager= gas_corrections.Gas_Correction_Manager()
gas_transmittance_table = manager.read_gas_transmittance_table("<path to sensor's transmittance table>")

def compute_co2_transmittance(l1b_solz, l1b_senz, sensor_wavelengths):
	print(’Compute CO2 transmittance ...’)
	
	l1_data = gas_corrections.L1_Data()
	l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
	l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
	l1_data.num_pixels = len(l1_data.cos_solar_zenith)
	l1_data.wavelengths = sensor_wavelengths
	l1_data.num_wavelengths = len(l1_data.wavelengths)
	
	gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)
	gas_transmittances = gas_corrections.co2_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)
	
	return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith

```

**Example: Ozone (Uses Ancillary Data instead of a lookup table)**

In this example it is assumed that the user has already read an array of ozone concentration values and the corresponding lat/lon grid on which they are defined from an ancillary data file (such as a MERRA2 Met file). In this example they are stored in variables `ozmap`, `oz_lat`, and `oz_lon`. In addition, an array (`ozone_absorption_cross_section`) of ozone absorption cross section values was created with cross section values at each wavelength at which transmittances are to be computed. 

The lat/lon grid on which the L1B data is defined has has been stored in arrays `l1b_lat` and `l1b_lon` (each element of these arrays has a lat or a lon value corresponding to each pixel of the L1B image). And finally just as in the CO2 example, the solar zenith and sensor zenith angles at each pixel in the L1B image are stored in `l1b_solz` and `l1b_senz` respectively. 

With all of these inputs available, they can be passed into a function that computes the ozone transmittances.

```python
def compute_ozone_transmittance(ozmap, oz_lat, oz_lon, ozone_absorption_cross_section, l1b_lat, l1b_lon, l1b_solz, l1b_senz):
	print(’Compute Ozone transmittance ...’)
	
	cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
	cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
	
	l1_data = gas_corrections.L1_Data()
	l1_data.cos_solar_zenith = cos_solar_zenith
	l1_data.cos_sensor_zenith = cos_sensor_zenith
	l1_data.num_pixels = len(cos_solar_zenith)
	l1_data.num_wavelengths = len(ozone_absorption_cross_section)
	
	# Interpolate ozone map to the L1B grid
	func = interpolate.RegularGridInterpolator( \
			(np.flip(oz_lat), oz_lon), \
			np.flip(ozmap, 0) \
		)
		
	ozone_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())
	self.ozone_concentration = ozone_concentration
	
	ancillary_data = gas_corrections.Ancillary_Data()
	ancillary_data.ozone_absorption_cross_section = ozone_absorption_cross_section
	ancillary_data.ozone_concentration = ozone_concentration
	
	gas_transmittances = gas_corrections.ozone_transmittance(l1_data=l1_data, ancillary_data=ancillary_data)
	
	return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith

```

The other six species follow the same pattern as the CO2 example above (`co_transmittance`, `ch4_transmittance`, `n2o_transmittance`, `o2_transmittance`, `h2o_transmittance`), and NO2 follows a similar pattern as Ozone.

> **Note:** `num_pixels` and `num_wavelengths` on `L1_Data` must match the dimensions of the arrays you assign to it — this is enforced with assertions, and mismatches will raise an error rather than silently producing incorrect results.

## Running the test suite

The test suite is not included in the published wheels. See the DEVELOPMENT.md file for instructions on how to install a version of the repo for development and testing.
