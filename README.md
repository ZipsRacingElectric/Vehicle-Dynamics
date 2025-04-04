![ZRE Logo](./images/Logo_with_zippy_subtext_white.png "Zips Racing Electric Logo")

## Getting Started
This repository contains a collection of scripts and models used for vehicle dynamics analysis. Refer to the documentation for each folder. Some scripts may make use of a vehicle parameter model that reflects the latest vehicle design.

### Using tire data with this repository
Commiting tire data to this repository may violate confidentiality agreements. On your cloned repository, use the following directory only to store this information:

```
/MATLAB/tire_data/
```

All contents of this directory get ignored when your commit your changes. For many of the vehicle dynamics scripts, you will need to grab the tire_data folder from the Zips Racing OneDrive and place the contents in this directory.

### Using common vehicle data across multiple scripts
The vehicle model script at /MATLAB/vehicle_data/vehicle.m allows you to use consistent vehicle parameters across all of the Simulink models by loading data from a spreadsheet containing vehicle parameters.

To write a MATLAB script that uses the vehicle parameters, first load the file in your MATLAB .m script:

```
% create vehicle object from vehicle data
zr25 = vehicle('../vehicle_data/zr25_data.xlsx');
```
This will create a vehicle object 'zr25' in the MATLAB workspace, where parameters can be accessed like this:

```
ratio = zr25.wheelbase / ((zr25.track_width_front + zr25.track_width_front) / 2);
```

If you're using Simulink, you need to generate Simulink parameters from the vehicle object:

```
zr25.create_simulink_parameters();
```

Once the Simulink parameters are loaded into the workspace, you can then load the Simulink model:

```
% open simulink model
open_system('my_sim');
```

Within Simulink, you can change block parameters to reference the Simulink parameters.

### Storing Simulation Data
Generally, it is not recommended to store simulation output data (.mat, .mldatx), especially from batch simulation runs. These files are large and can cause the repository to exceed storage limits. Instead, save them in a separate folder and update the .gitignore file to ignore the folder contents.

### Prerequisite Tools

- [MATLAB](https://www.mathworks.com/products/matlab.html)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [README-Template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)