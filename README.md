![ZRE Logo](./images/Logo_with_zippy_subtext_white.png "Zips Racing Electric Logo")

## Getting Started
This repository contains a collection of scripts and models used for vehidle dyanmics analysis. Refer to the documentation for each folder. Some scripts may make use of a vehicle parameter model which reflects the latest vehicle design.

### Using tire data with this repository
Commiting tire data to this repository may violate confidentiality agreements. On your cloned repository, use the following directory only to store this information:

/MATLAB/tire_model/

All contents of this directory get ignored when your commit your changes. For some of the MATLAB scripts you may need to grab the required tire model from the OneDrive and place it in this directory.

### Using vehicle data across multiple scripts
The vehicle model scripts under /MATLAB/vehicle_data/ allow you to use consistent vehicle parameters across all of the Simulink models.

To use the vehicle model, first load the file from your MATLAB .m script:

```
% load vehicle simulink parameters
run("../vehicle_data/zr25.m");
```

In Simulink, block parameters can be replaced with the variable names defined in the vehicle model. After you load the vehicle model, you can call your script to run the Simulink model:

```
% open simulink model
open_system('my_sim');
```

### Prerequisite Tools

- [MATLAB](https://www.mathworks.com/products/matlab.html)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [README-Template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)

