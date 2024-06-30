![ZRE Logo](./images/Logo_with_zippy_subtext_white.png "Zips Racing Electric Logo")
# Table of Contents

- [Table of Contents](#table-of-contents)
- [Project Name](#project-name)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Hardware](#hardware)
    - [OneDrive Documentation](#onedrive-documentation)
    - [Installing](#installing)
  - [Running unit tests](#running-unit-tests)
    - [Break down into end to end tests](#break-down-into-end-to-end-tests)
    - [And coding style tests](#and-coding-style-tests)
  - [Deployment](#deployment)
  - [Built With](#built-with)
  - [Versioning](#versioning)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)


# Project

1. This is a template repository to assist you with creating new embedded projects for Zips Racing Electric. Please fill out this README.md doc. Replace this text with a project description of your own.
2. Make sure to update the repository description.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

1. If you are new to the team, read the onboarding document on the OneDrive: **/Engineering/Subsystems/Software/Oboarding/Software New Member Onboarding.docx**
   
2. Please also review our embedded software standards: **/Engineering/Design Standards/Embedded Software Standards.docx**
   
3. If you haven't already, please install the pre-requesite software below.

### Prerequisites

- [MPLAB X IDE](https://www.microchip.com/en-us/tools-resources/develop/mplab-x-ide#tabs)

- [MPLAB XC Compiler](https://www.microchip.com/en-us/tools-resources/develop/mplab-xc-compilers/xc16)

### Hardware
To compile your software and deploy it on our hardware, you will need the following hardware:

- DT0003-AA AMS PCB
- Compatible AMS test harness
- [MPLAB ICD 4 / ICD 5 In-Circuit Debugger](https://www.microchip.com/en-us/development-tool/dv164045#)
- [Tag-Connect TC2030-MCP 6-Pin Cable](https://www.tag-connect.com/product/tc2030-mcp-6-pin-cable-with-rj12-modular-plug-for-microchip-icd)
- 12V Power supply

The project is setup to require external power to the board in order to program the microcontroller.

### OneDrive Documentation
- Additional hardware and system documentation for the AMS can be found on the OneDrive at **/Engineering/Subsystems/HV Systems/Documentation/AMS and High Voltage PCB**.
  
- Electrical documentation for the accumulator can be found at **/Engineering/Subsystems/HV Systems/Documentation/Accumulator**.

- Documentation for CAN bus system can be found at **/Engineering/Subsystems/GLV Systems/Documentation/ELA0002 Wiring Harness**.

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running unit tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Versioning

At the end of each vehicle season, the repository is copied and renamed for the new vehicle. The old repository is then archived.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [README-Template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)

