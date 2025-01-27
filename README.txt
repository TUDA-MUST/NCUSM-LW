# Code for Non-contact Ultrasonic Stress Measurement using Lamb Waves (Short: NCUSM-LW)

In this github repository you will find all the code necessary to evaluate the data and plot the figures of the manuscript entitled `Non-contact Ultrasonic Stress Measurement using Lamb Waves`. It is submitted to Journal of Nondestructive Testing and Evaluation. It uses Matlab 2023b (Mathworks) or higher.

Find the paper here: [TODO]
Find the data here: https://tudatalib.ulb.tu-darmstadt.de/handle/tudatalib/4444 DOI: https://doi.org/10.48328/tudatalib-1654

## Table of Contents
- [Abstract](#abstract)
- [Getting Started](#getting-started)
- [Features and File Descriptions](#features-and-file-descriptions)
- [License](#license)

## Abstract

Stress measurement is essential in many applications, such as aerospace or construction, and ultrasonic stress measurement systems are removable and non-destructive. Using air-coupled ultrasound has the advantage of being non-contact. 
In other works, an air-coupled ultrasonic phased array is used to adjust the coupling angle and then measure the stress using the conventional transit-time method. In this work, we investigate using the coupling angle of air-coupled Lamb waves directly to measure normal stress in the specimen.  
The coupling angle is dependent on the phase velocity, which, in turn, changes with the stress. We model that effect using numerical simulations with the semi-analytical finite element method. Ultrasonic measurements are conducted on foam-filled sandwich panels with two 0.5 mm steel face sheets during full-scale bending tests according to EN 14509:2013. The ultrasonic stress measurement setup consists of an air-coupled phased array for transmission, a MEMS-microphone array for reception, and a laser Doppler vibrometer for reference. We measure the stress via the coupling angle using either the transmit or the receive phased array, and for reference, we measure the stress via the transit-time using the group velocity. 
The measurement method with the coupling angle method works, both with the transmit or the receive array, with a repeatability of 5.3 MPa and 4.1 MPa, respectively. The transit-time measurement performs better than the coupling angle method with a repeatability of 2.1 MPa since time measurement is more accurate. However, the coupling angle methods measures phase velocity instead of group velocity. Therefore, both methods can be advantageously combined.

   
## Getting Started

1. Clone the repo.
2. Replace the empty folder `data` and its subfolders with the datafolder from here [TODO] (unzip).
3. Make sure Matlab 2023b or higher is installed.
4. Run whatever file you would like.

## Features and File Descriptions

 - The file `a_evaluate_data_of_hydraulic_press.m` will convert the data recorded by the hydraulic press into a more manageable format and saves the new format. The new converted file is already contained in the data, so there is no need to run this file unless there are changes to it.
 - The file `a_plot_data_of_hydraulic_press.m` will plot Fig. 5 of the manuscript and more figures which are not in the manuscript.
 - The file `b_evaluate_data_of_ultrasonic_system.m` reads all the raw ultrasonic data and calculates coupling angles, transit-times and more and saves them in a new file. The new converted file is also already contained in the data, so there is no need to run this file fully unless there are changes to it. However during runtime it plots images of each of the raw ultrasonic data, which is used as Fig. 8, Fig. 9 and Fig. 10 in the manuscript.
 - The file `b_plot_raw_data_of_ultrasonic_system.m` plots coupling angles and transit-times in a format that was not included in the manuscript.
 - The file `c_combine_all_measurements.m` uses the data from the simulation, the hydraulic press and the ultrasonic system and plots Fig. 11, Fig. 12 and Fig. 13 and outputs the values for Table 1 in the manuscript. 
 - In the folder `simulation` there are three scripts to plot Fig. 1, Fig. 2 and Fig. 3 in the manuscript.
 - In the folder `functions` there are some subfunctions that are required by some other script.

## License

This project is licensed under the MIT License. See the LICENSE file for details.