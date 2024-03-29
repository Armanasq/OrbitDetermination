# Orbital Mechanics Simulation README

## Introduction

This repository contains MATLAB code designed to simulate the orbital motion of a satellite around Earth. The simulation employs Keplerian orbital elements to calculate the satellite's position and velocity vectors over time. The primary objective of this code is to provide a visual representation and analysis of the satellite's trajectory in Earth-Centered Inertial (ECI) coordinates.

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Code Overview](#code-overview)
    - [Orbital Elements Calculation](#orbital-elements-calculation)
    - [Orbit Visualization](#orbit-visualization)
    - [Coordinate Transformations](#coordinate-transformations)
4. [Formulas Used](#formulas-used)
5. [Additional Information](#additional-information)
6. [Contact](#contact)

## Installation

1. Download and install MATLAB on your system.
2. Clone this repository to your local machine using the following command:

   ```bash
   git clone https://github.com/your-username/orbital-mechanics-simulation.git
   ```

## Usage

1. Open MATLAB and navigate to the cloned repository.
2. Run the `main.m` script to execute the simulation.
3. Examine the generated plots to visualize the satellite's orbit.

## Code Overview

### Orbital Elements Calculation

The initial conditions of the satellite, such as position and velocity, are defined in Cartesian coordinates. These values are then converted into classical orbital elements using the `rv_from_r0v0` and `coe_from_rv` functions.

### Orbit Visualization

The satellite's orbit is visualized using MATLAB's plotting capabilities. The script generates a 3D plot of the satellite's trajectory in ECI coordinates, overlaid on a sphere representing Earth.

### Coordinate Transformations

Coordinate transformations are applied to convert the satellite's position from perifocal coordinates to ECI coordinates. The matrices \( R_{z_{\Omega}} \), \( R_{x_i} \), and \( R_{z_{\omega}} \) are used for these transformations.

## Formulas Used

1. **Orbital Elements Calculation:**
    - Semi-Major Axis: ![Semi-Major Axis Formula](https://latex.codecogs.com/png.latex?\bg_white%20\dpi{250}a=\frac{h^2}{\mu(1-e^2)})

    - Eccentric Anomaly: ![Eccentric Anomaly Formula](https://latex.codecogs.com/png.latex?\bg_white%20\dpi{250}\tan\left(\frac{\theta}{2}\right)%20\sqrt{\frac{1-%20e}{1+e}}\tan\left(\frac{E}{2}\right))
    

2. **Eccentric Anomaly Iteration:**
    - Newton's Method: ![Newton's Method Formula](https://latex.codecogs.com/png.latex?\bg_white%20\dpi{250}E=E-\frac{E-e\sin(E)-M}{1-e\cos(E)})

3. **Orbit Positions Calculation:**
    - Perifocal Coordinates: ![Perifocal Coordinates Formula](https://latex.codecogs.com/png.latex?\bg_white%20\dpi{250}p=a(\cos(E)-e),\quad(q=a\sqrt{1-e^2}\sin(E)))

4. **Coordinate Transformations:**
    - ![Coordinate Transformations Formula](https://latex.codecogs.com/png.latex?\bg_white%20\dpi{250}r_{ECI}=R_{z_{\Omega}}^{-1}\cdot%20R_{x_i}^{-1}\cdot%20R_{z_{\omega}}^{-1}\cdot%20r_{pq})

## Additional Information

The simulation also provides information such as the apogee and perigee of the satellite's orbit, as well as its final position and velocity vectors. Eccentricity, inclination, and other orbital elements are calculated and displayed.

The visual representation includes an illustration of the Earth's surface, the satellite's orbit trajectory, and key orbital parameters.

## Contact

For any questions, feedback, or suggestions regarding this code, please contact:

**Author:** A. Asgharpoor  
**Email:** A.Asgharpoor@ut.ac.ir
