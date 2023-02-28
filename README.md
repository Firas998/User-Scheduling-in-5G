# IoT and 5G Wireless Networks

This project contains the implementation of algorithms for deriving optimal routes between sensors in an IoT network, as well as designing optimal packet schedulers for 5G wireless networks. The algorithms implemented in this project include:

- Preprocessing
- Greedy implementation
- LPsolver
- Two different implementations of dynamic programming
- An implementation of branch and bound algorithm

## Introduction

Internet of Things (IoT) is characterised by the deployment of a large amount of connected objects. These objects are autonomous devices that use sensors to monitor their environment and communicate through a wireless radio channel. When these sensors are not powerful, they cannot communicate over a long distance, so that additional intermediate relay nodes may be required to connect them together. Sensors can themselves route the messages from other sensors or relays. In this project, we derive the optimal routes between sensors by using the notion of Steiner tree in graph theory.

In 5G, an antenna transmits data packets to smartphones (or users) through a wireless medium, which is divided into a set of frequency channels. Figure 1 is an example of simultaneous transmission towards three users using three channels. The higher the power dedicated to a user, the higher the data rate it can experience. The exact dependence between power and data rate is however user and channel specific. With the same transmit power, a user close to the antenna will enjoy for example a higher data rate than a user far away. A wireless packet scheduler is thus responsible to allocate channels to users and to divide the total power budget of the antenna among the available channels. The goal of this project is to design optimal packet schedulers in this context.


## Usage

To use the algorithms implemented in this project, follow these steps:

1. Clone this repository to your local machine.
2. Install requirements : pip install -r requirements.txt
2. Compile the code using your preferred compiler.
3. Run the compiled executable with the appropriate input files.

## Contributors

This project was developed by Firas Ben Jedidia and Mattia Mauro.

