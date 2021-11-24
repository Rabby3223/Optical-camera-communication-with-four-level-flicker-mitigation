# Four-level-flicker-mitigation-for-OCC
This project includes signal generation and decoding process for a four-level optical camera communicaiton (OCC) system, with a four-level flicker mitigation scheme published in OECC 2021 and captured data. 
## Background
Optical camera communication (OCC) is an optical wireless communication technique that uses a camera as the receiver. Light sources such as light-emitting diode (LED) luminaires, can be used as the transmitters. 

In most of the complementary metal-oxide-semiconductor (CMOS) sensors, the rolling shutter is used to enable pixels to expose row by row or column by column. Therefore, each image can capture multiple data and significantly boost the data rate when the LED is modulated much faster than the frame rate. In the OCC system, one important issue is visual flicker, which may increase the chance of migraines, headaches, and even repetitive behavior among persons with autism. Therefore, mitigation schemes must be employed to achieve flicker-mitigation transmission in the practical OCC system.
## Install
This project uses Matlab. 

Download the code and data.zip file.

Unzip the compressed data and move the folders to the directory of code file. 

## Usage
### Transmitter side
tx1_PWM_signal_generation_proposed.m can generate PAM4 singals with proposed four level flicker mitigation scheme. 

tx1_PWN_signal_generation_manchester.m can generate PAM4 signals with manchester-like flicker mitigation scheme. 

Data packetization, flicker mitigation, PWM modulation, header insertion are included in this file. 

The generated PWM signals can be input to AWG (Channel2, Siglent SDG5162). In this OCC system, a MOSFET based driving circuit is used to drive LED luminaire.
### Receiver side
rx1_Decoding.m can decode the received OCC signals and estimate the BER performance. 

Scaling, synchronization, equalization, packet reconstruction, flicker mitigation decoding, and BER calculation are included. 

The attached data were output by a self-developed App. The smartphone (OnePlus 5T) captured the light reflected by a poster and generated an image in YUV format (960x1280). The Y components were used for averaging row by row. For each image frame, the generated 960 elements were coded by base64 and saved as a csv file in smartphone. 

<img height="200" width="200" src="https://github.com/Rabby3223/Four-level-flicker-mitigation-for-OCC/blob/main/imgForReadme/2.jpg">

The floder name was set as pam4_symbolRate_redundantBitLength_blockLength or pam4_symbolRate_manchester. Properly set the parameters (symbol_rate, block_bit_length, flicker_mitigation) when running the decoding file. 

<img height="576" width="364" src="https://github.com/Rabby3223/Four-level-flicker-mitigation-for-OCC/blob/main/imgForReadme/folders.jpg)">
## License
This project is licensed under the MIT License.

base64decode.m is from Plotly Graphing Library for MATLAB.

Cite as: L. Liu and L. Chen, "Four-level Flicker-mitigation Coding Scheme in the Non-line-of-sight Optical Camera Communication System," in 26th Optoelectronics and Communications Conference (OECC), paper M4B.5, Hong Kong, China, Jul. 2021, pp. 1-3.
