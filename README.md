# INO_CMT

Main header file:
Cosmic_trig.h

Main cpp file:
Cosmic_trig.cpp

Test Bench Codes:

Cosmic_trig_tb.cpp: This code just uses an initialized one-dimensional input data (can be used just to check if the code is running properly)

Cosmic_trig_tb_nofileread.cpp: This code just uses an initialized two-dimensional input data (can be used just to check if the code is running properly)

Cosmic_trig_tb_fileread.cpp: This code reads the data (two-dimensional) from an input file

GenerateHits.h and GenerateHits.cpp files are written to generate the hits randomly if the real data is not available

Input data (supplied by Yuvaraj): eve.bin

test_trial.cpp: Code which reads eve.bin file and writes the data into a convenient format in a text file, which is read by Cosmic_trig_tb_fileread.cpp test bench
