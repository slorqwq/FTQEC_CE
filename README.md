# Fault-Tolerant Quantum Error Correction for Constant-Excitation Stabilizer Codes under Coherent Noise

##  Development environment
- IDE: Dev-C++ 5.11
- Compiler: TDM-GCC 4.9.2 64-bit Release

## How to Run and Modify 1213_FT_Coherent_v2.1.cpp

### To modify the parameters related to errors, edit the following lines:
-   Located at lines 211–214:
    -   double alpha = 1; // CNOT error
    -   double beta = 1;  // measurement error
    -   double gamma = 0.01; // idle error
    -   double delta = 1;  // ancilla state preparation error

### Enabling CC Errors 
- By default, the simulation runs without CC errors.
- To simulate correlated coherent errors, make the following changes in the code:
    - Uncomment lines 233–234
    - Comment out lines 235–236
