# MILP Model
An OFDM-based MILP formulation to solve the routing, modulation and spectrum assignment for various VPNs inside the same physical transparent network.

An example of data file for a simple six-node network is contained inside the /data directory.

## How to use

1. You need a licensed version of [AMPL](https://ampl.com/) with CPLEX solver
2. Open the ampl.exe windows executable
3. Model the formulation "model "*/milp-vpn-rmsa/vpn_rmsa.mod;"
4. Load the data parameters "data "*/milp-vpn-rmsa/data/six-node.dat; (You can you your own data file)
5. Select the cplex solver "options solver cplex;"
6. Solve with "solve;"
