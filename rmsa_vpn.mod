# VPN RMSA complete model

# Defining Sets ---------------------------------------------------------------------------------------------------------

set N;                           # Set of nodes in the network      -	i                       / j / k / u / l / m / n
set E within (N cross N);        # Set of links in the network      -	(m,n)
set D within (N cross N);        # Set of demands in the network    -	(i,j)                   / (k,u)
set M;                           # Set of modulation formats        - 	(z)
set T;                           # Set of time periods              -     (t)
            
set DxE:= D cross E;                        # Set of D by E             -	( (i,j),(m,n) )               / (i,j,m,n)
set DxT:= D cross T;                        # Set of D by T             -	( (i,j),(t) )                 / (i,j,t)
set DxTxE:= D cross T cross E;              # Set of D by T by E        -	( (i,j),(t),(m,n) )           / (i,j,t,m,n)
set DxTxM:= D cross T cross M;              # Set of D by T by M        -	( (i,j),(t),(z) )             / (i,j,t,z)
set DxTxDxT:= D cross T cross D cross T;    # Set of D by T by D by T   -	( (i,j),(t),(k,l),(u) )       / (i,j,t,k,l,u)

# Defining Constants ----------------------------------------------------------------------------------------------------

param G:=1;         # Filter Guard Band between wavebands
param H:=1000;      # A large number H
param r:=12.5;      # Bandwidth per subcarrier (Ghz per subcarrier)

# Defining Input Parameters ---------------------------------------------------------------------------------------------

param v{(i,j,t) in DxT};    # Traffic demand from node i to node j at time t
param d{(m,n) in E};        # Length of each link from node m to node n
param d_{z in M};           # Maximum distance support by modulation format z
param ef{z in M};           # Efficiency of modulation format z

# Defining Formulation's variables --------------------------------------------------------------------------------------

var P{(i,j,t,m,n) in DxTxE} integer >=0;    # Traffic flow from demand (i,j,t)
var A{(i,j,t,m,n) in DxTxE} binary;         # Traffic flow indicator from demand (i,j,t)
var s{(i,j,t) in DxT} integer;              # Starting subcarrier index of demand (i,j,t)
var W{(i,j,t,k,l,u) in DxTxDxT} binary;     # Preceding slot indicator for demand (i,j,t) relative to (k,l,u)
var p{(i,j,t) in DxT} >=0;                  # Traffic from demand (i,j,t)
var e{(i,j,t,z) in DxTxM} binary;           # Modulation format for demand (i,j,t)
var f{(i,j,t,z) in DxTxM} >=0;              # Equation to calculate the load
var c{t in T} integer >=0;                  # Maximum consumed slot index

# Objective Function ----------------------------------------------------------------------------------------------------


minimize F_OBJETIVO: sum{t in T}c[t] + sum{(i,j,t) in DxT}p[i,j,t];


# Flow Balance Constraints ----------------------------------------------------------------------------------------------

subject to BALANCE_EQUATION { m in N, (i,j,t) in DxT }:
sum{(i,j,m,n) in DxE} P[i,j,t,m,n] - sum{(i,j,n,m) in DxE} P[i,j,t,n,m] =
(if i = m  			then p[i,j,t]
 else if  j = m 	then -p[i,j,t]
 else 				0) ;


subject to DEFINING_A { (i,j,t,m,n) in DxTxE }:
A[i,j,t,m,n] >= P[i,j,t,m,n]/H ;


subject to BALANCOovisnjg { (i,j,t) in DxT, (m,n) in E,(m,k) in E diff{(m,n)} }:
A[i,j,t,m,n] + A[i,j,t,m,k] <= 1 ;

# Modulation Level Constraints  ------------------------------------------------------------------------------------------

subject to DEFINING_p{ (i,j) in D, t in T, z in M }:
p[i,j,t] >= f[i,j,t,z] - (1-e[i,j,t,z])*H ;

subject to DEFINING_MOD_TRAFFIC{ (i,j) in D, z in M, t in T }:
f[i,j,t,z] >= v[i,j,t]/(r*ef[z]) ;

subject to ONLY_ONE_MODULATION { (i,j) in D, t in T }: 
sum{z in M} e[i,j,t,z] <= 1 ;

subject to DEFINING_e1{ (i,j) in D, t in T }:
sum{z in M} e[i,j,t,z] >= v[i,j,t]/H ;

subject to DISTANCE_TO_MODULATION{ (i,j) in D, t in T }:
sum{(m,n) in E} A[i,j,t,m,n] <= sum{z in M} d_[z]*e[i,j,t,z] ;

# Spectrum Continuity and Contiguity Constraints  ------------------------------------------------------------------------

subject to NON_OVERLAPPING_STARTING_SLOTS{ (i,j,t) in DxT, (k,l,u) in DxT diff{(i,j,t)} }:
W[i,j,t,k,l,u] + W[k,l,u,i,j,t] <= 1 ;


subject to NON_OVERLAPPING_STARTING_SLOTS2{ (k,l,u) in DxT, (i,j,t) in DxT diff {(k,l,u)}, (m,n) in E}:
W[i,j,t,k,l,u] + W[k,l,u,i,j,t] >= (A[i,j,t,m,n] + A[k,l,u,m,n])-1 ;


subject to NON_OVERLAPPING_FREQUENCIES{ (i,j,t) in DxT, (k,l,u) in DxT diff{(i,j,t)}}:
s[i,j,t] + p[i,j,t] + G <= s[k,l,u] + H*(1 - W[i,j,t,k,l,u]) ;


subject to NON_OVERLAPPING_FREQUENCIES2{ (i,j,t) in DxT, (k,l,u) in DxT diff{(i,j,t)}}:
s[k,l,u] + p[k,l,u] + G <= s[i,j,t] + H*(1 - W[k,l,u,i,j,t]) ;


# Additional Constraints  ------------------------------------------------------------------------------------------------

subject to FORCING_TRAFFIC{ (i,j) in D, t in T}:
sum{(m,n) in E}A[i,j,t,m,n] >= v[i,j,t]/H ;

subject to DEFINING_MAXSLOTNUMBER{ (i,j) in D, t in T }:
c[t] >= s[i,j,t] + p[i,j,t] ;

subject to STARTING_SLOT_POSITIVE { (i,j) in D, t in T }:
s[i,j,t] >= 0 ;