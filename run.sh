echo "Compiling Code"
ifort modules.f90 UTIL_PRE_SIM.f90 UTIL_POST_PROC.f90 TDMA.f90 UTIL_SOLVE.f90 UTIL_IMMERSED_BNDRY.f90 UTIL_BC.f90 UTIL_STATS.f90 MAIN.f90
echo "Done Compiling"
echo "Initiating Simulation"
./a.out
echo "Simulation Done"