# Orbit_Transfer_Direct_Shooting_Method

This project presents a optimal orbit transfer developed for a low-thrust planar transfer from a Low Earth Orbit (LOR) to a Geosynchronous Orbit (GEO). The optimal orbit transfer problem is modeled using a two-point boundary value problem (TPBVP). In order to minimize the fuel consumption, or what is the same the final time. The two point boundary value problem is solved using Legendre-Gauss-Radau collocation. IPOPT was used to solve highter order differential equaitons and ADiGator a differential equation solver was used to find the differential equation of highter order differential equations using MATLAB.
The problem was solved for 2 different objectives, (with thrust T and thrust angle ß as control inputs) 
  1) To minimize the final time or to find the fastest path.
  2) To maximize the final mass.
