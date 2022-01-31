#include <iostream>
#include "units.h"

int main()
{
    SIUnit u_Hz( "Hertz", "Hz", { {"time", -1} } );
    SIUnit u_pc( "parsec", "pc", { {"length", 1} }, 3.0857e16 );
    SIUnit u_G( "Gauss", "G",
              {
                  { "mass",     1 },  // kg
                  { "time",    -2 },  // s^-2
                  { "current", -1 }   // A^-1
              },
              1.0e-4 );
    SIUnit u_cm( "centimeter", "cm", { {"length", 1} }, 1e-2 );

    SIUnit u_pcG = u_pc / u_G;
    std::cout << u_pcG.name << " (" << u_pcG << ") = " << u_pcG.base_units_str() << std::endl;

    SIUnit u_DM = u_pc / (u_cm * u_cm * u_cm);
    std::cout << u_DM.name << " (" << u_DM << ") = " << u_DM.base_units_str() << std::endl;

    return 0;
}
