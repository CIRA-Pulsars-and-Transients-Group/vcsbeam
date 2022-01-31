#ifndef __UNITS_H__
#define __UNITS_H__

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <complex>

/** Class for base (SI) units */
class SIUnit
{
    public:
        std::string name;  /** The full name of the unit */
        std::string abbr;  /** The abbreviated name of the unit */

        /**
         * Constructor for the SIUnit class
         *
         * @param name The (full) name for the unit
         * @param abbr The abbreviated name for the unit
         *
         * This is a minimum constructor type, which must include a name and
         * an abbreviation for the unit. The unit is set to be dimensionless
         * (i.e. a scalar).
         */
        SIUnit( std::string name, std::string abbr, double scale_factor = 1.0 );


        /**
         * Constructor for the SIUnit class, including dimension specification
         *
         * @param name  The (full) name for the unit
         * @param abbr  The abbreviated name for the unit
         * @param dims  A formatted string specifying the base SI dimensions
         *
         * This constructor can be used to fully specify the base SI dimensions
         * of a newly defined unit. `dims` is a `map` whose keys are strings
         * drawn from the set ["length", "time", "mass", "current",
         * "temperature", "mole", "candela"], and whose values are integer-
         * double pairs, representing the dimension index and scaling factor
         * respectively.
         *
         * Strings that are not members of the above set are ignored, and
         * members of the above set which are not found in `dims` are assigned
         * values of <0,1.0>.
         *
         * For example, to generate the Gauss unit (for magnetic fields),
         * the call would be:
         *
         * ```
         * SIUnit unitGauss( "Gauss", "G",
         *     {
         *         { "mass",     1 },  // kg
         *         { "time",    -2 },  // s^-2
         *         { "current", -1 }   // A^-1
         *     },
         *     1.0e-4 );               // 1 G = 10^{-4} kg s^-2 A^-1
         *
         * ```
         */
        SIUnit( std::string name,
                std::string abbr,
                std::map<std::string, int> dims,
                double scale_factor = 1.0 );

        /**
         * Deconstructor for the SIUnit class
         */
        ~SIUnit() {}

        /**
         * Overloads the << operator
         *
         * The default behaviour is to output the abbreviated name.
         */
        friend std::ostream& operator<<(std::ostream& os, const SIUnit& u);

        /**
         * Overloads the (binary) * operator
         *
         * This "multiplies" the units together.
         */
        SIUnit operator*(const SIUnit& u);

        /**
         * Overloads the (binary) / operator
         *
         * This "divides" the units.
         */
        SIUnit operator/(const SIUnit& u);

        /**
         * Returns a string representation in base SI units
         */
        std::string base_units_str();

    private:

        std::map<std::string, int> dims; /** The dimensions and scales of each SI base unit */
        double scale_factor;             /** The scale factor for conversion to base SI units */

        /**
         * Initialises the dimensions
         *
         * Initialises the dimensions to a scalar type (where all SI base
         * units have dimension 0, and scale factor 1.0). This function is
         * intended only to be used by various constructors.
         */
        void init_dims();

};


#endif
