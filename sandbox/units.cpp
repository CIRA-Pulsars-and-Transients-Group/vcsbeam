#include "units.h"

SIUnit::SIUnit(
        std::string name,
        std::string abbr,
        double scale_factor ):
    name(name),
    abbr(abbr),
    scale_factor(scale_factor)
{
    this->init_dims();
}

SIUnit::SIUnit(
        std::string name,
        std::string abbr,
        std::map<std::string, int> dims,
        double scale_factor ):
    name(name),
    abbr(abbr),
    scale_factor(scale_factor)
{
    this->init_dims();

    // Iterate through the keys in the provided dims, and,
    // if one of the keys is also one of the standard ones
    // initialised in init_dims(), then replace its value with
    // the provided value.
    std::map<std::string, int>::iterator it, it_found;
    for (it = dims.begin(); it != dims.end(); it++)
    {
        it_found = this->dims.find( it->first );
        if (it_found != this->dims.end())
            it_found->second = it->second;
    }
}

void SIUnit::init_dims()
{
    this->dims = {
            { "length",      0 },
            { "time",        0 },
            { "mass",        0 },
            { "current",     0 },
            { "temperature", 0 },
            { "mole",        0 },
            { "candela",     0 }
        };
}


std::ostream& operator<<(std::ostream& os, const SIUnit& u)
{
    os << u.abbr;
    return os;
}

SIUnit SIUnit::operator/(const SIUnit& u)
{
    // Call the default copy constructor
    SIUnit uNew = u;
    uNew.name = this->name + " per " + u.name;
    uNew.abbr = this->abbr + " (" + u.abbr + ")^-1";

    // Subtract the dimension indices
    std::map<std::string, int>::iterator it, it_found;
    for (it = this->dims.begin(); it != this->dims.end(); it++)
    {
        it_found = uNew.dims.find( it->first );
        if (it_found != uNew.dims.end())
            it_found->second = it->second - it_found->second;
    }

    // Divide the scale_factors
    uNew.scale_factor = this->scale_factor / uNew.scale_factor;

    return uNew;
}
SIUnit SIUnit::operator*(const SIUnit& u)
{
    // Call the default copy constructor
    SIUnit uNew = u;
    uNew.name = this->name + " " + u.name;
    uNew.abbr = this->abbr + " " + u.abbr;

    // Add the dimension indices
    std::map<std::string, int>::iterator it, it_found;
    for (it = this->dims.begin(); it != this->dims.end(); it++)
    {
        it_found = uNew.dims.find( it->first );
        if (it_found != uNew.dims.end())
            it_found->second += it->second;
    }

    // Multiply the scale_factors
    uNew.scale_factor *= this->scale_factor;

    return uNew;
}


std::string SIUnit::base_units_str()
{
    char buffer[32];
    std::string out = "";

    if (this->scale_factor != 1.0)
    {
        snprintf(buffer, sizeof(buffer), "%g ", this->scale_factor);
        std::string s(buffer);
        out += s;
    }

    if (this->dims["mass"] != 0)
        out += "kg^"  + std::to_string( this->dims["mass"]        ) + " ";
    if (this->dims["length"] != 0)
        out += "m^"   + std::to_string( this->dims["length"]      ) + " ";
    if (this->dims["time"] != 0)
        out += "s^"   + std::to_string( this->dims["time"]        ) + " ";
    if (this->dims["current"] != 0)
        out += "A^"   + std::to_string( this->dims["current"]     ) + " ";
    if (this->dims["temperature"] != 0)
        out += "K^"   + std::to_string( this->dims["temperature"] ) + " ";
    if (this->dims["mole"] != 0)
        out += "mol^" + std::to_string( this->dims["mole"]        ) + " ";
    if (this->dims["candela"] != 0)
        out += "cd^"  + std::to_string( this->dims["candela"]     ) + " ";

    return out;
}
