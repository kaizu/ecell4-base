#ifndef __ECELL4_SPECIES_HPP
#define __ECELL4_SPECIES_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "config.h"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "UnitSpecies.hpp"


namespace ecell4
{

class Species
{
public:

    typedef UnitSpecies::serial_type serial_type; //XXX: std::string
    typedef std::vector<UnitSpecies> container_type;

protected:

    typedef utils::get_mapper_mf<std::string, std::string>::type
    attributes_container_type;

public:

    Species()
        : serial_("")
    {
        ; // do nothing
    }

    explicit Species(const serial_type& name)
        : serial_(name)
    {
        ;
    }

    Species(
        const serial_type& name, const std::string& radius, const std::string& D,
        const std::string location = "")
        : serial_(name)
    {
        set_attribute("radius", radius);
        set_attribute("D", D);
        set_attribute("location", location);
    }

    const serial_type serial() const
    {
        return serial_;
    }

    void add_unit(const UnitSpecies& usp);

    const std::vector<UnitSpecies> units() const
    {
        std::vector<std::string> unit_serials;
        boost::split(unit_serials, serial_, boost::is_any_of("."));

        std::vector<UnitSpecies> units_;
        for (std::vector<std::string>::const_iterator i(unit_serials.begin());
            i != unit_serials.end(); ++i)
        {
            UnitSpecies usp;
            usp.deserialize(*i);
            units_.insert(std::lower_bound(units_.begin(), units_.end(), usp), usp);
        }
        return units_;
    }

    const attributes_container_type& attributes() const
    {
        return attributes_;
    }

    std::vector<std::pair<std::string, std::string> > list_attributes();
    std::string get_attribute(const std::string& name_attr) const;
    void set_attribute(const std::string& name_attr, const std::string& value);
    void set_attributes(const Species& sp);
    void overwrite_attributes(const Species& sp);
    void remove_attribute(const std::string& name_attr);
    bool has_attribute(const std::string& name_attr) const;

    bool operator==(const Species& rhs) const;
    bool operator!=(const Species& rhs) const;
    bool operator<(const Species& rhs) const;
    bool operator>(const Species& rhs) const;

    Integer count(const Species& sp) const;

    /** for epdp
     */
    serial_type name() const
    {
        return serial();
    }

protected:

    serial_type serial_;
    attributes_container_type attributes_;
};

Species format_species(const Species& sp);

inline Species::serial_type unique_serial(const Species& sp)
{
    return format_species(sp).serial();
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm,
    const ecell4::Species& sp)
{
    strm << sp.serial();
    return strm;
}

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::Species>
{
    std::size_t operator()(const ecell4::Species& val) const
    {
        return hash<ecell4::Species::serial_type>()(val.serial());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_SPECIES_HPP */
