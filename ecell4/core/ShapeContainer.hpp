#ifndef __ECELL4_SHAPE_CONTAINER_HPP
#define __ECELL4_SHAPE_CONTAINER_HPP

#include "config.h"

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif


#include <utility>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "Real3.hpp"
#include "Identifier.hpp"
#include "Shape.hpp"
#include "Species.hpp"
#include "PlanarSurface.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <vector>
#include <map>
#include <algorithm>
#include <functional>

namespace ecell4
{

struct PlanarSurfaceID:
    public Identifier<PlanarSurfaceID, unsigned long long, int>
{
    typedef Identifier<PlanarSurfaceID, unsigned long long, int> base_type;
    PlanarSurfaceID(const value_type& value=value_type(0, 0)): base_type(value)
    {
        ;
    }
};

}


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
struct hash<ecell4::PlanarSurfaceID>
{
    std::size_t operator()(const ecell4::PlanarSurfaceID& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
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





namespace ecell4
{

template <typename T_>
struct compare_second {
    bool operator() (const T_ &lhs, const T_ &rhs) const
    {
        return lhs.second < rhs.second;
    }
};

class PlanarSurfaceContainer
{
public:
    typedef PlanarSurfaceID surface_id_type;
    typedef PlanarSurface    surface_type;
    typedef Species          species_type;

    typedef std::pair<species_type, surface_type> species_surface_pair;
    typedef std::vector<std::pair<surface_id_type, species_surface_pair> > surface_container_type;
    
    typedef std::pair<surface_id_type, Real> id_distance_pair;

protected:
    // ID -> index of surface
    typedef utils::get_mapper_mf<
        surface_id_type, surface_container_type::size_type>::type surface_map_type;


public:
    PlanarSurfaceContainer(void);
    Integer num_surfaces(void) const
    {
        return surfaces_.size();
    }
    const surface_container_type list_surfaces() const
    {
        return surfaces_;
    }
    const surface_container_type list_surfaces(species_type const &sp) const 
    {
        surface_container_type retval;
        for(surface_container_type::const_iterator it = surfaces_.begin(); it != surfaces_.end(); it++) 
        {
            if (it->second.first == sp)
            {
                retval.push_back(*it);
            }
        }
        return retval;
    }
    const species_surface_pair &get_surface(surface_id_type const &id) const
    {
        surface_map_type::const_iterator i(index_map_.find(id));
        if (i != index_map_.end()) 
        {
            return surfaces_[(*i).second].second;
        }
        else
        {
            throw NotFound("Surface not found");
        }
    }

    bool has_surface(surface_id_type const &id) const
    {
        surface_map_type::const_iterator i(index_map_.find(id));
        return (i != index_map_.end());
    }
    bool update_surface(const surface_id_type &id, const species_type &sp, const surface_type &surface)
    {
        surface_map_type::const_iterator i(index_map_.find(id));
        if (i == index_map_.end())
        {
            surface_container_type::size_type idx(surfaces_.size());
            index_map_[id] = idx;
            surfaces_.push_back(std::make_pair(id, std::make_pair(sp, surface)));
            return true;
        }
        else
        {
            surfaces_[(*i).second] = std::make_pair(id, std::make_pair(sp, surface));
            return false;
        }
    }

    std::vector<id_distance_pair> list_id_distance_pair(const Real3 &pos, bool do_sort = false) const
    {
        std::vector<id_distance_pair> ret;
        for(surface_container_type::const_iterator it = surfaces_.begin(); it != surfaces_.end(); it++)
        {
            Real dist = std::abs(it->second.second.is_inside(pos));
            ret.push_back( std::make_pair(it->first, dist) );
        }
        if (do_sort == true) {
            std::sort(ret.begin(), ret.end(), compare_second<id_distance_pair>() );
        }
        return ret;
    }

    Real3 apply_reflection(const Real3 &from, const Real3 &displacement) const;

    void disable_surface_reflection(void)
    {
       this->surface_reflection_enabled_ = false;
    }
    void enable_surface_reflection(void)
    {
       this->surface_reflection_enabled_ = true;
    }

protected:
    surface_container_type surfaces_;
    surface_map_type index_map_;

    bool surface_reflection_enabled_;
};

}   //ecell4


#endif
