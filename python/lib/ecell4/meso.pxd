from libcpp.string cimport string
from libcpp cimport bool

from core cimport *


## Cpp_MesoscopicWorld
#  ecell4::meso::MesoscopicWorld
cdef extern from "ecell4/meso/MesoscopicWorld.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_MesoscopicWorld "ecell4::meso::MesoscopicWorld":
        Cpp_MesoscopicWorld(Cpp_Position3&, Integer, Integer, Integer) except +
        Cpp_MesoscopicWorld(Cpp_Position3&, Integer, Integer, Integer, shared_ptr[Cpp_RandomNumberGenerator]) except +
        void set_t(Real)
        Real t()
        Real volume()
        Real subvolume()
        Integer num_subvolumes()
        void set_edge_lengths(Cpp_Position3&)
        Cpp_Position3 edge_lengths()
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        Integer num_molecules(Cpp_Species &, Integer)
        Integer num_molecules_exact(Cpp_Species &, Integer)
        Integer num_molecules(Cpp_Species &, Cpp_Global)
        Integer num_molecules_exact(Cpp_Species &, Cpp_Global)
        vector[Cpp_Species] list_species()
        void add_molecules(Cpp_Species &sp, Integer &num, Integer)
        void remove_molecules(Cpp_Species &sp, Integer &num, Integer)
        void add_molecules(Cpp_Species &sp, Integer &num, Cpp_Global)
        void remove_molecules(Cpp_Species &sp, Integer &num, Cpp_Global)
        # void save(string)
        # void load(string)
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)

## MesoscopicWorld
#  a python wrapper for Cpp_MesoscopicWorld
cdef class MesoscopicWorld:
    cdef shared_ptr[Cpp_MesoscopicWorld]* thisptr

cdef MesoscopicWorld MesoscopicWorld_from_Cpp_MesoscopicWorld(
    shared_ptr[Cpp_MesoscopicWorld] m)

## Cpp_MesoscopicSimulator
#  ecell4::meso::MesoscopicSimulator
cdef extern from "ecell4/meso/MesoscopicSimulator.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_MesoscopicSimulator "ecell4::meso::MesoscopicSimulator":
        Cpp_MesoscopicSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_MesoscopicWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real)
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        vector[Cpp_ReactionRule] last_reactions()
        void initialize()
        # Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_MesoscopicWorld] world()
        void run(Real)
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## MesoscopicSimulator
#  a python wrapper for Cpp_MesoscopicSimulator
cdef class MesoscopicSimulator:
    cdef Cpp_MesoscopicSimulator* thisptr
