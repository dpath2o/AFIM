ice_domain.o ice_domain.d : ice_domain.F90
ice_domain.o : ice_kinds_mod.o
ice_domain.o : ice_constants.o
ice_domain.o : ice_communicate.o
ice_domain.o : ice_broadcast.o
ice_domain.o : ice_blocks.o
ice_domain.o : ice_distribution.o
ice_domain.o : ice_boundary.o
ice_domain.o : ice_exit.o
ice_domain.o : ice_fileunits.o
ice_domain.o : icepack_intfc.o
ice_domain.o : ice_domain_size.o
