ice_read_write.o ice_read_write.d : ice_read_write.F90
ice_read_write.o : ice_kinds_mod.o
ice_read_write.o : ice_constants.o
ice_read_write.o : ice_communicate.o
ice_read_write.o : ice_broadcast.o
ice_read_write.o : ice_domain.o
ice_read_write.o : ice_domain_size.o
ice_read_write.o : ice_blocks.o
ice_read_write.o : ice_exit.o
ice_read_write.o : ice_fileunits.o
ice_read_write.o : ice_gather_scatter.o
