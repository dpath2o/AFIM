ice_broadcast.o ice_broadcast.d : ice_broadcast.F90
ice_broadcast.o : ice_kinds_mod.o
ice_broadcast.o : ice_communicate.o
ice_broadcast.o : ice_exit.o
ice_broadcast.o : icepack_intfc.o
