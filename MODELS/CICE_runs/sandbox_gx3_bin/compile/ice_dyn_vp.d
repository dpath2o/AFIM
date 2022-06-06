ice_dyn_vp.o ice_dyn_vp.d : ice_dyn_vp.F90
ice_dyn_vp.o : ice_kinds_mod.o
ice_dyn_vp.o : ice_blocks.o
ice_dyn_vp.o : ice_boundary.o
ice_dyn_vp.o : ice_communicate.o
ice_dyn_vp.o : ice_constants.o
ice_dyn_vp.o : ice_domain.o
ice_dyn_vp.o : ice_domain_size.o
ice_dyn_vp.o : ice_dyn_shared.o
ice_dyn_vp.o : ice_fileunits.o
ice_dyn_vp.o : ice_flux.o
ice_dyn_vp.o : ice_global_reductions.o
ice_dyn_vp.o : ice_grid.o
ice_dyn_vp.o : ice_exit.o
ice_dyn_vp.o : icepack_intfc.o
ice_dyn_vp.o : ice_arrays_column.o
ice_dyn_vp.o : ice_state.o
ice_dyn_vp.o : ice_timers.o
