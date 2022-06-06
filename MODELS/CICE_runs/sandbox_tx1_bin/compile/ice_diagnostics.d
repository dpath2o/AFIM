ice_diagnostics.o ice_diagnostics.d : ice_diagnostics.F90
ice_diagnostics.o : ice_kinds_mod.o
ice_diagnostics.o : ice_communicate.o
ice_diagnostics.o : ice_constants.o
ice_diagnostics.o : ice_calendar.o
ice_diagnostics.o : ice_domain_size.o
ice_diagnostics.o : ice_fileunits.o
ice_diagnostics.o : ice_exit.o
ice_diagnostics.o : icepack_intfc.o
ice_diagnostics.o : ice_arrays_column.o
ice_diagnostics.o : ice_blocks.o
ice_diagnostics.o : ice_broadcast.o
ice_diagnostics.o : ice_domain.o
ice_diagnostics.o : ice_flux.o
ice_diagnostics.o : ice_flux_bgc.o
ice_diagnostics.o : ice_global_reductions.o
ice_diagnostics.o : ice_grid.o
ice_diagnostics.o : ice_state.o
