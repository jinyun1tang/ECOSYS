configure_file(blk10.h  blk10.h  COPYONLY)
configure_file(blk11a.h blk11a.h COPYONLY)
configure_file(blk11b.h blk11b.h COPYONLY)
configure_file(blk12a.h  blk12a.h COPYONLY)
configure_file(blk12b.h  blk12b.h COPYONLY)
configure_file(blk13a.h  blk13a.h COPYONLY)
configure_file(blk13b.h  blk13b.h COPYONLY)
configure_file(blk13c.h  blk13c.h COPYONLY)
configure_file(blk13d.h  blk13d.h COPYONLY)
configure_file(blk14.h   blk14.h  COPYONLY)
configure_file(blk15a.h  blk15a.h COPYONLY)
configure_file(blk15b.h  blk15b.h COPYONLY)
configure_file(blk16.h   blk16.h  COPYONLY)
configure_file(blk17.h   blk17.h  COPYONLY)
configure_file(blk18a.h  blk18a.h COPYONLY)
configure_file(blk18b.h  blk18b.h COPYONLY)
configure_file(blk19a.h  blk19a.h COPYONLY)
configure_file(blk19b.h  blk19b.h COPYONLY)
configure_file(blk19c.h  blk19c.h COPYONLY)
configure_file(blk19d.h  blk19d.h COPYONLY)
configure_file(blk1cp.h  blk1cp.h COPYONLY)
configure_file(blk1cr.h  blk1cr.h COPYONLY)
configure_file(blk1g.h   blk1g.h COPYONLY)
configure_file(blk1n.h   blk1n.h COPYONLY)
configure_file(blk1p.h   blk1p.h COPYONLY)
configure_file(blk1s.h   blk1s.h COPYONLY)
configure_file(blk1u.h   blk1u.h COPYONLY)
configure_file(blk20a.h  blk20a.h COPYONLY)
configure_file(blk20b.h  blk20b.h COPYONLY)
configure_file(blk20c.h  blk20c.h COPYONLY)
configure_file(blk20d.h  blk20d.h COPYONLY)
configure_file(blk20e.h  blk20e.h COPYONLY)
configure_file(blk20f.h  blk20f.h COPYONLY)
configure_file(blk21a.h  blk21a.h COPYONLY)
configure_file(blk21b.h  blk21b.h COPYONLY)
configure_file(blk22a.h  blk22a.h COPYONLY)
configure_file(blk22b.h  blk22b.h COPYONLY)
configure_file(blk22c.h  blk22c.h COPYONLY)
configure_file(blk2a.h   blk2a.h COPYONLY)
configure_file(blk2b.h   blk2b.h COPYONLY)
configure_file(blk2c.h   blk2c.h COPYONLY)
configure_file(blk3.h    blk3.h COPYONLY)
configure_file(blk5.h    blk5.h COPYONLY)
configure_file(blk6.h    blk6.h COPYONLY)
configure_file(blk8a.h   blk8a.h COPYONLY)
configure_file(blk8b.h   blk8b.h COPYONLY)
configure_file(blk9a.h   blk9a.h COPYONLY)
configure_file(blk9b.h   blk9b.h COPYONLY)
configure_file(blk9c.h   blk9c.h COPYONLY)
configure_file(blkc.h    blkc.h COPYONLY)
configure_file(blktest.h  blktest.h COPYONLY)
configure_file(filec.h  filec.h COPYONLY)
configure_file(files.h  files.h COPYONLY)
configure_file(parameters.h parameters.h COPYONLY)
set(ECOSYS_SOURCES 
	BLOCKDATA001.f
	day.f
	endrun.f
	erosion.f
	exec.f
	extract.f
	foutp.f
	fouts.f
	grosub.f
	hfunc.f
	hour1.f
	nitro.f
	outpd.f
	outph.f
	outsd.f
	outsh.f
	omwriter.f
	readi.f
	readq.f
	reads.f
	redist.f
	routp.f
	routq.f
	routs.f
	solute.f
	split.f
	splitc.f
	starte.f
	startq.f
	starts.f
	stomate.f
	trnsfr.f
	trnsfrs.f
	uptake.f
	visual.f
	watsub.f
	woutp.f
	woutq.f
	wouts.f
	wthr.f
	splitp.c
	splits.c	
	)

	include(add_ecosys_library)
	add_ecosys_library(core ${ECOSYS_SOURCES})

	set(ECOSYS_LIBRARIES core;${ECOSYS_LIBRARIES} PARENT_SCOPE)

	if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
	  install(TARGETS core DESTINATION lib)
	  file(GLOB HEADERS *.h)
	  install(FILES ${HEADERS} DESTINATION include/ecosys/core)
	endif()
