SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (local_track_finder_sources
  ${src_dir}/LocalTrackFinderBR.cc
  ${src_dir}/LocalTrackFinderBR.h
  PARENT_SCOPE
)

SET (local_track_finder_data
  ${src_dir}/LocalTrackFinderBR.xml
  ${src_dir}/LocalTrackFinderBR.xsd
  PARENT_SCOPE
)
