SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (global_track_merger_sources
  ${src_dir}/GlobalTrackMergerBR.cc
  ${src_dir}/GlobalTrackMergerBR.h
  PARENT_SCOPE
)

SET (global_track_merger_data
  ${src_dir}/GlobalTrackMergerBR.xml
  ${src_dir}/GlobalTrackMergerBR.xsd
  PARENT_SCOPE
)
