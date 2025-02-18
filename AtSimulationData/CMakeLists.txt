# Create a library called "libAtSimulationData"

Set(LIBRARY_NAME AtSimulationData)
Set(LINKDEF ${LIBRARY_NAME}LinkDef.h)

set(DEPENDENCIES
  FairRoot::Base # FairRunIdGenerator
  FairRoot::FairTools # FairLogger
  Boost::headers
  ROOT::Core
  )

set(SRCS
  AtStack.cxx
  AtMCTrack.cxx
  AtVertexPropagator.cxx
  AtMCPoint.cxx
)

set(TEST_SRCS
  AtVertexPropagatorTest.cxx
)

attpcroot_generate_tests(${LIBRARY_NAME}Tests
  SRCS ${TEST_SRCS}
  DEPS ${LIBRARY_NAME}
  )

generate_target_and_root_library(${LIBRARY_NAME}
  LINKDEF ${LINKDEF}
  SRCS ${SRCS}
  DEPS_PUBLIC ${DEPENDENCIES}
  )

