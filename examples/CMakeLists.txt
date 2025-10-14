# Generate and install examples.

ADD_SUBDIRECTORY(cpp-examples)
ADD_SUBDIRECTORY(python-examples)
INSTALL(DIRECTORY benchmarks DESTINATION examples)
INSTALL(DIRECTORY extras DESTINATION examples)
