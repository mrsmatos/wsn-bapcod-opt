# The CARP application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP bapdock /CARP/src/run.jl /CARP/data/egl-e1-A.dat -u 3549 -o /CARP/sol/egl-e1-A.sol

# toy instance (optimal is 7.0)
docker run --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP bapdock /CARP/src/run.jl /CARP/data/toy.dat

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP bapdock /CARP/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP bapdock /CARP/src/run.jl -b /CARP/eglese.batch

# The application directory (/ABSOLUTE_PATH_TO_CARP_APP) was mounted with -v as /CARP inside the container. Also, it is possible to mount a different directory to read/write solutions:
docker run --rm -v /ABSOLUTE_PATH_TO_CARP_APP:/CARP -v /ABSOLUTE_PATH_TO_OUTPUT:/OUT bapdock /CARP/src/run.jl /CARP/data/egl-e1-A.dat -u 3549 -o /OUT/egl-e1-A.sol

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/egl-e1-A.dat -u 3549 -o sol/egl-e1-A.sol

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see eglese.batch before for adjustments):
./VRPSolver -b eglese.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
