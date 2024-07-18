# The HVRP application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP bapdock /HVRP/src/run.jl /HVRP/data/c50_13hvrp.txt -u 3185.2 -o /HVRP/sol/c50_13hvrp.sol

# toy instance (optimal is 165.86)
docker run --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP bapdock /HVRP/src/run.jl /HVRP/data/toy.txt

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP bapdock /HVRP/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP bapdock /HVRP/src/run.jl -b /HVRP/baldacci.batch

# The application directory (/ABSOLUTE_PATH_TO_HVRP_APP) was mounted with -v as /HVRP inside the container. Also, it is possible to mount a different directory to read/write solutions:
docker run --rm -v /ABSOLUTE_PATH_TO_HVRP_APP:/HVRP -v /ABSOLUTE_PATH_TO_OUTPUT:/OUT bapdock /HVRP/src/run.jl /HVRP/data/c50_13hvrp.txt -u 3185.2 -o /OUT/c50_13hvrp.sol

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/c50_13hvrp.txt -u 3185.2 -o sol/c50_13hvrp.sol

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see baldacci.batch before for adjustments):
./VRPSolver -b baldacci.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
