#define MAXBCREGIONS 50000
#define MAXBCBLOCKS 50000
#define MAXSTRINGLENGTH 512
#define MAXWORDLENGTH 32
#define MAXCARDLENGTH 4
#define MAXNAMELENGTH 20
#define MAXPERMREGIONS 35000
!#define MAXINITREGIONS 80000
#define MAXINITREGIONS 8400000
#define MAXSRC 10
#define MAXSRCTIMES 100
#define IUNIT1 15
#define IUNIT2 16
#define IUNIT3 17
#define IUNIT4 18
#define HHISTORY_LENGTH 1000
! HHISTORY_LENGTH is the length of the array used to store the differencing
! values h.

#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 3

! Macros that are used as 'dm_index' values.  --RTM
#define ONEDOF 1
#define NPHASEDOF 2
#define THREENPDOF 3
#define NDOF 4
#define NPHANCOMPDOF 5
#define NPHANSPECDOF 6
#define NPHANSPECNCOMPDOF 7
#define VARDOF 8

#define GLOBAL 1
#define LOCAL 2
#define NATURAL 3

! modes
#define NULL_MODE 0
#define RICHARDS_MODE 1
#define MPH_MODE 2
#define COND_MODE 3
#define TWOPH_MODE 4
#define VADOSE_MODE 5
#define LIQUID_MODE 6
#define OWG_MODE 7
#define FLASH_MODE 8
#define TH_MODE 9
#define THC_MODE 10

! grid types
#define STRUCTURED 1
#define UNSTRUCTURED 2
#define STRUCTURED_CARTESIAN 10
#define STRUCTURED_CYLINDRICAL 11
#define STRUCTURED_SPHERICAL 12

! condition types
#define DIRICHLET_BC 1
#define NEUMANN_BC 2
#define MASS_RATE 3
#define ZERO_GRADIENT_BC 5
#define HYDROSTATIC_BC 6

! coupler types
#define INITIAL_COUPLER_TYPE 1
#define BOUNDARY_COUPLER_TYPE 2
#define SRC_SINK_COUPLER_TYPE 3
#define COUPLER_IPHASE_INDEX 1

! connection types
#define INTERNAL_CONNECTION_TYPE 1
#define BOUNDARY_CONNECTION_TYPE 2
#define INITIAL_CONNECTION_TYPE 3
#define SRC_SINK_CONNECTION_TYPE 4

! dofs for each mode
#define RICHARDS_PRESSURE_DOF 1
#define RICHARDS_CONCENTRATION_DOF 3
#define RICHARDS_TEMPERATURE_DOF 2
#define RICHARDS_ENTHALPY_DOF 3

#define MPH_PRESSURE_DOF 1
#define MPH_CONCENTRATION_DOF 3
#define MPH_TEMPERATURE_DOF 2
#define MPH_ENTHALPY_DOF 3

#define THC_PRESSURE_DOF 1
#define THC_CONCENTRATION_DOF 3
#define THC_TEMPERATURE_DOF 2
#define THC_ENTHALPY_DOF 3
