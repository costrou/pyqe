"""
Contains all of the default key documentation.
Moving to module so that I can one day make this work
with the already done documentation (autoupdating then)
"""

def getPWDocForKey(name, key):
    """Adds doc strings to key info lists
    """
    doc_map = {
        "CONTROL": control_doc,
        "SYSTEM": system_doc,
        "ELECTRONS": electrons_doc,
        "IONS": ions_doc,
        "CELL": cell_doc
        }

    namelist_doc = doc_map.get(name)
    if not namelist_doc:
        error_str = "unknown namelist {0} unable to add docs"
        raise Exception(error_str.format(name))
    else:
        key_doc = namelist_doc.get(key)
        if not key_doc:
            error_str = "no documentation for key '{0}'"
            raise Exception(error_str.format(key))
        return key_doc

cell_doc = {
    'cell_dynamics': """
Specify the type of dynamics for the cell.
For different type of calculation different possibilities
are allowed and different default values apply:

CASE ( calculation = 'vc-relax' )
  'none':    no dynamics
  'sd':      steepest descent ( not implemented )
  'damp-pr': damped (Beeman) dynamics of the Parrinello-Rahman
             extended lagrangian
  'damp-w':  damped (Beeman) dynamics of the new Wentzcovitch
             extended lagrangian
  'bfgs':    BFGS quasi-newton algorithm (default)
             ion_dynamics must be 'bfgs' too
CASE ( calculation = 'vc-md' )
  'none':    no dynamics
  'pr':      (Beeman) molecular dynamics of the Parrinello-Rahman
             extended lagrangian
  'w':       (Beeman) molecular dynamics of the new Wentzcovitch
             extended lagrangian
""",
    'press': """
Target pressure [KBar] in a variable-cell md or relaxation run.
""",
    'wmass': """
Fictitious cell mass [amu] for variable-cell simulations
(both 'vc-md' and 'vc-relax')
""",
    'cell_factor': """
Used in the construction of the pseudopotential tables.
It should exceed the maximum linear contraction of the
cell during a simulation.
""",
    'press_conv_thr': """
Convergence threshold on the pressure for variable cell
relaxation ('vc-relax' : note that the other convergence
thresholds for ionic relaxation apply as well).
""",
    'cell_dofree': """
Select which of the cell parameters should be moved:

all     = all axis and angles are moved
x       = only the x component of axis 1 (v1_x) is moved
y       = only the y component of axis 2 (v2_y) is moved
z       = only the z component of axis 3 (v3_z) is moved
xy      = only v1_x and v2_y are moved
xz      = only v1_x and v3_z are moved
yz      = only v2_y and v3_z are moved
xyz     = only v1_x, v2_y, v3_z are moved
shape   = all axis and angles, keeping the volume fixed
volume  = the volume changes, keeping all angles fixed (i.e. only celldm(1) changes)
2Dxy    = only x and y components are allowed to change
2Dshape = as above, keeping the area in xy plane fixed

BEWARE: if axis are not orthogonal, some of these options do not
 work (symmetry is broken). If you are not happy with them,
 edit subroutine init_dofree in file Modules/cell_base.f90
"""
}

ions_doc = {
    'ion_dynamics': """
Specify the type of ionic dynamics.

For different type of calculation different possibilities are
allowed and different default values apply:

CASE ( calculation = 'relax' )
    'bfgs' :   (default)   use BFGS quasi-newton algorithm,
                           based on the trust radius procedure,
                           for structural relaxation
    'damp' :               use damped (quick-min Verlet)
                           dynamics for structural relaxation
                           Can be used for constrained
                           optimisation: see CONSTRAINTS card

CASE ( calculation = 'md' )
    'verlet' : (default)   use Verlet algorithm to integrate
                           Newton's equation. For constrained
                           dynamics, see CONSTRAINTS card
    'langevin'             ion dynamics is over-damped Langevin
    'langevin-smc'         over-damped Langevin with Smart Monte Carlo:
                           see R.J.Rossky, JCP, 69, 4628(1978)


CASE ( calculation = 'vc-relax' )
    'bfgs' :   (default)   use BFGS quasi-newton algorithm;
                           cell_dynamics must be 'bfgs' too
    'damp' :               use damped (Beeman) dynamics for
                           structural relaxation
CASE ( calculation = 'vc-md' )
    'beeman' : (default)   use Beeman algorithm to integrate
                           Newton's equation
""",
    'ion_positions': """
'default '  : if restarting, use atomic positions read from the
              restart file; in all other cases, use atomic
              positions from standard input.

'from_input' : restart the simulation with atomic positions read
              from standard input, even if restarting.
""",
    'pot_extrapolation': """
 Used to extrapolate the potential from preceding ionic steps.

   'none'        :  no extrapolation

   'atomic'      :  extrapolate the potential as if it was a sum of
                    atomic-like orbitals

   'first_order' :  extrapolate the potential with first-order
                    formula

   'second_order':  as above, with second order formula

Note: 'first_order' and 'second-order' extrapolation make sense
only for molecular dynamics calculations
""",
    'wfc_extrapolation': """
Used to extrapolate the wavefunctions from preceding ionic steps.

   'none'        :  no extrapolation

   'first_order' :  extrapolate the wave-functions with first-order
                    formula.

   'second_order':  as above, with second order formula.

Note: 'first_order' and 'second-order' extrapolation make sense
only for molecular dynamics calculations
""",
    'remove_rigid_rot': """
This keyword is useful when simulating the dynamics and/or the
thermodynamics of an isolated system. If set to true the total
torque of the internal forces is set to zero by adding new forces
that compensate the spurious interaction with the periodic
images. This allows for the use of smaller supercells.

BEWARE: since the potential energy is no longer consistent with
the forces (it still contains the spurious interaction with the
repeated images), the total energy is not conserved anymore.
However the dynamical and thermodynamical properties should be
in closer agreement with those of an isolated system.
Also the final energy of a structural relaxation will be higher,
but the relaxation itself should be faster.
""",
    'ion_temperature': """
'rescaling'   control ionic temperature via velocity rescaling
              (first method) see parameters "tempw", "tolp", and
              "nraise" (for VC-MD only). This rescaling method
              is the only one currently implemented in VC-MD

'rescale-v'   control ionic temperature via velocity rescaling
              (second method) see parameters "tempw" and "nraise"

'rescale-T'   control ionic temperature via velocity rescaling
              (third method) see parameter "delta_t"

'reduce-T'    reduce ionic temperature every "nraise" steps
              by the (negative) value "delta_t"

'berendsen'   control ionic temperature using "soft" velocity
              rescaling - see parameters "tempw" and "nraise"

'andersen'    control ionic temperature using Andersen thermostat
              see parameters "tempw" and "nraise"

'initial'     initialize ion velocities to temperature "tempw"
              and leave uncontrolled further on

'not_controlled' (default) ionic temperature is not controlled
""",
    'tempw': """
Starting temperature (Kelvin) in MD runs
target temperature for most thermostats.
""",
    'tolp': """
Tolerance for velocity rescaling. Velocities are rescaled if
the run-averaged and target temperature differ more than tolp.
""",
    'delta_t': """
if ion_temperature='rescale-T':
       at each step the instantaneous temperature is multiplied
       by delta_t; this is done rescaling all the velocities.

if ion_temperature='reduce-T':
       every 'nraise' steps the instantaneous temperature is
       reduced by -delta_T (i.e. delta_t < 0 is added to T)

The instantaneous temperature is calculated at the end of
every ionic move and BEFORE rescaling. This is the temperature
reported in the main output.

For delta_t < 0, the actual average rate of heating or cooling
should be roughly C*delta_t/(nraise*dt) (C=1 for an
ideal gas, C=0.5 for a harmonic solid, theorem of energy
equipartition between all quadratic degrees of freedom).
""",
    'nraise': """
if ion_temperature='reduce-T':
       every 'nraise' steps the instantaneous temperature is
       reduced by -delta_T (.e. delta_t is added to the temperature)

if ion_temperature='rescale-v':
       every 'nraise' steps the average temperature, computed from
       the last nraise steps, is rescaled to tempw

if ion_temperature='rescaling' and calculation='vc-md':
       every 'nraise' steps the instantaneous temperature
       is rescaled to tempw

if ion_temperature='berendsen':
       the "rise time" parameter is given in units of the time step:
       tau = nraise*dt, so dt/tau = 1/nraise

if ion_temperature='andersen':
       the "collision frequency" parameter is given as nu=1/tau
       defined above, so nu*dt = 1/nraise
""",
    'refold_pos': """
This keyword applies only in the case of molecular dynamics or
damped dynamics. If true the ions are refolded at each step into
the supercell.
""",
    'upscale': """
Max reduction factor for conv_thr during structural optimization
conv_thr is automatically reduced when the relaxation
approaches convergence so that forces are still accurate,
but conv_thr will not be reduced to less that
conv_thr / upscale.
""",
    'bfgs_ndim': """
Number of old forces and displacements vectors used in the
PULAY mixing of the residual vectors obtained on the basis
of the inverse hessian matrix given by the BFGS algorithm.
When bfgs_ndim = 1, the standard quasi-Newton BFGS method is
used.
(bfgs only)
""",
    'trust_radius_max': """
Maximum ionic displacement in the structural relaxation.
(bfgs only)
""",
    'trust_radius_min': """
Minimum ionic displacement in the structural relaxation
BFGS is reset when trust_radius < trust_radius_min.
(bfgs only)
""",
    'trust_radius_ini': """
Initial ionic displacement in the structural relaxation.
(bfgs only)
""",
    'w_1': """
Parameters used in line search based on the Wolfe conditions.
(bfgs only)
""",
    'w_2': """
Parameters used in line search based on the Wolfe conditions.
(bfgs only)
"""
}

electrons_doc = {
    'electron_maxstep': """
maximum number of iterations in a scf step
""",
    'scf_must_converge': """
If .false. do not stop molecular dynamics or ionic relaxation
when electron_maxstep is reached. Use with care.
""",
    'conv_thr': """
Convergence threshold for selfconsistency:
   estimated energy error < conv_thr
(note that conv_thr is extensive, like the total energy).
For non-self-consistent calculations, conv_thr is used
to set the default value of the threshold (ethr) for
iterative diagonalizazion: see diago_thr_init
""",
    'adaptive_thr': """
If .TRUE. this turns on the use of an adaptive conv_thr for
the inner scf loops when using EXX.

""", 
    'conv_thr_init': """
When adaptive_thr = .TRUE. the convergence threshold for
each scf cycle is given by:
max( conv_thr, conv_thr_multi * dexx )
""", 
'conv_thr_multi': """
When adaptive_thr = .TRUE. the convergence threshold for
each scf cycle is given by:
max( conv_thr, conv_thr_multi * dexx )
""", 
    'mixing_mode': """
'plain' :    charge density Broyden mixing

'TF' :       as above, with simple Thomas-Fermi screening
            (for highly homogeneous systems)

'local-TF':  as above, with local-density-dependent TF screening
             (for highly inhomogeneous systems)
""", 
    'mixing_beta': """
mixing factor for self-consistency
""", 
    'mixing_ndim': """
number of iterations used in mixing scheme.
If you are tight with memory, you may reduce it to 4 or so.
""", 
    'mixing_fixed_ns': """
For DFT+U : number of iterations with fixed ns ( ns is the
  atomic density appearing in the Hubbard term ).
""", 
    'diagonalization': """
'david':  Davidson iterative diagonalization with overlap matrix
          (default). Fast, may in some rare cases fail.

'cg' :    conjugate-gradient-like band-by-band diagonalization
          Typically slower than 'david' but it uses less memory
          and is more robust (it seldom fails)

'cg-serial', 'david-serial': obsolete, use "-ndiag 1 instead"
          The subspace diagonalization in Davidson is performed
          by a fully distributed-memory parallel algorithm on
          4 or more processors, by default. The allocated memory
          scales down with the number of procs. Procs involved
          in diagonalization can be changed with command-line
          option "-ndiag N". On multicore CPUs it is often
          convenient to let just one core per CPU to work
          on linear algebra.
""", 
    'ortho_para': """
OBSOLETE: use command-line option " -ndiag XX" instead 
""", 
    'diago_thr_init': """
Convergence threshold (ethr) for iterative diagonalization
(the check is on eigenvalue convergence).
For scf calculations: default is 1.D-2 if starting from a
superposition of atomic orbitals; 1.D-5 if starting from a
charge density. During self consistency the threshold
is automatically reduced (but never below 1.D-13) when
approaching convergence.
For non-scf calculations: default is (conv_thr/N elec)/10.
""", 
    'diago_cg_maxiter': """
For conjugate gradient diagonalization:
max number of iterations
""", 
    'diago_david_ndim': """
For Davidson diagonalization: dimension of workspace
(number of wavefunction packets, at least 2 needed).
A larger value may yield a somewhat faster algorithm
but uses more memory. The opposite holds for smaller values.
Try diago_david_ndim=2 if you are tight on memory or if
your job is large: the speed penalty is often negligible
""", 
    'diago_full_acc': """
If .TRUE. all the empty states are diagonalized at the same level
of accuracy of the occupied ones. Otherwise the empty states are
diagonalized using a larger threshold (this should not affect
total energy, forces, and other ground-state properties).
""", 
    'efield': """
Amplitude of the finite electric field (in Ry a.u.;
1 a.u. = 36.3609*10^10 V/m). Used only if lelfield=.TRUE.
and if k-points (K_POINTS card) are not automatic.
""", 
    'efield_cart': """
Finite electric field (in Ry a.u.=36.3609*10^10 V/m) in
cartesian axis. Used only if lelfield=.TRUE. and if
k-points (K_POINTS card) are automatic.
""", 
    'startingpot': """
'atomic': starting potential from atomic charge superposition
          ( default for scf, *relax, *md )

'file'  : start from existing "charge-density.xml" file in the
          directory specified by variables "prefix" and "outdir"
          For nscf and bands calculation this is the default
          and the only sensible possibility.
""", 
    'startingwfc': """
'atomic': start from superposition of atomic orbitals
          If not enough atomic orbitals are available,
          fill with random numbers the remaining wfcs
          The scf typically starts better with this option,
          but in some high-symmetry cases one can "loose"
          valence states, ending up in the wrong ground state.

'atomic+random': as above, plus a superimposed "randomization"
          of atomic orbitals. Prevents the "loss" of states
          mentioned above.

'random': start from random wfcs. Slower start of scf but safe.
          It may also reduce memory usage in conjunction with
          diagonalization='cg'

'file':   start from an existing wavefunction file in the
          directory specified by variables "prefix" and "outdir"
""",
    'tqr': """
If .true., use the real-space algorithm for augmentation
charges in ultrasoft pseudopotentials.
Must faster execution of ultrasoft-related calculations,
but numerically less accurate than the default algorithm.
Use with care and after testing!
"""
}

system_doc = {
    'ibrav': """
 Bravais-lattice index. If ibrav /= 0, specify EITHER
  [ celldm(1)-celldm(6) ] OR [ A,B,C,cosAB,cosAC,cosBC ]
  but NOT both. The lattice parameter "alat" is set to
  alat = celldm(1) (in a.u.) or alat = A (in Angstrom);
  see below for the other parameters.
  For ibrav=0 specify the lattice vectors in CELL_PARAMETER,
  optionally the lattice parameter alat = celldm(1) (in a.u.)
  or = A (in Angstrom), or else it is taken from CELL_PARAMETERS

ibrav      structure                   celldm(2)-celldm(6)
                                     or: b,c,cosab,cosac,cosbc
  0          free
      crystal axis provided in input: see card CELL_PARAMETERS

  1          cubic P (sc)
      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)

  2          cubic F (fcc)
      v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)

  3          cubic I (bcc)
      v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)

  4          Hexagonal and Trigonal P        celldm(3)=c/a
      v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)

  5          Trigonal R, 3fold axis c        celldm(4)=cos(alpha)
      The crystallographic vectors form a three-fold star around
      the z-axis, the primitive cell is a simple rhombohedron:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
      where c=cos(alpha) is the cosine of the angle alpha between
      any pair of crystallographic vectors, tx, ty, tz are:
        tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
 -5          Trigonal R, 3fold axis <111>    celldm(4)=cos(alpha)
      The crystallographic vectors form a three-fold star around
      <111>. Defining a' = a/sqrt(3) :
      v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
      where u and v are defined as
        u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
      and tx, ty, tz as for case ibrav=5
      Note: if you prefer x,y,z as axis in the cubic limit,
            set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
            See also the note in flib/latgen.f90

  6          Tetragonal P (st)               celldm(3)=c/a
      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)

  7          Tetragonal I (bct)              celldm(3)=c/a
      v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)

  8          Orthorhombic P                  celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)

  9          Orthorhombic base-centered(bco) celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
 -9          as 9, alternate description
      v1 = (a/2,-b/2,0),  v2 = (a/2,-b/2,0),  v3 = (0,0,c)

 10          Orthorhombic face-centered      celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)

 11          Orthorhombic body-centered      celldm(2)=b/a
                                             celldm(3)=c/a
      v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)

 12          Monoclinic P, unique axis c     celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(4)=cos(ab)
      v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
      where gamma is the angle between axis a and b.
-12          Monoclinic P, unique axis b     celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(5)=cos(ac)
      v1 = (a,0,0), v2 = (0,b,0), v3 = (c*cos(beta),0,c*sin(beta))
      where beta is the angle between axis a and c

 13          Monoclinic base-centered        celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(4)=cos(ab)
      v1 = (  a/2,         0,                -c/2),
      v2 = (b*cos(gamma), b*sin(gamma), 0),
      v3 = (  a/2,         0,                  c/2),
      where gamma is the angle between axis a and b

 14          Triclinic                       celldm(2)= b/a,
                                             celldm(3)= c/a,
                                             celldm(4)= cos(bc),
                                             celldm(5)= cos(ac),
                                             celldm(6)= cos(ab)
      v1 = (a, 0, 0),
      v2 = (b*cos(gamma), b*sin(gamma), 0)
      v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
           c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
                     - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
  where alpha is the angle between axis b and c
         beta is the angle between axis a and c
        gamma is the angle between axis a and b
""",
    'celldm': """
Crystallographic constants - see the "ibrav" variable.
Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.
Only needed values (depending on "ibrav") must be specified
alat = celldm(1) is the lattice parameter "a" (in BOHR)
If ibrav=0, only celldm(1) is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'A': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'B': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'C': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'cosAB': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'cosAC': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'cosBC': """
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
""",
    'nat': """
number of atoms in the unit cell
""",
    'ntyp': """
number of types of atoms in the unit cell
""",
    'nbnd': """
number of electronic states (bands) to be calculated.
Note that in spin-polarized calculations the number of
k-point, not the number of bands per k-point, is doubled
""",
    'tot_charge': """
total charge of the system. Useful for simulations with charged cells.
By default the unit cell is assumed to be neutral (tot_charge=0).
tot_charge=+1 means one electron missing from the system,
tot_charge=-1 means one additional electron, and so on.

In a periodic calculation a compensating jellium background is
inserted to remove divergences if the cell is not neutral.
""",
    'tot_magnetization': """
total majority spin charge - minority spin charge.
Used to impose a specific total electronic magnetization.
If unspecified then tot_magnetization variable is ignored and
the amount of electronic magnetization is determined during
the self-consistent cycle.
""",
    'starting_magnetization': """
starting spin polarization on atomic type 'i' in a spin
polarized calculation. Values range between -1 (all spins
down for the valence electrons of atom type 'i') to 1
(all spins up). Breaks the symmetry and provides a starting
point for self-consistency. The default value is zero, BUT a
value MUST be specified for AT LEAST one atomic type in spin
polarized calculations, unless you constrain the magnetization
(see "tot_magnetization" and "constrained_magnetization").
Note that if you start from zero initial magnetization, you
will invariably end up in a nonmagnetic (zero magnetization)
state. If you want to start from an antiferromagnetic state,
you may need to define two different atomic species
corresponding to sublattices of the same atomic type.
starting_magnetization is ignored if you are performing a
non-scf calculation, if you are restarting from a previous
run, or restarting from an interrupted run.
If you fix the magnetization with "tot_magnetization",
you should not specify starting_magnetization.
In the spin-orbit case starting with zero
starting_magnetization on all atoms imposes time reversal
symmetry. The magnetization is never calculated and
kept zero (the internal variable domag is .FALSE.).
""",
    'ecutwfc': """
kinetic energy cutoff (Ry) for wavefunctions
""",
    'ecutrho': """
kinetic energy cutoff (Ry) for charge density and potential
For norm-conserving pseudopotential you should stick to the
default value, you can reduce it by a little but it will
introduce noise especially on forces and stress.
If there are ultrasoft PP, a larger value than the default is
often desirable (ecutrho = 8 to 12 times ecutwfc, typically).
PAW datasets can often be used at 4*ecutwfc, but it depends
on the shape of augmentation charge: testing is mandatory.
The use of gradient-corrected functional, especially in cells
with vacuum, or for pseudopotential without non-linear core
correction, usually requires an higher values of ecutrho
to be accurately converged.
""",
    'ecutfock': """
kinetic energy cutoff (Ry) for the exact exchange operator in
EXX type calculations. By default this is the same as ecutrho
but in some EXX calculations significant speed-up can be found
by reducing ecutfock, at the expense of some loss in accuracy.
Currently only implemented for the optimized gamma point only
calculations.
""",
    'nr1': """
three-dimensional FFT mesh (hard grid) for charge
density (and scf potential). If not specified
the grid is calculated based on the cutoff for
charge density (see also "ecutrho")
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nr2': """
three-dimensional FFT mesh (hard grid) for charge
density (and scf potential). If not specified
the grid is calculated based on the cutoff for
charge density (see also "ecutrho")
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nr3': """
three-dimensional FFT mesh (hard grid) for charge
density (and scf potential). If not specified
the grid is calculated based on the cutoff for
charge density (see also "ecutrho")
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nr1s': """
three-dimensional mesh for wavefunction FFT and for the smooth
part of charge density ( smooth grid ).
Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nr2s': """
three-dimensional mesh for wavefunction FFT and for the smooth
part of charge density ( smooth grid ).
Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nr3s': """
three-dimensional mesh for wavefunction FFT and for the smooth
part of charge density ( smooth grid ).
Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
Note: you must specify all three dimensions for this setting to
be used.
""",
    'nosym': """
if (.TRUE.) symmetry is not used. Note that
- if the k-point grid is provided in input, it is used "as is"
  and symmetry-inequivalent k-points are not generated;
- if the k-point grid is automatically generated, it will
  contain only points in the irreducible BZ for the bravais
  lattice, irrespective of the actual crystal symmetry.
A careful usage of this option can be advantageous
- in low-symmetry large cells, if you cannot afford a k-point
  grid with the correct symmetry
- in MD simulations
- in calculations for isolated atoms
""",
    'nosym_evc': """
if(.TRUE.) symmetry is not used but the k-points are
forced to have the symmetry of the Bravais lattice;
an automatically generated k-point grid will contain
all the k-points of the grid and the points rotated by
the symmetries of the Bravais lattice which are not in the
original grid. If available, time reversal is
used to reduce the k-points (and the q => -q symmetry
is used in the phonon code). To disable also this symmetry set
noinv=.TRUE..
""",
    'noinv': """
if (.TRUE.) disable the usage of k => -k symmetry
(time reversal) in k-point generation
""",
    'no_t_rev': """
if (.TRUE.) disable the usage of magnetic symmetry operations
that consist in a rotation + time reversal.
""",
    'force_symmorphic': """
if (.TRUE.) force the symmetry group to be symmorphic by disabling
symmetry operations having an associated fractionary translation
""",
    'use_all_frac': """
if (.TRUE.) do not discard symmetry operations with an
associated fractionary translation that does not send the
real-space FFT grid into itself. These operations are
incompatible with real-space symmetrization but not with the
new G-space symmetrization. BEWARE: do not use for phonons!
The phonon code still uses real-space symmetrization.
""",
    'occupations': """
'smearing':     gaussian smearing for metals
                see variables 'smearing' and 'degauss'

'tetrahedra' :  especially suited for calculation of DOS
                (see P.E. Bloechl, PRB49, 16223 (1994))
                Requires uniform grid of k-points,
                automatically generated (see below)
                Not suitable (because not variational) for
                force/optimization/dynamics calculations

'fixed' :       for insulators with a gap

'from_input' :  The occupation are read from input file,
                card OCCUPATIONS. Option valid only for a
                single k-point, requires "nbnd" to be set
                in input. Occupations should be consistent
                with the value of "tot_charge".
""",
    'one_atom_occupations': """
This flag is used for isolated atoms (nat=1) together with
occupations='from_input'. If it is .TRUE., the wavefunctions
are ordered as the atomic starting wavefunctions, independently
from their eigenvalue. The occupations indicate which atomic
states are filled.
The order of the states is written inside the UPF
pseudopotential file.
In the scalar relativistic case:
S -> l=0, m=0
P -> l=1, z, x, y
D -> l=2, r^2-3z^2, xz, yz, xy, x^2-y^2
In the noncollinear magnetic case (with or without spin-orbit),
each group of states is doubled. For instance:
P -> l=1, z, x, y for spin up, l=1, z, x, y for spin down.
Up and down is relative to the direction of the starting
magnetization.
In the case with spin-orbit and time-reversal
(starting_magnetization=0.0) the atomic wavefunctions are
radial functions multiplied by spin-angle functions.
For instance:
P -> l=1, j=1/2, m_j=-1/2,1/2. l=1, j=3/2,
     m_j=-3/2, -1/2, 1/2, 3/2.
In the magnetic case with spin-orbit the atomic wavefunctions
can be forced to be spin-angle functions by setting
starting_spin_angle to .TRUE..
""",
    'starting_spin_angle': """
In the spin-orbit case when domag=.TRUE., by default,
the starting wavefunctions are initialized as in scalar
relativistic noncollinear case without spin-orbit.
By setting starting_spin_angle=.TRUE. this behaviour can
be changed and the initial wavefunctions are radial
functions multiplied by spin-angle functions.
When domag=.FALSE. the initial wavefunctions are always
radial functions multiplied by spin-angle functions
independently from this flag.
When lspinorb is .FALSE. this flag is not used.
""",
    'degauss': """
value of the gaussian spreading (Ry) for brillouin-zone
integration in metals.
""",
    'smearing': """
'gaussian', 'gauss':
    ordinary Gaussian spreading (Default)

'methfessel-paxton', 'm-p', 'mp':
    Methfessel-Paxton first-order spreading
    (see PRB 40, 3616 (1989)).

'marzari-vanderbilt', 'cold', 'm-v', 'mv':
    Marzari-Vanderbilt cold smearing
    (see PRL 82, 3296 (1999))

'fermi-dirac', 'f-d', 'fd':
    smearing with Fermi-Dirac function
""",
    'nspin': """
nspin = 1 :  non-polarized calculation (default)

nspin = 2 :  spin-polarized calculation, LSDA
             (magnetization along z axis)

nspin = 4 :  spin-polarized calculation, noncollinear
             (magnetization in generic direction)
             DO NOT specify nspin in this case;
             specify "noncolin=.TRUE." instead
""",
    'noncolin': """
if .true. the program will perform a noncollinear calculation.
""",
    'ecfixed': """
ecfixed, qcutz, q2sigma:  parameters for modified functional to be
used in variable-cell molecular dynamics (or in stress calculation).
"ecfixed" is the value (in Rydberg) of the constant-cutoff;
"qcutz" and "q2sigma" are the height and the width (in Rydberg)
of the energy step for reciprocal vectors whose square modulus
is greater than "ecfixed". In the kinetic energy, G^2 is
replaced by G^2 + qcutz * (1 + erf ( (G^2 - ecfixed)/q2sigma) )
See: M. Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995)
""",
    'qcutz': """
ecfixed, qcutz, q2sigma:  parameters for modified functional to be
used in variable-cell molecular dynamics (or in stress calculation).
"ecfixed" is the value (in Rydberg) of the constant-cutoff;
"qcutz" and "q2sigma" are the height and the width (in Rydberg)
of the energy step for reciprocal vectors whose square modulus
is greater than "ecfixed". In the kinetic energy, G^2 is
replaced by G^2 + qcutz * (1 + erf ( (G^2 - ecfixed)/q2sigma) )
See: M. Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995)
""",
    'q2sigma': """
ecfixed, qcutz, q2sigma:  parameters for modified functional to be
used in variable-cell molecular dynamics (or in stress calculation).
"ecfixed" is the value (in Rydberg) of the constant-cutoff;
"qcutz" and "q2sigma" are the height and the width (in Rydberg)
of the energy step for reciprocal vectors whose square modulus
is greater than "ecfixed". In the kinetic energy, G^2 is
replaced by G^2 + qcutz * (1 + erf ( (G^2 - ecfixed)/q2sigma) )
See: M. Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995)
""",
    'input_dft': """
Exchange-correlation functional: eg 'PBE', 'BLYP' etc
See Modules/functionals.f90 for allowed values.
Overrides the value read from pseudopotential files.
Use with care and if you know what you are doing!
""",
    'exx_fraction': """
Fraction of EXX for hybrid functional calculations. In the case of
input_dft='PBE0', the default value is 0.25, while for input_dft='B3LYP'
the exx_fraction default value is 0.20.
""",
    'screening_parameter': """
screening_parameter for HSE like hybrid functionals.
See J. Chem. Phys. 118, 8207 (2003)
and J. Chem. Phys. 124, 219906 (2006) for more informations.
""",
    'exxdiv_treatment': """
Specific for EXX. It selects the kind of approach to be used
for treating the Coulomb potential divergencies at small q vectors.

gygi-baldereschi : appropriate for cubic and quasi-cubic supercells
vcut_spherical : appropriate for cubic and quasi-cubic supercells
vcut_ws : appropriate for strongly anisotropic supercells, see also
          ecutvcut.
none : sets Coulomb potential at G,q=0 to 0.0 (required for GAU-PBE)
""",
    'x_gamma_extrapolation': """
Specific for EXX. If true, extrapolate the G=0 term of the
potential (see README in examples/EXX_example for more)
Set this to .false. for GAU-PBE.
""",
    'ecutvcut': """
Reciprocal space cutoff for correcting
Coulomb potential divergencies at small q vectors.
""",
    'nqx1': """
three-dimensional mesh for q (k1-k2) sampling of
the Fock operator (EXX). Can be smaller than
the number of k-points.

Currently this defaults to the size of the k-point mesh used.
 In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.
""",
    'nqx2': """
three-dimensional mesh for q (k1-k2) sampling of
the Fock operator (EXX). Can be smaller than
the number of k-points.

Currently this defaults to the size of the k-point mesh used.
 In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.
""",
    'nqx3': """
three-dimensional mesh for q (k1-k2) sampling of
the Fock operator (EXX). Can be smaller than
the number of k-points.

Currently this defaults to the size of the k-point mesh used.
 In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.
""",
    'lda_plus_u': """
Specify lda_plus_u = .TRUE. to enable DFT+U calculations
See: Anisimov, Zaanen, and Andersen, PRB 44, 943 (1991);
     Anisimov et al., PRB 48, 16929 (1993);
     Liechtenstein, Anisimov, and Zaanen, PRB 52, R5467 (1994).
You must specify, for each species with a U term, the value of
U and (optionally) alpha, J of the Hubbard model (all in eV):
see lda_plus_u_kind, Hubbard_U, Hubbard_alpha, Hubbard_J
""",
    'lda_plus_u_kind': """
Specifies the type of DFT+U calculation:
                  0   simplified version of Cococcioni and de Gironcoli,
                      PRB 71, 035105 (2005), using Hubbard_U
                  1   rotationally invariant scheme of Liechtenstein et al.,
                      using Hubbard_U and Hubbard_J
""",
    'Hubbard_U': """
Hubbard_U(i): U parameter (eV) for species i, DFT+U calculation
""",
    'Hubbard_J0': """
Hubbard_J0(i): J0 parameter (eV) for species i, DFT+U+J calculation,
see PRB 84, 115108 (2011) for details.
""",
    'Hubbard_alpha': """
Hubbard_alpha(i) is the perturbation (on atom i, in eV)
used to compute U with the linear-response method of
Cococcioni and de Gironcoli, PRB 71, 35105 (2005)
(only for lda_plus_u_kind=0)
""",
    'Hubbard_beta': """
Hubbard_beta(i) is the perturbation (on atom i, in eV)
used to compute J0 with the linear-response method of
Cococcioni and de Gironcoli, PRB 71, 35105 (2005)
(only for lda_plus_u_kind=0). See also
PRB 84, 115108 (2011).
""",
    'Hubbard_J': """
Hubbard_J(i,ityp): J parameters (eV) for species ityp,
used in DFT+U calculations (only for lda_plus_u_kind=1)
For p orbitals:  J = Hubbard_J(1,ityp);
For d orbitals:  J = Hubbard_J(1,ityp), B = Hubbard_J(2,ityp);
For f orbitals:  J = Hubbard_J(1,ityp), E2 = Hubbard_J(2,ityp),
                 E3= Hubbard_J(3,ityp).
If B or E2 or E3 are not specified or set to 0 they will be
calculated from J using atomic ratios.
""",
    'starting_ns_eigenvalue': """
In the first iteration of an DFT+U run it overwrites
the m-th eigenvalue of the ns occupation matrix for the
ispin component of atomic species I. Leave unchanged
eigenvalues that are not set. This is useful to suggest
the desired orbital occupations when the default choice
takes another path.
""",
    'U_projection_type': """
Only active when lda_plus_U is .true., specifies the type
of projector on localized orbital to be used in the DFT+U
scheme.

Currently available choices:
'atomic': use atomic wfc's (as they are) to build the projector

'ortho-atomic': use Lowdin orthogonalized atomic wfc's

'norm-atomic':  Lowdin normalization of atomic wfc. Keep in mind:
                atomic wfc are not orthogonalized in this case.
                This is a "quick and dirty" trick to be used when
                atomic wfc from the pseudopotential are not
                normalized (and thus produce occupation whose
                value exceeds unity). If orthogonalized wfc are
                not needed always try 'atomic' first.

'file':         use the information from file "prefix".atwfc that must
                have been generated previously, for instance by pmw.x
                (see PP/src/poormanwannier.f90 for details).

'pseudo':       use the pseudopotential projectors. The charge density
                outside the atomic core radii is excluded.
                N.B.: for atoms with +U, a pseudopotential with the
                all-electron atomic wavefunctions is required (i.e.,
                as generated by ld1.x with lsave_wfc flag).

NB: forces and stress currently implemented only for the
'atomic' and 'pseudo' choice.
""",
    'edir': """
The direction of the electric field or dipole correction is
parallel to the bg(:,edir) reciprocal lattice vector, so the
potential is constant in planes defined by FFT grid points;
edir = 1, 2 or 3. Used only if tefield is .TRUE.
""",
    'emaxpos': """
Position of the maximum of the saw-like potential along crystal
axis "edir", within the  unit cell (see below), 0 < emaxpos < 1
Used only if tefield is .TRUE.
""",
    'eopreg': """
Zone in the unit cell where the saw-like potential decreases.
( see below, 0 < eopreg < 1 ). Used only if tefield is .TRUE.
""",
    'eamp': """
Amplitude of the electric field, in ***Hartree*** a.u.;
1 a.u. = 51.4220632*10^10 V/m). Used only if tefield=.TRUE.
The saw-like potential increases with slope "eamp" in the
region from (emaxpos+eopreg-1) to (emaxpos), then decreases
to 0 until (emaxpos+eopreg), in units of the crystal
vector "edir". Important: the change of slope of this
potential must be located in the empty region, or else
unphysical forces will result.
""",
    'angle1': """
The angle expressed in degrees between the initial
magnetization and the z-axis. For noncollinear calculations
only; index i runs over the atom types.
""",
    'angle2': """
The angle expressed in degrees between the projection
of the initial magnetization on x-y plane and the x-axis.
For noncollinear calculations only.
""",
    'constrained_magnetization': """
Used to perform constrained calculations in magnetic systems.
Currently available choices:

'none':
         no constraint

'total':
         total magnetization is constrained by
         adding a penalty functional to the total energy:

         LAMBDA * SUM_{i} ( magnetization(i) - fixed_magnetization(i) )**2

         where the sum over i runs over the three components of
         the magnetization. Lambda is a real number (see below).
         Noncolinear case only. Use "tot_magnetization" for LSDA

'atomic':
         atomic magnetization are constrained to the defined
         starting magnetization adding a penalty:

         LAMBDA * SUM_{i,itype} ( magnetic_moment(i,itype) - mcons(i,itype) )**2

         where i runs over the cartesian components (or just z
         in the collinear case) and itype over the types (1-ntype).
         mcons(:,:) array is defined from starting_magnetization,
         (and angle1, angle2 in the non-collinear case). lambda is
         a real number

'total direction':
          the angle theta of the total magnetization
          with the z axis (theta = fixed_magnetization(3))
          is constrained:

          LAMBDA * ( arccos(magnetization(3)/mag_tot) - theta )**2

          where mag_tot is the modulus of the total magnetization.

'atomic direction':
          not all the components of the atomic
          magnetic moment are constrained but only the cosine
          of angle1, and the penalty functional is:

          LAMBDA * SUM_{itype} ( mag_mom(3,itype)/mag_mom_tot - cos(angle1(ityp)) )**2

N.B.: symmetrization may prevent to reach the desired orientation
      of the magnetization. Try not to start with very highly symmetric
      configurations or use the nosym flag (only as a last remedy)
""",
    'fixed_magnetization': """
total magnetization vector (x,y,z components) to be kept
fixed when constrained_magnetization='total'
""",
    'lambda': """
parameter used for constrained_magnetization calculations
N.B.: if the scf calculation does not converge, try to reduce lambda
      to obtain convergence, then restart the run with a larger lambda
""",
    'report': """
It is the number of iterations after which the program
write all the atomic magnetic moments.
""",
    'lspinorb': """
if .TRUE. the noncollinear code can use a pseudopotential with
spin-orbit.
""",
    'assume_isolated': """
Used to perform calculation assuming the system to be
isolated (a molecule or a cluster in a 3D supercell).

Currently available choices:

'none' (default): regular periodic calculation w/o any correction.

'makov-payne', 'm-p', 'mp' : the Makov-Payne correction to the
         total energy is computed. An estimate of the vacuum
         level is also calculated so that eigenvalues can be
         properly aligned. ONLY FOR CUBIC SYSTEMS (ibrav=1,2,3)
         Theory:
         G.Makov, and M.C.Payne,
         "Periodic boundary conditions in ab initio
         calculations" , Phys.Rev.B 51, 4014 (1995)

'martyna-tuckerman', 'm-t', 'mt' : Martyna-Tuckerman correction
         to both total energy and scf potential. Adapted from:
         G.J. Martyna, and M.E. Tuckerman,
         "A reciprocal space based method for treating long
         range interactions in ab-initio and force-field-based
         calculation in clusters", J.Chem.Phys. 110, 2810 (1999)

'esm' :  Effective Screening Medium Method.
         For polarized or charged slab calculation, embeds
         the simulation cell within an effective semi-
         infinite medium in the perpendicular direction
         (along z). Embedding regions can be vacuum or
         semi-infinite metal electrodes (use 'esm_bc' to
         choose boundary conditions). If between two
         electrodes, an optional electric field
         ('esm_efield') may be applied. Method described in
         M. Otani and O. Sugino, "First-principles
         calculations of charged surfaces and interfaces:
         A plane-wave nonrepeated slab approach," PRB 73,
         115407 (2006).
         NB: Requires cell with a_3 lattice vector along z,
         normal to the xy plane, with the slab centered
         around z=0. Also requires symmetry checking to be
         disabled along z, either by setting 'nosym' = .TRUE.
         or by very slight displacement (i.e., 5e-4 a.u.)
         of the slab along z.
         See 'esm_bc', 'esm_efield', 'esm_w', 'esm_nfit'.
""",
    'esm_bc': """
If assume_isolated = 'esm', determines the boundary
conditions used for either side of the slab.

Currently available choices:

'pbc' (default): regular periodic calculation (no ESM).

'bc1' : Vacuum-slab-vacuum (open boundary conditions)

'bc2' : Metal-slab-metal (dual electrode configuration).
        See also 'esm_efield'.

'bc3' : Vacuum-slab-metal
""",
    'esm_w': """
If assume_isolated = 'esm', determines the position offset
[in a.u.] of the start of the effective screening region,
measured relative to the cell edge. (ESM region begins at
z = +/- [L_z/2 + esm_w] ).
""",
    'esm_efield': """
If assume_isolated = 'esm' and esm_bc = 'bc2', gives the
magnitude of the electric field [Ry/a.u.] to be applied
between semi-infinite ESM electrodes.
""",
    'esm_nfit': """
If assume_isolated = 'esm', gives the number of z-grid points
for the polynomial fit along the cell edge.
""",
    'vdw_corr': """
Type of Van der Waals correction. Allowed values:

   'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d': semiempirical Grimme's DFT-D2.
    Optional variables: "london_s6", "london_rcut"
    S. Grimme, J. Comp. Chem. 27, 1787 (2006),
    V. Barone et al., J. Comp. Chem. 30, 934 (2009).

    'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler': Tkatchenko-Scheffler
     dispersion corrections with first-principle derived C6 coefficients
     (implemented in CP only). Optional variables: "ts_vdw_econv_thr", "ts_vdw_isolated"
     See A. Tkatchenko and M. Scheffler, Phys. Rev. Lett. 102, 073005 (2009)

    'XDM', 'xdm': Exchange-hole dipole-moment model. Optional variables: "xdm_a1", "xdm_a2"
     A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007)
         A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 136, 174109 (2012)

Note that non-local functionals (eg vdw-DF) are NOT specified here but in "input_dft"
""",
    'london': """
OBSOLESCENT, same as vdw_corr='DFT-D'
""",
    'london_s6': """
global scaling parameter for DFT-D. Default is good for PBE.
""",
    'london_rcut': """
cutoff radius (a.u.) for dispersion interactions
""",
    'xdm': """
OBSOLESCENT, same as vdw_corr='xdm'
""",
    'xdm_a1': """
Damping function parameter a1 (adimensional). This value should change
with the exchange-correlation functional. The default corresponds to
PW86PBE.
For other functionals, see:
   http://gatsby.ucmerced.edu/wiki/XDM_damping_function_parameters
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013)
""",
    'xdm_a2': """
Damping function parameter a2 (angstrom). This value should change
with the exchange-correlation functional. The default corresponds to
PW86PBE.
For other functionals, see:
   http://gatsby.ucmerced.edu/wiki/XDM_damping_function_parameters
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013)
""",
    'space_group': """
The number of the space group of the crystal, as given
                in the International Tables of Crystallography A (ITA).
                This allows to give in input only the inequivalent atomic
                positions. The positions of all the symmetry equivalent atoms
                are calculated by the code. Used only when the atomic positions
                are of type crystal_sg.
""",
    'uniqueb': """
Used only for monoclinic lattices. If .TRUE. the b
                 unique ibrav (-12 or -13) are used, and symmetry
                 equivalent positions are chosen assuming that the
                 two fold axis or the mirror normal is parallel to the
                 b axis. If .FALSE. it is parallel to the c axis.
""",
    'origin_choice': """
Used only for space groups that in the ITA allow
                 the use of two different origins. origin_choice=1,
                 means the first origin, while origin_choice=2 is the
                 second origin.
""",
    'rhombohedral': """
Used only for rhombohedral space groups.
                 When .TRUE. the coordinates of the inequivalent atoms are
                 given with respect to the rhombohedral axes, when .FALSE.
                 the coordinates of the inequivalent atoms are given with
                 respect to the hexagonal axes. They are converted internally
                 to the rhombohedral axes and ibrav=5 is used in both cases.
"""
}


control_doc = {
    'calculation': """
A string describing the task to be performed
vc = variable cell
""",
    'title': """
reprinted on output
""",
    'verbosity': """
Currently two verbosity levels are implemented:
'high' and 'low'. 'debug' and 'medium' have the same
effect as 'high'; 'default' and 'minimal', as 'low'
""",
    'restart_mode': """
'from_scratch'  : from scratch. This is the normal way
                  to perform a PWscf calculation
'restart'       : from previous interrupted run. Use this
                  switch only if you want to continue an
                  interrupted calculation, not to start a
                  new one, or to perform non-scf calculations.
                  Works only if the calculation was cleanly
                  stopped using variable "max_seconds", or
                  by user request with an "exit file" (i.e.:
                  create a file "prefix".EXIT, in directory
                  "outdir"; see variables "prefix", "outdir")
""",
    'wf_collect': """
This flag controls the way wavefunctions are stored to disk :

True    collect wavefunctions from all processors, store them
        into the output data directory "outdir"/"prefix".save,
        one wavefunction per k-point in subdirs K000001/,
        K000001/, etc.. Use this if you want wavefunctions
        to be readable on a different number of processors.

False   do not collect wavefunctions, leave them in temporary
        local files (one per processor). The resulting format
        will be readable only by jobs running on the same
        number of processors and pools. Requires less I/O
        than the previous case.

Note that this flag has no effect on reading, only on writing.
""",
    'nstep': """
number of ionic + electronic steps
""",
    'iprint': """
band energies are written every "iprint" iterations
""",
    'tstress': """
calculate stress. It is set to True automatically if
calculation='vc-md' or 'vc-relax'
""",
    'tprnfor': """
calculate forces. It is set to .TRUE. automatically if
calculation='relax','md','vc-md'
""",
    'dt': """
time step for molecular dynamics, in Rydberg atomic units
(1 a.u.=4.8378 * 10^-17 s : beware, the CP code uses
 Hartree atomic units, half that much!!!)
""",
    'outdir': """
input, temporary, output files are found in this directory,
see also "wfcdir"
""",
    'wfcdir': """
this directory specifies where to store files generated by
each processor (*.wfc{N}, *.igk{N}, etc.). Useful for
machines without a parallel file system: set "wfcdir" to
a local file system, while "outdir" should be a parallel
or networkfile system, visible to all processors. Beware:
in order to restart from interrupted runs, or to perform
further calculations using the produced data files, you
may need to copy files to "outdir". Works only for pw.x.
""",
    'prefix': """
prepended to input/output filenames:
prefix.wfc, prefix.rho, etc.
""",
    'lkpoint_dir': """
If .false. a subdirectory for each k_point is not opened
in the "prefix".save directory; Kohn-Sham eigenvalues are
stored instead in a single file for all k-points. Currently
doesn't work together with "wf_collect"
""",
    'max_seconds': """
jobs stops after "max_seconds" CPU time. Use this option
in conjunction with option "restart_mode" if you need to
split a job too long to complete into shorter jobs that
fit into your batch queues. Default is 150 days.
""",
    'etot_conv_thr': """
convergence threshold on total energy (a.u) for ionic
minimization: the convergence criterion is satisfied
when the total energy changes less than "etot_conv_thr"
between two consecutive scf steps. Note that "etot_conv_thr"
is extensive, like the total energy.
See also "forc_conv_thr" - both criteria must be satisfied
""",
    'forc_conv_thr': """
convergence threshold on forces (a.u) for ionic minimization:
the convergence criterion is satisfied when all components of
all forces are smaller than "forc_conv_thr".
See also "etot_conv_thr" - both criteria must be satisfied
""",
    'disk_io': """
Specifies the amount of disk I/O activity
'high':   save all data to disk at each SCF step

'medium': save wavefunctions at each SCF step unless
          there is a single k-point per process (in which
          case the behavior is the same as 'low')

'low' :   store wfc in memory, save only at the end

'none':   do not save anything, not even at the end
          ('scf', 'nscf', 'bands' calculations; some data
           may be written anyway for other calculations)

Default is 'low' for the scf case, 'medium' otherwise.
Note that the needed RAM increases as disk I/O decreases!
It is no longer needed to specify 'high' in order to be able
to restart from an interrupted calculation (see "restart_mode")
but you cannot restart in disk_io='none'
""",
    'pseudo_dir': """
directory containing pseudopotential files
""",
    'tefield': """
If True a saw-like potential simulating an electric field
is added to the bare ionic potential. See variables "edir",
"eamp", "emaxpos", "eopreg" for the form and size of
the added potential.
""",
    'dipfield': """
If True and tefield=True a dipole correction is also
added to the bare ionic potential - implements the recipe
of L. Bengtsson, PRB 59, 12301 (1999). See variables "edir",
"emaxpos", "eopreg" for the form of the correction. Must
be used ONLY in a slab geometry, for surface calculations,
with the discontinuity FALLING IN THE EMPTY SPACE.
""",
    'lelfield': """
If .TRUE. a homogeneous finite electric field described
through the modern theory of the polarization is applied.
This is different from "tefield=.true." !
""",
    'nberrycyc': """
In the case of a finite electric field  ( lelfield == .TRUE. )
it defines the number of iterations for converging the
wavefunctions in the electric field Hamiltonian, for each
external iteration on the charge density
""",
    'lorbm': """
If .TRUE. perform orbital magnetization calculation.
If finite electric field is applied (lelfield=.true.)
only Kubo terms are computed
[for details see New J. Phys. 12, 053032 (2010)].
The type of calculation is 'nscf' and should be performed
on an automatically generated uniform grid of k points.
Works ONLY with norm-conserving pseudopotentials.
""",
    'lberry': """
If .TRUE. perform a Berry phase calculation
See the header of PW/src/bp_c_phase.f90 for documentation
""",
    'gdir': """
For Berry phase calculation: direction of the k-point
strings in reciprocal space. Allowed values: 1, 2, 3
1=first, 2=second, 3=third reciprocal lattice vector
For calculations with finite electric fields
(lelfield==.true.) "gdir" is the direction of the field
""",
    'nppstr': """
For Berry phase calculation: number of k-points to be
calculated along each symmetry-reduced string
The same for calculation with finite electric fields
(lelfield=.true.)
"""}

