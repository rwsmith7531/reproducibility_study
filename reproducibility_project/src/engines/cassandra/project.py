"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys
import flow
import foyer
import mbuild as mb
from flow import environments
from mbuild import box
import mosdef_cassandra as mc
import unyt as u


def get_molecule(job):
    from reproducibility_project.src.molecules.methane_ua import MethaneUA
    from reproducibility_project.src.molecules.pentane_ua import PentaneUA
    from mbuild.lib.molecules.water import WaterSPC

    molecule = job.sp.molecule

    molecule_dict = {
        "methaneUA": MethaneUA(),
        "pentaneUA": PentaneUA(),
        "benzeneUA": None,
        "waterSPC/E": WaterSPC(),
        "ethanolAA": None,
    }

    return molecule_dict[molecule]


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.operation
@Project.pre(lambda job: job.sp.engine == "cassandra")
def initialize(job):

    sys.path.append(Project().root_directory() + "/..")
    from reproducibility_project.src.utils.forcefields import load_ff

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )

    systems = construct_system(job.sp)
    #    ff = load_ff(job.sp.forcefield_name)

    #    for compound in compounds:
    #        if compound is None:
    #            continue
    #
    #        param_compound = ff.apply(compound)
    #
    if job.sp.N_vap is None:
        mols_to_add = [[job.sp.N_liquid]]
    else:
        mols_to_add = [[job.sp.N_liquid], [job.sp.N_vap]]

    ff = load_ff(job.sp.forcefield_name)
    molecule = get_molecule(job)
    param_molecule = ff.apply(molecule)

    box_liq_length = job.sp.box_L_liq * u.nm
    box_liq = mb.Box(lengths=3 * [box_liq_length.to_value("angstrom")])
    boxes = [box_liq]

    if job.sp.box_L_vap is not None:
        box_vap_length = job.sp.box_L_vap * u.nm
        box_vap = mb.Box(lengths=3 * [box_vap_length])
        boxes.append(box_vap)

    species_list = [param_molecule]

    mc_system = mc.System(boxes, species_list, mols_to_add=mols_to_add)

    ensemble = job.sp.ensemble

    ensemble_mapping = {"GEMC-NVT": "gemc", "NPT": "npt", "NVT": "nvt"}
    ensemble = ensemble_mapping[ensemble]

    moveset = mc.MoveSet(ensemble, species_list)

    thermo_props = [
        "energy_total",
        "mass_density",
        "volume",
        "pressure",
        "nmols",
    ]

    if job.sp.cutoff_style == "hard":
        cutoff_style = "cut"

        # "rcut_min": 1.0 * u.angstrom,

    custom_args = {
        "run_name": "equil",
        "charge_style": "none",
        "cutoff_style": cutoff_style,
        "vdw_cutoff": (job.sp.r_cut * u.nm).to_value("angstrom") * u.angstrom,
        "units": "sweeps",
        "steps_per_sweep": job.sp.N_liquid,
        "coord_freq": 500,
        "prop_freq": 10,
        "properties": thermo_props,
    }

    if ensemble == "npt":
        # custom_args["pressure"] = (job.sp.pressure * u.kPa).to_value("bar")
        custom_args["pressure"] = job.sp.pressure * u.kPa
    #        custom_args["vdw_cutoff"] = len(boxes) * [
    #            (job.sp.r_cut * u.nm).to_value("angstrom") * u.angstrom
    #        ]
    # "vdw_cutoff": (job.sp.r_cut * u.nm).to_value("angstrom") * u.angstrom,
    mc.run(
        system=mc_system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=job.sp.temperature * u.K,
        **custom_args
    )

    import pdb

    pdb.set_trace()


# for idx, compound in enumerate(compounds):
#    if compound is not None:
#        import pdb; pdb.set_trace()


if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
