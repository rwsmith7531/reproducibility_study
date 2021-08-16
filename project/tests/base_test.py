import os

import foyer
import pytest
import signac

from project.src import xmls

from .utils import create_gsd


class BaseTest:
    @pytest.fixture
    def spceff(self, name="spce.xml"):
        abs_path = os.path.dirname(os.path.abspath(xmls.__file__))
        return foyer.Forcefield(forcefield_files=str(abs_path) + "/" + name)

    @pytest.fixture
    def job_gsdfile(self, tmp_job):
        filename = tmp_job.fn("trajectory.gsd")
        create_gsd(filename)
        return tmp_job

    @pytest.fixture
    def tmp_project(self):
        with signac.TemporaryProject(name="test") as p:
            return p

    @pytest.fixture
    def tmp_job(self, tmp_project, statepoint={"a": 0}):
        tmp_project.open_job(statepoint).init()
        for job in tmp_project:
            return job

    @pytest.fixture
    def trappe_ua(self):
        return foyer.forcefields.load_TRAPPE_UA()

    @pytest.fixture
    def gsdfile_random(self, tmp_path):
        filename = tmp_path / "traj_random.gsd"
        create_gsd(filename)
        return filename

    @pytest.fixture
    def gsdfile_xstal(self, tmp_path):
        filename = tmp_path / "traj_xstal.gsd"
        create_gsd(filename, system="xstal")
        return filename
