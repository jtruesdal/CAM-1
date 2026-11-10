"""
Implementation of the CAM single column test.

This is a CAM specific test:
Verifies that single column is working by checking that T and Q
only have round off differences.

(1) 3D
(2) scam

"""

import os
import re
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.test_status import *
from CIME.status import append_testlog


logger = logging.getLogger(__name__)

class SCT(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'default',
                                       run_one_description = '3D CAM run',
                                       run_two_description = 'SCAM run')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist = 'CAMIOP'")
        if self._case.get_value("CAM_DYCORE") == "se":
            append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scale_dry_air_mass     = 0.0D0")

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")
        self._case.set_value("BFBFLAG","TRUE")

        # A testmod may choose the SCAM column by setting PTS_LAT/PTS_LON.
        # A valid PTS_LAT/PTS_LON would put the 3D run into single-column
        # mode (scol_lat/scol_lon in nuopc.runconfig), so save the values
        # for case two and reset case one to global.
        pts_lat = float(self._case.get_value("PTS_LAT"))
        pts_lon = float(self._case.get_value("PTS_LON"))
        if pts_lat > -999. and pts_lon > -999.:
            with open(self._scm_column_file(), "w") as f:
                f.write("{!r} {!r}\n".format(pts_lat, pts_lon))
            self._case.set_value("PTS_LAT", -999.)
            self._case.set_value("PTS_LON", -999.)

    def _case_one_ncdata(self):
        """The initial file case one actually uses (from its generated atm_in), or None.

        The SCAM case must start from the same state as the 3-D run.  Its namelist
        defaults can differ from the 3-D case's (single-column configurations have
        their own ncdata defaults at some resolutions, e.g. ne3np4), so the test
        passes the 3-D case's file explicitly instead of relying on the defaults.
        """
        atm_in = os.path.join(self._case1.get_value("CASEROOT"), "CaseDocs", "atm_in")
        if not os.path.isfile(atm_in):
            return None
        with open(atm_in) as f:
            for line in f:
                m = re.match(r"\s*ncdata\s*=\s*['\"]([^'\"]+)['\"]", line)
                if m:
                    return m.group(1)
        return None

    def _scm_column_file(self):
        return os.path.join(self._case1.get_value("CASEROOT"), "sct_scm_column")

    def _scm_column(self, default_lat, default_lon):
        """SCAM column from the testmod (saved in case one) or the default."""
        if os.path.isfile(self._scm_column_file()):
            with open(self._scm_column_file()) as f:
                lat, lon = (float(x) for x in f.read().split())
            return lat, lon
        return default_lat, default_lon



    def _case_two_setup(self):
        case_name = self._case.get_value("CASE")
        RUN_STARTDATE = self._case1.get_value("RUN_STARTDATE")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "NDENS    = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "MFILT    = 1,10000,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "nhtfrq   = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "fincl2='T','Q','TDIFF','QDIFF','LANDFRAC'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "iopfile = '../"+case_name+".cam.h1i."+RUN_STARTDATE+"-00000.nc'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist = 'YEARLY'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_cambfb_mode                = .true.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_use_obs_uv         = .true.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_relaxation         = .false.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_use_3dfrc          = .true.")
        for comp in self._case.get_values("COMP_CLASSES"):
            self._case.set_value("NTASKS_{}".format(comp), 1)
            self._case.set_value("NTHRDS_{}".format(comp), 1)
            self._case.set_value("ROOTPE_{}".format(comp), 0)

        if self._case.get_value("COMP_INTERFACE") == "mct":
            self._case.set_value("PTS_MODE","TRUE")

        self._case.set_value("BFBFLAG","TRUE")

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS").replace('-camiop','')
        self._case.set_value("CAM_CONFIG_OPTS","{} -scam camfrc ".format(CAM_CONFIG_OPTS))
        if self._case.get_value("CAM_DYCORE") == "se":
            lat, lon = self._scm_column(44.80320177421346, 276.7082039324993)
            self._case.set_value("PTS_LAT",lat)
            self._case.set_value("PTS_LON",lon)
            append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scale_dry_air_mass     = 0.0D0")
        else:
            append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "ncdata = '../"+case_name+".cam.i."+RUN_STARTDATE+"-00000.nc'")
            lat, lon = self._scm_column(-20.0, 140.0)
            self._case.set_value("PTS_LAT",lat)
            self._case.set_value("PTS_LON",lon)

        self._case.set_value("STOP_N",5)
        self._case.case_setup(test_mode=True, reset=True)

    def _component_compare_test(self, suffix1, suffix2,
                                success_change=False,
                                ignore_fieldlist_diffs=False):
        with self._test_status:
            stat,netcdf_filename,err=run_cmd('ls ./run/case2run/*.cam.h1*0000.nc ')
            stat,DIFFs,err=run_cmd('ncdump -ff -p 9,17 -v QDIFF,TDIFF '+netcdf_filename+' | egrep //\.\*DIFF | sed s/^\ \*// | sed s/^0,/0.0,/ | sed s/^0\;/0.0\;/ | sed s/\[,\;\].\*// |  uniq')
            array_of_DIFFs=DIFFs.split("\n")
            answer=max([abs(float(x)) for x in array_of_DIFFs])
            comments="Checking QDIFF,TDIFF in SCAM run."
            append_testlog(comments, self._orig_caseroot)
            # Test for greater than round off changes.
            if answer < 1e-10:
                self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, self._run_one_suffix, self._run_two_suffix), TEST_PASS_STATUS)
                comments="QDIFF,TDIFF: PASS"
            else:
                self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, self._run_one_suffix, self._run_two_suffix), TEST_FAIL_STATUS)
                comments="QDIFF,TDIFF: Difference greater than round off."
            append_testlog(comments, self._orig_caseroot)

    def _case_two_custom_prerun_action(self):
        self._set_case_two_ncdata()
        self._case_two_derecho_select()

    def _set_case_two_ncdata(self):
        """SE: the SCAM case starts from case one's initial file.  Done just before case two
        runs: case one has run by then (its atm_in exists), and case two's namelists are
        generated afterwards, at the start of its run."""
        if self._case1.get_value("CAM_DYCORE") != "se":
            return
        ncdata = self._case_one_ncdata()
        if not ncdata:
            logger.warning("SCT: case one's ncdata not found; the SCAM case uses its own default")
            return
        caseroot = self._case2.get_value("CASEROOT")
        line = "ncdata = '"+ncdata+"'"
        user_nl = os.path.join(caseroot, "user_nl_cam")
        if os.path.isfile(user_nl) and line in open(user_nl).read():
            return
        append_to_user_nl_files(caseroot = caseroot, component = "cam", contents = line)

    def _case_two_derecho_select(self):
        """ On NCAR derecho system the mpibind script causes ESMF in the second job to think it is using 128 tasks when it should only use 1
        changing the env variable PBS_SELECT solves this issue
        """
        machine = self._case2.get_value("MACH")
        if "derecho" in machine:
            os.environ["PBS_SELECT"] = "1:ncpus=1:mpiprocs=1:ompthreads=1:mem=230gb:Qlist=cpu:ngpus=0"
            
                      
