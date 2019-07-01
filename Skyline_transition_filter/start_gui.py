#!/usr/bin/python2.7


import sys
from gui import *
import skyline_filter


def GUIMain(gui):
    sys.stderr = gui
    sys.stdout = gui
    print("=" * 80)
    for parameter in gui.parameters.values():
        exec("{}".format(parameter))
        exec("skyline_filter.{}".format(parameter))
    skyline_filter.main()
    print("=" * 80)


if __name__ == '__main__':
    gui = GUI(
        "Parse Skyline File",
        [
            GUIParameter(
                name="BASE_PATH",
                info="Test/",
                widget_type=Tkinter.Entry,
            ),
            GUIParameter(
                name="CURATE_PTMS",
                info="1",
                widget_type=Tkinter.Checkbutton,
                default=1
            ),
            GUIParameter(
                name="USE_EXPLICIT_SNAPSHOT_VARIANTS",
                info="1",
                widget_type=Tkinter.Checkbutton,
                default=1
            ),
            GUIParameter(
                name="REMOVE_PRECURSOR_TRANSITIONS",
                info="1",
                widget_type=Tkinter.Checkbutton,
                default=1
            ),
            GUIParameter(
                name="TRANSITION_FILE_NAME",
                info='BASE_PATH + "Peptide PTM Report_2000_irt_30W.csv"',
                widget_type=Tkinter.Entry,
                default='BASE_PATH + '
            ),
            GUIParameter(
                name="DATABASE_FILE_NAME",
                info='BASE_PATH + "db.peff"',
                widget_type=Tkinter.Entry,
                default='BASE_PATH + '
            ),
            GUIParameter(
                name="WINDOW_FILE_NAME",
                info='BASE_PATH + "Var_Wind_30W.csv"',
                widget_type=Tkinter.Entry,
                default='BASE_PATH + '
            ),
            GUIParameter(
                name="DELIMITER",
                info=",",
                widget_type=Tkinter.Entry,
                default=","
            ),
            GUIParameter(
                name="OUTPUT_PATH",
                info='BASE_PATH + "Results/"',
                widget_type=Tkinter.Entry,
                default='BASE_PATH + '
            ),
            GUIParameter(
                name="FILTERED_TRANSITION_FILE_NAME",
                info='OUTPUT_PATH + "filtered_transitions.csv"',
                widget_type=Tkinter.Entry,
                default='OUTPUT_PATH + '
            ),
            GUIParameter(
                name="UNIQUE_TRANSITION_FILE_NAME",
                info='OUTPUT_PATH + "unique_transitions.csv"',
                widget_type=Tkinter.Entry,
                default='OUTPUT_PATH + '
            ),
            GUIParameter(
                name="PRECURSOR_PPM_ERROR",
                info="10",
                widget_type=Tkinter.Entry,
                default=10
            ),
            GUIParameter(
                name="PRODUCT_PPM_ERROR",
                info="10",
                widget_type=Tkinter.Entry,
                default=10
            ),
            GUIParameter(
                name="CHEMICAL_PTMS",
                info='{"N-term": {"Pr": "Fixed"},"K": {"Pr": "Fixed"}}',
                widget_type=Tkinter.Entry,
                default={"N-term": {"Pr": "Fixed"},"K": {"Pr": "Fixed"}}
            ),
            GUIParameter(
                name="MASS_TRANSLATIONS",
                info='{42.0106: "Ac", 98.0368: "Ac", 70.0419: "Bu", 126.0681: "Bu", 68.0262: "Cr", 124.0524: "Cr", 27.9949: "Fo", 84.0211: "Fo", 86.0368: "Hib", 142.063: "Hib", 14.0156: "Me", 28.0313: "Me2", 84.0575: "Me2", 42.047: "Me3", 98.0732: "Me3", 56.0262: "Pr", 112.0524: "Pr", 100.016: "Su", 156.0423: "Su", 114.0429: "Ub", 170.1: "Ub"}',
                widget_type=Tkinter.Entry,
                default={42.0106: "Ac", 98.0368: "Ac", 70.0419: "Bu", 126.0681: "Bu", 68.0262: "Cr", 124.0524: "Cr", 27.9949: "Fo", 84.0211: "Fo", 86.0368: "Hib", 142.063: "Hib", 14.0156: "Me", 28.0313: "Me2", 84.0575: "Me2", 42.047: "Me3", 98.0732: "Me3", 56.0262: "Pr", 112.0524: "Pr", 100.016: "Su", 156.0423: "Su", 114.0429: "Ub", 170.1: "Ub"}
            ),
            GUIParameter(
                name="TRANSITION_FILE_COLUMN_NAMES",
                info='{"pep_mod_seq": "Peptide Modified Sequence Monoisotopic Masses", "pep_mz": "Precursor Mz", "prod_mz": "Product Mz", "prod_ion": "Fragment Ion", "prod_ion_type": "Fragment Ion Type", "protein": "Protein Name"}',
                widget_type=Tkinter.Entry,
                default={"pep_mod_seq": "Peptide Modified Sequence Monoisotopic Masses", "pep_mz": "Precursor Mz", "prod_mz": "Product Mz", "prod_ion": "Fragment Ion", "prod_ion_type": "Fragment Ion Type", "protein": "Protein Name"}
            ),
        ],
        GUIMain
    )
