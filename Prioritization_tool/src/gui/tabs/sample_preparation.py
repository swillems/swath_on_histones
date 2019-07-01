#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk


class SamplePreparationFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        simple_entries = [
            "CLEAVAGE_AAS",
            "MAX_MISCLEAVAGES",
            "MAX_PEPTIDE_LENGTH"
        ]
        for row, variable_name in enumerate(simple_entries):
            gui.setParameter(
                self,
                variable_name_tuple=("SAMPLE_PREP", variable_name),
                row=row,
            )
        browse_entries = [
            ("AA_FILE_NAME", "Select amino acid file to open"),
            ("PTM_FILE_NAME", "Select PTM file to open"),
            ("PROTEIN_FILE_NAME", "Select protein file to open"),
        ]
        for row, (variable_name, dialog_title) in enumerate(browse_entries):
            gui.setBrowsableParameter(
                self,
                variable_name_tuple=("SAMPLE_PREP", variable_name),
                row=row + len(simple_entries),
                dialog_title=dialog_title,
            )
        gui.varfix_ptm_frame = gui.ptmVarFixFrame(self)
        gui.varfix_ptm_frame.grid(row=len(simple_entries + browse_entries), column=0, columnspan=3)
        tk.Button(
            self,
            text="Digest",
            command=lambda: gui.actOnFunctions(
                (
                    gui.loadAADict,
                ), (
                    gui.loadPTMDict,
                ), (
                    gui.loadProteinDict,
                ), (
                    gui.prioritization.generateIsoforms,
                )
            )
        ).grid(row=len(simple_entries + browse_entries) + 1, column=0, sticky='w')
