#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk


class PrioritizationFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        for row, variable_name in enumerate(
            [
                "MAX_PTMC_SIZE",
                "SEQUENTIAL_SEARCHES",
            ]
        ):
            gui.setParameter(
                self,
                variable_name_tuple=("PRIORITIZATION", variable_name),
                row=row,
            )
        gui.setBrowsableParameter(
            self,
            variable_name_tuple=("PRIORITIZATION", "PRIORITIES_FILE_NAME"),
            row=row + 1,
            dialog_title="Save as priorities file name",
            saveButton=True
        )
        tk.Button(
            self,
            text="Prioritize",
            command=lambda: gui.actOnFunctions(
                (
                    gui.prioritization.calculatePTMCWeights,
                ), (
                    gui.prioritization.prioritizeAllPTMs,
                ), (
                    gui.prioritization.writePrioritizedPTMsToCSV,
                )
            )
        ).grid(row=row + 2, column=0, sticky='w')
