#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk


class AnnotationFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        for row, variable_name in enumerate(
            [
                "MAX_C13",
                "PPM",
            ]
        ):
            gui.setParameter(
                self,
                variable_name_tuple=("MS", variable_name),
                row=row,
            )
        gui.setBrowsableParameter(
            self,
            variable_name_tuple=("MS", "PRECURSOR_FILE_NAME"),
            row=row + 1,
            dialog_title="Select MGF file to open",
        )
        tk.Button(
            self,
            text="Match candidates to precursors",
            command=lambda: gui.actOnFunctions(
                (
                    gui.loadPrecursorList,
                ), (
                    gui.prioritization.annotatePrecursors,
                )
            )
        ).grid(row=row + 2, column=0, sticky='w')
