#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk
from src.gui.generic.scrolling import ScrollFrame
from src.gui.generic.file_handling import BrowseFrame


class AminoAcidFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        self.scrollFrame = ScrollFrame(self)
        self.scrollFrame.pack(side="left", fill="both")
        self.frame = tk.Frame(self)
        self.frame.pack(side="left", fill="both", expand=1)
        self.scrollFrame.addReloadButton(
            lambda: gui.actOnFunction(
                self.reload,
                gui
            )
        )
        gui.setBrowsableParameter(
            self.frame,
            variable_name_tuple=("SAMPLE_PREP", "AA_FILE_NAME"),
            row=0,
            dialog_title="Select AA file to open",
        )
        self.scrollFrame.setSelectFunction(self.selectFunction)

    def reload(self, gui):
        if gui.prioritization.PARAMETERS["SAMPLE_PREP"]["AA_FILE_NAME"]:
            self.scrollFrame.clear()
            gui.prioritization.loadAADict()
            self.scrollFrame.addMultiple(gui.prioritization.AA_DICT)

    def selectFunction(self, value):
        print(value)
