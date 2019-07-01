#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk


class ProteinFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        gui.protein_scroller, gui.protein_modification_frame = gui.splitScrollFrame(self)
        def onselect(event):
            w = event.widget
            index = int(w.curselection()[0])
            value = w.get(index)
            print(gui.prioritization.PROTEIN_DICT[value])
        gui.protein_scroller.children["listbox"].bind('<<ListboxSelect>>', onselect)
