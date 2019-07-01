#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk


class PTMFrame(tk.Frame):

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        gui.ptm_scroller, gui.PTM_modification_frame = gui.splitScrollFrame(self)
        gui.setBrowsableParameter(
            gui.PTM_modification_frame,
            variable_name_tuple=("SAMPLE_PREP", "PTM_FILE_NAME"),
            row=0,
            dialog_title="Select PTM library file to open",
        )
        tk.Button(
            gui.PTM_modification_frame,
            text="Reload",
            command=gui.updatePTMDict
        ).grid(row=1, column=0, sticky='w')
        tk.Button(
            gui.PTM_modification_frame,
            text="Save",
            command=lambda: gui.actOnFunction(
                gui.prioritization.save,
                gui.prioritization.PARAMETERS["SAMPLE_PREP"]["PTM_FILE_NAME"],
                item_to_save=gui.prioritization.PTM_DICT
            )
        ).grid(row=1, column=1, sticky='w')
        gui.manual_PTM_frame = gui.manualPTMFrame(gui.PTM_modification_frame)
        gui.manual_PTM_frame.grid(row=2, column=0, sticky='w')
        def onselect(event):
            # Note here that Tkinter passes an event object to onselect()
            w = event.widget
            index = int(w.curselection()[0])
            ptm_name = w.get(index)
            gui.current_manual_ptm_fields["Name"].set(ptm_name)
            values = gui.prioritization.PTM_DICT[ptm_name]
            gui.current_manual_ptm_fields["Chemical Formula"].set(values["CF"])
            gui.current_manual_ptm_fields["References"].set(values["DR"])
            gui.current_manual_ptm_fields["Keywords"].set(values["KW"])
            gui.current_manual_ptm_fields["Monoisotopic Mass"].set(values["MM"])
            gui.current_manual_ptm_fields["Target"].set(values["TG"])
            print(gui.prioritization.PTM_DICT[ptm_name])
        gui.ptm_scroller.children["listbox"].bind('<<ListboxSelect>>', onselect)
