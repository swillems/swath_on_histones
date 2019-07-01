#!/usr/bin/python2.7


from __future__ import print_function
import Tkinter
import tkFileDialog
from collections import OrderedDict
from time import asctime
import sys


class GUIParameter(object):

    def __init__(
        self,
        name,
        widget_type,
        info=None,
        default=None
    ):
        self.name = name
        self.info = info
        self.widget_type = widget_type
        self.default = default
        self.widget = None
        self.value = None

    def __str__(self):
        return "{} = {}".format(
            self.name,
            self.value.get()
        )


class GUI(object):

    def __init__(self, title, parameters, externalRunCommand):
        sys.stderr = self
        sys.stdout = self
        self.root = Tkinter.Tk()
        self.root.winfo_toplevel().title(title)
        self.parameters = OrderedDict(
            (parameter.name, parameter) for parameter in parameters
        )
        self.externalRunCommand = externalRunCommand
        self.createForm()
        self.root.mainloop()

    def createForm(self):
        self.setCustomWidgets()
        self.setDefaultButtons()
        self.setDisplay()

    def setCustomWidgets(self):
        for i, parameter in enumerate(self.parameters.values()):
            self.setCustomWidget(i, parameter)

    def setDefaultButtons(self):
        for i, (name, command) in enumerate(
            [
                # ('Check Parameters', self.checkParameterCommand),
                ('Load Parameters', self.loadParameterCommand),
                ('Save Parameters', self.saveParameterCommand),
                ('Run', self.runCommand),
            ]
        ):
            self.setDefaultButton(name, command, len(self.parameters) + i)

    def setCustomWidget(self, row, parameter):
        Tkinter.Label(
            self.root,
            text=parameter.name
        ).grid(row=row, sticky='w')
        if parameter.widget_type == Tkinter.Entry:
            value = Tkinter.StringVar()
            widget = parameter.widget_type(
                self.root,
                textvariable=value
            )
        elif parameter.widget_type == Tkinter.Checkbutton:
            value = Tkinter.IntVar()
            widget = parameter.widget_type(
                self.root,
                variable=value
            )
        if parameter.default is not None:
            value.set(parameter.default)
        widget.grid(row=row, column=1, sticky='w', )
        info_button = Tkinter.Button(
            self.root,
            text='Info',
            command=lambda info=parameter.info: print(info)
        )
        info_button.grid(row=row, column=2)
        parameter.widget = widget
        parameter.value = value

    def setDefaultButton(self, name, command, row):
        default_button = Tkinter.Button(
            self.root,
            text=name,
            command=command
        )
        default_button.grid(
            row=row,
            column=0,
            columnspan=3,
            sticky='we'
        )

    # def checkParameterCommand(self):
    #     self.triggerActive(False)
    #     # TODO
    #     print("Not implemented")
    #     self.triggerActive(True)

    def loadParameterCommand(self):
        self.triggerActive(False)
        file_name = tkFileDialog.askopenfilename(
            parent=self.root,
            title='Select parameter file'
        )
        if file_name:
            with open(file_name, "rb") as infile:
                for line in infile:
                    if not line.startswith("#") and line.strip():
                        name, value = line.split(" = ")
                        self.parameters[name].value.set(value.rstrip())
            print(
                "Succesfully loaded parameters from {}".format(
                    file_name
                )
            )
        self.triggerActive(True)

    def saveParameterCommand(self):
        self.triggerActive(False)
        file_name = tkFileDialog.asksaveasfilename(
            parent=self.root,
            title='Select parameter file'
        )
        if file_name:
            with open(file_name, "wb") as outfile:
                for parameter in self.parameters.values():
                    outfile.write("{}\n".format(parameter))
            print(
                "Succesfully saved parameters to {}".format(
                    file_name
                )
            )
        self.triggerActive(True)

    def runCommand(self):
        self.triggerActive(False)
        self.externalRunCommand(self)
        self.triggerActive(True)

    def triggerActive(self, active):
        if active:
            state = "normal"
        else:
            state = "disabled"
        for widget in self.root.winfo_children():
            widget.config(state=state)

    def setDisplay(self):
        self.display_indent = 0
        self.display = Tkinter.Text(self.root)
        self.display.grid(
            row=0,
            column=3,
            # rowspan=len(self.parameters) + 4
            rowspan=len(self.parameters) + 3,
            sticky="ns"
        )
        self.display.config(state="disabled")

    def write(self, text):
        self.display.config(state="normal")
        # text = text.rstrip()
        # time = asctime()
        # space = 80 - len(time) - 3
        # chunks = len(text) / space + 1
        # if len(text) % space == 0:
        #     chunks -= 1
        # for i in xrange(chunks):
        #     self.display.insert(
        #         Tkinter.END,
        #         "{} > {}\n".format(
        #             time,
        #             text[space * i: space * (i + 1)]
        #         )
        #     )
        self.display.insert(Tkinter.END, text)
        self.display.update()
        self.display.config(state="disabled")
        # indent_space = "  " * self.display_indent
        # self.display.config(state="normal")
        # self.display.insert(
        #     Tkinter.END,
        #     "{} > {}{}\n".format(asctime(), indent_space, text)
        # )
        # self.display.update()
        # self.display.config(state="disabled")
        # self.display_indent += indent
