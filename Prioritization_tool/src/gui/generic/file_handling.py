#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk
import bisect


class BrowseFrame(tk.Frame):

    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        def browseFile():
            self.toggleActiveStatus(True)
            if saveButton:
                browseFileName = tkFileDialog.asksaveasfilename
            else:
                browseFileName = tkFileDialog.askopenfilename
            file_name = browseFileName(
                parent=parameter_frame,
                title=dialog_title
            )
            if file_name:
                variable.set(file_name)
                # if tkMessageBox.askyesno(
                #     "Load file",
                #     "Do you want to immediately load this file?"
                # ):
                #     self.actOnFunction(loadFileFunction)
            self.toggleActiveStatus(False)
        # def loadFile():
        #     self.actOnFunction(loadFileFunction)
        # def saveFile():
        #     self.actOnFunction(saveFileFunction)
        tk.Label(
            parameter_frame,
            text=variable_name_tuple[-1]
        ).grid(row=row, column=0, sticky='w')
        variable = self.getVariable(variable_name_tuple)
        tk.Entry(
            parameter_frame,
            textvariable=variable,
        ).grid(row=row, column=1, sticky='w')
        options = [("Browse", browseFile)]
        if not browseOnly:
            options.append(("Load", loadFile))
            options.append(("Save", saveFile))
        for column, (text, command) in enumerate(options):
            tk.Button(
                parameter_frame,
                text=text,
                command=command
            ).grid(row=row, column=column + 2, sticky='w')





    def allowDuplicates(self):
        self.allow_duplicates = True

    def getIndex(self):
        return self.elements.curselection()[0]

    def all(self):
        return self.elements.get(0, "end")

    def clear(self):
        self.elements.delete(0, "end")

    def get(self):
        return self.elements.get(self.getIndex)

    def remove(self):
        self.elements.delete(self.getIndex)

    def add(self, element):
        if not self.allow_duplicates and element in self:
            return
        new_index = bisect.bisect_left(self.all(), element)
        self.elements.insert(new_index, element)

    def add_multiple(self, elements):
        for element in elements:
            self.add(element)

    def moveTo(self, other):
        element = self.get()
        if not other.allow_duplices and element in other:
            return
        other.add(element)
        self.remove(element)

    def __contains__(self, element):
        if element in self.all():
            return True
        return False


if __name__ == "__main__":
    import Tkinter as tk
    root = tk.Tk()
    x = ScrollFrame(root)
    x.add(4)
    x.pack(fill="both", expand=True)
    x.add_multiple(range(2, 8, 2))
    x.add_multiple(range(0, 10, 3))
