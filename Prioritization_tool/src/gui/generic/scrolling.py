#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk
import bisect


class ScrollFrame(tk.Frame):

    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        self.elements = tk.Listbox(self, name="elements")
        self.scrollbar = tk.Scrollbar(self, orient="vertical")
        self.scrollbar.config(command=self.elements.yview)
        self.elements.config(yscrollcommand=self.scrollbar.set)
        self.allow_duplicates = False
        self.elements.pack(side="left", fill="both")
        self.scrollbar.pack(side="left", fill="y")

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

    def addMultiple(self, elements):
        for element in elements:
            self.add(element)

    def moveTo(self, other):
        element = self.get()
        if not other.allow_duplices and element in other:
            return
        other.add(element)
        self.remove(element)

    def setSelectFunction(self, function):
        def onselect(*args):
            function(self.get())
        self.elements.bind('<<ListboxSelect>>', function)

    def addReloadButton(self, function):
        self.reload_button = tk.Button(
            self,
            text="Reload",
            command=function
        )
        for widget in self.children.values():
            widget.pack_forget()
        self.reload_button.pack(side="bottom", fill="x")
        self.elements.pack(side="left", fill="both")
        self.scrollbar.pack(side="left", fill="y")

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
